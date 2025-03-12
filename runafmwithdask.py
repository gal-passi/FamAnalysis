from pathlib import Path
import os
import glob
import json
import time

from dask.distributed import Client
import dask.dataframe as dd
from dask import delayed, compute


from Protein import Protein

AFM_HEADER = 3
AFM_COL_NAMES = ["uniprot_id", "protein_variant", "am_pathogenicity"]
AFM_LOC = "./DB/AFM/AlphaMissense_aa_substitutions.tsv"
DB_PATH = Path("./DB/")
PROTEINS = "proteins_30_fams"
PROTEIN_PATH = DB_PATH / PROTEINS


def all_proteins():
    """
    :return: iterable of all Protein objects in DB
    """
    # test_cases = ["A3GALT2", "A4GALT", "AADACL4", "A1CF", "A2M"]
    # for prot_name in test_cases:
    for path in glob.glob(os.path.join(PROTEIN_PATH, "*")):
        prot_name = os.path.basename(path)
        yield Protein(ref_name=prot_name)


def all_mutations(computed=[], use_alias=False):
    """Yields only necessary mutation info to avoid large memory usage."""
    for prot in all_proteins():
        for mutation in prot.generate_mutations():
            if mutation.name not in computed:
                uniprot_id = mutation.protein.Uid
                uniprot_ids = [uniprot_id]
                if use_alias:
                    reviewed_uids = mutation.protein.all_uids()["reviewed"]
                    if isinstance(reviewed_uids, str):
                        reviewed_uids = [reviewed_uids]
                    for rid in reviewed_uids:
                        if rid not in uniprot_ids:
                            uniprot_ids.append(rid)
                uniprot_ids = sorted(uniprot_ids)
                yield {
                    "protein_name": mutation.protein_name,
                    "name": mutation.name,
                    "variant": f"{mutation.origAA}{mutation.loc}{mutation.changeAA}",
                    "uniprot_ids": list(uniprot_ids),
                }


def partition_and_save():
    """
    Partitions the AlphaMissense dataset by the first letter of 'protein_variant' and saves them separately.
    This step should be done once as preprocessing.
    """
    df = dd.read_csv(AFM_LOC, sep="\t", header=AFM_HEADER, usecols=AFM_COL_NAMES)

    df.to_parquet("afm_partitions", partition_on=["uniprot_id"], engine="pyarrow")


def afm_score_batch(uniprot_id, mutations):
    """
    Searches AlphaMissense scores for a batch of mutations using a single partition load.

    :param uniprot_id: str, the first letter of mutations
    :param mutations: list of mutation objects
    :param results: dict, shared results dictionary to store scores
    :return: None (modifies results in-place)
    """
    results = {}
    partition_path = f"afm_partitions/uniprot_id={uniprot_id}"
    try:
        data = dd.read_parquet(partition_path, engine="pyarrow")
    except FileNotFoundError:
        return results

    # Compute the dataset to load it into memory
    data = data.compute()

    for mut in mutations:
        if mut["name"] in results:
            continue  # Skip if entry already exists

        filtered_data = data[data["protein_variant"] == mut["variant"]]
        score = (
            float(filtered_data["am_pathogenicity"].iloc[0])
            if not filtered_data.empty
            else ""
        )

        if score:
            results[mut["name"]] = {
                "uniprot_id": uniprot_id,
                "protein_name": mut["protein_name"],
                "afm_score": score,
            }
    return results


def process_mutations(tasks):
    """
    Processes a batch of mutations using Dask.
    """
    grouped_mutations = {}

    for mutation in tasks:
        for uid in mutation["uniprot_ids"]:
            uniprot_id = uid
            if uniprot_id not in grouped_mutations:
                grouped_mutations[uniprot_id] = []
            grouped_mutations[uniprot_id].append(mutation)

    computations = [
        afm_score_batch(uniprot_id, mutations)
        for uniprot_id, mutations in grouped_mutations.items()
    ]
    results_list = compute(*computations)
    # Example: ({}, {}, {'Q127H': {'uniprot_id': 'U3KPV4', 'protein_name': 'A3GALT2', 'afm_score': 0.2044}})
    results = {}

    for batch_result in results_list:
        for key, value in batch_result.items():
            if key not in results:
                results[key] = value
    return results


if __name__ == "__main__":

    # start_time_partition = time.time()

    client = Client(n_workers=2)
    print(client.dashboard_link)

    # partition_and_save()
    
    # end_time_partition = time.time()

    # print(
    #     f"Time to generate partitions: {end_time_partition - start_time_partition:.4f} seconds"
    # )


    start_time_mutations = time.time()
    tasks = list(all_mutations())
    end_time_mutations = time.time()

    print(f"Total mutations: {len(tasks)}")
    print(
        f"Time to generate mutations: {end_time_mutations - start_time_mutations:.4f} seconds"
    )

    start_time_dask = time.time()

    results = process_mutations(tasks)

    with open("resultsafmwithdask.csv", "w") as f:
        f.write("uids,protein_name,name,afm_score\n")
        for name, result in results.items():
            f.write(
                f"{result['uniprot_id']},{result['protein_name']},{name},{result['afm_score']}\n"
            )

    end_time_dask = time.time()
    print(
        f"Time to run all Dask computations(without alias): {end_time_dask - start_time_dask:.4f} seconds"
    )

    tasks = list(all_mutations(computed=results.keys(), use_alias=True))
    start_time_dask = time.time()

    results = process_mutations(tasks)

    with open("resultsafmwithdask.csv", "a") as f:
        for name, result in results.items():
            f.write(
                f"{result['uniprot_id']},{result['protein_name']},{name},{result['afm_score']}\n"
            )
    end_time_dask = time.time()
    print(
        f"Time to run all Dask computations (with alias): {end_time_dask - start_time_dask:.4f} seconds"
    )

    client.shutdown()
