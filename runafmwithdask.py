from pathlib import Path
import os
import glob
import json

from dask.distributed import Client
import dask.dataframe as dd
from dask import delayed, compute


from Protein import Protein

AFM_HEADER = 3
AFM_COL_NAMES = ["uniprot_id", "protein_variant", "am_pathogenicity"]
AFM_LOC = "./DB/AFM/AlphaMissense_aa_substitutions.tsv"
DB_PATH = Path("./DB/")
PROTEINS = 'proteins_30_fams'
PROTEIN_PATH = DB_PATH / PROTEINS

def all_proteins():
    """
    :return: iterable of all Protein objects in DB
    """
    for path in glob.glob(os.path.join(PROTEIN_PATH, '*')):
        prot_name = os.path.basename(path)
        yield Protein(ref_name=prot_name)


# def all_mutations():
#     for prot in all_proteins():
#         for mutation in prot.generate_mutations():
#             yield mutation

def all_mutations():
    """Yields only necessary mutation info to avoid large memory usage."""
    for prot in all_proteins():
        for mutation in prot.generate_mutations():
            yield {
                "protein_name": mutation.protein_name,
                "name": mutation.name,
                "variant": f"{mutation.origAA}{mutation.loc}{mutation.changeAA}",
            }

def partition_and_save():
    """
    Partitions the AlphaMissense dataset by the first letter of 'protein_variant' and saves them separately.
    This step should be done once as preprocessing.
    """
    df = dd.read_csv(
        AFM_LOC, sep='\t', header=AFM_HEADER,
        usecols=AFM_COL_NAMES
    )

    df['first_letter'] = df['protein_variant'].str[0]

    df.to_parquet('afm_partitions', partition_on=['first_letter'], engine='pyarrow')


def afm_score(protein_name, name, variant):
    """
    Searches for AlphaMissense score using pre-partitioned dataset (by first letter of variant).
    
    :param variant: str, mutation in [AA][location][AA] format
    :return: float AlphaMissense score if found, else None
    """
    first_letter = variant[0]  # Get first letter
    partition_path = f"afm_partitions/first_letter={first_letter}"  # Path to partitioned data

    # Load only the relevant partition
    try:
        data = dd.read_parquet(partition_path, engine='pyarrow')
    except FileNotFoundError:
        return None  # If no partition exists, return None

    # Perform filtering
    score = data[data['protein_variant'] == variant]['am_pathogenicity'].compute()

    # Return first match or None if not found
    return {'protein_name': protein_name, 'name': name, 'afm_score': float(score.iloc[0]) if not score.empty else ''}

# if __name__ == "__main__":
#     tasks = list(all_mutations())
#     client = Client(n_workers=4)
#     print(client.dashboard_link)
#     # partition_and_save() 
#     computaion = [delayed(afm_score(variant=f"{mutation.origAA}{mutation.loc}{mutation.changeAA}", protein_name=mutation.protein_name, name=mutation.name)) for mutation in tasks]
#     results = compute(*computaion)
#     with open('resultsdaskfull.csv', 'w') as f:
#         for result in results:
#             f.write(f"{result['protein_name']},{result['name']},{result['afm_score']}\n")
#     client.shutdown()

def afm_score_batch(first_letter, mutations):
    """
    Searches AlphaMissense scores for a batch of mutations using a single partition load.

    :param first_letter: str, the first letter of mutations
    :param mutations: list of mutation objects
    :return: list of dicts with {protein_name, name, afm_score}
    """
    partition_path = f"afm_partitions/first_letter={first_letter}"  # Path to partitioned data

    # Load only the relevant partition
    try:
        data = dd.read_parquet(partition_path, engine='pyarrow')
    except FileNotFoundError:
        return [{'protein_name': mut['protein_name'], 'name': mut['name'], 'afm_score': ''} for mut in mutations]  # If no partition exists

    # Create a list of all variants we need to search
    variant_list = [mut['variant'] for mut in mutations]

    # Filter only the needed variants
    filtered_data = data[data['protein_variant'].isin(variant_list)].compute()

    # Create a lookup dictionary for fast retrieval
    score_dict = dict(zip(filtered_data['protein_variant'], filtered_data['am_pathogenicity']))

    # Return results, looking up the score if it exists
    results = []
    for mut in mutations:
        score = score_dict.get(mut['variant'], '')
        results.append({'protein_name': mut['protein_name'], 'name': mut['name'], 'afm_score': score})

    return results

if __name__ == "__main__":
    tasks = list(all_mutations())
    print(f"Total mutations: {len(tasks)}")
    client = Client(n_workers=4)
    print(client.dashboard_link)
    # partition_and_save() 
    # Group mutations by the first letter of their variant
    grouped_mutations = {}
    for mutation in tasks:
        first_letter = mutation['variant'][0]  # Get the first letter
        if first_letter not in grouped_mutations:
            grouped_mutations[first_letter] = []
        grouped_mutations[first_letter].append(mutation)

    # Create Dask computations for each group
    computations = [
        delayed(afm_score_batch)(first_letter, mutations)
        for first_letter, mutations in grouped_mutations.items()
    ]

    # Compute all results in parallel
    results_list = compute(*computations)

    # Flatten the results list (since we return a list per batch)
    results = [item for sublist in results_list for item in sublist]

    # Save results
    with open('resultsafmwithdask.csv', 'w') as f:
        for result in results:
            f.write(f"{result['protein_name']},{result['name']},{result['afm_score']}\n")

    profile_data = client.profile()
    with open("full_dask_profile.json", "w") as f:
        json.dump(profile_data, f, indent=4)

    client.shutdown()