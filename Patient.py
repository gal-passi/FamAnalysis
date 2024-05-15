import Analyze
import os
import pandas as pd
import utils
import Mutation
import Protein
import math
import pickle


class Patient:
    Patient_db = "DB/Patients"

    def __init__(self, id):
        """
        constructor for patient object
        :param id: str identification of patient
        """

        self.id = id
        self._analyzer = Analyze.ProteinAnalyzer()
        if not os.path.exists(f"{self.Patient_db}/{id}.csv"):
            raise NameError("make sure a csv file for the given patient was made\n"
                            "format: protein,variant,bert,firm,eve_score,eve_prediction,interface")
        self._data = pd.read_csv(f"{self.Patient_db}/{id}.csv")
        self._data['firm'] = (pd.to_numeric(self._data['firm'], errors='coerce').fillna(0))
        self.n_records = len(self._data)

    def proteins(self):
        for p_name in self._data['protein'].to_list():
            yield Protein.Protein(ref_name=p_name)

    def mutations(self):
        for m_name, p_name in zip(self._data['variant'].to_list(), self._data['protein'].to_list()):
            yield Mutation.Mutation(m_name, p_name)

    def top_bert(self, n=None, in_literature=False, return_raw=False):
        """
        returns most harmful mutations by bert score
        :param n: number of mutations to return
        :in_literature: will return only mutation of proteins mentioned in litriture
        :return_raw: return as [(protein_name, mutation_name, bert_score),]
        :return: [Mutations]
        """
        n = math.ceil(self.n_records * 0.1) if not n else n
        df = self._data.nsmallest(n, 'bert')
        if return_raw and not in_literature:
            return list(zip(df['variant'].to_list(), df['protein'].to_list(), df['bert'].to_list()))
        if return_raw and in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            mentioned_mutations = [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]
            return [(mutation.extended_description, mutation._protein.name, mutation.bert_score)
                    for mutation in mentioned_mutations]
        if not return_raw and not in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return list(utils.generate_mutations(raw))
        else:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]

    def top_eve(self, n=5, in_literature=False, return_raw=False):
        """
        returns most harmful mutations by eve score will filter by non-Benign mutation by eve_prediction
        :param n: number of mutations to return
        :in_literature: will return only mutation of proteins mentioned in litriture
        :return_raw: [(protein_name, mutation_name, eve_score, eve_prediction),]
        :return: [Mutations]
        """
        n = math.ceil(self.n_records * 0.1) if not n else n
        df = self._data[self._data['eve_score'] > -1]
        #df = df[self._data['eve_prediction'] != 'Benign']
        df = df.nlargest(n, 'eve_score')
        if return_raw and not in_literature:
            return list(zip(df['variant'].to_list(), df['protein'].to_list(), df['eve_score'].to_list(), df['eve_prediction'].to_list()))
        if return_raw and in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            mentioned_mutations = [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]
            return [(mutation.extended_description, mutation._protein.name, mutation.eve_score, mutation.eve_prediction)
                    for mutation in mentioned_mutations]
        if not return_raw and not in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return list(utils.generate_mutations(raw))
        else:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]

    def top_firm(self, n=None, in_literature=False, return_raw=False):
        """
        returns most harmful mutations by firm score
        :param n: number of mutations to return
        :param in_literature: will return only mutation of proteins mentioned in litriture
        :param return_raw: [(protein_name, mutation_name, firm_score),]
        :return: [Mutations]
        """
        n = math.ceil(self.n_records * 0.1) if not n else n
        df = self._data[self._data['firm'] > -1]
        df = df.nlargest(n, 'firm')
        if return_raw and not in_literature:
            return list(zip(df['variant'].to_list(), df['protein'].to_list(), df['firm'].to_list()))
        if return_raw and in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            mentioned_mutations = [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]
            return [(mutation.extended_description, mutation._protein.name, mutation.firm_score)
                    for mutation in mentioned_mutations]
        if not return_raw and not in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return list(utils.generate_mutations(raw))
        else:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]

    def top_interface(self, mode='bound', n=4, in_literature=False, return_raw=False):
        """
        returns most harmful mutations by intersection score. in bound mode return all results < n.
        in top_n mode return top n mutations
        :param mode: 'bound'|'top_n'
        :param n: number of mutations to return in top_n mode
        :param in_literature: will return only mutation of proteins mentioned in litriture
        :param return_raw: [(protein_name, mutation_name, interface_score),]
        :return: [Mutation]
        """
        df = self._data[self._data['interface'] != 'inf']
        if mode == 'top_n':
            df = df.nsmallest(n, 'interface')
        else:
            df = df[self._data['interface'] <= n]

        if return_raw and not in_literature:
            return list(zip(df['variant'].to_list(), df['protein'].to_list(), df['interface'].to_list()))
        if return_raw and in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            mentioned_mutations = [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]
            return [(mutation.extended_description, mutation._protein.name, mutation.min_interface())
                    for mutation in mentioned_mutations]
        if not return_raw and not in_literature:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return list(utils.generate_mutations(raw))
        else:
            raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
            return [mutation for mutation in utils.generate_mutations(raw) if mutation.mentioned]

    def in_literature(self):
        return [mutation for mutation in self.mutations() if mutation.mentioned]


    def top_all(self, mode='union', n=None, interface=4, interface_mode='bound', in_literature=False, eve_n=5):
        """
        returns a set of mutation from all categories.
        :param mode: 'union'|'intersect'
        :param interface: cutoff for interface
        :param interface_mode: 'bound'|'top_n' default is bound
        :param in_literature: return only entities in literature
        :param n: number of mutations to return from each catagory if none take top 10%
        """
        n = math.ceil(self.n_records * 0.1) if not n else n
        if mode == 'union':
            return set(self.top_firm(n=n, in_literature=in_literature)) | \
                   set(self.top_bert(n=n, in_literature=in_literature)) | \
                   set(self.top_eve(n=eve_n, in_literature=in_literature)) |  \
                   set(self.top_interface(n=interface, mode=interface_mode, in_literature=in_literature))
        if mode == 'intersect':
            return set(self.top_firm(n=n, in_literature=in_literature)) & \
                   set(self.top_bert(n=n, in_literature=in_literature)) & \
                   set(self.top_eve(n=eve_n, in_literature=in_literature)) &  \
                   set(self.top_interface(n=interface, mode=interface_mode, in_literature=in_literature))

    def in_n_params(self, n_params, top_n=None, interface=4, interface_mode='bound', in_literature=False, eve_n=5):
        """
        returns a list of mutation found in n_params out of 4
        :param n_params: number of parameters in intersection
        :param top_n: number of mutations in each parameter default is 10%
        :param interface: cutoff for interface function default is 4
        :param interface_mode: 'bound'|'top_n' default is bound
        :param in_literature: will only return mutations found in litriture
        :return: [Mutation]
        """

        if n_params == 1:
            return self.top_all(n=top_n, interface=interface, interface_mode=interface_mode, in_literature=in_literature)
        if n_params == 4:
            return self.top_all(n=top_n, mode='intersect', interface=interface, interface_mode=interface_mode, in_literature=in_literature)
        firm, bert, eve, inter = set(self.top_firm(n=top_n, in_literature=in_literature)), \
                                 set(self.top_bert(n=top_n, in_literature=in_literature)), \
                                 set(self.top_eve(n=eve_n, in_literature=in_literature)), \
                                 set(self.top_interface(n=interface, mode=interface_mode, in_literature=in_literature))
        if n_params == 2:
            return (firm & (bert | eve | inter)) | (bert & (firm | eve | inter)) | \
                   (eve & (firm | bert | inter)) | (inter & (firm | bert | eve))
        if n_params == 3:
            return (firm & bert & eve) | (firm & bert & inter) | (firm & eve & inter) | (bert & eve & inter)

    def get_pairs(self):
        """
        return pairs of damaged mutations
        :return: {(prot1, prot2)}
        """
        with open(os.path.join(self.Patient_db, 'patients_pairs.txt'), 'rb') as file:
            pairs = pickle.load(file)
        if self.id in pairs:
            return set(utils.generate_pairs(pairs[self.id]))
        return self._analyzer.find_mutation_pairs(self)  # TODO NOT IMPLEMENTED YET

    def top_rank(self, thr = 3):
        """
        returns top ranking mutations as set of Mutation obj
        Notice does not preform rank calculations
        """
        df = self._data[self._data['rank'] >= thr]
        raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
        return set(utils.generate_mutations(raw))