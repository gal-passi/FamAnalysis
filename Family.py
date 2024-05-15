import pandas as pd
import math
import os
import ast
import Patient
import Protein
import Mutation
import utils


class Family:
    DATABASE = "DB/Family"

    def __init__(self, family_id=None, patients=None, use_n1=False):
        """
        constructor for Family class
        :param family_id: if included will ignore patients parameter and load from DB
        :param patients: list of patients ids in the family
        :param use_n1: use mutations only in n-1 members
        """
        if not family_id and not patients:
            raise ValueError("Usage: must include either family id or list of patients")
        if family_id:
            with open(os.path.join(self.DATABASE, family_id + '.txt')) as file:
                patients = file.read().split('\n')[:-1]
        with open(os.path.join(self.DATABASE, family_id + '_recurring' + '.txt')) as file:
            recurring = {utils.create_mutation(m) for m in ast.literal_eval(file.read())}
        self.id = family_id
        self._patients = {patient: Patient.Patient(patient) for patient in patients}
        self.size = len(patients)
        self._data = pd.read_csv(os.path.join(self.DATABASE, family_id + '_data' + '.csv')) if not use_n1 else \
            pd.read_csv(os.path.join(self.DATABASE, family_id + '_data_n1' + '.csv'))
        self.n_records = len(self._data)
        self.recurring = recurring


    def proteins(self):
        for p_name in self._data['protein'].to_list():
            yield Protein.Protein(ref_name=p_name)

    def mutations(self):
        for m_name, p_name in zip(self._data['variant'].to_list(), self._data['protein'].to_list()):
            yield Mutation.Mutation(m_name, p_name, load_only=True)

    def members(self):
        for p in self._patients.values():
            yield p

    def multiple_recurrence(self):
        """
        returns mutations in at least one other family (all members)
        """
        with open(os.path.join(self.DATABASE, 'n_2_recurring' + '.txt')) as file:
            tmp = {utils.create_mutation(m) for m in ast.literal_eval(file.read())}
        return tmp & self.recurring

    def recurring_pairs(self):
        first, res = True, set()
        for p in self.members():
            if first:
                res |= p.get_pairs()
                first = False
            else:
                res &= p.get_pairs()
        return res

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

    def top_rank(self, thr=3):
        """
        returns top ranking mutations as set of Mutation obj
        Notice does not preform rank calculations
        """
        df = self._data[self._data['rank'] >= thr]
        raw = list(zip(df['variant'].to_list(), df['protein'].to_list()))
        return set(utils.generate_mutations(raw))

'''
    # TODO very simple implemetation can add many features
    def intersect_all(self, top_n=None, inteface_cutoff=4, interface_mode='bound', n_members=None, in_literature=False):
        """
        return recurring top mutations over family members usimg Patients.top_all()
        :param top_n: number of top mutation to quary from each patient default is 5
        :param inteface_cutoff: cutoff for Patient.top_interface
        :param interface_mode: mode for Patient.top_interface
        :param in_literature: for Patient.top_all() function
        :param n_members: number of patients in family with a common mutation to return currently only the family size
        """
        first, res = True, set()
        for p in self.members():
            if first:
                res |= p.top_all(mode='union', n=top_n, interface=inteface_cutoff,
                                 interface_mode=interface_mode, in_literature=in_literature)
                first = False
            else:
                res &= p.top_all(mode='union', n=top_n, interface=inteface_cutoff,
                                 interface_mode=interface_mode, in_literature=in_literature)
        return res
'''