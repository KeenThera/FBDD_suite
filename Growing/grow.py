author__ = "dalong"
__version__ = "2020.12.23"
# imports
import json
from itertools import combinations
from itertools import islice
import random
import sqlite3
import time

# import chemistry
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdMolDescriptors

# import draws
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG = True

class mutation:
    # output limitation
    # _out_product_smiles=[]

    def __init__(self):
        self.load_reaction()
        self.load_buildingblock("fragments.db")

    def load_reaction(self):
        # read reaction json file
        reaction_dict_1 = {}
        reaction_dict_2 = {}

        with open('reactions.json') as f:
            data = json.load(f)
            for name in data:
                # print(i)
                sma = data[name]['inverse_smarts']
                rxns = rdChemReactions.ReactionFromSmarts(sma)
                # print(rxns.GetNumReactantTemplates())
                if rxns.GetNumReactantTemplates() == 1:
                    reaction_dict_1[name] = rxns
                if rxns.GetNumReactantTemplates() == 2:
                    reaction_dict_2[name] = rxns

        # initialize the reaction dict
        self.reaction_dict_1 = reaction_dict_1
        self.reaction_dict_2 = reaction_dict_2

    def load_transformations(self):
        # read the transformation json into a dict
        transform_dict = {}
        with open('transformations.json', 'r') as f:
            data = json.load(f)
            for item in data:
                name = item['name']
                sma = item['smarts']
                transform_dict[name] = rdChemReactions.ReactionFromSmarts(sma)
        # initialize the transformation dict
        self.transform_dict = transform_dict

    def load_buildingblock(self, fragment_db, mw_max=300, logp_min=-3, logp_max=5, ringcount_max=2, rotB_max=6, num=1000):
        # building block list
        now = time.strftime("%Y-%m-%d-%H_%M_%S",time.localtime(time.time()))
        save = open("building_block_list_" + now + "_log.csv", 'w')
        # building block selection rules
        # molecular weight, logp, ring count, rotatable bond number
        '''
        mw_max = 300
        logp_min = -3
        logp_max = 5
        ringcount_max = 2
        rotB_max = 6
        '''

        # initialize the bb database
        mol_list = []
        conn = sqlite3.connect(fragment_db)
        c = conn.cursor()
        sql_string = '''
        SELECT * FROM fragments where source = {0}
        and MW < {1}
        and LOGP > {2} 
        and LOGP < {3}
        and RINGCOUNT < {4}
        and ROTATABLEBOND < {5}
        ORDER BY RANDOM() limit {6}
        '''

        # source selection, Enamine, Enamine_NP, FCH, PB
        for source in ['Enamine', 'Enamine_NP', 'FCH', 'PB']:
            sql = sql_string.format("'"+source+"'", mw_max, logp_min, logp_max, ringcount_max, rotB_max, num)
            c.execute(sql)
            rs = c.fetchall()
            for row in rs:
                save.write(str(row)+"\n")
                smi = row[1]
                product_id = row[2]
                mol = Chem.MolFromSmiles(smi)
                mol.SetProp("ID", product_id)
                mol.SetProp("MW",str(row[3]))
                mol.SetProp("logP", str(row[4]))
                mol.SetProp("Ring_count", str(row[5]))
                mol.SetProp("Rotatable_bond_count", str(row[6]))
                mol_list.append(mol)

        # initialize the mol_list
        self.mol_list = mol_list
        save.close()

    # set smiles
    def load_mol(self, input_smiles):
        self.input_smiles = input_smiles
        self.mol = Chem.MolFromSmiles(input_smiles)

    # output smiles list
    def output_setting(self, output_number_limit):

        self._out_product_smiles = []
        self.output_number_limit = output_number_limit
        self.cutoff = 1.5

    def remove_dump(self):
        self._out_product_smiles = list(set(self._out_product_smiles))

        return self._out_product_smiles

    # find the matched molecules for specific reaction
    def find_match_mol(self, pattern):
        matched_mol = []

        for mol in self.mol_list:
            try:
                if mol.HasSubstructMatch(pattern):
                    matched_mol.append(mol)
            except Exception as e:
                print("read error!")

        return matched_mol

    def reaction(self, rxn, react):
        products = rxn.RunReactants(react)
        try:
            uniq = set(
                [Chem.MolToSmiles(Chem.RemoveHs(x[0]), isomericSmiles=True, kekuleSmiles=True) for x in products])
            for smi in uniq:
                self._out_product_smiles.append(smi)
        except Exception as e:
            # print(e)
            pass

    # first cycle mutation using one reactant reactions
    def reaction_1(self):
        mol = self.mol
        for name in self.reaction_dict_1:
            rxn = self.reaction_dict_1[name]
            if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                self.reaction(rxn, (mol,))
                self.remove_dump()
                if len(self._out_product_smiles) >= self.output_number_limit * self.cutoff:
                    break
        return self.remove_dump()

    # second cycle mutation using two reactants reactions
    def reaction_2(self):
        mol = self.mol

        for name in self.reaction_dict_2:
            rxn = self.reaction_dict_2[name]

            if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                reactant_list = self.find_match_mol(rxn.GetReactantTemplate(1))
                if len(reactant_list) > 0:
                    for mol2 in reactant_list:
                        reacts = (mol, mol2)
                        self.reaction(rxn, reacts)
                        self.remove_dump()
                        if len(self._out_product_smiles) >= self.output_number_limit * self.cutoff:
                            break
            elif mol.HasSubstructMatch(rxn.GetReactantTemplate(1)):
                reactant_list = self.find_match_mol(rxn.GetReactantTemplate(0))
                if len(reactant_list) > 0:
                    for mol2 in reactant_list:
                        reacts = (mol2, mol)
                        self.reaction(rxn, reacts)
                        self.remove_dump()
                        if len(self._out_product_smiles) >= self.output_number_limit * self.cutoff:
                            break
            else:
                continue
        return self.remove_dump()

    def transformation(self):
        mol = self.mol
        for name in self.transform_dict:
            rxn = self.transform_dict[name]
            if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                self.reaction(rxn, (mol,))
                self.remove_dump()
                if len(self._out_product_smiles) >= self.output_number_limit * self.cutoff:
                    break
        return self.remove_dump()

    def nitrogen_walk(self, num_N=1):
        """
        Perform positional analogue scanning to sequentially replace aromatic cH wth n
        :param mol_in: input molecule
        :param num_N: number of nitrogens to replace in each analogue
        :return: list of analogue molecules
        """
        used = set()
        aromatic_cH = Chem.MolFromSmarts("[cH]")
        match_atms = [x[0] for x in self.mol.GetSubstructMatches(aromatic_cH)]
        n_combos = combinations(match_atms, num_N)
        for combo in n_combos:
            new_mol = Chem.RWMol(self.mol)
            for idx in combo:
                atm = new_mol.GetAtomWithIdx(idx)
                atm.SetAtomicNum(7)
            smi = Chem.MolToSmiles(new_mol)
            if smi not in used:
                used.add(smi)
                self._out_product_smiles.append(smi)

        return self.remove_dump()

    def attach_atom(self, atomic_symbol="F", smarts="[cH]", num_sub=1):
        """
        Perform positional analogue scanning to sequentially add a single atom substituent
        :param mol_in: input molecule
        :param atomic_symbol: symbol for atom to be attached
        :param smarts: smarts defining the position to be substituted
        :param num_sub: number of groups to substitute at each iteration
        :return:
        """
        pt = Chem.GetPeriodicTable()
        atomic_num = pt.GetAtomicNumber(atomic_symbol)
        used = set()
        query = Chem.MolFromSmarts(smarts)
        match_atms = [x[0] for x in self.mol.GetSubstructMatches(query)]
        n_combos = combinations(match_atms, num_sub)
        for combo in n_combos:
            new_mol = Chem.RWMol(self.mol)
            for idx in combo:
                new_idx = new_mol.AddAtom(Chem.Atom(atomic_num))
                new_mol.AddBond(idx, new_idx, order=Chem.rdchem.BondType.SINGLE)
            Chem.SanitizeMol(new_mol)
            smi = Chem.MolToSmiles(new_mol)
            if smi not in used:
                used.add(smi)
                self._out_product_smiles.append(smi)

        return self.remove_dump()

    def add_six_elem_ring(self):
        match = self.mol.GetSubstructMatches(Chem.MolFromSmarts('[*R1]@[*R1]!@[*R1]@[*R1]'))
        # newAtom = random.choice([Chem.rdchem.Atom(x) for x in {'C', 'I', 'N', 'O', 'P', 'S', 'Se', 'Si'}])
        for item in match:
            em = Chem.RWMol(self.mol)
            beg = item[0]
            end = item[-1]
            aid1 = em.AddAtom(random.choice([Chem.rdchem.Atom(x) for x in {'C', 'N', 'O', 'S'}]))
            em.AddBond(beg, aid1, order=Chem.rdchem.BondType.SINGLE)
            aid2 = em.AddAtom(random.choice([Chem.rdchem.Atom(x) for x in {'C', 'N', 'O', 'S'}]))
            em.AddBond(aid1, aid2, order=Chem.rdchem.BondType.SINGLE)
            em.AddBond(aid2, end, order=Chem.rdchem.BondType.SINGLE)
            newMol = em.GetMol()
            try:
                Chem.SanitizeMol(newMol)
                self._out_product_smiles.append(Chem.MolToSmiles(newMol, isomericSmiles=True, kekuleSmiles=True))
            except:
                problems = Chem.DetectChemistryProblems(newMol)
                print(len(problems))

        self.remove_dump()

    def result(self):
        # remove duplicate smiles
        return self.remove_dump()

    def wash_mol(self):
        import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "smi")
        ob_mol = openbabel.OBMol()
        _out_new_product_smiles = []
        for smi in self._out_product_smiles:
            rdkit_mol = Chem.MolFromSmiles(smi)
            if rdkit_mol:
                _out_new_product_smiles.append(smi)
            else:
                obConversion.ReadString(ob_mol, smi)
                obConversion.Convert()
                _out_new_product_smiles.append(obConversion.WriteString(ob_mol).strip())

        self._out_product_smiles = _out_new_product_smiles

    def clean(self):
        self.input_smiles = ''
        self._out_product_smiles = []
        self.output_number_limit = 0
