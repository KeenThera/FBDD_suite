# imports
import os
from itertools import islice
import json
import copy

# import chemistry
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS


class crossover:

    def __init__(self):
        self.result = []

    def loadSmiles(self, smi_1, smi_2):
        self.mol_1 = Chem.MolFromSmiles(smi_1)
        self.mol_2 = Chem.MolFromSmiles(smi_2)

    def reNumber(self, mol, query):
        new_mol = copy.deepcopy(mol)
        matched = new_mol.GetSubstructMatch(query)
        new_index = list(matched)
        no_matched = []
        for i in range(len(new_mol.GetAtoms())):
            if i in new_index:
                pass
            else:
                new_index.append(i)
                no_matched.append(i)
        return Chem.RenumberAtoms(new_mol, new_index)

    def addIsolabel(self, mol, query, base):
        mol = self.reNumber(mol, query)
        matched_num = len(query.GetAtoms())
        atom_num = len(mol.GetAtoms())
        atom_index = [i for i in range(atom_num)]
        matched = [i for i in range(matched_num)]
        no_matched = [i for i in atom_index if i not in matched]

        pairs = []
        break_atoms = set()
        for i in matched:
            atom = mol.GetAtomWithIdx(i)
            neighbors = atom.GetNeighbors()
            for n_atom in neighbors:
                if n_atom.GetIdx() in no_matched:
                    pairs.append((i, n_atom.GetIdx()))
                    break_atoms.add(i)
                    break_atoms.add(n_atom.GetIdx())
                    print(i, n_atom.GetIdx())

        for i in break_atoms:
            atom = mol.GetAtomWithIdx(i)
            atom.SetIsotope(atom.GetIsotope() + base)
            atom.SetProp('atomLabel', str(atom.GetIsotope() + atom.GetIdx()))

        return mol, break_atoms

    def MCSGen(self):
        result = rdFMCS.FindMCS([self.mol_1, self.mol_2], atomCompare=rdFMCS.AtomCompare.CompareElements,
                                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                matchValences=False, ringMatchesRingOnly=True, completeRingsOnly=False)
        return result

    def cutoff(self):
        results = self.MCSGen()
        self.query = results.queryMol

        mol_1_atom_num = len(self.mol_1.GetAtoms())
        mol_2_atom_num = len(self.mol_2.GetAtoms())
        query_atom_num = len(self.query.GetAtoms())

        sub_atom_fraction_1 = query_atom_num / mol_1_atom_num
        sub_atom_fraction_2 = query_atom_num / mol_2_atom_num

        flag_1 = 0.3 < sub_atom_fraction_1 < 0.8
        flag_2 = 0.3 < sub_atom_fraction_2 < 0.8

        if flag_1 and flag_2:
            return True
        else:
            return False

    def updateIsotope(self):
        mol_1 = copy.deepcopy(self.mol_1)
        mol_2 = copy.deepcopy(self.mol_2)

        my_mol_1, break_atom_1 = self.addIsolabel(self.mol_1, self.query, 1000)
        my_mol_2, break_atom_2 = self.addIsolabel(self.mol_2, self.query, 2000)
        break_atom_1.update(break_atom_2)
        new_break_atoms = break_atom_1

        for i in new_break_atoms:
            for mol in [my_mol_1, my_mol_2]:
                a = mol.GetAtomWithIdx(i)
                if a.GetIsotope() < 1000:
                    a.SetIsotope(3000 + i)
        self.mol_1 = my_mol_1
        self.mol_2 = my_mol_2

    def check_mol(self, mol):
        # check the problem
        problems = Chem.DetectChemistryProblems(mol)
        if len(problems) >= 1:
            return False
        # remove one mol object with two molecules
        frag_num = len(Chem.GetMolFrags(mol))
        if frag_num > 1:
            return False

        if mol.GetNumAtoms() < 10:
            return False

        return True

    def cross(self):
        if not self.cutoff():
            return None
        self.updateIsotope()
        mol_1 = copy.deepcopy(self.mol_1)
        mol_2 = copy.deepcopy(self.mol_2)

        mod_mol = Chem.ReplaceSubstructs(mol_1,
                                         self.query,
                                         mol_2,
                                         replaceAll=True)
        for mol in mod_mol:
            modi = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if len(modi) < 2:
                if self.check_mol(mol):
                    self.result.append(mol)
                continue

            ed = Chem.RWMol(mol)
            pair = []
            for i in ed.GetAtoms():
                if 900 < i.GetIsotope() < 1500:
                    pair.append(i.GetIdx())
                if i.GetIsotope() > 2500:
                    pair.append(i.GetIdx())
            print(pair)
            if len(pair) == 2:
                ed.AddBond(pair[0], pair[1], order=Chem.rdchem.BondType.SINGLE)
                res = ed.GetMol()
                if self.check_mol(res):
                    self.result.append(res)
            else:
                pass

    def cleanIsotope(self):
        if len(self.result) <= 0:
            return None
        else:
            for mol in self.result:
                for a in mol.GetAtoms():
                    a.SetIsotope(0)
                    a.SetProp('atomLabel', str(a.GetIdx()))

            return True

    def clean(self):
        self.mol_1 = Chem.MolFromSmiles("")
        self.mol_2 = Chem.MolFromSmiles("")
        self.query = Chem.MolFromSmiles("")
        self.result = []
