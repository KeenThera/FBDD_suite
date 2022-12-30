import __future__
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem.Lipinski import NumRotatableBonds
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem import Crippen


class Filter:

    def __init__(self, choice, input_smiles):
        self.choice = choice
        self.input_smiles = input_smiles
        self.mol = Chem.MolFromSmiles(input_smiles)

    def LogP_filter(self):
        mol = Chem.MolFromSmiles(self.input_smiles)
        mol_log_p = Crippen.MolLogP(mol)
        if 5 >= mol_log_p >= -3:
            return True
        else:
            return False

    def Element_filter(self):
        F_count = 0
        Br_count = 0
        Cl_count = 0
        I_count = 0
        S_count = 0

        for atom in self.mol.GetAtoms():
            symb = atom.GetSymbol()
            if symb == "F":
                F_count += 1
            elif symb == "Br":
                Br_count += 1
            elif symb == "Cl":
                Cl_count += 1
            elif symb == "I":
                I_count += 1
            elif symb == "S":
                S_count += 1
            else:
                pass

        return all([F_count <= 3, Br_count < 3, Cl_count <= 3, I_count <= 1, S_count <= 1])

    def bad_substructure_filter(self):
        mol = self.mol
        bad1 = Chem.MolFromSmarts("[cR2][cR2][cR2]")
        bad2 = Chem.MolFromSmarts("[A][cR2][a][cR2][A]")

        number_of_bad1 = len(mol.GetSubstructMatches(bad1, maxMatches=3))
        if number_of_bad1 > 1:
            return False

        number_of_bad2 = len(mol.GetSubstructMatches(bad2, maxMatches=3))
        if number_of_bad2 > 1:
            return False

        return True

    def Mozziconacci_filter(self):
        '''
        To pass the filter a molecule should be:
            >> of Rotatable bonds: Max 15
            >> of Rings: Max 6
            >> of Oxygens: Min 1
            >> of Nitrogens: Min 1
            >> of Halogens: Max 7
        '''
        mol = self.mol
        halogen = Chem.MolFromSmarts("[*;#9,#17,#35,#53,#85]")
        number_of_halogens = len(mol.GetSubstructMatches(halogen, maxMatches=8))
        if number_of_halogens > 7:
            return False

        oxygen = Chem.MolFromSmarts("[#8]")
        number_of_oxygen = len(mol.GetSubstructMatches(oxygen, maxMatches=2))
        if number_of_oxygen < 1:
            return False

        nitrogen = Chem.MolFromSmarts("[#7]")
        number_of_nitrogen = len(mol.GetSubstructMatches(nitrogen, maxMatches=2))
        if number_of_nitrogen < 1:
            return False

        num_rotatable_bonds = NumRotatableBonds(mol)
        if num_rotatable_bonds > 15:
            return False

        ring_count = Chem.rdmolops.GetSSSR(mol)
        if ring_count > 6:
            return False

        return True

    def Lipinski_filter(self):
        mol = Chem.MolFromSmiles(self.input_smiles)
        violation_counter = 0

        exact_mwt = CalcExactMolWt(mol)
        if exact_mwt >= 500:
            violation_counter = violation_counter + 1

        num_hydrogen_bond_donors = CalcNumHBD(mol)
        if num_hydrogen_bond_donors >= 5:
            violation_counter = violation_counter + 1

        num_hydrogen_bond_acceptors = CalcNumHBA(mol)
        if num_hydrogen_bond_acceptors >= 10:
            violation_counter = violation_counter + 1

        mol_log_p = Crippen.MolLogP(mol)
        if mol_log_p > 5:
            violation_counter = violation_counter + 1

        if violation_counter < 2:
            return True
        else:
            return False

    def Load_fragment_filter(self):

