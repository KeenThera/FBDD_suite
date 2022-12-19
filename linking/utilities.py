from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdMolTransforms
from copy import deepcopy
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import math

def geometry_center(mol):
    XYZ = Point3D(0,0,0)
    for i, atom in enumerate(mol.GetAtoms()):
        XYZ += mol.GetConformer().GetAtomPosition(i)
    num = i+1
    out= XYZ/num
    return out

# calculate the distance between one atom to the plan perpendicular to the two center points
def distance(p4, c1, c2):
    v1 = Point3D.DirectionVector(c1,c2)
    c3=(c1+c2)/2
    d = c3.DotProduct(v1)
    n2 = v1/v1.Length()
    n3=v1*d/v1.LengthSq()
    k = Point3D.DotProduct(n2, p4-n3)
    
    return k

# calculate the distance of every atom from mol to the plan perpendicular to the two center points
def distance_list(mol, c1, c2):
    out = []
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        k=distance(positions, c1, c2)
        out.append(abs(k))
    
    return out

# calculate the index of smallest number of the atom distance list
def top_min_idx(distance_list):  
    Lst = distance_list[:]
    k = 1
    index_k = []
    for i in range(k):
        index_i = Lst.index(min(Lst))
        index_k.append(index_i)
        Lst[index_i] = float('inf')

    #print(index_k)
    for i in range(k):
        print(distance_list[index_k[i]])
    
    return index_k


def reaction(sma, react):
    rxn = Chem.rdChemReactions.ReactionFromSmarts(sma)
    out_product = []
    try:
        products = rxn.RunReactants(react)
        for mol_tuple in products:
            Chem.SanitizeMol(mol_tuple[0])
            out_product.append(mol_tuple[0])
    except Exception as e:
        pass
    
    return out_product

def setIsotype(mol, idx, num):
    atom = mol.GetAtomWithIdx(idx)
    atom.SetIsotope(num)
    
    return mol

def getIsotype(mol):
    c=0
    d=0
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetIsotope() == 3:
            c=i
        if atom.GetIsotope() == 4:
            d=i
    if d == 0:
        print(Chem.MolToSmiles(mol))
    return c,d

def reNumber(mol, query):
    new_mol = deepcopy(mol)
    matched = new_mol.GetSubstructMatch(query)
    new_index = list(matched)
    no_matched = []
    for i in range(len(new_mol.GetAtoms())):
        if i in new_index:
            pass
        else:
            new_index.append(i)
            no_matched.append(i)
    out = Chem.RenumberAtoms(new_mol, new_index)
    return out

def reorgmol(mol):
    newmol=deepcopy(mol)
    newmol=Chem.RemoveHs(newmol)
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp('atom_idx', idx)
    cx_smi = Chem.MolToCXSmiles(newmol)
    newmol = Chem.MolFromSmiles(cx_smi)
    return newmol


def set_protected_atom_label(mol, lst):
    for idx in lst:
        mol.GetAtomWithIdx(idx).SetProp('_protected', '1')
    return mol

def clean_isotype(mol):
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)
    return mol

def minimize(mol):
    mol = Chem.AddHs(mol)
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    ff.Minimize(maxIts=10000)
    
    return mol
