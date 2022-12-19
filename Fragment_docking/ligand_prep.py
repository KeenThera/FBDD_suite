from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import subprocess

def protonate_smiles(smiles: str, pH: float) -> str:

    obabel = os.path.join(r"C:\Program Files\OpenBabel-3.1.1", 'obabel.exe')
    # cmd list format raises errors, therefore one string
    cmd = f' "{obabel}" -:"{smiles}" -ismi -ocan -p{pH}'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode('utf-8')
    if cmd_return.returncode != 0:
        raise ValueError('Ligand protonation with OpenBabel failed')

    return output.strip()
  
  
  
def ligPrep(smi):
    protonated = protonate_smiles(smi, 7.4)
    
    mol = Chem.MolFromSmiles(protonated)
    mol_3D=Chem.AddHs(mol, addCoords=True)
    result = AllChem.EmbedMolecule(mol_3D,randomSeed=0xf00d)
    
    preparator = MoleculePreparation()
    preparator.prepare(mol_3D)
    preparator.show_setup() # optional
    pdbqt_string = preparator.write_pdbqt_string()
    
    return pdbqt_string
  
  
  
