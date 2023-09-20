import argparse
import os
import subprocess
import pathlib
from pathlib import Path
import time
# molecule manipulation
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem



def protonate_smiles(smiles: str, pH: float) -> str:
    obabel = os.path.join(r"C:\Program Files\OpenBabel-3.1.1", 'obabel.exe')
    # cmd list format raises errors, therefore one string
    cmd = f' "{obabel}" -:"{smiles}" -ismi -ocan -p{pH}'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode('utf-8')
    if cmd_return.returncode != 0:
        raise ValueError('Ligand protonation with OpenBabel failed')
    return output.strip()


def ligPrep(smi, name, workdir):
    protonated = protonate_smiles(smi, 7.4)

    mol = Chem.MolFromSmiles(protonated)
    mol_3D = Chem.AddHs(mol, addCoords=True)
    result = AllChem.EmbedMolecule(mol_3D, randomSeed=0xf00d)

    preparation = MoleculePreparation()
    preparation.prepare(mol_3D)
    #preparation.show_setup()  # optional
    pdbqt_string = preparation.write_pdbqt_string()
    with open(Path(workdir, name+".pdbqt"), 'w') as out:
        out.write(pdbqt_string)


def process_sdf(sdf_file, workdir):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        name = mol.GetProp("_Name")
        smi = Chem.MolToSmiles(mol)
        ligPrep(smi, name, workdir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version',
                    version='%(prog)s version : v 1.0', help='show the FBDDsuit version')
    parser.add_argument("--sdffile", "-sdf", type=str)
    parser.add_argument("--outputdir", "-outdir", type=str)
    args = parser.parse_args()
    #print(args.smiles)
    #print(args.output)
    workdir=Path(args.outputdir)
    workdir.mkdir(parents=True, exist_ok=True)
    print("Starting to prepare the ligand file!")
    start_time = time.time()
    process_sdf(args.sdffile, args.outputdir)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Finished after {execution_time:.2f}s")

if __name__ == '__main__':
    main()
