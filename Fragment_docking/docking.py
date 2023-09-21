import subprocess
from pathlib import Path
import argparse
import time


class Dock:
    def __init__(self, vina_exe, receptor_file, ligand_file, config_filename):
        self.vina_exe = vina_exe
        self.receptor_file = receptor_file
        self.ligand_file = ligand_file
        self.config_filename = config_filename

    def do_dock(self, out_filename):
        cmd = "{0} --ligand {1} --receptor {2} --config {3} --out {4}".format(self.vina_exe, self.ligand_file,
           self.receptor_file, self.config_filename, out_filename)
        #print(cmd)
        cmd_return = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output: str = cmd_return.stdout.decode('utf-8')

        return output


def process(workdir):
    p = Path(workdir, "ligands")
    for ligand_filename in p.glob('*.pdbqt'):
        ligand_filename = str(ligand_filename)
        #print(ligand_filename)
        dock = Dock("vina.exe", "receptor.pdbqt", ligand_filename, "config.txt")
        ligand_filebase = ligand_filename.split('.')[0]
        dock.do_dock(ligand_filebase + "_out.pdbqt")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version',
                    version='%(prog)s version : v 1.0', help='show the FBDDsuit version')
    parser.add_argument("--workdir", "-dir", type=str)
    args = parser.parse_args()
    print("starting docking!")
    start_time = time.time()
    process(args.workdir)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Finished after {execution_time:.2f}s")

if __name__ == '__main__':
    main()
