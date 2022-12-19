def protein_preparation(pdbfile: str, outpdbqtfile: str):
    python_exe = r"C:\Users\shien\miniconda3\python.exe"

    protein_prep = os.path.join(r"C:\Users\shien\miniconda3\Lib\site-packages\AutoDockTools_py3-1.5.7.post1-py3.9.egg\AutoDockTools\Utilities24", 'prepare_receptor4.py')

    # cmd list format raises errors, therefore one string
    cmd = f' "{python_exe}" "{protein_prep}" -r "{pdbfile}" -A bonds_hydrogens -e -o "{outpdbqtfile}" '
    cmd_return = subprocess.run(cmd, capture_output=False, shell=True)
    #output = cmd_return.stdout.decode('utf-8')
    print(cmd)

    if cmd_return.returncode != 0:
        raise ValueError('Protein prep  failed')

    return True
 

def get_center(fname):
    from statistics import mean
    X=[]
    Y=[]
    Z=[]
    with open(fname) as inf:
        for line in inf:
            if line.startswith('HETATM'):
                clean = line.strip()
                #print(clean.split())
                X.append(float(clean.split()[6]))
                Y.append(float(clean.split()[7]))
                Z.append(float(clean.split()[8]))
    
    x=mean(X)
    y=mean(Y)
    z=mean(Z)
    
    return x,y,z
 
def config_gen():
  pass
