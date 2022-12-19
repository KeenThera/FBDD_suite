from utilities import *

class link:
    """
    class to do fragment linking
    """
    
    def __init__(self, mol1, mol2):
        self.mol1 = mol1
        self.mol2 = mol2
        
        mol1_NOHs = Chem.RemoveAllHs(mol1)
        mol2_NOHs = Chem.RemoveAllHs(mol2)
        self.mol1_NOHs = mol1_NOHs
        self.mol2_NOHs = mol2_NOHs
        
        smi1 = Chem.MolToSmiles(mol1_NOHs)
        smi2 = Chem.MolToSmiles(mol2_NOHs)
        self.smi1 = smi1
        self.smi2 = smi2
        
        merged = Chem.CombineMols(mol1, mol2)
        merged_3D = Chem.RemoveAllHs(merged)
        merged_2D = deepcopy(merged_3D)
        AllChem.Compute2DCoords(merged_2D)
        merged_smi = Chem.MolToSmiles(merged_2D)
        self.merged = merged
        self.merged_3D = merged_3D
        self.merged_2D = merged_2D
        self.merged_smi = merged_smi
        
        # the geometry center
        self.c1 = geometry_center(mol1)
        self.c2 = geometry_center(mol2)
        
        #results
        self.out3D_list = []
    
    def set(self, mol1, mol2):
        self.mol1 = mol1
        self.mol2 = mol2
    
    def anchor_analysis(self):
        mol1 = deepcopy(self.mol1_NOHs)
        mol2 = deepcopy(self.mol2_NOHs)
    
        mol1_dist_list = distance_list(mol1, self.c1, self.c2)
        mol2_dist_list = distance_list(mol2, self.c1, self.c2)
    
        mol1_anchor_index = top_min_idx(mol1_dist_list)
        mol2_anchor_index = top_min_idx(mol2_dist_list)
    
        mol1_anchor_pos = mol1.GetConformer().GetAtomPosition(mol1_anchor_index[0])
        mol2_anchor_pos = mol2.GetConformer().GetAtomPosition(mol2_anchor_index[0])
    
        dist = mol1_anchor_pos.Distance(mol2_anchor_pos)
        self.mol1_anchor_index = mol1_anchor_index[0]
        self.mol2_anchor_index = mol2_anchor_index[0]
        print("mol1 anchor", mol1_anchor_index[0])
        print("mol2 anchor", mol2_anchor_index[0])
        self.dist = dist
        
        
    def distance_check(self, cutoff = 5):
        if self.dist >= 5:
            return True
        else:
            return False
    
    def linker_smiles(self, batch=2.45 ):
        dist = self.dist
        length = math.ceil(dist/batch)
        linker_smiles = length*"C"
        linker_mol = Chem.MolFromSmiles(linker_smiles)
        self.linker_length = length
        self.linker_smiles = linker_smiles
        self.linker_mol = linker_mol
        
        
    def do_link(self):
        comp = Chem.MolFromSmiles(self.merged_smi)
        comp = reNumber(comp, self.merged_2D)
        a = self.mol1_anchor_index
        b = self.mol2_anchor_index + self.mol1.GetNumHeavyAtoms()
        print(a, b)
        
        mol3 = Chem.MolFromSmiles(self.linker_smiles)
        
        test = Chem.CombineMols(comp,mol3)
        test=Chem.RWMol(test)
        base = self.mol1.GetNumHeavyAtoms() + self.mol2.GetNumHeavyAtoms()
        c =  base
        d =  self.linker_length - 1 + base
        test.AddBond(a,c,Chem.BondType.SINGLE)
        test.AddBond(b,d,Chem.BondType.SINGLE)
        
        test.GetAtomWithIdx(c).SetNumExplicitHs( 0 )
        test.GetAtomWithIdx(d).SetNumExplicitHs( 0 )
        test.UpdatePropertyCache(strict=False)
        my_test=test.GetMol()
        Chem.SanitizeMol(my_test)
        
        lst=my_test.GetSubstructMatch(self.merged_2D)
        #print(lst)
        set_protected_atom_label(my_test, lst)
        self.init_mol = my_test
        
    
    def mutation(self):
        mol_new = self.init_mol
        lst = mol_new.GetSubstructMatch(self.merged_2D)
        mol_new = set_protected_atom_label(mol_new, lst)
        
        #Chem.SanitizeMol(self.linker_mol)
        transforms = ["[CH2:1][C:2]>>[C:1]([N:2])=O", "[C!H0:1]C[C!H0:2]>>[C:1]1C[C:2]C1",
                     '[C!H0:1]>>[N:1]','[CH2:1]>>[O:1]', '[C:1][C:2]>>[C:1]#[C:2]', 
                      "[C:1][C:2][C:3]>>[C:1]1=[C:2][N:3]N=N1",
                     '[C:1][C:2][C:3][C:4]>>[N:1]1[C:2][C:3][N:4]CC1',
                     '[C:1][C:2][C:3][C:4]>>[N:1]1[C:2][C:3][C:4]CC1']
        generate_mols = []
        for sma in transforms:
            tmp = reaction(sma, (mol_new,))
            generate_mols+=tmp
        
        self.generate_mols = generate_mols
        
    def align(self):
        for generate_mol in self.generate_mols:
            generate_3D = Chem.AddHs(generate_mol, addCoords=True)
            result = AllChem.EmbedMolecule(generate_3D, randomSeed=0xf00d)
            c = generate_3D.GetConformer()
            lst=generate_3D.GetSubstructMatch(self.merged_3D)
            amap=[(lst[i],i)for i in range(self.merged_3D.GetNumAtoms())]
            Chem.rdMolAlign.AlignMol(generate_3D, self.merged_3D, atomMap=amap)
            self.out3D_list.append(generate_3D)
    
    def optimize(self):
        for generate_mol in self.generate_mols:
            generate_3D = Chem.AddHs(generate_mol, addCoords=True)
            result = AllChem.EmbedMolecule(generate_3D, randomSeed=0xf00d)
            c = generate_3D.GetConformer()
            lst=generate_3D.GetSubstructMatch(self.merged_3D)
            print(lst)
            
            for i in range(self.merged_3D.GetNumHeavyAtoms()):
                p = self.merged_3D.GetConformer().GetAtomPosition(i)
                c.SetAtomPosition(lst[i],p)
                
            mp = AllChem.MMFFGetMoleculeProperties(generate_3D)
            ff = AllChem.MMFFGetMoleculeForceField(generate_3D, mp)
            for i in lst:
                ff.MMFFAddPositionConstraint(i, 0, 1.e4)
            ff.Minimize(maxIts=100)
            generate_3D = Chem.RemoveHs(generate_3D)
            generate_3D = Chem.AddHs(generate_3D, addCoords=True)
            
            self.out3D_list.append(generate_3D)
    
    def process(self, method):
        self.anchor_analysis()
        if self.distance_check():
            self.linker_smiles()
            self.do_link()
            self.mutation()

            if method == "align":
                self.align()
            elif method == "opt":
                self.optimize()
            else:
                raise ValueError("Only align and opt supported")
        else:
            raise ValueError("the distance is {} we need anchor distance >= 5".format(self.distance))
            
    def show(self):
        return self.out3D_list
    
    def remove_duplicates(self):
        uniq = set()
        clean_mols = []
        for mol in self.generate_mols:
            smi = Chem.MolToSmiles(Chem.RemoveHs(mol_tuple[0]), isomericSmiles=True, kekuleSmiles=False)
            if smi in uniq:
                pass
            else:
                uniq.add(smi)
                clean_mols.append(mol)
        self.clean_mols = clean_mols

    def save(self, fname):
        w = Chem.SDWriter(fname)
        for out in self.clean_mols:
            w.write(out)
        w.close()

