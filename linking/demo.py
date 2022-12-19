from utilities import *
from linking import *

suppl = Chem.SDMolSupplier('source.sdf')
mylist= [mol for mol in suppl]
mol1 = mylist[0]
mol2 = mylist[1]

my = link(mol1,mol2)
my.process("align")
#my.show()
my.save("demo_linking.sdf")
