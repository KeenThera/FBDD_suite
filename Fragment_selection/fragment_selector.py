# This script is designed for chemical library filtration.
# The resulting output can be utilized as a virtural fragment library.

from rdkit import Chem
from filter import Filter

# Function to load molecules from an SDF file
def load_sdf(sdf_filename):
    suppl = Chem.SDMolSupplier(sdf_filename)
    molecules = [mol for mol in suppl if mol is not None]
    return molecules

# Molecular selector function to filter molecules based on criteria
def selector(mol_list):
    selected_list = []
    myfilter = Filter("mol")  # Create a Filter object with "mol" as its initialization parameter
    for mol in mol_list:
        myfilter.load_mol(mol)  # Load the molecule into the filter object
        if myfilter.organic_filter() and myfilter.fragment_filter():
            selected_list.append(mol)  # If it passes the organic and fragment filters, add it to the selected list
    return selected_list

# Function to write selected molecules to an SDF file
def write_sdf(selected_list, sdf_out_filename):
    writer = Chem.SDWriter(sdf_out_filename)
    for mol in selected_list:
        writer.write(mol)

# Main program
if __name__ == "__main__":
    input_sdf_filename = 'input.sdf'  # Input SDF file name
    output_sdf_filename = 'output.sdf'  # Output SDF file name

    # Load molecules from the input SDF file
    molecules = load_sdf(input_sdf_filename)

    # Use the molecular selector to filter molecules
    selected_molecules = selector(molecules)

    # Write the selected molecules to a new SDF file
    write_sdf(selected_molecules, output_sdf_filename)

    print("Chemical library filtering completed. Result saved to", output_sdf_filename)

