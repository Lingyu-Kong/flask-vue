from ase import Atoms
from pymatgen.core import Lattice, Structure

def Atoms2Structure(atoms: Atoms):
    lattice = atoms.get_cell().tolist()  
    species = atoms.get_chemical_symbols()  
    positions = atoms.get_positions() 
    struc = Structure(Lattice(lattice), species, positions)
    return struc