from ase import Atoms
from pymatgen.core import Lattice, Structure

def Atoms2Structure(atoms: Atoms):
    lattice = atoms.get_cell().tolist()  
    species = atoms.get_chemical_symbols()  
    positions = atoms.get_positions() 
    struc = Structure(Lattice(lattice), species, positions, coords_are_cartesian=True)
    return struc

def Structure2Atoms(structure: Structure) -> Atoms:  
    symbols = [site.specie.symbol for site in structure]  
    positions = [site.coords for site in structure]  
    cell = structure.lattice.matrix  
    ase_atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)  
    return ase_atoms