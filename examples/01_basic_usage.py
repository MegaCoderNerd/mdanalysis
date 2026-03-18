#!/usr/bin/env python
"""
Example 1: Basic MDAnalysis Usage

This script demonstrates:
- Loading a molecular system
- Creating a Universe object
- Basic atom selections
- Accessing atom properties
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD

def main():
    print("=" * 60)
    print("Example 1: Basic MDAnalysis Usage")
    print("=" * 60)

    # Create a Universe from topology and trajectory files
    # Here we use built-in test data (adenylate kinase)
    print("\n1. Loading molecular system...")
    u = mda.Universe(PSF, DCD)
    print(f"   ✓ Universe created successfully")
    print(f"   - Topology file: {PSF}")
    print(f"   - Trajectory file: {DCD}")

    # Inspect the system
    print("\n2. System information:")
    print(f"   - Total atoms: {len(u.atoms)}")
    print(f"   - Total residues: {len(u.residues)}")
    print(f"   - Total segments: {len(u.segments)}")
    print(f"   - Number of trajectory frames: {len(u.trajectory)}")
    print(f"   - Current timestep: {u.trajectory.ts.frame}")

    # Access box dimensions
    print("\n3. Simulation box:")
    dims = u.dimensions
    print(f"   - Dimensions: {dims[:3]} Å")
    print(f"   - Angles: {dims[3:]} degrees")

    # Basic atom selections
    print("\n4. Basic atom selections:")

    # Select all atoms (default)
    all_atoms = u.atoms
    print(f"   - All atoms: {len(all_atoms)}")

    # Select protein atoms
    protein = u.select_atoms('protein')
    print(f"   - Protein atoms: {len(protein)}")

    # Select C-alpha atoms
    ca_atoms = u.select_atoms('name CA')
    print(f"   - C-alpha atoms: {len(ca_atoms)}")

    # Select backbone atoms
    backbone = u.select_atoms('backbone')
    print(f"   - Backbone atoms: {len(backbone)}")

    # Access atom properties
    print("\n5. Atom properties (first C-alpha):")
    first_ca = ca_atoms[0]
    print(f"   - Atom index: {first_ca.index}")
    print(f"   - Atom name: {first_ca.name}")
    print(f"   - Residue: {first_ca.resname}{first_ca.resid}")
    print(f"   - Position: {first_ca.position}")
    print(f"   - Mass: {first_ca.mass} u")

    # Access group properties (as NumPy arrays)
    print("\n6. AtomGroup properties (C-alphas):")
    print(f"   - Positions shape: {ca_atoms.positions.shape}")
    print(f"   - Masses shape: {ca_atoms.masses.shape}")
    print(f"   - Total mass: {ca_atoms.total_mass():.2f} u")

    # Calculate group properties
    print("\n7. Group calculations:")
    com = ca_atoms.center_of_mass()
    cog = ca_atoms.center_of_geometry()
    rgyr = ca_atoms.radius_of_gyration()

    print(f"   - Center of mass: {com}")
    print(f"   - Center of geometry: {cog}")
    print(f"   - Radius of gyration: {rgyr:.2f} Å")

    # Residue-level information
    print("\n8. Residue information:")
    print(f"   - Total residues: {len(u.residues)}")
    print(f"   - First residue: {u.residues[0].resname}{u.residues[0].resid}")
    print(f"   - Last residue: {u.residues[-1].resname}{u.residues[-1].resid}")

    # Count residue types
    residue_names = [r.resname for r in u.residues]
    unique_residues = set(residue_names)
    print(f"   - Unique residue types: {len(unique_residues)}")
    print(f"   - Residue types: {', '.join(sorted(unique_residues))}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)

if __name__ == '__main__':
    main()
