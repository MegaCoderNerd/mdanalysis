#!/usr/bin/env python
"""
Example 3: Advanced Atom Selections

This script demonstrates:
- Various selection syntaxes
- Boolean operations on selections
- Geometric selections
- Dynamic selections
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD

def main():
    print("=" * 60)
    print("Example 3: Advanced Atom Selections")
    print("=" * 60)

    # Load universe
    u = mda.Universe(PSF, DCD)

    # 1. Selection by atom name
    print("\n1. Selection by atom name:")
    ca_atoms = u.select_atoms('name CA')
    print(f"   name CA: {len(ca_atoms)} atoms")

    cb_atoms = u.select_atoms('name CB')
    print(f"   name CB: {len(cb_atoms)} atoms")

    # 2. Selection by residue name
    print("\n2. Selection by residue name:")
    ala = u.select_atoms('resname ALA')
    print(f"   resname ALA: {len(ala)} atoms")

    # Multiple residue names
    charged = u.select_atoms('resname ARG LYS ASP GLU')
    print(f"   charged residues (ARG, LYS, ASP, GLU): {len(charged)} atoms")

    # 3. Selection by residue ID
    print("\n3. Selection by residue ID:")
    res42 = u.select_atoms('resid 42')
    print(f"   resid 42: {len(res42)} atoms")

    # Residue range
    res_range = u.select_atoms('resid 1-10')
    print(f"   resid 1-10: {len(res_range)} atoms")

    # Specific residues
    specific = u.select_atoms('resid 1 5 10 15')
    print(f"   resid 1 5 10 15: {len(specific)} atoms")

    # 4. Boolean operations
    print("\n4. Boolean operations:")

    # AND
    ca_in_helix = u.select_atoms('name CA and resid 1-20')
    print(f"   name CA AND resid 1-20: {len(ca_in_helix)} atoms")

    # OR
    ca_or_cb = u.select_atoms('name CA or name CB')
    print(f"   name CA OR name CB: {len(ca_or_cb)} atoms")

    # NOT
    non_protein = u.select_atoms('not protein')
    print(f"   NOT protein: {len(non_protein)} atoms")

    # Complex
    complex_sel = u.select_atoms('(protein and not backbone) and resid 1-50')
    print(f"   (protein AND NOT backbone) AND resid 1-50: {len(complex_sel)} atoms")

    # 5. Predefined selections
    print("\n5. Predefined selections:")
    protein = u.select_atoms('protein')
    print(f"   protein: {len(protein)} atoms")

    backbone = u.select_atoms('backbone')
    print(f"   backbone: {len(backbone)} atoms")

    nucleic = u.select_atoms('nucleic')
    print(f"   nucleic: {len(nucleic)} atoms")

    # 6. Geometric selections
    print("\n6. Geometric selections:")

    # Atoms around a point
    ref_residue = u.select_atoms('resid 42')
    around = u.select_atoms('around 5.0 resid 42')
    print(f"   around 5.0 Å of resid 42: {len(around)} atoms")

    # Spherical zone
    sphzone = u.select_atoms('sphzone 6.0 resid 42')
    print(f"   sphzone 6.0 Å of resid 42: {len(sphzone)} atoms")

    # Protein atoms near a specific residue
    near_42 = u.select_atoms('protein and around 8.0 resid 42')
    print(f"   protein within 8.0 Å of resid 42: {len(near_42)} atoms")

    # 7. Property-based selections
    print("\n7. Property-based selections:")

    # By mass
    heavy = u.select_atoms('mass > 15.0')
    print(f"   mass > 15.0: {len(heavy)} atoms")

    # 8. Selection arithmetic
    print("\n8. Selection arithmetic (combining AtomGroups):")

    # Union
    ca = u.select_atoms('name CA')
    cb = u.select_atoms('name CB')
    union = ca | cb  # or ca + cb
    print(f"   CA | CB (union): {len(union)} atoms")

    # Intersection
    helix_ca = u.select_atoms('resid 1-20')
    intersection = ca & helix_ca
    print(f"   CA & helix (intersection): {len(intersection)} atoms")

    # Difference
    not_ca = protein - ca
    print(f"   protein - CA (difference): {len(not_ca)} atoms")

    # 9. Dynamic vs Static selections
    print("\n9. Dynamic selections (update each frame):")

    # Jump to frame 0
    u.trajectory[0]

    # Dynamic selection that updates each frame
    # (useful for geometric criteria that change)
    dynamic_sel = u.select_atoms('around 5.0 resid 42', updating=True)
    initial_count = len(dynamic_sel)

    # Jump to another frame
    u.trajectory[50]
    final_count = len(dynamic_sel)

    print(f"   Frame 0: {initial_count} atoms around resid 42")
    print(f"   Frame 50: {final_count} atoms around resid 42")
    print(f"   Difference: {final_count - initial_count} atoms")

    # 10. Accessing selection information
    print("\n10. Selection information:")
    selection = u.select_atoms('resid 1-5 and name CA')
    print(f"   Selection: 'resid 1-5 and name CA'")
    print(f"   Number of atoms: {len(selection)}")
    print(f"   Atom indices: {selection.indices[:10]}")  # First 10
    print(f"   Atom names: {selection.names[:10]}")
    print(f"   Residue IDs: {[a.resid for a in selection]}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)

if __name__ == '__main__':
    main()
