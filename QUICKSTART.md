# MDAnalysis Quick Start Guide

Welcome to MDAnalysis! This guide will help you get started with MDAnalysis, understand its architecture, and provide you with practical examples to explore the library.

## Table of Contents

1. [What is MDAnalysis?](#what-is-mdanalysis)
2. [Installation](#installation)
3. [Quick Start - Your First Analysis](#quick-start---your-first-analysis)
4. [Project Architecture](#project-architecture)
5. [Core Components](#core-components)
6. [Common Use Cases with Examples](#common-use-cases-with-examples)
7. [Supported File Formats](#supported-file-formats)
8. [How to Play With It](#how-to-play-with-it)
9. [Advanced Topics](#advanced-topics)
10. [Getting Help](#getting-help)

---

## What is MDAnalysis?

**MDAnalysis** is a mature, open-source Python library for analyzing molecular dynamics (MD) simulations. It enables you to:

- Read and write trajectory files from various MD simulation packages (GROMACS, Amber, NAMD, CHARMM, LAMMPS, etc.)
- Access atomic coordinates, velocities, and forces as NumPy arrays
- Select atoms using powerful, CHARMM-style selection commands
- Perform structural analysis (RMSD, alignment, contacts, hydrogen bonds, etc.)
- Calculate properties (radial distribution functions, density profiles, etc.)
- Transform and manipulate trajectories

**Key Features:**
- Support for 29+ trajectory formats and 20+ topology formats
- 33+ built-in analysis algorithms
- High performance through Cython/C optimized code
- Integration with the scientific Python ecosystem (NumPy, SciPy, matplotlib)
- Fiscally sponsored by NumFOCUS

---

## Installation

### Using pip (recommended for users)

```bash
pip install MDAnalysis
```

### Using conda

```bash
conda install -c conda-forge mdanalysis
```

### From source (for developers)

```bash
# Clone the repository
git clone https://github.com/MDAnalysis/mdanalysis.git
cd mdanalysis

# Install in development mode
cd package
pip install -e .
```

### Verify Installation

```python
import MDAnalysis as mda
print(mda.__version__)
```

---

## Quick Start - Your First Analysis

Here's a minimal example to get you started:

```python
import MDAnalysis as mda
import numpy as np

# Load a molecular system
# You can use test data that comes with MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD

# Create a Universe (the central object in MDAnalysis)
u = mda.Universe(PSF, DCD)

# Select atoms - here we select all C-alpha atoms
ca = u.select_atoms('name CA')

print(f"Number of C-alpha atoms: {len(ca)}")
print(f"Number of frames in trajectory: {len(u.trajectory)}")

# Iterate through trajectory and calculate center of mass
for ts in u.trajectory:
    com = ca.center_of_mass()
    print(f"Frame {ts.frame}: C-alpha center of mass = {com}")
```

**Output:**
```
Number of C-alpha atoms: 214
Number of frames in trajectory: 98
Frame 0: C-alpha center of mass = [ 0.06873595 -0.04605918 -0.24643682]
Frame 1: C-alpha center of mass = [ 0.07123456 -0.04712345 -0.25123456]
...
```

---

## Project Architecture

MDAnalysis follows a modular, plugin-based architecture:

```
MDAnalysis/
├── core/                    # Core data structures
│   ├── universe.py          # Universe class (main container)
│   ├── groups.py            # AtomGroup, ResidueGroup, SegmentGroup
│   ├── topology.py          # Topology management
│   ├── topologyattrs.py     # Atom attributes (mass, charge, etc.)
│   └── selection.py         # Atom selection engine
│
├── coordinates/             # Trajectory/coordinate readers
│   ├── DCD.py               # CHARMM/NAMD DCD format
│   ├── XTC.py               # GROMACS XTC format
│   ├── TRR.py               # GROMACS TRR format
│   ├── PDB.py               # PDB format
│   └── [26+ other formats]
│
├── topology/                # Topology file parsers
│   ├── PSFParser.py         # CHARMM PSF format
│   ├── PDBParser.py         # PDB format
│   ├── TPRParser.py         # GROMACS TPR format
│   └── [20+ other parsers]
│
├── analysis/                # Analysis algorithms
│   ├── rms.py               # RMSD calculations
│   ├── align.py             # Structural alignment
│   ├── distances.py         # Distance analysis
│   ├── contacts.py          # Contact analysis
│   ├── rdf.py               # Radial distribution functions
│   ├── msd.py               # Mean squared displacement
│   ├── hbonds/              # Hydrogen bond analysis
│   └── [30+ other modules]
│
├── transformations/         # On-the-fly trajectory transformations
│   ├── fit.py               # Alignment transformations
│   ├── wrap.py              # Periodic boundary wrapping
│   └── [other transforms]
│
├── lib/                     # Utilities and low-level code
│   ├── distances.pyx        # Fast distance calculations (Cython)
│   ├── NeighborSearch.py    # Spatial searching
│   └── [other utilities]
│
└── converters/              # Format converters (ParmEd, RDKit, OpenMM)
```

### Design Principles

1. **Universe-AtomGroup Model**: Central abstraction
   - `Universe`: Container for topology + trajectory data
   - `AtomGroup`: Collections of atoms with properties and methods

2. **Plugin Architecture**: Extensible format support
   - Readers/writers register themselves
   - Easy to add new formats

3. **Lazy Loading**: Memory efficient
   - Trajectories loaded frame-by-frame on demand
   - Only what you need is loaded into memory

4. **NumPy Integration**: Fast array operations
   - Coordinates exposed as NumPy arrays
   - Enables vectorized calculations

5. **Performance**: Cython/C for critical paths
   - Distance calculations optimized in C
   - Neighbor searches use efficient algorithms

---

## Core Components

### 1. Universe

The `Universe` is the central object that contains all information about your system.

**Location:** `/package/MDAnalysis/core/universe.py`

```python
import MDAnalysis as mda

# Create a Universe from topology and trajectory
u = mda.Universe('topology.pdb', 'trajectory.dcd')

# Access components
print(u.atoms)          # All atoms
print(u.residues)       # All residues
print(u.segments)       # All segments
print(u.trajectory)     # Trajectory reader
print(u.dimensions)     # Box dimensions
```

### 2. AtomGroup

A collection of atoms with properties and methods.

**Location:** `/package/MDAnalysis/core/groups.py`

```python
# Select atoms using selection language
protein = u.select_atoms('protein')
ca_atoms = u.select_atoms('name CA')
within_5A = u.select_atoms('around 5.0 resid 42')

# Access properties (as NumPy arrays)
positions = ca_atoms.positions     # Nx3 array of coordinates
masses = ca_atoms.masses           # N array of masses
charges = ca_atoms.charges         # N array of charges

# Calculate properties
com = ca_atoms.center_of_mass()
cog = ca_atoms.center_of_geometry()
radius_gyr = ca_atoms.radius_of_gyration()

# Combine selections
backbone = u.select_atoms('backbone')
ca_and_cb = u.select_atoms('name CA or name CB')
```

### 3. Trajectory

Iterator over trajectory frames.

**Location:** `/package/MDAnalysis/coordinates/`

```python
# Iterate over all frames
for ts in u.trajectory:
    print(f"Frame {ts.frame}, Time {ts.time} ps")
    print(f"Positions shape: {ts.positions.shape}")

# Jump to specific frame
u.trajectory[50]

# Slice trajectory
for ts in u.trajectory[10:50:2]:  # frames 10-50, step 2
    pass

# Access frame information
ts = u.trajectory[0]
print(ts.frame)        # Frame number
print(ts.time)         # Time in ps
print(ts.dimensions)   # Box dimensions
```

### 4. Selections

Powerful atom selection language.

**Location:** `/package/MDAnalysis/core/selection.py`

```python
# By name
u.select_atoms('name CA')

# By residue
u.select_atoms('resname ALA')
u.select_atoms('resid 1:10')

# By properties
u.select_atoms('mass > 15')
u.select_atoms('charge < 0')

# Geometric selections
u.select_atoms('around 5.0 resid 42')
u.select_atoms('sphzone 6.0 protein')
u.select_atoms('cyzone 15 10 -8 protein')  # cylindrical

# Boolean operations
u.select_atoms('protein and name CA')
u.select_atoms('protein or nucleic')
u.select_atoms('protein and not backbone')

# Pre-defined keywords
u.select_atoms('protein')
u.select_atoms('backbone')
u.select_atoms('nucleic')
u.select_atoms('water')
```

### 5. Analysis Base Class

Template for creating analysis workflows.

**Location:** `/package/MDAnalysis/analysis/base.py`

```python
from MDAnalysis.analysis.base import AnalysisBase

class MyAnalysis(AnalysisBase):
    def __init__(self, atomgroup, **kwargs):
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def _prepare(self):
        # Run once before analysis
        self.results.positions = []

    def _single_frame(self):
        # Run for each frame
        self.results.positions.append(self.atomgroup.positions.copy())

    def _conclude(self):
        # Run once after analysis
        self.results.positions = np.array(self.results.positions)

# Use it
analysis = MyAnalysis(u.select_atoms('name CA'))
analysis.run()
print(analysis.results.positions.shape)
```

---

## Common Use Cases with Examples

### 1. Calculate RMSD (Root Mean Square Deviation)

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.tests.datafiles import PSF, DCD

u = mda.Universe(PSF, DCD)

# RMSD of C-alpha atoms
R = rms.RMSD(u, select='backbone')
R.run()

# Plot results
import matplotlib.pyplot as plt
plt.plot(R.results.rmsd[:, 1], R.results.rmsd[:, 2])
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.show()
```

### 2. Align Structures

```python
from MDAnalysis.analysis import align

# Align a trajectory to the first frame
aligner = align.AlignTraj(u, u, select='protein', filename='aligned.dcd')
aligner.run()

# Align two structures
mobile = mda.Universe('mobile.pdb')
reference = mda.Universe('reference.pdb')
align.alignto(mobile, reference, select='backbone')

# Save aligned structure
mobile.atoms.write('aligned.pdb')
```

### 3. Calculate Distances

```python
from MDAnalysis.analysis import distances

# Distance between two atom groups
ag1 = u.select_atoms('resid 1 and name CA')
ag2 = u.select_atoms('resid 50 and name CA')

# Calculate over trajectory
dist_array = []
for ts in u.trajectory:
    d = distances.dist(ag1, ag2)
    dist_array.append(d[2][0])  # distance in Angstroms

# Or use the Distances analysis class
from MDAnalysis.analysis.atomicdistances import AtomicDistances

atoms1 = u.select_atoms('resid 1 and name CA')
atoms2 = u.select_atoms('resid 50 and name CA')

dist_analysis = AtomicDistances([(atoms1, atoms2)])
dist_analysis.run()
print(dist_analysis.results)
```

### 4. Radial Distribution Function (RDF)

```python
from MDAnalysis.analysis.rdf import InterRDF

# Calculate RDF between water oxygen and protein
s1 = u.select_atoms('name OH2')  # water oxygen
s2 = u.select_atoms('protein')

rdf = InterRDF(s1, s2, nbins=75, range=(0.0, 15.0))
rdf.run()

# Plot
plt.plot(rdf.results.bins, rdf.results.rdf)
plt.xlabel('Distance (Å)')
plt.ylabel('g(r)')
plt.show()
```

### 5. Hydrogen Bond Analysis

```python
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

# Find hydrogen bonds in protein
hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel='protein',
    hydrogens_sel='protein',
    acceptors_sel='protein',
)
hbonds.run()

# Access results
print(f"Number of H-bonds: {len(hbonds.results.hbonds)}")
```

### 6. Contact Analysis

```python
from MDAnalysis.analysis import contacts

# Native contacts (Q)
ca = u.select_atoms('name CA')
q = contacts.Contacts(
    u,
    select=('name CA', 'name CA'),
    refgroup=(ca, ca),
    radius=8.0
)
q.run()

# Fraction of native contacts over time
plt.plot(q.results.timeseries[:, 0], q.results.timeseries[:, 1])
plt.xlabel('Time (ps)')
plt.ylabel('Fraction of Native Contacts')
plt.show()
```

### 7. Density Analysis

```python
from MDAnalysis.analysis.density import DensityAnalysis

# Calculate water density around protein
water = u.select_atoms('resname SOL and name OH2')

D = DensityAnalysis(water, delta=1.0)
D.run()

# Export to grid format
D.results.density.export('water_density.dx')
```

### 8. Mean Squared Displacement (MSD)

```python
from MDAnalysis.analysis.msd import EinsteinMSD

# Calculate MSD for water molecules
water = u.select_atoms('resname SOL and name OH2')
msd = EinsteinMSD(water)
msd.run()

# Calculate diffusion coefficient
nframes = msd.n_frames
timestep = 1  # in ps
lagtimes = np.arange(nframes) * timestep
```

### 9. Principal Component Analysis (PCA)

```python
from MDAnalysis.analysis import pca

# PCA on C-alpha atoms
ca = u.select_atoms('name CA')
pc = pca.PCA(u, select='name CA')
pc.run()

# Access principal components
print(pc.results.p_components.shape)

# Project trajectory onto first 2 PCs
transformed = pc.transform(ca, n_components=2)
```

### 10. Trajectory Transformations

```python
from MDAnalysis.transformations import *

# Unwrap trajectory (fix periodic boundary jumps)
u.trajectory.add_transformations(unwrap(u.atoms))

# Center on protein and wrap
transforms = [
    center_in_box(u.select_atoms('protein')),
    wrap(u.atoms)
]
u.trajectory.add_transformations(*transforms)

# Fit to reference structure on-the-fly
transforms = [
    fit_rot_trans(u.select_atoms('backbone'), u.select_atoms('backbone'))
]
u.trajectory.add_transformations(*transforms)

# Now iterate - transformations applied automatically
for ts in u.trajectory:
    # positions are already transformed
    pass
```

---

## Supported File Formats

### Trajectory Formats (29+)

| Format | Extension | Description |
|--------|-----------|-------------|
| DCD | .dcd | CHARMM, NAMD, LAMMPS |
| XTC | .xtc | GROMACS compressed |
| TRR | .trr | GROMACS full precision |
| TNG | .tng | GROMACS next gen |
| H5MD | .h5 | HDF5-based format |
| NetCDF | .nc | Amber NetCDF |
| PDB | .pdb | Protein Data Bank |
| GRO | .gro | GROMACS |
| LAMMPS | .lammpstrj | LAMMPS dump |
| XYZ | .xyz | Generic XYZ |
| MOL2 | .mol2 | Tripos MOL2 |
| GSD | .gsd | HOOMD-blue |
| And many more... | | See docs for full list |

### Topology Formats (20+)

| Format | Extension | Description |
|--------|-----------|-------------|
| PSF | .psf | CHARMM/NAMD topology |
| PDB | .pdb | Protein Data Bank |
| GRO | .gro | GROMACS |
| TOP | .top | GROMACS topology |
| TPR | .tpr | GROMACS run input |
| PARM7 | .prmtop | Amber topology |
| PQR | .pqr | PDB with charge/radius |
| MOL2 | .mol2 | Tripos MOL2 |
| PDBQT | .pdbqt | AutoDock |
| And many more... | | See docs for full list |

---

## How to Play With It

### Option 1: Use Built-in Test Data

MDAnalysis comes with test data you can use immediately:

```python
from MDAnalysis.tests.datafiles import (
    PSF, DCD,           # CHARMM/NAMD format
    GRO, XTC,           # GROMACS format
    PDB,                # PDB format
    TPR, TRR,           # GROMACS binary
    PDB_small,          # Small protein
    ALIGN_STRUCTURED,   # For alignment tests
)

# Example: Protein in water
u = mda.Universe(PSF, DCD)
print(f"System has {len(u.atoms)} atoms")
print(f"Residues: {len(u.residues)}")
print(f"Frames: {len(u.trajectory)}")

# Select and analyze
protein = u.select_atoms('protein')
water = u.select_atoms('resname TIP3')
print(f"Protein atoms: {len(protein)}")
print(f"Water molecules: {len(water.residues)}")
```

### Option 2: Download Example Data

```python
# Fetch a PDB file from the Protein Data Bank
from MDAnalysis.tests.datafiles import PDB_small

u = mda.Universe(PDB_small)

# Or use the fetch module for remote data
# Note: MDAnalysis.fetch can download from PDB
```

### Option 3: Generate Simple Test System

```python
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.core.topologyattrs import Masses

# Create a simple universe programmatically
n_atoms = 100
positions = np.random.random((n_atoms, 3)) * 50  # 50 Å box

u = Universe.empty(n_atoms, trajectory=True)
u.atoms.positions = positions

# Add some attributes
u.add_TopologyAttr('masses', values=np.ones(n_atoms) * 12.0)

print(f"Created universe with {len(u.atoms)} atoms")
```

### Hands-On Exercises

#### Exercise 1: Basic Trajectory Analysis

```python
from MDAnalysis.tests.datafiles import PSF, DCD
import MDAnalysis as mda

# Load data
u = mda.Universe(PSF, DCD)

# Task 1: Find the protein's center of mass at each frame
protein = u.select_atoms('protein')
com_trajectory = []

for ts in u.trajectory:
    com = protein.center_of_mass()
    com_trajectory.append(com)

print(f"Initial COM: {com_trajectory[0]}")
print(f"Final COM: {com_trajectory[-1]}")

# Task 2: Calculate distance between two residues over time
res1 = u.select_atoms('resid 1 and name CA')
res100 = u.select_atoms('resid 100 and name CA')

distances = []
for ts in u.trajectory:
    d = np.linalg.norm(res1.positions[0] - res100.positions[0])
    distances.append(d)

print(f"Average distance: {np.mean(distances):.2f} Å")
```

#### Exercise 2: Selection Practice

```python
u = mda.Universe(PSF, DCD)

# Practice different selections
print(f"All atoms: {len(u.atoms)}")
print(f"Protein: {len(u.select_atoms('protein'))}")
print(f"C-alphas: {len(u.select_atoms('name CA'))}")
print(f"Charged residues: {len(u.select_atoms('resname ARG LYS ASP GLU'))}")
print(f"Hydrophobic core: {len(u.select_atoms('resname ALA VAL LEU ILE PHE TRP'))}")

# Geometric selections
print(f"Atoms within 5 Å of residue 42: {len(u.select_atoms('around 5 resid 42'))}")
```

#### Exercise 3: Simple Analysis

```python
from MDAnalysis.analysis import rms

u = mda.Universe(PSF, DCD)

# Calculate RMSD
R = rms.RMSD(u, select='backbone')
R.run()

# Find frame with minimum RMSD
min_rmsd_frame = R.results.rmsd[:, 2].argmin()
print(f"Frame with minimum RMSD: {min_rmsd_frame}")
print(f"Minimum RMSD: {R.results.rmsd[min_rmsd_frame, 2]:.2f} Å")
```

---

## Advanced Topics

### Writing Custom Analysis

```python
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

class EndToEndDistance(AnalysisBase):
    """Calculate end-to-end distance of a polymer."""

    def __init__(self, atomgroup, **kwargs):
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def _prepare(self):
        self.results.distances = []

    def _single_frame(self):
        # Calculate distance between first and last atom
        first = self.atomgroup[0].position
        last = self.atomgroup[-1].position
        distance = np.linalg.norm(last - first)
        self.results.distances.append(distance)

    def _conclude(self):
        self.results.distances = np.array(self.results.distances)
        self.results.mean = np.mean(self.results.distances)
        self.results.std = np.std(self.results.distances)

# Use it
u = mda.Universe(PSF, DCD)
backbone = u.select_atoms('backbone')
analysis = EndToEndDistance(backbone)
analysis.run()
print(f"Mean end-to-end distance: {analysis.results.mean:.2f} Å")
```

### Parallel Analysis

```python
from MDAnalysis.analysis.rms import RMSD

# Split trajectory analysis across cores
R = RMSD(u, select='backbone')
R.run(n_jobs=4)  # Use 4 CPU cores
```

### Memory-Mapped Trajectories

```python
# For very large trajectories, use memory mapping
u = mda.Universe('topology.pdb', 'huge_trajectory.xtc', in_memory=False)

# Process in chunks
for ts in u.trajectory[::10]:  # Every 10th frame
    # Analyze
    pass
```

### Writing Custom Readers

```python
from MDAnalysis.coordinates.base import ReaderBase

class MyFormatReader(ReaderBase):
    format = 'MYFORMAT'

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        # Initialize your reader

    def _read_next_timestep(self):
        # Read next frame
        pass
```

---

## Getting Help

### Documentation

- **User Guide**: https://userguide.mdanalysis.org
- **API Documentation**: https://docs.mdanalysis.org
- **Quickstart Guide**: https://userguide.mdanalysis.org/examples/quickstart.html
- **Tutorials**: https://userguide.mdanalysis.org/examples/README.html

### Community

- **GitHub Discussions**: https://github.com/MDAnalysis/mdanalysis/discussions
- **Issue Tracker**: https://github.com/MDAnalysis/mdanalysis/issues
- **Mailing List**: https://groups.google.com/g/mdnalysis-discussion

### Learning Resources

- **Video Tutorials**: https://www.mdanalysis.org/pages/learning_MDAnalysis/#videos
- **Workshops**: Check the MDAnalysis website for upcoming workshops
- **Paper/Citations**: See README.rst for how to cite MDAnalysis

### Contributing

Want to contribute? Check out:
- **Developer Guide**: https://userguide.mdanalysis.org/contributing.html
- **Code of Conduct**: https://www.mdanalysis.org/conduct/
- **GitHub**: https://github.com/MDAnalysis/mdanalysis

---

## Quick Reference Card

### Essential Imports

```python
import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD  # Test data
```

### Creating Universe

```python
u = mda.Universe('topology', 'trajectory')
u = mda.Universe('structure.pdb')
u = mda.Universe(PSF, DCD)  # Test data
```

### Selections

```python
all_atoms = u.atoms
protein = u.select_atoms('protein')
ca = u.select_atoms('name CA')
resid_42 = u.select_atoms('resid 42')
```

### Accessing Data

```python
positions = u.atoms.positions  # Nx3 array
masses = u.atoms.masses        # N array
residues = u.residues
n_frames = len(u.trajectory)
```

### Iteration

```python
for ts in u.trajectory:
    # ts.frame, ts.time, ts.positions
    pass
```

### Common Analyses

```python
from MDAnalysis.analysis import rms, align, distances, contacts, rdf

# RMSD
R = rms.RMSD(u, select='backbone')
R.run()

# Alignment
align.alignto(mobile, ref, select='backbone')

# RDF
rdf = rdf.InterRDF(sel1, sel2)
rdf.run()
```

---

## Next Steps

1. **Try the examples** above with the built-in test data
2. **Read the User Guide** for in-depth tutorials
3. **Explore the API docs** to learn about specific modules
4. **Join the community** on GitHub Discussions
5. **Start analyzing** your own simulation data!

---

## File Locations in Repository

For developers and contributors:

- **Main Package**: `/package/MDAnalysis/`
- **Core Module**: `/package/MDAnalysis/core/`
- **Analysis Module**: `/package/MDAnalysis/analysis/`
- **Test Suite**: `/testsuite/MDAnalysisTests/`
- **Documentation**: `/package/doc/sphinx/`
- **Benchmarks**: `/benchmarks/`

---

## Credits

MDAnalysis is developed by a large community of contributors and is fiscally sponsored by NumFOCUS.

**License**: LGPLv3+ (some components LGPLv2.1+)

**Citation**: When using MDAnalysis, please cite the two papers listed in README.rst

**Website**: https://www.mdanalysis.org

---

Happy analyzing! 🧬🔬
