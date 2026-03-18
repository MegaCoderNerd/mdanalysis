# MDAnalysis Quick Reference Card

A one-page reference for common MDAnalysis operations.

## Installation

```bash
pip install MDAnalysis
# or
conda install -c conda-forge mdanalysis
```

## Basic Imports

```python
import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD  # Test data
```

## Creating a Universe

```python
# From files
u = mda.Universe('topology.pdb', 'trajectory.dcd')
u = mda.Universe('topology.gro', 'traj1.xtc', 'traj2.xtc')  # Multiple trajectories

# Using test data
u = mda.Universe(PSF, DCD)
```

## Atom Selections

```python
# By name/type
u.select_atoms('name CA')
u.select_atoms('type C')

# By residue
u.select_atoms('resname ALA')
u.select_atoms('resid 1-10')
u.select_atoms('resid 1 5 10 15')

# Predefined
u.select_atoms('protein')
u.select_atoms('backbone')
u.select_atoms('nucleic')
u.select_atoms('water')

# Geometric
u.select_atoms('around 5.0 resid 42')
u.select_atoms('sphzone 6.0 protein')

# Boolean
u.select_atoms('protein and name CA')
u.select_atoms('resname ALA or resname GLY')
u.select_atoms('protein and not backbone')

# Dynamic (updates each frame)
u.select_atoms('around 5.0 resid 42', updating=True)
```

## AtomGroup Properties

```python
ag = u.select_atoms('protein')

# Coordinates (NumPy arrays)
ag.positions       # Nx3 array
ag.velocities      # Nx3 array (if available)
ag.forces          # Nx3 array (if available)

# Topology attributes
ag.names           # Atom names
ag.masses          # Masses
ag.charges         # Charges
ag.resnames        # Residue names
ag.resids          # Residue IDs

# Properties
ag.n_atoms         # Number of atoms
ag.n_residues      # Number of residues
ag.total_mass()    # Sum of masses
ag.total_charge()  # Sum of charges
```

## AtomGroup Calculations

```python
ag.center_of_mass()        # Center of mass
ag.center_of_geometry()    # Geometric center
ag.radius_of_gyration()    # Radius of gyration
ag.principal_axes()        # Principal axes
ag.moment_of_inertia()     # Moment of inertia tensor
```

## Trajectory Iteration

```python
# Iterate all frames
for ts in u.trajectory:
    print(ts.frame, ts.time)

# Slice trajectory
for ts in u.trajectory[0:100:10]:  # Frames 0-100, step 10
    pass

# Access specific frame
u.trajectory[50]

# Frame information
ts = u.trajectory.ts
ts.frame           # Frame number
ts.time            # Time (ps)
ts.dt              # Timestep (ps)
ts.dimensions      # Box dimensions [a, b, c, α, β, γ]
ts.positions       # All positions
```

## Common Analyses

### RMSD
```python
from MDAnalysis.analysis import rms

R = rms.RMSD(u, select='backbone')
R.run()
rmsd_data = R.results.rmsd  # [frame, time, RMSD]
```

### Alignment
```python
from MDAnalysis.analysis import align

# Align two structures
align.alignto(mobile, reference, select='backbone')

# Align trajectory
aligner = align.AlignTraj(u, u, select='protein', filename='aligned.dcd')
aligner.run()
```

### RDF
```python
from MDAnalysis.analysis.rdf import InterRDF

s1 = u.select_atoms('name O')
s2 = u.select_atoms('name H')
rdf = InterRDF(s1, s2, nbins=75, range=(0.0, 15.0))
rdf.run()
# Results: rdf.results.bins, rdf.results.rdf
```

### Hydrogen Bonds
```python
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel='protein',
    hydrogens_sel='protein',
    acceptors_sel='protein'
)
hbonds.run()
```

### Contacts
```python
from MDAnalysis.analysis import contacts

ca = u.select_atoms('name CA')
q = contacts.Contacts(u, select=('name CA', 'name CA'),
                      refgroup=(ca, ca), radius=8.0)
q.run()
```

### Distances
```python
from MDAnalysis.analysis.atomicdistances import AtomicDistances

ag1 = u.select_atoms('resid 1 and name CA')
ag2 = u.select_atoms('resid 50 and name CA')
dist = AtomicDistances([(ag1, ag2)])
dist.run()
```

## Transformations

```python
from MDAnalysis.transformations import *

# Unwrap coordinates
u.trajectory.add_transformations(unwrap(u.atoms))

# Center and wrap
u.trajectory.add_transformations(
    center_in_box(u.select_atoms('protein')),
    wrap(u.atoms)
)

# Fit to reference
u.trajectory.add_transformations(
    fit_rot_trans(u.select_atoms('backbone'),
                  u.select_atoms('backbone'))
)
```

## Writing Output

```python
# Write single structure
u.atoms.write('output.pdb')
protein.write('protein.gro')

# Write trajectory
with mda.Writer('output.dcd', u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)

# Write aligned trajectory
align.AlignTraj(u, reference, select='backbone',
                filename='aligned.dcd').run()
```

## Custom Analysis Template

```python
from MDAnalysis.analysis.base import AnalysisBase

class MyAnalysis(AnalysisBase):
    def __init__(self, atomgroup, **kwargs):
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def _prepare(self):
        self.results.data = []

    def _single_frame(self):
        # Analysis for current frame
        value = calculate_something(self.atomgroup)
        self.results.data.append(value)

    def _conclude(self):
        self.results.data = np.array(self.results.data)

# Usage
analysis = MyAnalysis(u.select_atoms('protein'))
analysis.run()
```

## Performance Tips

```python
# Use views, not copies
positions = ag.positions  # View into timestep

# Vectorize with NumPy
distances = np.linalg.norm(positions - ref, axis=1)

# Process subset of frames
for ts in u.trajectory[::10]:  # Every 10th frame
    pass

# Parallel analysis (where supported)
analysis.run(n_jobs=4)
```

## Supported Formats

**Trajectories**: DCD, XTC, TRR, TNG, NetCDF, H5MD, PDB, GRO, LAMMPS, XYZ, MOL2, and 20+ more

**Topologies**: PSF, PDB, GRO, TOP, TPR, PARM7, MOL2, PQR, PDBQT, and 15+ more

## Getting Help

- **QUICKSTART.md**: Comprehensive guide in this repo
- **examples/**: Ready-to-run examples in this repo
- **ARCHITECTURE.md**: Technical architecture docs in this repo
- **User Guide**: https://userguide.mdanalysis.org
- **API Docs**: https://docs.mdanalysis.org
- **Discussions**: https://github.com/MDAnalysis/mdanalysis/discussions
- **Issues**: https://github.com/MDAnalysis/mdanalysis/issues

## Common Patterns

### Calculate property over trajectory
```python
values = []
for ts in u.trajectory:
    value = calculate(u.select_atoms('protein'))
    values.append(value)
values = np.array(values)
```

### Find frames matching criteria
```python
matching_frames = []
for ts in u.trajectory:
    if condition(ts):
        matching_frames.append(ts.frame)
```

### Compare two atom groups
```python
ag1 = u.select_atoms('protein')
ag2 = u.select_atoms('resid 1-50')
common = ag1 & ag2        # Intersection
either = ag1 | ag2        # Union
diff = ag1 - ag2          # Difference
```

### Track multiple properties
```python
results = {'com': [], 'rgyr': [], 'time': []}
for ts in u.trajectory:
    results['com'].append(protein.center_of_mass())
    results['rgyr'].append(protein.radius_of_gyration())
    results['time'].append(ts.time)
```

---

**Version**: MDAnalysis 2.11.0-dev0
**License**: LGPLv3+
**Website**: https://www.mdanalysis.org
