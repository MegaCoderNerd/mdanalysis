# MDAnalysis Architecture

This document provides a comprehensive overview of the MDAnalysis architecture, design patterns, and internal structure.

## Table of Contents

1. [Overview](#overview)
2. [Core Design Principles](#core-design-principles)
3. [Main Components](#main-components)
4. [Data Flow](#data-flow)
5. [Plugin Architecture](#plugin-architecture)
6. [Performance Optimizations](#performance-optimizations)
7. [Testing Infrastructure](#testing-infrastructure)
8. [Build and Distribution](#build-and-distribution)

---

## Overview

MDAnalysis is designed as a **modular, extensible Python library** for molecular dynamics analysis. The architecture emphasizes:

- **Separation of concerns**: Clear boundaries between topology, coordinates, and analysis
- **Extensibility**: Plugin-based architecture for file formats
- **Performance**: Strategic use of compiled code (Cython/C) for critical paths
- **Ease of use**: High-level Python API with NumPy integration

### Project Structure

```
mdanalysis/
├── package/                  # Main library code
│   └── MDAnalysis/          # Python package
│       ├── core/            # Core data structures
│       ├── coordinates/     # Trajectory readers/writers
│       ├── topology/        # Topology parsers
│       ├── analysis/        # Analysis algorithms
│       ├── lib/             # Utilities and Cython code
│       ├── transformations/ # Trajectory transformations
│       ├── converters/      # External format converters
│       ├── auxiliary/       # Auxiliary data readers
│       └── ...
│
├── testsuite/               # Comprehensive test suite
│   └── MDAnalysisTests/     # Mirror of package structure
│
├── benchmarks/              # Performance benchmarks (ASV)
├── maintainer/              # Maintenance scripts
└── .github/                 # CI/CD workflows
```

---

## Core Design Principles

### 1. Universe-AtomGroup Model

The central abstraction in MDAnalysis is the **Universe-AtomGroup** model:

```
Universe
  ├── Topology (static structure information)
  │   ├── Atoms
  │   ├── Residues
  │   ├── Segments
  │   └── TopologyAttrs (masses, charges, etc.)
  │
  ├── Trajectory (dynamic coordinate information)
  │   └── Timesteps (positions, velocities, forces)
  │
  └── AtomGroups (selections and analysis)
      ├── Properties (positions, masses, etc.)
      └── Methods (center_of_mass, radius_of_gyration, etc.)
```

**Key insight**: Topology and trajectory are **separate but linked**. This allows:
- Multiple trajectories for the same topology
- Memory-efficient trajectory streaming
- Flexible analysis workflows

### 2. Lazy Loading and Memory Efficiency

**Principle**: Only load what you need, when you need it.

- **Trajectories**: Read frame-by-frame on demand (not loaded into memory)
- **Topology**: Parse only when accessed
- **Arrays**: NumPy views where possible (avoid copies)

```python
# This does NOT load entire trajectory into memory
u = mda.Universe('topology.pdb', 'huge_trajectory.xtc')

# Only current frame is in memory
for ts in u.trajectory:
    # ts.positions is a view, not a copy
    analyze(ts.positions)
```

### 3. NumPy Integration

**Principle**: All numerical data exposed as NumPy arrays.

Benefits:
- Native vectorization for performance
- Integration with scientific Python ecosystem
- Familiar API for users
- Zero-copy operations where possible

```python
# All these return NumPy arrays
positions = atomgroup.positions  # shape: (n_atoms, 3)
masses = atomgroup.masses        # shape: (n_atoms,)
velocities = timestep.velocities # shape: (n_atoms, 3)

# Enables vectorized operations
com = np.average(positions, weights=masses, axis=0)
```

### 4. Plugin Architecture

**Principle**: Format readers/writers are discoverable plugins.

- **Automatic registration**: Readers register themselves on import
- **Format detection**: Automatic based on file extension
- **Extensibility**: Easy to add new formats without modifying core

```python
# Format automatically detected from extension
u = mda.Universe('file.pdb')   # Uses PDBReader
u = mda.Universe('file.gro')   # Uses GROReader
u = mda.Universe('file.dcd')   # Uses DCDReader
```

### 5. Performance Through Compilation

**Principle**: Python for API, Cython/C for performance-critical code.

**Python layer**:
- High-level API
- Business logic
- User-facing interfaces

**Cython/C layer**:
- Distance calculations
- Neighbor searches
- Coordinate transformations
- Inner loops

File locations:
- `/package/MDAnalysis/lib/*.pyx` - Cython modules
- `/package/MDAnalysis/lib/src/` - C source code
- `/package/MDAnalysis/lib/include/` - C headers

---

## Main Components

### 1. Core Module (`/package/MDAnalysis/core/`)

The heart of MDAnalysis. Contains fundamental data structures.

#### Universe (`universe.py`)

**Purpose**: Main container for molecular system.

**Responsibilities**:
- Load topology and trajectory files
- Manage topology-trajectory relationship
- Provide atom selection interface
- Coordinate different components

**Key methods**:
```python
__init__(topology, *trajectories)  # Initialize from files
select_atoms(selection_string)     # Create AtomGroup
add_TopologyAttr(attr_name)        # Add topology attribute
load_new(trajectory)               # Load new trajectory
```

**Internal structure**:
```python
Universe
  ._topology      # Topology object
  ._trajectory    # TrajectoryReader object
  .atoms          # All atoms (AtomGroup)
  .residues       # All residues (ResidueGroup)
  .segments       # All segments (SegmentGroup)
```

#### AtomGroup (`groups.py`)

**Purpose**: Represent collections of atoms with properties and methods.

**Key features**:
- Array-like interface (`ag[0]`, `ag[1:10]`)
- Boolean operations (`ag1 & ag2`, `ag1 | ag2`, `ag1 - ag2`)
- Property access via TopologyAttrs
- Method dispatch for calculations

**Hierarchy**:
```
GroupBase (abstract base)
  ├── AtomGroup
  ├── ResidueGroup
  └── SegmentGroup
```

**Properties** (via TopologyAttrs):
- `positions`, `velocities`, `forces` (dynamic, from trajectory)
- `masses`, `charges`, `radii` (static, from topology)
- `names`, `types`, `resnames`, etc.

**Methods**:
- `center_of_mass()`, `center_of_geometry()`
- `radius_of_gyration()`
- `total_mass()`, `total_charge()`
- `write(filename)` - Write to file

#### Topology (`topology.py`, `topologyattrs.py`)

**Purpose**: Store and manage structural information.

**TopologyAttrs**:
Attributes are implemented as separate classes:
- `Atomnames`, `Atomids`, `Masses`, `Charges`
- `Resnames`, `Resids`, `Resnums`
- `Segids`, `Bonds`, `Angles`, `Dihedrals`

**Advantages**:
- Only store what's available in input file
- Uniform interface for access
- Easy to extend with new attributes

**Guessing**:
If attributes missing, MDAnalysis can guess:
- Masses (from atom types)
- Bonds (from distances)
- Atom types (from atom names)

#### Selection Engine (`selection.py`)

**Purpose**: Parse and execute atom selection commands.

**Architecture**:
1. **Parser**: Tokenize selection string
2. **Selection objects**: Represent criteria (ByName, ByResid, etc.)
3. **Evaluator**: Apply selections to atom indices
4. **Boolean logic**: Combine selections (AND, OR, NOT)

**Selection types**:
- **Static**: Evaluated once, cached
- **Dynamic**: Re-evaluated each frame (for geometric selections)

**Example flow**:
```
"name CA and resid 1-10"
    ↓ tokenize
["name", "CA", "and", "resid", "1-10"]
    ↓ parse
[ByName("CA"), And, ByResid(range(1,11))]
    ↓ evaluate
array([0, 15, 32, ...])  # atom indices
    ↓ create
AtomGroup(indices)
```

### 2. Coordinates Module (`/package/MDAnalysis/coordinates/`)

Handles trajectory and structure file I/O.

#### Base Classes

**`ReaderBase`** (`base.py`):
- Abstract base for all readers
- Defines reader interface
- Provides iteration protocol

**Key methods**:
```python
__init__(filename)          # Open file
_read_next_timestep()       # Read next frame (abstract)
__iter__()                  # Iteration protocol
__getitem__(frame)          # Random access
close()                     # Close file handle
```

**`WriterBase`** (`base.py`):
- Abstract base for trajectory writers
- Frame-by-frame writing interface

**`Timestep`** (`timestep.pyx`):
- Contains single frame data
- Implemented in Cython for performance
- Attributes: `positions`, `velocities`, `forces`, `dimensions`

#### Format Readers

Each format has its own reader module:

**Example: DCD Reader** (`DCD.py`):
```python
class DCDReader(ReaderBase):
    format = 'DCD'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self, filename):
        # Open file, read header
        # Determine n_frames, dt, etc.

    def _read_next_timestep(self):
        # Read binary data
        # Populate self.ts (Timestep object)
        return self.ts
```

**Format registration**:
Readers register themselves via `format` attribute and are automatically discovered.

#### Chain Reader

**Purpose**: Treat multiple trajectory files as one continuous trajectory.

```python
u = mda.Universe('topology.pdb',
                 'traj1.dcd', 'traj2.dcd', 'traj3.dcd')
# Seamlessly iterates across all files
```

### 3. Topology Module (`/package/MDAnalysis/topology/`)

Parses topology files to extract structural information.

#### Parser Architecture

**`TopologyReaderBase`** (`base.py`):
```python
class TopologyReaderBase:
    def parse(self):
        # Returns Topology object
        pass
```

**Parser responsibilities**:
1. Read file
2. Extract atoms, residues, segments
3. Parse topology attributes (masses, bonds, etc.)
4. Return `Topology` object

**Example: PSF Parser** (`PSFParser.py`):
```python
class PSFParser(TopologyReaderBase):
    format = 'PSF'

    def parse(self):
        # Parse PSF sections
        # Create Atoms, Residues, Segments
        # Parse Bonds, Angles, Dihedrals
        # Return Topology
```

#### Topology Guessing

**Purpose**: Fill in missing information.

**Guessers** (`topology/guessers.py`):
- `guess_masses()` - From atom names/types
- `guess_bonds()` - From distances
- `guess_types()` - From atom names
- `guess_angles()`, `guess_dihedrals()` - From bonds

### 4. Analysis Module (`/package/MDAnalysis/analysis/`)

Pre-built analysis algorithms.

#### AnalysisBase Framework

**Purpose**: Template for analysis workflows.

**Lifecycle**:
```
analysis.run()
    ↓
1. _prepare()           # Setup (once)
    ↓
2. for frame in trajectory:
       _single_frame()   # Process frame
    ↓
3. _conclude()          # Finalize (once)
    ↓
results ready
```

**Benefits**:
- Consistent API across analyses
- Automatic frame slicing (`start`, `stop`, `step`)
- Parallel execution support
- Progress bar integration

**Example structure**:
```python
class MyAnalysis(AnalysisBase):
    def _prepare(self):
        self.results.data = []

    def _single_frame(self):
        # Access current frame via self._ts
        value = calculate_something(self._ts)
        self.results.data.append(value)

    def _conclude(self):
        self.results.data = np.array(self.results.data)
        self.results.mean = self.results.data.mean()
```

#### Key Analysis Modules

**RMSD** (`rms.py`):
- `RMSD`: RMSD over trajectory
- `rmsd()`: Single RMSD calculation

**Alignment** (`align.py`):
- `alignto()`: Align two structures
- `AlignTraj`: Align trajectory to reference

**Distances** (`distances.py`):
- Distance calculations
- Distance arrays
- Self-distance calculations

**RDF** (`rdf.py`):
- `InterRDF`: RDF between two groups
- Binning and normalization

**Contacts** (`contacts.py`):
- Native contacts (Q)
- Contact matrices

### 5. Transformations Module (`/package/MDAnalysis/transformations/`)

**Purpose**: On-the-fly trajectory modifications.

**Architecture**:
Transformations are callables that modify `Timestep` objects.

**Usage**:
```python
u.trajectory.add_transformations(
    unwrap(u.atoms),
    center_in_box(u.select_atoms('protein')),
    wrap(u.atoms)
)

# Now all iteration applies transformations automatically
for ts in u.trajectory:
    # ts is already transformed
    pass
```

**Available transformations**:
- `unwrap()` - Fix periodic boundary jumps
- `wrap()` - Apply periodic boundaries
- `center_in_box()` - Center selection
- `fit_rot_trans()` - Align to reference
- `translate()` - Translate coordinates
- `rotate()` - Rotate coordinates

**Implementation**:
```python
class TransformationBase:
    def __call__(self, ts):
        # Modify ts in-place
        # Must return ts
        return ts
```

### 6. Library Module (`/package/MDAnalysis/lib/`)

Utilities and performance-critical code.

#### Compiled Extensions (Cython)

**`c_distances.pyx`**:
- Fast distance calculations
- Distance arrays
- Self-distances
- Uses vectorized operations

**`qcprot.pyx`**:
- QCP rotation algorithm
- Fast RMSD with optimal rotation
- Used by alignment code

**`nsgrid.pyx`**:
- Grid-based neighbor search
- Fast spatial queries
- O(N) performance for neighbor finding

**`transformations.py`**:
- Coordinate transformation matrices
- Euler angles
- Quaternions
- Rotation/translation utilities

#### Pure Python Utilities

**`log.py`**:
- Logging configuration
- Progress bars (tqdm integration)

**`util.py`**:
- General utilities
- File handling
- Path utilities

**`mdamath.py`**:
- Mathematical functions
- Angle calculations
- Vector operations

---

## Data Flow

### Typical Analysis Workflow

```
1. User creates Universe
   ↓
2. Parse topology file
   → TopologyParser reads file
   → Creates Topology object
   → Creates Atoms/Residues/Segments
   ↓
3. Open trajectory file
   → TrajectoryReader opens file
   → Reads header
   → Positions at frame 0
   ↓
4. User selects atoms
   → Selection parser tokenizes string
   → Creates selection objects
   → Evaluates to indices
   → Returns AtomGroup
   ↓
5. User iterates trajectory
   → For each frame:
     - Reader loads coordinates
     - Updates Timestep
     - Applies transformations (if any)
     - Yields to user
   ↓
6. User accesses properties
   → AtomGroup.positions
   → Gets indices from AtomGroup
   → Looks up in Timestep.positions
   → Returns NumPy array (view)
   ↓
7. User performs calculations
   → NumPy operations
   → Compiled functions (if needed)
   → Store results
```

### Memory Layout

**During iteration**:
```
Memory:
  ├── Topology (small, persistent)
  │   ├── Atom metadata
  │   ├── Residue metadata
  │   └── Bonds/attributes
  │
  ├── Timestep (medium, reused each frame)
  │   ├── positions[n_atoms, 3]
  │   ├── velocities[n_atoms, 3] (if present)
  │   └── forces[n_atoms, 3] (if present)
  │
  └── User analysis data (variable)

Disk:
  └── Trajectory file (streamed, not in memory)
```

**Key point**: Only one frame in memory at a time (unless user explicitly copies).

---

## Plugin Architecture

### Format Reader Registration

**Automatic discovery**:
1. Reader class defines `format` attribute
2. On import, registers with format registry
3. `Universe` looks up reader by file extension

**Example**:
```python
# In DCD.py
class DCDReader(ReaderBase):
    format = 'DCD'
    # ...

# In coordinates/__init__.py
_READERS = {}  # format → reader class

def register_reader(format_name, reader_class):
    _READERS[format_name.upper()] = reader_class

# Readers auto-register on import
```

### Adding New Formats

To add a new format:

1. Create reader class:
```python
# myformat.py
from MDAnalysis.coordinates.base import ReaderBase

class MYFORMATReader(ReaderBase):
    format = 'MYFORMAT'

    def __init__(self, filename):
        super().__init__(filename)
        # Open file, parse header

    def _read_next_timestep(self):
        # Read frame data
        # Populate self.ts
        return self.ts
```

2. Import in `coordinates/__init__.py`:
```python
from .myformat import MYFORMATReader
```

3. Done! Now `Universe('file.myformat')` works.

---

## Performance Optimizations

### 1. Compiled Extensions

**Strategy**: Identify bottlenecks, rewrite in Cython/C.

**Hot paths**:
- Distance calculations → `c_distances.pyx`
- RMSD calculations → `qcprot.pyx`
- Neighbor searches → `nsgrid.pyx`
- Timestep operations → `timestep.pyx`

**Typical speedup**: 10-100x over pure Python.

### 2. NumPy Vectorization

**Strategy**: Use NumPy array operations instead of Python loops.

```python
# Slow (Python loop)
for i in range(len(positions)):
    distances[i] = np.linalg.norm(positions[i] - ref)

# Fast (vectorized)
distances = np.linalg.norm(positions - ref, axis=1)
```

### 3. Memory Views

**Strategy**: Avoid copying data when possible.

```python
# AtomGroup.positions returns a view into Timestep._pos
# No copy made!
positions = atomgroup.positions

# Slicing also returns views
ca_positions = atomgroup[:10].positions  # Still a view
```

### 4. Lazy Evaluation

**Strategy**: Defer computation until needed.

**Examples**:
- Topology attributes computed on demand
- Selections cached until invalidated
- Bond guessing only if bonds accessed

### 5. Parallel Processing

**Available in some analyses**:
```python
analysis.run(n_jobs=4)  # Use 4 CPU cores
```

**Implementation**: Uses `joblib` for parallel frame processing.

---

## Testing Infrastructure

### Organization

```
testsuite/
├── MDAnalysisTests/         # Main test package
│   ├── analysis/            # Analysis module tests
│   ├── coordinates/         # Reader/writer tests
│   ├── core/                # Core component tests
│   ├── topology/            # Topology parser tests
│   └── ...
│
└── data/                    # Test data files
```

### Test Data

**Location**: `/testsuite/MDAnalysisTests/data/`

**Access**: Via `datafiles.py` module
```python
from MDAnalysis.tests.datafiles import PSF, DCD
```

**Types**:
- Small trajectories for fast tests
- Edge cases (empty files, corrupted data)
- Format examples from various packages

### Test Types

1. **Unit tests**: Individual functions/methods
2. **Integration tests**: Component interactions
3. **Regression tests**: Prevent bug reintroduction
4. **Format tests**: File I/O correctness

### Running Tests

```bash
# All tests
pytest testsuite/

# Specific module
pytest testsuite/MDAnalysisTests/analysis/test_rms.py

# With coverage
pytest --cov=MDAnalysis testsuite/
```

### Continuous Integration

**Platforms tested**:
- Linux (Ubuntu latest)
- macOS
- Windows

**Python versions**: 3.11, 3.12, 3.13, 3.14

**Checks**:
- Unit tests
- Code coverage
- Linting (pylint, flake8)
- Type checking (mypy)
- Documentation builds

---

## Build and Distribution

### Build Process

1. **Cython compilation**:
   - `.pyx` files → `.c` files (Cython)
   - `.c` files → `.so` shared libraries (C compiler)

2. **Python package**:
   - Pure Python modules copied
   - Compiled extensions included

### Build Tools

- **setuptools**: Package building
- **Cython**: Compile `.pyx` files
- **NumPy**: Build-time dependency (for C headers)

### Build Commands

```bash
# Development install (editable)
cd package
pip install -e .

# Build wheel
python -m build

# Install from wheel
pip install dist/MDAnalysis-*.whl
```

### Dependencies

**Runtime**:
- NumPy >= 1.26.0
- SciPy >= 1.5.0
- matplotlib >= 1.5.1
- GridDataFormats
- mmtf-python
- joblib
- tqdm

**Build-time**:
- Cython >= 0.28
- NumPy >= 2.0 (for building)
- C compiler

**Optional**:
- netCDF4, h5py (advanced formats)
- chemfiles (format support)
- parmed (conversions)

### Distribution

**PyPI**: `pip install MDAnalysis`

**conda-forge**: `conda install -c conda-forge mdanalysis`

---

## Code Organization Best Practices

### File Naming Conventions

- **Modules**: lowercase with underscores (`atom_group.py`)
- **Classes**: PascalCase (`AtomGroup`)
- **Functions**: lowercase with underscores (`center_of_mass()`)
- **Constants**: UPPERCASE with underscores (`TRAJECTORY_FORMAT`)

### Import Organization

```python
# Standard library
import os
import sys

# Third-party
import numpy as np
from scipy import spatial

# MDAnalysis
from ..core.groups import AtomGroup
from . import util
```

### Documentation

**Docstring format**: NumPy style

```python
def function(arg1, arg2):
    """
    Brief description.

    Longer description if needed.

    Parameters
    ----------
    arg1 : type
        Description of arg1
    arg2 : type
        Description of arg2

    Returns
    -------
    type
        Description of return value
    """
```

---

## Summary

MDAnalysis architecture emphasizes:

1. **Modularity**: Clean separation of components
2. **Extensibility**: Easy to add new formats and analyses
3. **Performance**: Strategic use of compiled code
4. **Usability**: High-level Python API with NumPy integration
5. **Robustness**: Comprehensive testing and CI

The **Universe-AtomGroup** model provides a powerful and intuitive abstraction for MD analysis, while the **plugin architecture** ensures flexibility and extensibility.

For more details, see:
- **User Guide**: https://userguide.mdanalysis.org
- **API Documentation**: https://docs.mdanalysis.org
- **Developer Guide**: https://userguide.mdanalysis.org/contributing.html
