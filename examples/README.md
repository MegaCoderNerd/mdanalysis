# MDAnalysis Examples

This directory contains example scripts demonstrating how to use MDAnalysis for common molecular dynamics analysis tasks.

## Getting Started

Before running these examples, make sure you have MDAnalysis installed:

```bash
pip install MDAnalysis
```

## Available Examples

### Basic Examples

1. **01_basic_usage.py** - Introduction to Universe, AtomGroup, and basic selections
2. **02_trajectory_iteration.py** - Iterating through trajectory frames
3. **03_selections.py** - Advanced atom selection examples
4. **04_accessing_data.py** - Accessing coordinates, velocities, and properties

### Analysis Examples

5. **05_rmsd_calculation.py** - Calculate RMSD over trajectory
6. **06_structure_alignment.py** - Align structures and trajectories
7. **07_distance_analysis.py** - Calculate distances between atom groups
8. **08_rdf_calculation.py** - Radial distribution function
9. **09_hydrogen_bonds.py** - Hydrogen bond analysis
10. **10_contact_analysis.py** - Native contacts and contact maps

### Advanced Examples

11. **11_custom_analysis.py** - Writing custom analysis classes
12. **12_transformations.py** - On-the-fly trajectory transformations
13. **13_density_analysis.py** - Density calculations and grid export
14. **14_pca_analysis.py** - Principal component analysis
15. **15_parallel_analysis.py** - Parallel analysis with multiple cores

## Running Examples

Each example is self-contained and can be run directly:

```bash
python 01_basic_usage.py
```

Most examples use built-in test data from MDAnalysis, so you don't need to provide your own files to get started.

## Using Your Own Data

To analyze your own simulation data, simply replace the test data paths with your file paths:

```python
# Instead of:
from MDAnalysis.tests.datafiles import PSF, DCD
u = mda.Universe(PSF, DCD)

# Use:
u = mda.Universe('your_topology.pdb', 'your_trajectory.dcd')
```

## Example Template

Here's a template for creating your own analysis script:

```python
#!/usr/bin/env python
"""
Brief description of what this script does.
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD

def main():
    # Load universe
    u = mda.Universe(PSF, DCD)

    # Select atoms
    protein = u.select_atoms('protein')

    # Perform analysis
    for ts in u.trajectory:
        # Your analysis here
        pass

    # Save/display results
    print("Analysis complete!")

if __name__ == '__main__':
    main()
```

## Further Resources

- **Quick Start Guide**: See `../QUICKSTART.md` for comprehensive overview
- **User Guide**: https://userguide.mdanalysis.org
- **API Docs**: https://docs.mdanalysis.org
- **Tutorials**: https://userguide.mdanalysis.org/examples/README.html

## Getting Help

- **GitHub Discussions**: https://github.com/MDAnalysis/mdanalysis/discussions
- **Issue Tracker**: https://github.com/MDAnalysis/mdanalysis/issues
