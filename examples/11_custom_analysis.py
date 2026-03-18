#!/usr/bin/env python
"""
Example 11: Writing Custom Analysis Classes

This script demonstrates:
- Creating a custom analysis class
- Using the AnalysisBase framework
- Storing and accessing results
- Adding custom options and parameters
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.tests.datafiles import PSF, DCD

# Example 1: Simple custom analysis
class EndToEndDistance(AnalysisBase):
    """
    Calculate end-to-end distance of a polymer or protein.

    The end-to-end distance is calculated as the distance between
    the first and last atom in the selection.
    """

    def __init__(self, atomgroup, **kwargs):
        """
        Parameters
        ----------
        atomgroup : AtomGroup
            Atoms for which to calculate end-to-end distance
        """
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def _prepare(self):
        """Set up data structures before iteration."""
        # This runs once before iterating through frames
        self.results.distances = []
        self.results.times = []

    def _single_frame(self):
        """Calculate data for a single frame."""
        # This runs for each frame in the trajectory
        first_atom = self.atomgroup[0].position
        last_atom = self.atomgroup[-1].position
        distance = np.linalg.norm(last_atom - first_atom)

        self.results.distances.append(distance)
        self.results.times.append(self._ts.time)

    def _conclude(self):
        """Finalize results after iteration."""
        # This runs once after all frames are processed
        self.results.distances = np.array(self.results.distances)
        self.results.times = np.array(self.results.times)
        self.results.mean = np.mean(self.results.distances)
        self.results.std = np.std(self.results.distances)


# Example 2: More complex custom analysis
class BackboneDistanceMap(AnalysisBase):
    """
    Calculate a distance map between C-alpha atoms.

    For each frame, calculates the pairwise distances between
    all C-alpha atoms and stores them.
    """

    def __init__(self, universe, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.ca_atoms = universe.select_atoms('name CA')
        self.n_atoms = len(self.ca_atoms)

    def _prepare(self):
        """Prepare data structures."""
        # Store distance maps for each frame
        self.results.distance_maps = []

    def _single_frame(self):
        """Calculate distance map for current frame."""
        positions = self.ca_atoms.positions

        # Calculate pairwise distance matrix
        diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
        distances = np.sqrt((diff ** 2).sum(axis=-1))

        self.results.distance_maps.append(distances)

    def _conclude(self):
        """Calculate average and standard deviation."""
        self.results.distance_maps = np.array(self.results.distance_maps)

        # Calculate average distance map over all frames
        self.results.average_map = np.mean(self.results.distance_maps, axis=0)
        self.results.std_map = np.std(self.results.distance_maps, axis=0)


# Example 3: Analysis with customizable parameters
class RadiusOfGyrationByResidue(AnalysisBase):
    """
    Calculate radius of gyration for specific residue types.

    Parameters can be customized to select different residue types.
    """

    def __init__(self, universe, residue_names=None, **kwargs):
        """
        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe
        residue_names : list of str, optional
            List of residue names to analyze (default: all)
        """
        super().__init__(universe.trajectory, **kwargs)
        self.universe = universe

        if residue_names is None:
            self.selection = universe.select_atoms('protein')
            self.label = 'all protein'
        else:
            resname_str = ' '.join(residue_names)
            self.selection = universe.select_atoms(f'resname {resname_str}')
            self.label = resname_str

    def _prepare(self):
        """Initialize result arrays."""
        self.results.rgyr = []
        self.results.times = []

    def _single_frame(self):
        """Calculate Rgyr for current frame."""
        rgyr = self.selection.radius_of_gyration()
        self.results.rgyr.append(rgyr)
        self.results.times.append(self._ts.time)

    def _conclude(self):
        """Calculate statistics."""
        self.results.rgyr = np.array(self.results.rgyr)
        self.results.times = np.array(self.results.times)
        self.results.mean = np.mean(self.results.rgyr)
        self.results.std = np.std(self.results.rgyr)


def main():
    print("=" * 60)
    print("Example 11: Writing Custom Analysis Classes")
    print("=" * 60)

    # Load universe
    u = mda.Universe(PSF, DCD)

    # Example 1: End-to-end distance
    print("\n1. End-to-End Distance Analysis:")
    print("   Analyzing protein backbone...")

    backbone = u.select_atoms('backbone')
    e2e = EndToEndDistance(backbone)
    e2e.run()

    print(f"   Mean end-to-end distance: {e2e.results.mean:.2f} ± {e2e.results.std:.2f} Å")
    print(f"   Min distance: {e2e.results.distances.min():.2f} Å")
    print(f"   Max distance: {e2e.results.distances.max():.2f} Å")

    # Example 2: Run with frame slicing
    print("\n2. Running analysis on subset of frames:")
    e2e_subset = EndToEndDistance(backbone)
    e2e_subset.run(start=0, stop=50, step=5)  # Frames 0-50, every 5th frame

    print(f"   Analyzed {len(e2e_subset.results.distances)} frames")
    print(f"   Mean distance: {e2e_subset.results.mean:.2f} Å")

    # Example 3: Distance map analysis
    print("\n3. Backbone Distance Map Analysis:")
    print("   Calculating C-alpha distance maps...")

    dist_map = BackboneDistanceMap(u)
    dist_map.run(stop=10, verbose=True)  # Just first 10 frames for speed

    print(f"   Number of C-alpha atoms: {dist_map.n_atoms}")
    print(f"   Average distance map shape: {dist_map.results.average_map.shape}")
    print(f"   Average distance (all pairs): {dist_map.results.average_map.mean():.2f} Å")

    # Example 4: Analysis with custom parameters
    print("\n4. Radius of Gyration by Residue Type:")

    # Analyze hydrophobic residues
    hydrophobic = RadiusOfGyrationByResidue(
        u,
        residue_names=['ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET']
    )
    hydrophobic.run()

    print(f"   Hydrophobic residues:")
    print(f"   Mean Rgyr: {hydrophobic.results.mean:.2f} ± {hydrophobic.results.std:.2f} Å")

    # Analyze charged residues
    charged = RadiusOfGyrationByResidue(
        u,
        residue_names=['ARG', 'LYS', 'ASP', 'GLU']
    )
    charged.run()

    print(f"\n   Charged residues:")
    print(f"   Mean Rgyr: {charged.results.mean:.2f} ± {charged.results.std:.2f} Å")

    # Example 5: Accessing and saving results
    print("\n5. Accessing and Saving Results:")

    # Results are stored in the results attribute
    print(f"   Results attributes: {dir(e2e.results)}")
    print(f"   Distance array shape: {e2e.results.distances.shape}")
    print(f"   Time array shape: {e2e.results.times.shape}")

    # Save to file
    output_file = '/tmp/e2e_distances.txt'
    data = np.column_stack([e2e.results.times, e2e.results.distances])
    np.savetxt(output_file, data,
               header='Time(ps) Distance(Angstrom)',
               fmt='%10.2f %10.4f')
    print(f"   Data saved to {output_file}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("\nKey points for custom analysis classes:")
    print("1. Inherit from AnalysisBase")
    print("2. Implement _prepare(), _single_frame(), and _conclude()")
    print("3. Store results in self.results")
    print("4. Use run() method with start, stop, step parameters")
    print("5. Access frame data via self._ts")
    print("=" * 60)

if __name__ == '__main__':
    main()
