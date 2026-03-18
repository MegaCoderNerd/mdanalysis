#!/usr/bin/env python
"""
Example 2: Trajectory Iteration

This script demonstrates:
- Iterating through trajectory frames
- Accessing frame information
- Slicing trajectories
- Tracking properties over time
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD

def main():
    print("=" * 60)
    print("Example 2: Trajectory Iteration")
    print("=" * 60)

    # Load universe
    u = mda.Universe(PSF, DCD)
    ca_atoms = u.select_atoms('name CA')

    print(f"\nLoaded trajectory with {len(u.trajectory)} frames")

    # Example 1: Iterate through all frames
    print("\n1. Iterating through first 5 frames:")
    for i, ts in enumerate(u.trajectory[:5]):
        print(f"   Frame {ts.frame}: Time = {ts.time:.2f} ps, "
              f"Timestep = {ts.dt:.2f} ps")

    # Example 2: Access specific frame
    print("\n2. Jump to specific frame (frame 50):")
    u.trajectory[50]
    ts = u.trajectory.ts
    print(f"   Current frame: {ts.frame}")
    print(f"   Current time: {ts.time:.2f} ps")
    print(f"   Box dimensions: {ts.dimensions[:3]}")

    # Example 3: Slice trajectory
    print("\n3. Slicing trajectory (every 10th frame):")
    frame_numbers = []
    for ts in u.trajectory[::10]:
        frame_numbers.append(ts.frame)
    print(f"   Processed frames: {frame_numbers}")

    # Example 4: Track property over time
    print("\n4. Tracking center of mass over trajectory:")
    com_trajectory = []
    times = []

    for ts in u.trajectory:
        com = ca_atoms.center_of_mass()
        com_trajectory.append(com)
        times.append(ts.time)

    com_trajectory = np.array(com_trajectory)
    times = np.array(times)

    print(f"   Initial COM (t={times[0]:.1f} ps): {com_trajectory[0]}")
    print(f"   Final COM (t={times[-1]:.1f} ps): {com_trajectory[-1]}")

    # Calculate COM displacement
    displacement = np.linalg.norm(com_trajectory[-1] - com_trajectory[0])
    print(f"   Total COM displacement: {displacement:.2f} Å")

    # Example 5: Calculate average property
    print("\n5. Calculate time-averaged radius of gyration:")
    rgyr_values = []

    for ts in u.trajectory:
        rgyr = ca_atoms.radius_of_gyration()
        rgyr_values.append(rgyr)

    rgyr_values = np.array(rgyr_values)
    print(f"   Mean Rgyr: {rgyr_values.mean():.2f} ± {rgyr_values.std():.2f} Å")
    print(f"   Min Rgyr: {rgyr_values.min():.2f} Å (frame {rgyr_values.argmin()})")
    print(f"   Max Rgyr: {rgyr_values.max():.2f} Å (frame {rgyr_values.argmax()})")

    # Example 6: Conditional frame processing
    print("\n6. Find frames where Rgyr > mean + std:")
    threshold = rgyr_values.mean() + rgyr_values.std()
    expanded_frames = np.where(rgyr_values > threshold)[0]
    print(f"   Threshold: {threshold:.2f} Å")
    print(f"   Number of expanded frames: {len(expanded_frames)}")
    if len(expanded_frames) > 0:
        print(f"   First few expanded frames: {expanded_frames[:5]}")

    # Example 7: Trajectory metadata
    print("\n7. Trajectory metadata:")
    print(f"   Total frames: {u.trajectory.n_frames}")
    print(f"   Timestep (dt): {u.trajectory.dt:.2f} ps")
    print(f"   Total time: {u.trajectory.totaltime:.2f} ps")
    print(f"   Has velocities: {u.trajectory.ts.has_velocities}")
    print(f"   Has forces: {u.trajectory.ts.has_forces}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)

if __name__ == '__main__':
    main()
