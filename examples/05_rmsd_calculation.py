#!/usr/bin/env python
"""
Example 5: RMSD Calculation

This script demonstrates:
- Calculating RMSD over a trajectory
- Using the RMSD analysis class
- RMSD relative to different references
- Visualizing RMSD data
"""

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms
from MDAnalysis.tests.datafiles import PSF, DCD

def main():
    print("=" * 60)
    print("Example 5: RMSD Calculation")
    print("=" * 60)

    # Load universe
    u = mda.Universe(PSF, DCD)
    print(f"\nLoaded system with {len(u.atoms)} atoms")
    print(f"Trajectory has {len(u.trajectory)} frames")

    # 1. Basic RMSD calculation (backbone atoms)
    print("\n1. RMSD of backbone atoms relative to first frame:")
    R = rms.RMSD(u, select='backbone', ref_frame=0)
    R.run()

    # Results stored in R.results.rmsd
    # Format: [frame, time (ps), RMSD (Å)]
    rmsd_data = R.results.rmsd
    print(f"   Shape of results: {rmsd_data.shape}")
    print(f"   Columns: [frame, time, RMSD]")
    print(f"\n   First 5 frames:")
    for i in range(5):
        frame, time, rmsd_value = rmsd_data[i]
        print(f"   Frame {int(frame):3d}: Time {time:6.1f} ps, RMSD {rmsd_value:6.3f} Å")

    # 2. RMSD statistics
    print("\n2. RMSD Statistics (entire trajectory):")
    rmsd_values = rmsd_data[:, 2]  # Extract RMSD column
    print(f"   Mean RMSD: {rmsd_values.mean():.3f} Å")
    print(f"   Std RMSD:  {rmsd_values.std():.3f} Å")
    print(f"   Min RMSD:  {rmsd_values.min():.3f} Å (frame {rmsd_values.argmin()})")
    print(f"   Max RMSD:  {rmsd_values.max():.3f} Å (frame {rmsd_values.argmax()})")

    # 3. RMSD for different atom selections
    print("\n3. RMSD for different atom selections:")

    selections = {
        'C-alpha': 'name CA',
        'Backbone': 'backbone',
        'Protein': 'protein',
        'Core residues': 'resid 10-100 and name CA'
    }

    for name, selection in selections.items():
        R = rms.RMSD(u, select=selection, ref_frame=0)
        R.run()
        final_rmsd = R.results.rmsd[-1, 2]
        mean_rmsd = R.results.rmsd[:, 2].mean()
        print(f"   {name:15s}: Mean = {mean_rmsd:.3f} Å, Final = {final_rmsd:.3f} Å")

    # 4. RMSD relative to average structure
    print("\n4. RMSD relative to average structure:")

    # First, calculate average structure
    u.trajectory[0]  # Reset to first frame
    avg_positions = np.zeros_like(u.atoms.positions)

    for ts in u.trajectory:
        avg_positions += u.atoms.positions
    avg_positions /= len(u.trajectory)

    # Create a reference universe with average structure
    ref = mda.Universe(PSF)
    ref.atoms.positions = avg_positions

    # Calculate RMSD relative to average
    R_avg = rms.RMSD(u, reference=ref, select='name CA')
    R_avg.run()

    print(f"   Mean RMSD to average: {R_avg.results.rmsd[:, 2].mean():.3f} Å")

    # 5. Frame-to-frame RMSD
    print("\n5. Frame-to-frame RMSD (consecutive frames):")

    frame_to_frame = []
    ca = u.select_atoms('name CA')

    for i, ts in enumerate(u.trajectory[:-1]):
        pos1 = ca.positions.copy()
        u.trajectory[i+1]
        pos2 = ca.positions

        # Simple RMSD calculation (after centering)
        pos1_centered = pos1 - pos1.mean(axis=0)
        pos2_centered = pos2 - pos2.mean(axis=0)
        rmsd = np.sqrt(((pos1_centered - pos2_centered) ** 2).sum() / len(ca))
        frame_to_frame.append(rmsd)

    frame_to_frame = np.array(frame_to_frame)
    print(f"   Mean frame-to-frame RMSD: {frame_to_frame.mean():.3f} Å")
    print(f"   Max frame-to-frame RMSD: {frame_to_frame.max():.3f} Å")

    # 6. Find most similar frames
    print("\n6. Find frames most similar to first frame:")

    u.trajectory[0]
    R = rms.RMSD(u, select='name CA', ref_frame=0)
    R.run()

    rmsd_values = R.results.rmsd[:, 2]
    most_similar_indices = np.argsort(rmsd_values)[:5]

    print(f"   5 most similar frames to frame 0:")
    for idx in most_similar_indices:
        print(f"   Frame {idx:3d}: RMSD = {rmsd_values[idx]:.3f} Å")

    # 7. RMSD convergence check
    print("\n7. RMSD convergence (last 10% of trajectory):")

    cutoff = int(0.9 * len(rmsd_values))
    last_10_percent = rmsd_values[cutoff:]

    print(f"   Mean RMSD (last 10%): {last_10_percent.mean():.3f} Å")
    print(f"   Std RMSD (last 10%):  {last_10_percent.std():.3f} Å")

    if last_10_percent.std() < 0.5:
        print(f"   ✓ Structure appears converged (std < 0.5 Å)")
    else:
        print(f"   ✗ Structure may still be equilibrating")

    # 8. Save RMSD data
    print("\n8. Saving RMSD data:")
    output_file = '/tmp/rmsd_data.txt'
    np.savetxt(output_file, R.results.rmsd,
               header='Frame Time(ps) RMSD(Angstrom)',
               fmt='%10.1f %10.2f %10.4f')
    print(f"   RMSD data saved to {output_file}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("Tip: You can plot RMSD data with matplotlib:")
    print("  import matplotlib.pyplot as plt")
    print("  plt.plot(R.results.rmsd[:, 1], R.results.rmsd[:, 2])")
    print("  plt.xlabel('Time (ps)')")
    print("  plt.ylabel('RMSD (Å)')")
    print("  plt.show()")
    print("=" * 60)

if __name__ == '__main__':
    main()
