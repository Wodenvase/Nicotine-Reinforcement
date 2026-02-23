"""
Example: Single cigarette dose response.

Demonstrates:
- Fast nicotine kinetics
- Transient receptor activation
- Dopamine surge
- Return to baseline
"""

import sys
sys.path.append('../src')

from model import NicotineRewardModel, bolus_input
from visualization import Visualizer
import matplotlib.pyplot as plt


def main():
    print("=" * 60)
    print("EXAMPLE: SINGLE CIGARETTE DOSE")
    print("=" * 60)
    
    # Create model
    model = NicotineRewardModel()
    
    # Single cigarette: rapid bolus delivery
    # Typical cigarette delivers ~1-2 mg nicotine, peak brain conc ~0.2-0.5 μM
    I = bolus_input(
        amplitude=5.0,  # Fast input rate during smoking
        t_start=10,     # Starts at 10 min
        duration=5      # Smoking duration ~5 min
    )
    
    # Simulate for 5 hours
    print("\nSimulating 5-hour period...")
    result = model.simulate(
        t_span=(0, 300),
        input_func=I
    )
    
    # Analysis
    print("\nPeak responses:")
    print(f"  Nicotine:  {result['N'].max():.3f} μM")
    print(f"  Active receptors: {result['R_a'].max():.3f}")
    print(f"  Dopamine: {result['D'].max():.3f} AU")
    
    # Time to peak
    t_peak_D = result['t'][result['D'].argmax()]
    print(f"\nTime to peak dopamine: {t_peak_D:.1f} min")
    
    # Half-life (nicotine clearance)
    N_peak = result['N'].max()
    idx_peak = result['N'].argmax()
    idx_half = idx_peak + (result['N'][idx_peak:] < N_peak/2).argmax()
    t_half = result['t'][idx_half] - result['t'][idx_peak]
    print(f"Nicotine half-life: {t_half:.1f} min")
    
    # Visualization
    print("\nGenerating plots...")
    vis = Visualizer()
    
    # Time series
    fig1 = vis.plot_timeseries(
        result,
        variables=['N', 'R_a', 'R_d', 'D'],
        title='Single Cigarette: Acute Response',
        save_path='single_dose_timeseries.png'
    )
    
    # Phase portrait
    fig2 = vis.plot_phase_portrait_2d(
        result,
        var1='R_a',
        var2='D',
        title='Phase Portrait: Receptor Activation → Dopamine',
        save_path='single_dose_phase.png'
    )
    
    print("\n✓ Plots saved:")
    print("  - single_dose_timeseries.png")
    print("  - single_dose_phase.png")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("Key observations:")
    print("=" * 60)
    print("1. Fast nicotine kinetics (~2h half-life)")
    print("2. Rapid dopamine surge (reinforcement window)")
    print("3. Receptor desensitisation develops transiently")
    print("4. Full recovery to baseline (no long-term changes)")
    print("=" * 60)


if __name__ == "__main__":
    main()
