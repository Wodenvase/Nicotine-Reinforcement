"""
Example: Chronic nicotine exposure and steady-state regime.

Demonstrates:
- Receptor upregulation over time
- Elevated baseline dopamine (allostasis)
- Tolerance development
- Qualitative shift in system dynamics
"""

import sys
sys.path.append('../src')

from model import NicotineRewardModel, continuous_input
from visualization import Visualizer
from stability import StabilityAnalyzer
import matplotlib.pyplot as plt
import numpy as np


def main():
    print("=" * 60)
    print("EXAMPLE: CHRONIC EXPOSURE (30 DAYS)")
    print("=" * 60)
    
    # Create model
    model = NicotineRewardModel()
    
    # Continuous moderate input (like regular smoking or patch)
    I_chronic = 0.3  # μM/min (moderate constant exposure)
    I = continuous_input(I_chronic)
    
    # Simulate for 30 days (43,200 min)
    print(f"\nSimulating 30 days of constant exposure (I = {I_chronic} μM/min)...")
    print("This may take a moment...")
    
    result = model.simulate(
        t_span=(0, 43200),
        input_func=I,
        method='LSODA'  # Efficient for stiff systems
    )
    
    # Convert time to days for plotting
    t_days = result['t'] / (60 * 24)
    
    # Analysis
    print("\n" + "-" * 60)
    print("INITIAL STATE (Day 0)")
    print("-" * 60)
    print(f"  R_T: {result['R_T'][0]:.3f}")
    print(f"  D:   {result['D'][0]:.3f} AU")
    
    print("\n" + "-" * 60)
    print("FINAL STATE (Day 30)")
    print("-" * 60)
    print(f"  R_T: {result['R_T'][-1]:.3f}")
    print(f"  D:   {result['D'][-1]:.3f} AU")
    
    # Receptor upregulation
    upregulation_pct = (result['R_T'][-1] - result['R_T'][0]) / result['R_T'][0] * 100
    print(f"\nReceptor upregulation: {upregulation_pct:.1f}%")
    
    # Dopamine elevation
    dopamine_increase_pct = (result['D'][-1] - result['D'][0]) / result['D'][0] * 100
    print(f"Dopamine elevation: {dopamine_increase_pct:.1f}%")
    
    # Desensitisation
    print(f"\nDesensitised receptor fraction: {result['R_d'][-1]:.3f}")
    print(f"Active receptor fraction: {result['R_a'][-1]:.3f}")
    
    # Stability analysis at steady state
    print("\n" + "-" * 60)
    print("STABILITY ANALYSIS AT STEADY STATE")
    print("-" * 60)
    
    y_ss = np.array([result['N'][-1], result['R_a'][-1], 
                     result['R_d'][-1], result['R_T'][-1], result['D'][-1]])
    
    analyzer = StabilityAnalyzer(model)
    eigen_results = analyzer.eigenanalysis(y_ss, I_chronic)
    ts_results = analyzer.timescale_separation(y_ss, I_chronic)
    
    print(f"System stability: {'STABLE' if eigen_results['is_stable'] else 'UNSTABLE'}")
    print(f"Timescale separation ratio: {ts_results['separation_ratio']:.0f}×")
    
    # Visualization
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATIONS")
    print("=" * 60)
    
    # Create figure with all timescales
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    
    # Nicotine (fast)
    axes[0].plot(t_days, result['N'], color='#2E86AB', linewidth=2)
    axes[0].set_ylabel('Nicotine [μM]')
    axes[0].set_title('FAST: Nicotine Concentration', fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    # Dopamine (fast)
    axes[1].plot(t_days, result['D'], color='#9D4EDD', linewidth=2)
    axes[1].axhline(model.params['D_0'], color='gray', linestyle='--', 
                   alpha=0.5, label='Baseline')
    axes[1].set_ylabel('Dopamine [AU]')
    axes[1].set_title('FAST: Dopamine Response', fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Receptors (medium)
    axes[2].plot(t_days, result['R_a'], color='#06A77D', 
                linewidth=2, label='Active')
    axes[2].plot(t_days, result['R_d'], color='#D62828', 
                linewidth=2, linestyle='--', label='Desensitised')
    axes[2].set_ylabel('Receptor Fraction')
    axes[2].set_title('MEDIUM: Receptor Dynamics', fontweight='bold')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    # Total receptors (slow)
    axes[3].plot(t_days, result['R_T'], color='#F77F00', linewidth=2)
    axes[3].axhline(model.params['R_T_0'], color='gray', linestyle='--',
                   alpha=0.5, label='Initial')
    axes[3].set_ylabel('Total Receptors')
    axes[3].set_xlabel('Time [days]')
    axes[3].set_title('SLOW: Receptor Upregulation (Allostatic Shift)', fontweight='bold')
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.suptitle('Chronic Nicotine Exposure: Multi-Timescale Dynamics', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('chronic_exposure.png', dpi=300, bbox_inches='tight')
    print("\n✓ Saved: chronic_exposure.png")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("KEY OBSERVATIONS")
    print("=" * 60)
    print("1. Receptor upregulation occurs over weeks (allostasis)")
    print("2. Chronic dopamine elevation (new 'baseline')")
    print("3. Persistent receptor desensitisation (tolerance)")
    print("4. System achieves new steady state (dependence regime)")
    print("=" * 60)


if __name__ == "__main__":
    main()
