"""
Example: Cigarettes vs nicotine patch - delivery kinetics matter.

Demonstrates:
- Spike (cigarette) vs sustained (patch) delivery
- Why spikes are more reinforcing (dopamine dynamics)
- Control theory perspective: input shaping affects outcomes
"""

import sys
sys.path.append('../src')

from model import NicotineRewardModel, bolus_input, continuous_input
from visualization import Visualizer
import matplotlib.pyplot as plt
import numpy as np


def main():
    print("=" * 60)
    print("EXAMPLE: CIGARETTES VS PATCH")
    print("Delivery Profile Matters for Reinforcement")
    print("=" * 60)
    
    model = NicotineRewardModel()
    
    # Scenario 1: Cigarette (rapid bolus)
    print("\nScenario 1: Cigarette (rapid spike)...")
    
    def I_cigarette(t):
        """One cigarette per hour for 8 hours"""
        for i in range(8):
            t_start = i * 60
            if t_start <= t < t_start + 5:
                return 5.0  # High rate, short duration
        return 0.0
    
    result_cig = model.simulate(
        t_span=(0, 480),  # 8 hours
        input_func=I_cigarette,
        method='LSODA'
    )
    
    # Scenario 2: Nicotine patch (sustained delivery)
    print("Scenario 2: Nicotine patch (sustained release)...")
    
    # Match total nicotine delivered (area under curve)
    # 8 cigarettes × 5.0 × 5 min = 200 μM·min
    # Patch: 200 / 480 min = 0.417 μM/min
    
    I_patch_rate = 0.417
    result_patch = model.simulate(
        t_span=(0, 480),
        input_func=continuous_input(I_patch_rate),
        method='LSODA'
    )
    
    # Analysis: Peak dopamine (reinforcement proxy)
    print("\n" + "-" * 60)
    print("REINFORCEMENT METRICS")
    print("-" * 60)
    
    D_peak_cig = result_cig['D'].max()
    D_peak_patch = result_patch['D'].max()
    
    print(f"\nPeak dopamine:")
    print(f"  Cigarettes: {D_peak_cig:.3f} AU")
    print(f"  Patch:      {D_peak_patch:.3f} AU")
    print(f"  Ratio:      {D_peak_cig/D_peak_patch:.2f}×")
    
    # Dopamine variance (spikiness)
    D_var_cig = np.var(result_cig['D'])
    D_var_patch = np.var(result_patch['D'])
    
    print(f"\nDopamine variability (variance):")
    print(f"  Cigarettes: {D_var_cig:.4f}")
    print(f"  Patch:      {D_var_patch:.4f}")
    print(f"  Ratio:      {D_var_cig/D_var_patch:.1f}×")
    
    # Average levels
    D_mean_cig = np.mean(result_cig['D'])
    D_mean_patch = np.mean(result_patch['D'])
    
    print(f"\nAverage dopamine:")
    print(f"  Cigarettes: {D_mean_cig:.3f} AU")
    print(f"  Patch:      {D_mean_patch:.3f} AU")
    
    # Desensitisation
    print(f"\nFinal desensitised fraction:")
    print(f"  Cigarettes: {result_cig['R_d'][-1]:.3f}")
    print(f"  Patch:      {result_patch['R_d'][-1]:.3f}")
    
    # Visualization
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATIONS")
    print("=" * 60)
    
    # Create comparison figure
    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    
    t_hours = result_cig['t'] / 60
    
    # Nicotine profiles
    axes[0].plot(t_hours, result_cig['N'], 
                color='#D62828', linewidth=2.5, label='Cigarettes (spikes)')
    axes[0].plot(t_hours, result_patch['N'], 
                color='#2E86AB', linewidth=2.5, linestyle='--', label='Patch (sustained)')
    axes[0].set_ylabel('Nicotine [μM]')
    axes[0].set_title('INPUT: Delivery Profile', fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Dopamine response (KEY for reinforcement)
    axes[1].plot(t_hours, result_cig['D'], 
                color='#D62828', linewidth=2.5, label='Cigarettes', alpha=0.8)
    axes[1].plot(t_hours, result_patch['D'], 
                color='#2E86AB', linewidth=2.5, linestyle='--', label='Patch', alpha=0.8)
    axes[1].axhline(model.params['D_0'], color='gray', linestyle=':', alpha=0.5)
    axes[1].fill_between(t_hours, model.params['D_0'], result_cig['D'],
                        where=(result_cig['D'] > model.params['D_0']),
                        alpha=0.2, color='red', label='Cig reinforcement')
    axes[1].set_ylabel('Dopamine [AU]')
    axes[1].set_title('OUTPUT: Dopamine Response → Reinforcement Strength', 
                     fontweight='bold', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Receptor activation
    axes[2].plot(t_hours, result_cig['R_a'], 
                color='#D62828', linewidth=2, label='Cigarettes')
    axes[2].plot(t_hours, result_patch['R_a'], 
                color='#2E86AB', linewidth=2, linestyle='--', label='Patch')
    axes[2].set_ylabel('Active Receptors')
    axes[2].set_title('Receptor Activation Dynamics', fontweight='bold')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    # Desensitisation
    axes[3].plot(t_hours, result_cig['R_d'], 
                color='#D62828', linewidth=2, label='Cigarettes')
    axes[3].plot(t_hours, result_patch['R_d'], 
                color='#2E86AB', linewidth=2, linestyle='--', label='Patch')
    axes[3].set_ylabel('Desensitised\nReceptors')
    axes[3].set_xlabel('Time [hours]')
    axes[3].set_title('Desensitisation (Tolerance)', fontweight='bold')
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.suptitle('Cigarettes vs Patch: Why Delivery Kinetics Matter', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('patch_vs_cigarette.png', dpi=300, bbox_inches='tight')
    print("\n✓ Saved: patch_vs_cigarette.png")
    
    # Dopamine comparison (detailed)
    fig2, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(t_hours, result_cig['D'], 
           color='#D62828', linewidth=3, label='Cigarettes', alpha=0.9)
    ax.plot(t_hours, result_patch['D'], 
           color='#2E86AB', linewidth=3, linestyle='--', label='Patch', alpha=0.9)
    ax.axhline(model.params['D_0'], color='gray', linestyle=':', 
              linewidth=2, label='Baseline')
    
    # Annotate peaks
    for i in range(8):
        t_peak_idx = np.argmax(result_cig['D'][i*60:(i+1)*60]) + i*60
        if t_peak_idx < len(result_cig['D']):
            ax.annotate('', xy=(result_cig['t'][t_peak_idx]/60, result_cig['D'][t_peak_idx]),
                       xytext=(result_cig['t'][t_peak_idx]/60, model.params['D_0']),
                       arrowprops=dict(arrowstyle='<->', color='red', lw=1.5, alpha=0.5))
    
    ax.set_xlabel('Time [hours]')
    ax.set_ylabel('Dopamine [AU]')
    ax.set_title('Dopamine Spikes → Reinforcement', fontweight='bold', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('dopamine_comparison.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: dopamine_comparison.png")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("KEY OBSERVATIONS")
    print("=" * 60)
    print(f"1. Cigarettes produce {D_peak_cig/D_peak_patch:.1f}× higher dopamine peaks")
    print(f"2. Spiky dopamine = strong reinforcement (learning signal)")
    print(f"3. Patch has {D_var_cig/D_var_patch:.0f}× less variability → weaker reinforcement")
    print("4. Same total nicotine, VERY different behavioural outcomes")
    print("=" * 60)
    
    print("\n💡 Control theory insight:")
    print("Input shaping is a control strategy. Patches smooth the input,")
    print("reducing peak responses and breaking the reinforcement cycle.")
    print("\nThis is why patches aid cessation but don't replicate the 'hit'.")


if __name__ == "__main__":
    main()
