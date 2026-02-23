"""
Example: Withdrawal dynamics after chronic exposure.

Demonstrates:
- Dopamine crash below baseline (negative affect)
- Slow receptor downregulation
- Extended recovery period
- Why relapse is likely (system instability)
"""

import sys
sys.path.append('../src')

from model import NicotineRewardModel, continuous_input
from visualization import Visualizer
import matplotlib.pyplot as plt
import numpy as np


def main():
    print("=" * 60)
    print("EXAMPLE: WITHDRAWAL DYNAMICS")
    print("=" * 60)
    
    # Create model
    model = NicotineRewardModel()
    
    # Phase 1: Chronic exposure (30 days)
    print("\nPhase 1: Establishing dependence (30 days chronic exposure)...")
    I_chronic = 0.3  # μM/min
    
    result_chronic = model.simulate(
        t_span=(0, 43200),  # 30 days
        input_func=continuous_input(I_chronic),
        method='LSODA'
    )
    
    # Get final state (dependent state)
    y_dependent = np.array([
        result_chronic['N'][-1],
        result_chronic['R_a'][-1],
        result_chronic['R_d'][-1],
        result_chronic['R_T'][-1],
        result_chronic['D'][-1]
    ])
    
    print(f"  Dependent state achieved:")
    print(f"    R_T = {y_dependent[3]:.3f}")
    print(f"    D = {y_dependent[4]:.3f} AU")
    
    # Phase 2: Abrupt cessation (withdrawal)
    print("\nPhase 2: Abrupt cessation (withdrawal period: 14 days)...")
    
    # No nicotine input
    def I_withdrawal(t):
        return 0.0
    
    result_withdrawal = model.simulate(
        t_span=(0, 20160),  # 14 days
        input_func=I_withdrawal,
        y0=y_dependent,
        method='LSODA'
    )
    
    # Convert to days
    t_withdrawal_days = result_withdrawal['t'] / (60 * 24)
    
    # Analysis
    print("\n" + "-" * 60)
    print("WITHDRAWAL METRICS")
    print("-" * 60)
    
    # Minimum dopamine (nadir)
    D_nadir = result_withdrawal['D'].min()
    t_nadir_days = t_withdrawal_days[result_withdrawal['D'].argmin()]
    D_baseline = model.params['D_0']
    
    print(f"\nDopamine nadir: {D_nadir:.3f} AU (baseline: {D_baseline:.3f})")
    print(f"Time to nadir: {t_nadir_days:.2f} days")
    print(f"Deficit: {(D_nadir - D_baseline)/D_baseline * 100:.1f}%")
    
    # Recovery time (back to baseline)
    baseline_threshold = D_baseline * 0.95
    recovery_indices = np.where(result_withdrawal['D'] > baseline_threshold)[0]
    
    if len(recovery_indices) > 0:
        t_recovery_days = t_withdrawal_days[recovery_indices[0]]
        print(f"\nRecovery to baseline: {t_recovery_days:.2f} days")
    else:
        print(f"\nRecovery incomplete within 14 days")
    
    # Receptor downregulation
    R_T_final = result_withdrawal['R_T'][-1]
    R_T_initial = y_dependent[3]
    downregulation_pct = (R_T_initial - R_T_final) / R_T_initial * 100
    
    print(f"\nReceptor downregulation: {downregulation_pct:.1f}%")
    print(f"  Initial (dependent): {R_T_initial:.3f}")
    print(f"  Final (14 days): {R_T_final:.3f}")
    
    # Visualization
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATIONS")
    print("=" * 60)
    
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    
    # Nicotine
    axes[0].plot(t_withdrawal_days, result_withdrawal['N'], 
                color='#2E86AB', linewidth=2)
    axes[0].set_ylabel('Nicotine [μM]')
    axes[0].set_title('Nicotine Clearance', fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    # Dopamine (key withdrawal symptom)
    axes[1].plot(t_withdrawal_days, result_withdrawal['D'], 
                color='#9D4EDD', linewidth=2.5, label='Dopamine')
    axes[1].axhline(D_baseline, color='gray', linestyle='--', 
                   linewidth=2, label='Baseline', alpha=0.7)
    axes[1].axhline(D_nadir, color='red', linestyle=':', 
                   alpha=0.5, label=f'Nadir ({D_nadir:.2f})')
    axes[1].fill_between(t_withdrawal_days, D_baseline, result_withdrawal['D'],
                        where=(result_withdrawal['D'] < D_baseline),
                        alpha=0.3, color='red', label='Deficit region')
    axes[1].set_ylabel('Dopamine [AU]')
    axes[1].set_title('WITHDRAWAL: Dopamine Deficit → Negative Affect', 
                     fontweight='bold')
    axes[1].legend(loc='best')
    axes[1].grid(True, alpha=0.3)
    
    # Receptors
    axes[2].plot(t_withdrawal_days, result_withdrawal['R_a'], 
                color='#06A77D', linewidth=2, label='Active')
    axes[2].plot(t_withdrawal_days, result_withdrawal['R_d'], 
                color='#D62828', linewidth=2, linestyle='--', label='Desensitised')
    axes[2].set_ylabel('Receptor Fraction')
    axes[2].set_title('Receptor Re-sensitisation', fontweight='bold')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    # Total receptors (slow downregulation)
    axes[3].plot(t_withdrawal_days, result_withdrawal['R_T'], 
                color='#F77F00', linewidth=2.5, label='R_T')
    axes[3].axhline(model.params['R_T_0'], color='gray', linestyle='--',
                   linewidth=2, alpha=0.7, label='Pre-exposure level')
    axes[3].fill_between(t_withdrawal_days, model.params['R_T_0'], 
                        result_withdrawal['R_T'],
                        where=(result_withdrawal['R_T'] > model.params['R_T_0']),
                        alpha=0.3, color='orange', label='Upregulated region')
    axes[3].set_ylabel('Total Receptors')
    axes[3].set_xlabel('Time Since Cessation [days]')
    axes[3].set_title('SLOW: Receptor Downregulation', fontweight='bold')
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)
    
    plt.suptitle('Withdrawal: Multi-Timescale Recovery Dynamics', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('withdrawal_dynamics.png', dpi=300, bbox_inches='tight')
    print("\n✓ Saved: withdrawal_dynamics.png")
    
    # Comparison plot: Dependent vs Withdrawal
    fig2, ax = plt.subplots(figsize=(10, 6))
    
    ax.axhline(D_baseline, color='gray', linestyle='--', linewidth=2, 
              label='Baseline (never-user)', zorder=1)
    ax.axhline(y_dependent[4], color='blue', linestyle='--', linewidth=2,
              label='Dependent state', zorder=1)
    ax.plot(t_withdrawal_days, result_withdrawal['D'], 
           color='#9D4EDD', linewidth=3, label='Withdrawal trajectory')
    
    ax.set_xlabel('Days Since Cessation')
    ax.set_ylabel('Dopamine [AU]')
    ax.set_title('Withdrawal: Dopamine Falls Below Baseline', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('withdrawal_dopamine.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: withdrawal_dopamine.png")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("KEY OBSERVATIONS")
    print("=" * 60)
    print("1. Dopamine falls BELOW baseline → negative affect, anhedonia")
    print("2. Nadir occurs within first few days (acute withdrawal)")
    print("3. Receptor downregulation is SLOW (weeks to normalize)")
    print("4. Extended vulnerability period → high relapse risk")
    print("5. System is unstable: small nicotine input → large response")
    print("=" * 60)
    
    print("\n💡 Clinical insight:")
    print("The mismatch between fast symptoms (days) and slow recovery (weeks)")
    print("explains why relapse is common even after acute withdrawal ends.")


if __name__ == "__main__":
    main()
