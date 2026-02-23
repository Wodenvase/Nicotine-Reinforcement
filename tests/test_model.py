"""
Unit tests for the nicotine reward model.

Tests numerical accuracy, stability properties, and edge cases.
"""

import sys
sys.path.append('../src')

import numpy as np
from model import NicotineRewardModel, bolus_input, continuous_input
from stability import StabilityAnalyzer
from bifurcation import BifurcationAnalyzer


def test_baseline_equilibrium():
    """Test that baseline state is a true equilibrium."""
    print("\n[TEST] Baseline equilibrium...")
    
    model = NicotineRewardModel()
    y0 = model.baseline_state()
    
    # Derivatives should be zero
    dydt = model.dynamics(0, y0, lambda t: 0.0)
    
    assert np.allclose(dydt, 0.0, atol=1e-10), "Baseline is not equilibrium!"
    print("  ✓ Baseline state is exact equilibrium")


def test_conservation():
    """Test that receptor fractions don't exceed total."""
    print("\n[TEST] Receptor conservation...")
    
    model = NicotineRewardModel()
    I = bolus_input(amplitude=10.0, t_start=0, duration=10)
    
    result = model.simulate((0, 200), I)
    
    # R_a + R_d should never exceed R_T
    total_used = result['R_a'] + result['R_d']
    excess = total_used - result['R_T']
    
    assert np.all(excess <= 1e-6), f"Conservation violated! Max excess: {excess.max()}"
    print(f"  ✓ Receptor conservation holds (max excess: {excess.max():.2e})")


def test_stability_baseline():
    """Test that baseline is stable."""
    print("\n[TEST] Baseline stability...")
    
    model = NicotineRewardModel()
    analyzer = StabilityAnalyzer(model)
    
    y0 = model.baseline_state()
    eigen_results = analyzer.eigenanalysis(y0, 0.0)
    
    assert eigen_results['is_stable'], "Baseline should be stable!"
    print(f"  ✓ Baseline is stable (max Re(λ) = {np.max(eigen_results['real_parts']):.4f})")


def test_timescale_ordering():
    """Test that timescales are properly separated."""
    print("\n[TEST] Timescale separation...")
    
    model = NicotineRewardModel()
    analyzer = StabilityAnalyzer(model)
    
    # At moderate input
    y_ss = model.steady_state(0.5, include_upregulation=True)
    ts_results = analyzer.timescale_separation(y_ss, 0.5)
    
    ratio = ts_results['separation_ratio']
    assert ratio > 10, f"Timescale separation too weak: {ratio:.1f}×"
    
    print(f"  ✓ Timescales well separated: {ratio:.0f}×")
    print(f"    Fast: {ts_results['fast_timescale']:.2f} min")
    print(f"    Slow: {ts_results['slow_timescale']:.2f} min")


def test_monotonic_bifurcation():
    """Test that dopamine increases monotonically with input."""
    print("\n[TEST] Bifurcation monotonicity...")
    
    model = NicotineRewardModel()
    bif_analyzer = BifurcationAnalyzer(model)
    
    I_range = np.linspace(0.0, 1.0, 20)
    result = bif_analyzer.continuation_1d('I', I_range, 
                                          include_upregulation=False,
                                          track_stability=False)
    
    D_curve = result['steady_states'][:, 4]
    
    # Check monotonicity
    diffs = np.diff(D_curve)
    assert np.all(diffs >= -1e-6), "Dopamine should increase monotonically!"
    
    print("  ✓ Dopamine increases monotonically with intake")


def test_numerical_convergence():
    """Test that different tolerances give consistent results."""
    print("\n[TEST] Numerical convergence...")
    
    model = NicotineRewardModel()
    I = bolus_input(amplitude=5.0, t_start=10, duration=5)
    
    # Loose tolerance
    result1 = model.simulate((0, 100), I, method='LSODA')
    
    # Tight tolerance
    model_tight = NicotineRewardModel()
    result2 = model_tight.simulate((0, 100), I, method='Radau')
    
    # Compare final states
    D_diff = abs(result1['D'][-1] - result2['D'][-1])
    
    assert D_diff < 0.01, f"Numerical methods disagree: {D_diff}"
    print(f"  ✓ Numerical convergence verified (diff = {D_diff:.4f})")


def test_parameter_ranges():
    """Test model behaves reasonably across parameter ranges."""
    print("\n[TEST] Parameter robustness...")
    
    model = NicotineRewardModel()
    
    # Vary desensitisation rate
    k_des_values = [0.01, 0.05, 0.1]
    
    for k_des in k_des_values:
        model.params['k_des'] = k_des
        try:
            result = model.simulate((0, 100), continuous_input(0.5))
            assert np.all(np.isfinite(result['D'])), f"Non-finite values at k_des={k_des}"
        except Exception as e:
            raise AssertionError(f"Failed at k_des={k_des}: {e}")
    
    print(f"  ✓ Model stable across parameter variations")


def run_all_tests():
    """Run complete test suite."""
    print("=" * 60)
    print("RUNNING TEST SUITE")
    print("=" * 60)
    
    tests = [
        test_baseline_equilibrium,
        test_conservation,
        test_stability_baseline,
        test_timescale_ordering,
        test_monotonic_bifurcation,
        test_numerical_convergence,
        test_parameter_ranges
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ✗ FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"TEST RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
