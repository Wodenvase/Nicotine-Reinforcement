"""
Bifurcation analysis for the nicotine reward model.

This module implements:
- Parameter continuation methods
- Bifurcation detection (saddle-node, Hopf, etc.)
- Hysteresis analysis
- Transition to dependence regime

Critical for understanding:
- How intake rate changes system behaviour
- Sudden transitions vs gradual adaptation
- Why cessation is difficult (hysteresis)

References:
    Kuznetsov (2004) Elements of Applied Bifurcation Theory
    Koob & Le Moal (2008) Phil Trans R Soc B 363:3113-3123
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable
from scipy.optimize import fsolve, brentq
import warnings


class BifurcationAnalyzer:
    """
    Bifurcation analysis for the nicotine reward system.
    
    Studies qualitative changes in dynamics as parameters vary.
    """
    
    def __init__(self, model):
        """
        Args:
            model: Instance of NicotineRewardModel
        """
        self.model = model
    
    def continuation_1d(self, 
                       param_name: str,
                       param_range: np.ndarray,
                       include_upregulation: bool = True,
                       track_stability: bool = True) -> Dict:
        """
        One-parameter continuation: track equilibria as parameter varies.
        
        This is the core bifurcation tool. Shows how steady states
        change and potentially disappear/appear.
        
        Args:
            param_name: Parameter to vary (e.g., 'I' for input rate)
            param_range: Array of parameter values
            include_upregulation: Include slow R_T dynamics
            track_stability: Compute stability at each point (slower)
        
        Returns:
            Dictionary with:
                'param_values': Parameter values
                'steady_states': Array of steady states (n_params × 5)
                'stability': Boolean array if track_stability=True
                'eigenvalues': Max real eigenvalue if track_stability=True
        """
        n_params = len(param_range)
        steady_states = np.zeros((n_params, 5))
        is_stable = np.zeros(n_params, dtype=bool)
        max_eigs = np.zeros(n_params)
        
        # Special handling for input rate (not a model parameter)
        if param_name == 'I':
            varying_input = True
            original_param = None
        else:
            varying_input = False
            original_param = self.model.params.get(param_name)
        
        for i, p_val in enumerate(param_range):
            try:
                if varying_input:
                    # Use constant input value
                    I_const = p_val
                else:
                    # Set model parameter
                    self.model.params[param_name] = p_val
                    I_const = 0.0  # Or could be another parameter
                
                # Find steady state
                y_ss = self.model.steady_state(I_const, include_upregulation)
                steady_states[i] = y_ss
                
                # Stability analysis
                if track_stability:
                    from stability import StabilityAnalyzer
                    analyzer = StabilityAnalyzer(self.model)
                    eigen_results = analyzer.eigenanalysis(y_ss, I_const)
                    is_stable[i] = eigen_results['is_stable']
                    max_eigs[i] = np.max(eigen_results['real_parts'])
            
            except Exception as e:
                # Continuation failed (steady state doesn't exist?)
                steady_states[i] = np.nan
                is_stable[i] = False
                max_eigs[i] = np.nan
        
        # Restore original parameter
        if not varying_input and original_param is not None:
            self.model.params[param_name] = original_param
        
        result = {
            'param_values': param_range,
            'steady_states': steady_states,
        }
        
        if track_stability:
            result['stability'] = is_stable
            result['max_eigenvalue'] = max_eigs
        
        return result
    
    def detect_bifurcations(self, continuation_result: Dict) -> Dict:
        """
        Detect bifurcation points from continuation data.
        
        Looks for:
        - Saddle-node: fold in equilibrium curve
        - Stability changes: eigenvalue crossing zero
        
        Args:
            continuation_result: Output from continuation_1d()
        
        Returns:
            Dictionary with:
                'saddle_node_indices': Where folds occur
                'stability_change_indices': Where stability changes
        """
        param_vals = continuation_result['param_values']
        steady_states = continuation_result['steady_states']
        
        saddle_node_indices = []
        stability_change_indices = []
        
        # Detect folds (saddle-node) by checking for turning points
        # Look at any state variable (e.g., dopamine D)
        D_curve = steady_states[:, 4]
        
        for i in range(1, len(D_curve) - 1):
            if np.isnan(D_curve[i]):
                continue
            
            # Check for local extremum in D vs parameter
            if ((D_curve[i] > D_curve[i-1] and D_curve[i] > D_curve[i+1]) or
                (D_curve[i] < D_curve[i-1] and D_curve[i] < D_curve[i+1])):
                saddle_node_indices.append(i)
        
        # Detect stability changes
        if 'stability' in continuation_result:
            stability = continuation_result['stability']
            for i in range(1, len(stability)):
                if stability[i] != stability[i-1]:
                    stability_change_indices.append(i)
        
        return {
            'saddle_node_indices': saddle_node_indices,
            'stability_change_indices': stability_change_indices
        }
    
    def hysteresis_loop(self,
                       param_range_up: np.ndarray,
                       param_range_down: np.ndarray,
                       param_name: str = 'I') -> Dict:
        """
        Compute hysteresis loop: continuation up and down.
        
        Shows path-dependence: increasing vs decreasing parameter
        gives different states. Critical for understanding why
        cessation is difficult.
        
        Args:
            param_range_up: Parameter values for upward sweep
            param_range_down: Parameter values for downward sweep
            param_name: Parameter to vary
        
        Returns:
            Dictionary with:
                'up_sweep': Continuation result for increasing parameter
                'down_sweep': Continuation result for decreasing parameter
                'hysteresis_width': Width of hysteresis region
        """
        # Upward sweep
        up_sweep = self.continuation_1d(
            param_name, 
            param_range_up,
            include_upregulation=True,
            track_stability=True
        )
        
        # Downward sweep
        down_sweep = self.continuation_1d(
            param_name,
            param_range_down,
            include_upregulation=True,
            track_stability=True
        )
        
        # Estimate hysteresis width (crude approximation)
        # Find where stable branches differ
        D_up = up_sweep['steady_states'][:, 4]
        D_down = down_sweep['steady_states'][:, 4]
        
        # Find overlap region
        p_min = max(param_range_up[0], param_range_down[-1])
        p_max = min(param_range_up[-1], param_range_down[0])
        
        if p_max > p_min:
            # Interpolate and compare
            overlap_indices_up = (param_range_up >= p_min) & (param_range_up <= p_max)
            overlap_indices_down = (param_range_down >= p_min) & (param_range_down <= p_max)
            
            if np.any(overlap_indices_up) and np.any(overlap_indices_down):
                max_diff = np.nanmax(np.abs(
                    D_up[overlap_indices_up] - D_down[overlap_indices_down]
                ))
                hysteresis_width = max_diff
            else:
                hysteresis_width = 0.0
        else:
            hysteresis_width = 0.0
        
        return {
            'up_sweep': up_sweep,
            'down_sweep': down_sweep,
            'hysteresis_width': hysteresis_width
        }
    
    def critical_transition_threshold(self,
                                     param_name: str = 'I',
                                     param_range: Tuple[float, float] = (0.0, 2.0),
                                     threshold_var: str = 'D',
                                     threshold_value: float = 1.0) -> Optional[float]:
        """
        Find critical parameter value where system crosses threshold.
        
        Example: At what input rate does dopamine exceed a "dependence" level?
        
        Args:
            param_name: Parameter to vary
            param_range: (min, max) search range
            threshold_var: Which variable to threshold ('D', 'R_T', etc.)
            threshold_value: Threshold value
        
        Returns:
            Critical parameter value, or None if not found
        """
        var_indices = {'N': 0, 'R_a': 1, 'R_d': 2, 'R_T': 3, 'D': 4}
        var_idx = var_indices[threshold_var]
        
        def threshold_func(p_val):
            """Returns difference from threshold"""
            if param_name == 'I':
                I_const = p_val
            else:
                self.model.params[param_name] = p_val
                I_const = 0.0
            
            try:
                y_ss = self.model.steady_state(I_const, include_upregulation=True)
                return y_ss[var_idx] - threshold_value
            except:
                return np.nan
        
        # Check if threshold is crossed in range
        val_min = threshold_func(param_range[0])
        val_max = threshold_func(param_range[1])
        
        if np.isnan(val_min) or np.isnan(val_max):
            return None
        
        if val_min * val_max > 0:
            # Same sign - no crossing
            return None
        
        # Find crossing using bisection
        try:
            critical_val = brentq(threshold_func, param_range[0], param_range[1])
            return critical_val
        except:
            return None
    
    def intake_regime_classifier(self, I_values: np.ndarray) -> Dict:
        """
        Classify intake regimes based on steady-state response.
        
        Categories:
        - "Casual": Low dopamine response, no upregulation
        - "Regular": Moderate dopamine, beginning upregulation
        - "Dependent": High dopamine, significant upregulation
        
        Args:
            I_values: Array of intake rates to classify
        
        Returns:
            Dictionary with:
                'intake_rates': Input array
                'dopamine_levels': Steady-state dopamine
                'receptor_upregulation': R_T - R_T_0
                'regime': Classification string array
        """
        n = len(I_values)
        D_steady = np.zeros(n)
        R_T_steady = np.zeros(n)
        regimes = []
        
        R_T_0 = self.model.params['R_T_0']
        
        for i, I in enumerate(I_values):
            try:
                y_ss = self.model.steady_state(I, include_upregulation=True)
                D_steady[i] = y_ss[4]
                R_T_steady[i] = y_ss[3]
                
                # Classification heuristic
                D_increase = (D_steady[i] - self.model.params['D_0']) / self.model.params['D_0']
                R_T_increase = (R_T_steady[i] - R_T_0) / R_T_0
                
                if D_increase < 0.2 and R_T_increase < 0.05:
                    regime = "Casual"
                elif D_increase < 0.5 and R_T_increase < 0.15:
                    regime = "Regular"
                else:
                    regime = "Dependent"
                
                regimes.append(regime)
            
            except:
                D_steady[i] = np.nan
                R_T_steady[i] = np.nan
                regimes.append("Undefined")
        
        return {
            'intake_rates': I_values,
            'dopamine_levels': D_steady,
            'receptor_upregulation': R_T_steady - R_T_0,
            'regime': np.array(regimes)
        }


if __name__ == "__main__":
    """Demonstration of bifurcation analysis"""
    from model import NicotineRewardModel
    import matplotlib.pyplot as plt
    
    print("=" * 60)
    print("BIFURCATION ANALYSIS: INTAKE RATE")
    print("=" * 60)
    
    # Create model
    model = NicotineRewardModel()
    analyzer = BifurcationAnalyzer(model)
    
    # Vary intake rate from low to high
    I_range = np.linspace(0.0, 1.5, 50)
    
    print("\nRunning continuation analysis...")
    result = analyzer.continuation_1d(
        param_name='I',
        param_range=I_range,
        include_upregulation=True,
        track_stability=True
    )
    
    # Detect bifurcations
    bif_points = analyzer.detect_bifurcations(result)
    print(f"Stability changes detected at indices: {bif_points['stability_change_indices']}")
    
    # Classify regimes
    regime_analysis = analyzer.intake_regime_classifier(I_range)
    
    # Plot bifurcation diagram
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Dopamine vs intake rate
    D_curve = result['steady_states'][:, 4]
    stable = result['stability']
    
    axes[0].plot(I_range[stable], D_curve[stable], 'b-', linewidth=2, label='Stable')
    axes[0].plot(I_range[~stable], D_curve[~stable], 'r--', linewidth=2, label='Unstable')
    axes[0].axhline(model.params['D_0'], color='gray', linestyle=':', label='Baseline')
    axes[0].set_ylabel('Dopamine [AU]')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Bifurcation Diagram: Dopamine vs Intake Rate')
    
    # Receptor upregulation
    R_T_curve = result['steady_states'][:, 3]
    axes[1].plot(I_range, R_T_curve - model.params['R_T_0'], 'purple', linewidth=2)
    axes[1].axhline(0, color='gray', linestyle=':')
    axes[1].set_ylabel('Receptor Upregulation\n(ΔR_T)')
    axes[1].set_xlabel('Intake Rate [μM/min]')
    axes[1].grid(True, alpha=0.3)
    
    # Mark regimes
    for regime_name in ['Casual', 'Regular', 'Dependent']:
        mask = regime_analysis['regime'] == regime_name
        if np.any(mask):
            i_start = np.where(mask)[0][0]
            axes[1].axvline(I_range[i_start], color='orange', alpha=0.3, linestyle='--')
            axes[1].text(I_range[i_start], axes[1].get_ylim()[1] * 0.9, 
                        regime_name, rotation=90, va='top', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('bifurcation_diagram.png', dpi=150, bbox_inches='tight')
    print("\nBifurcation diagram saved to bifurcation_diagram.png")
    plt.show()
