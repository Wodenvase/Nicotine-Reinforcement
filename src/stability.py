"""
Stability analysis tools for the nicotine reward model.

This module provides:
- Linearisation and Jacobian computation
- Eigenvalue analysis for timescale separation
- Stability classification of fixed points
- Phase portrait generation

References:
    Strogatz (2015) Nonlinear Dynamics and Chaos
    Izhikevich (2007) Dynamical Systems in Neuroscience
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.linalg import eig
import warnings


class StabilityAnalyzer:
    """
    Stability analysis for the nicotine reward dynamical system.
    """
    
    def __init__(self, model):
        """
        Args:
            model: Instance of NicotineRewardModel
        """
        self.model = model
    
    def jacobian(self, y: np.ndarray, I_const: float) -> np.ndarray:
        """
        Compute Jacobian matrix at a given state.
        
        The Jacobian J_ij = ∂f_i/∂y_j represents the linearised dynamics
        near an equilibrium or trajectory point.
        
        Args:
            y: State vector [N, R_a, R_d, R_T, D]
            I_const: Constant input rate
        
        Returns:
            5×5 Jacobian matrix
        """
        N, R_a, R_d, R_T, D = y
        p = self.model.params
        
        R_available = R_T - R_a - R_d
        
        # Initialize Jacobian
        J = np.zeros((5, 5))
        
        # Row 0: dN/dt derivatives
        J[0, 0] = -p['k_N']  # ∂(dN/dt)/∂N
        
        # Row 1: dR_a/dt derivatives
        J[1, 0] = p['k_on'] * R_available  # ∂(dR_a/dt)/∂N
        J[1, 1] = -p['k_on'] * N - p['k_off'] - p['k_des']  # ∂(dR_a/dt)/∂R_a
        J[1, 2] = -p['k_on'] * N  # ∂(dR_a/dt)/∂R_d
        J[1, 3] = p['k_on'] * N   # ∂(dR_a/dt)/∂R_T
        
        # Row 2: dR_d/dt derivatives
        J[2, 1] = p['k_des']   # ∂(dR_d/dt)/∂R_a
        J[2, 2] = -p['k_res']  # ∂(dR_d/dt)/∂R_d
        
        # Row 3: dR_T/dt derivatives
        J[3, 4] = p['epsilon']  # ∂(dR_T/dt)/∂D
        
        # Row 4: dD/dt derivatives
        J[4, 1] = p['k_D']      # ∂(dD/dt)/∂R_a
        J[4, 4] = -p['k_clear']  # ∂(dD/dt)/∂D
        
        return J
    
    def eigenanalysis(self, y: np.ndarray, I_const: float) -> Dict:
        """
        Perform eigenvalue/eigenvector analysis at a fixed point.
        
        This reveals:
        - Stability (all Re(λ) < 0 → stable)
        - Timescales (τ = 1/|Re(λ)|)
        - Direction of fastest/slowest modes
        
        Args:
            y: State vector (typically a fixed point)
            I_const: Constant input rate
        
        Returns:
            Dictionary containing:
                'eigenvalues': Complex eigenvalues λ
                'eigenvectors': Corresponding eigenvectors
                'timescales': Characteristic timescales [min]
                'is_stable': Boolean stability flag
                'real_parts': Real parts of eigenvalues
                'imag_parts': Imaginary parts (for oscillations)
        """
        J = self.jacobian(y, I_const)
        
        # Compute eigenvalues and eigenvectors
        eigenvalues, eigenvectors = eig(J)
        
        # Extract real and imaginary parts
        real_parts = eigenvalues.real
        imag_parts = eigenvalues.imag
        
        # Stability: all real parts < 0
        is_stable = np.all(real_parts < 0)
        
        # Timescales: τ = 1/|Re(λ)|
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            timescales = 1.0 / np.abs(real_parts)
            timescales[np.isinf(timescales)] = np.nan
        
        return {
            'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors,
            'timescales': timescales,
            'is_stable': is_stable,
            'real_parts': real_parts,
            'imag_parts': imag_parts
        }
    
    def timescale_separation(self, y: np.ndarray, I_const: float) -> Dict:
        """
        Quantify timescale separation in the system.
        
        This is the mathematical foundation for the "multiple timescales"
        claim. Large ratios indicate clear separation.
        
        Args:
            y: State vector
            I_const: Constant input rate
        
        Returns:
            Dictionary with:
                'fast_timescale': Fastest mode [min]
                'slow_timescale': Slowest mode [min]
                'separation_ratio': τ_slow / τ_fast
                'all_timescales': Sorted array of timescales
        """
        eigen_results = self.eigenanalysis(y, I_const)
        timescales = eigen_results['timescales']
        
        # Sort and filter valid timescales
        valid_timescales = timescales[~np.isnan(timescales)]
        valid_timescales = np.sort(valid_timescales)
        
        if len(valid_timescales) == 0:
            return {
                'fast_timescale': np.nan,
                'slow_timescale': np.nan,
                'separation_ratio': np.nan,
                'all_timescales': np.array([])
            }
        
        fast = valid_timescales[0]
        slow = valid_timescales[-1]
        
        return {
            'fast_timescale': fast,
            'slow_timescale': slow,
            'separation_ratio': slow / fast if fast > 0 else np.inf,
            'all_timescales': valid_timescales
        }
    
    def classify_fixed_point(self, y: np.ndarray, I_const: float) -> str:
        """
        Classify fixed point topology (2D projection).
        
        For 2D subsystems:
        - Stable node: both λ < 0, real
        - Unstable node: both λ > 0, real
        - Saddle: λ have opposite signs
        - Stable spiral: Re(λ) < 0, Im(λ) ≠ 0
        - Unstable spiral: Re(λ) > 0, Im(λ) ≠ 0
        
        Args:
            y: State vector
            I_const: Constant input rate
        
        Returns:
            Classification string
        """
        eigen_results = self.eigenanalysis(y, I_const)
        eigenvalues = eigen_results['eigenvalues']
        real_parts = eigen_results['real_parts']
        imag_parts = eigen_results['imag_parts']
        
        # Full 5D system - give general stability
        if np.all(real_parts < 0):
            if np.any(imag_parts != 0):
                return "Stable focus (with oscillatory modes)"
            else:
                return "Stable node"
        elif np.all(real_parts > 0):
            return "Unstable node"
        else:
            return "Saddle (mixed stability)"
    
    def parameter_sensitivity(self, param_name: str, 
                             param_range: np.ndarray,
                             I_const: float,
                             include_upregulation: bool = True) -> Dict:
        """
        Compute stability as function of a parameter.
        
        Useful for understanding robustness and parameter dependencies.
        
        Args:
            param_name: Name of parameter to vary
            param_range: Array of parameter values to test
            I_const: Constant input rate
            include_upregulation: Include slow R_T dynamics
        
        Returns:
            Dictionary with arrays:
                'param_values': Parameter values tested
                'max_real_eigenvalue': Largest Re(λ) (determines stability)
                'is_stable': Boolean array
        """
        original_value = self.model.params[param_name]
        
        max_real_eigs = []
        is_stable = []
        
        for p_val in param_range:
            # Set parameter
            self.model.params[param_name] = p_val
            
            # Find steady state
            try:
                y_ss = self.model.steady_state(I_const, include_upregulation)
                eigen_results = self.eigenanalysis(y_ss, I_const)
                
                max_real_eigs.append(np.max(eigen_results['real_parts']))
                is_stable.append(eigen_results['is_stable'])
            except:
                max_real_eigs.append(np.nan)
                is_stable.append(False)
        
        # Restore original parameter
        self.model.params[param_name] = original_value
        
        return {
            'param_values': param_range,
            'max_real_eigenvalue': np.array(max_real_eigs),
            'is_stable': np.array(is_stable)
        }


def phase_portrait_2d(model, var1_idx: int, var2_idx: int,
                      var1_range: Tuple[float, float],
                      var2_range: Tuple[float, float],
                      I_const: float,
                      n_points: int = 20) -> Dict:
    """
    Generate 2D phase portrait (vector field).
    
    Projects the 5D system onto a 2D plane for visualization.
    
    Args:
        model: NicotineRewardModel instance
        var1_idx: Index of first variable (0=N, 1=R_a, 2=R_d, 3=R_T, 4=D)
        var2_idx: Index of second variable
        var1_range: (min, max) for first variable
        var2_range: (min, max) for second variable
        I_const: Constant input rate
        n_points: Grid resolution
    
    Returns:
        Dictionary with:
            'var1_grid': Meshgrid for first variable
            'var2_grid': Meshgrid for second variable
            'dvar1': Velocity field for first variable
            'dvar2': Velocity field for second variable
    """
    # Create meshgrid
    var1 = np.linspace(var1_range[0], var1_range[1], n_points)
    var2 = np.linspace(var2_range[0], var2_range[1], n_points)
    V1, V2 = np.meshgrid(var1, var2)
    
    # Compute derivatives at each point
    dV1 = np.zeros_like(V1)
    dV2 = np.zeros_like(V2)
    
    # Get a reference steady state for other variables
    y_ref = model.steady_state(I_const, include_upregulation=False)
    
    for i in range(n_points):
        for j in range(n_points):
            # Set state with current grid values
            y = y_ref.copy()
            y[var1_idx] = V1[i, j]
            y[var2_idx] = V2[i, j]
            
            # Compute derivatives
            dydt = model.dynamics(0, y, lambda t: I_const)
            dV1[i, j] = dydt[var1_idx]
            dV2[i, j] = dydt[var2_idx]
    
    return {
        'var1_grid': V1,
        'var2_grid': V2,
        'dvar1': dV1,
        'dvar2': dV2
    }


if __name__ == "__main__":
    """Demonstration of stability analysis"""
    from model import NicotineRewardModel
    
    # Create model
    model = NicotineRewardModel()
    
    # Find steady state under constant input
    I_const = 0.5  # μM/min
    y_ss = model.steady_state(I_const, include_upregulation=False)
    
    print("=" * 60)
    print("STABILITY ANALYSIS AT STEADY STATE")
    print("=" * 60)
    print(f"\nConstant input: {I_const} μM/min")
    print(f"\nSteady state:")
    print(f"  N   = {y_ss[0]:.4f} μM")
    print(f"  R_a = {y_ss[1]:.4f}")
    print(f"  R_d = {y_ss[2]:.4f}")
    print(f"  R_T = {y_ss[3]:.4f}")
    print(f"  D   = {y_ss[4]:.4f} AU")
    
    # Stability analysis
    analyzer = StabilityAnalyzer(model)
    eigen_results = analyzer.eigenanalysis(y_ss, I_const)
    
    print(f"\nStability: {'STABLE' if eigen_results['is_stable'] else 'UNSTABLE'}")
    print(f"\nEigenvalues (real parts):")
    for i, (λ, τ) in enumerate(zip(eigen_results['real_parts'], 
                                     eigen_results['timescales'])):
        print(f"  λ_{i+1} = {λ:.6f}  →  τ = {τ:.2f} min")
    
    # Timescale separation
    ts_results = analyzer.timescale_separation(y_ss, I_const)
    print(f"\nTimescale separation:")
    print(f"  Fastest mode: {ts_results['fast_timescale']:.2f} min")
    print(f"  Slowest mode: {ts_results['slow_timescale']:.2f} min")
    print(f"  Ratio: {ts_results['separation_ratio']:.1f}×")
    
    # Classification
    classification = analyzer.classify_fixed_point(y_ss, I_const)
    print(f"\nFixed point type: {classification}")
