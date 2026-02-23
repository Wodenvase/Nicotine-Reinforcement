"""
Core ODE model for nicotine feedback dynamics in the brain reward system.

This implements a multi-timescale system capturing:
- Fast: nicotine pharmacokinetics and dopamine release
- Medium: receptor desensitisation
- Slow: receptor upregulation (allostatic adaptation)

References:
    Gutkin et al. (2006) PNAS 103(4):1106-1111
    Quick & Lester (2002) J Neurobiol 53(4):457-478
    Benowitz et al. (2009) Handb Exp Pharmacol 192:29-60
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Callable, Dict, Optional, Tuple


class NicotineRewardModel:
    """
    Multi-timescale ODE model of nicotine action on the mesolimbic dopamine system.
    
    State variables:
        N(t)   - brain nicotine concentration [μM]
        R_a(t) - active α4β2 nAChRs [fraction of total]
        R_d(t) - desensitised receptors [fraction of total]
        R_T(t) - total receptor pool [dimensionless]
        D(t)   - dopamine level [arbitrary units]
    """
    
    def __init__(self, params: Optional[Dict[str, float]] = None):
        """
        Initialize model with physiological parameters.
        
        Args:
            params: Dictionary of parameter values. Uses defaults if None.
        """
        self.params = self._default_parameters()
        if params is not None:
            self.params.update(params)
    
    @staticmethod
    def _default_parameters() -> Dict[str, float]:
        """
        Physiologically-based parameter set.
        
        Timescales:
            Fast (minutes): k_N, k_clear
            Medium (hours): k_des, k_res
            Slow (days): epsilon
        
        Returns:
            Dictionary of parameter values with units in comments
        """
        return {
            # Nicotine pharmacokinetics (Benowitz et al. 2009)
            'k_N': 0.06,        # clearance rate [1/min] (t_1/2 ~ 120 min)
            
            # Receptor binding kinetics (Quick & Lester 2002)
            'k_on': 0.5,        # binding rate [1/(μM·min)]
            'k_off': 0.1,       # unbinding rate [1/min]
            'k_des': 0.05,      # desensitisation rate [1/min] 
            'k_res': 0.01,      # resensitisation rate [1/min]
            
            # Dopamine dynamics (Schultz 2002)
            'k_D': 1.0,         # dopamine release rate [AU/min]
            'k_clear': 0.5,     # dopamine clearance [1/min]
            
            # Allostatic adaptation (Koob & Le Moal 2008)
            'epsilon': 0.0001,  # upregulation rate [1/min] (very slow)
            'D_0': 0.5,         # baseline dopamine [AU]
            
            # Initial conditions
            'R_T_0': 1.0,       # initial receptor density (normalized)
        }
    
    def dynamics(self, t: float, y: np.ndarray, 
                 input_func: Callable[[float], float]) -> np.ndarray:
        """
        ODE system: dy/dt = f(t, y).
        
        Args:
            t: Time [min]
            y: State vector [N, R_a, R_d, R_T, D]
            input_func: Nicotine input function I(t) [μM/min]
        
        Returns:
            Time derivatives [dN/dt, dR_a/dt, dR_d/dt, dR_T/dt, dD/dt]
        """
        N, R_a, R_d, R_T, D = y
        p = self.params
        
        # Current nicotine input
        I = input_func(t)
        
        # Available receptor fraction (not active or desensitised)
        R_available = R_T - R_a - R_d
        R_available = max(0, R_available)  # Ensure non-negative
        
        # === Multi-timescale dynamics ===
        
        # Nicotine pharmacokinetics (FAST)
        dN_dt = I - p['k_N'] * N
        
        # Receptor activation & desensitisation (MEDIUM)
        dR_a_dt = (p['k_on'] * N * R_available 
                   - p['k_off'] * R_a 
                   - p['k_des'] * R_a)
        
        dR_d_dt = p['k_des'] * R_a - p['k_res'] * R_d
        
        # Dopamine release (FAST)
        dD_dt = p['k_D'] * R_a - p['k_clear'] * D
        
        # Receptor upregulation (ULTRA-SLOW)
        # Allostatic shift driven by chronic dopamine elevation
        dR_T_dt = p['epsilon'] * (D - p['D_0'])
        
        return np.array([dN_dt, dR_a_dt, dR_d_dt, dR_T_dt, dD_dt])
    
    def simulate(self, 
                 t_span: Tuple[float, float],
                 input_func: Callable[[float], float],
                 y0: Optional[np.ndarray] = None,
                 t_eval: Optional[np.ndarray] = None,
                 method: str = 'LSODA') -> Dict:
        """
        Simulate the model over a time interval.
        
        Args:
            t_span: Integration interval (t_start, t_end) [min]
            input_func: Nicotine input function I(t) [μM/min]
            y0: Initial state [N, R_a, R_d, R_T, D]. Uses baseline if None.
            t_eval: Time points for output. Auto-generated if None.
            method: Integration method ('LSODA', 'BDF', 'Radau')
        
        Returns:
            Dictionary with keys:
                't': time array
                'N': nicotine concentration
                'R_a': active receptors
                'R_d': desensitised receptors
                'R_T': total receptors
                'D': dopamine level
                'sol': raw scipy solution object
        """
        if y0 is None:
            y0 = self.baseline_state()
        
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 1000)
        
        # Solve ODE system with stiff solver (handles multiple timescales)
        sol = solve_ivp(
            fun=lambda t, y: self.dynamics(t, y, input_func),
            t_span=t_span,
            y0=y0,
            method=method,
            t_eval=t_eval,
            rtol=1e-8,
            atol=1e-10
        )
        
        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")
        
        return {
            't': sol.t,
            'N': sol.y[0],
            'R_a': sol.y[1],
            'R_d': sol.y[2],
            'R_T': sol.y[3],
            'D': sol.y[4],
            'sol': sol
        }
    
    def baseline_state(self) -> np.ndarray:
        """
        Baseline state (no nicotine, equilibrium).
        
        Returns:
            State vector [N=0, R_a=0, R_d=0, R_T=R_T_0, D=D_0]
        """
        return np.array([
            0.0,                    # N
            0.0,                    # R_a
            0.0,                    # R_d
            self.params['R_T_0'],   # R_T
            self.params['D_0']      # D
        ])
    
    def steady_state(self, I_const: float, 
                     include_upregulation: bool = True) -> np.ndarray:
        """
        Find steady state under constant nicotine input.
        
        This is useful for stability analysis and bifurcation studies.
        
        Args:
            I_const: Constant nicotine input rate [μM/min]
            include_upregulation: Whether to include slow R_T dynamics
        
        Returns:
            Steady state vector [N_ss, R_a_ss, R_d_ss, R_T_ss, D_ss]
        """
        from scipy.optimize import fsolve
        
        def steady_state_eqs(y):
            """Steady state: dy/dt = 0"""
            derivatives = self.dynamics(0, y, lambda t: I_const)
            if not include_upregulation:
                derivatives[3] = 0  # Fix R_T
            return derivatives
        
        # Initial guess: simulate to quasi-steady state
        y0 = self.baseline_state()
        if not include_upregulation:
            # Fast equilibrium only
            sol = self.simulate((0, 500), lambda t: I_const, y0=y0)
            y_guess = np.array([sol['N'][-1], sol['R_a'][-1], 
                               sol['R_d'][-1], sol['R_T'][-1], sol['D'][-1]])
        else:
            # Full equilibrium (very long time)
            sol = self.simulate((0, 100000), lambda t: I_const, y0=y0)
            y_guess = np.array([sol['N'][-1], sol['R_a'][-1], 
                               sol['R_d'][-1], sol['R_T'][-1], sol['D'][-1]])
        
        # Solve for exact steady state
        y_ss = fsolve(steady_state_eqs, y_guess)
        
        return y_ss


# === Input functions (nicotine delivery profiles) ===

def bolus_input(amplitude: float, t_start: float, duration: float) -> Callable:
    """
    Single bolus dose (e.g., cigarette).
    
    Args:
        amplitude: Peak input rate [μM/min]
        t_start: Onset time [min]
        duration: Delivery duration [min]
    
    Returns:
        Input function I(t)
    """
    def I(t):
        if t_start <= t < t_start + duration:
            return amplitude
        return 0.0
    return I


def continuous_input(rate: float) -> Callable:
    """
    Constant infusion (e.g., nicotine patch).
    
    Args:
        rate: Constant input rate [μM/min]
    
    Returns:
        Input function I(t)
    """
    return lambda t: rate


def repeated_bolus(amplitude: float, interval: float, 
                   n_doses: int, duration: float = 5.0) -> Callable:
    """
    Repeated dosing (e.g., regular smoking pattern).
    
    Args:
        amplitude: Per-dose input rate [μM/min]
        interval: Time between doses [min]
        n_doses: Number of doses
        duration: Delivery duration per dose [min]
    
    Returns:
        Input function I(t)
    """
    def I(t):
        for i in range(n_doses):
            t_start = i * interval
            if t_start <= t < t_start + duration:
                return amplitude
        return 0.0
    return I


def exponential_taper(initial_rate: float, half_life: float) -> Callable:
    """
    Exponentially decaying input (e.g., cessation with residual).
    
    Args:
        initial_rate: Initial input rate [μM/min]
        half_life: Decay half-life [min]
    
    Returns:
        Input function I(t)
    """
    k = np.log(2) / half_life
    return lambda t: initial_rate * np.exp(-k * t)


if __name__ == "__main__":
    """Quick validation test"""
    import matplotlib.pyplot as plt
    
    # Create model
    model = NicotineRewardModel()
    
    # Single cigarette bolus
    I = bolus_input(amplitude=5.0, t_start=10, duration=5)
    
    # Simulate
    result = model.simulate(
        t_span=(0, 300),
        input_func=I
    )
    
    # Plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    axes[0].plot(result['t'], result['N'], 'b-', linewidth=2)
    axes[0].set_ylabel('Nicotine [μM]')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(result['t'], result['R_a'], 'g-', label='Active', linewidth=2)
    axes[1].plot(result['t'], result['R_d'], 'r--', label='Desensitised', linewidth=2)
    axes[1].set_ylabel('Receptor Fraction')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(result['t'], result['D'], 'purple', linewidth=2)
    axes[2].set_ylabel('Dopamine [AU]')
    axes[2].set_xlabel('Time [min]')
    axes[2].grid(True, alpha=0.3)
    
    plt.suptitle('Single Cigarette Response', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
