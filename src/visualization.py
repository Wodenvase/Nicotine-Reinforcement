"""
Visualization utilities for the nicotine reward model.

Provides publication-quality plots for:
- Time series
- Phase portraits
- Bifurcation diagrams
- Parameter sensitivity
- Multi-panel comparative figures

Uses consistent styling for scientific presentation.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Dict, List, Optional, Tuple
import warnings


# Publication style settings
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


class Visualizer:
    """
    Visualization suite for nicotine reward dynamics.
    """
    
    STATE_LABELS = {
        'N': 'Nicotine [μM]',
        'R_a': 'Active Receptors',
        'R_d': 'Desensitised Receptors',
        'R_T': 'Total Receptors',
        'D': 'Dopamine [AU]'
    }
    
    STATE_COLORS = {
        'N': '#2E86AB',    # Blue
        'R_a': '#06A77D',  # Green
        'R_d': '#D62828',  # Red
        'R_T': '#F77F00',  # Orange
        'D': '#9D4EDD'     # Purple
    }
    
    @staticmethod
    def plot_timeseries(result: Dict, 
                       variables: List[str] = ['N', 'R_a', 'R_d', 'D'],
                       title: Optional[str] = None,
                       save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot time series of state variables.
        
        Args:
            result: Simulation result dictionary from model.simulate()
            variables: List of variables to plot
            title: Optional figure title
            save_path: Path to save figure
        
        Returns:
            Figure object
        """
        n_vars = len(variables)
        fig, axes = plt.subplots(n_vars, 1, figsize=(10, 2.5*n_vars), 
                                sharex=True)
        
        if n_vars == 1:
            axes = [axes]
        
        for ax, var in zip(axes, variables):
            ax.plot(result['t'], result[var], 
                   color=Visualizer.STATE_COLORS[var],
                   linewidth=2, label=var)
            ax.set_ylabel(Visualizer.STATE_LABELS[var])
            ax.grid(True, alpha=0.3)
            ax.set_xlim([result['t'][0], result['t'][-1]])
        
        axes[-1].set_xlabel('Time [min]')
        
        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_phase_portrait_2d(result: Dict,
                              var1: str, var2: str,
                              nullclines: bool = False,
                              title: Optional[str] = None,
                              save_path: Optional[str] = None) -> plt.Figure:
        """
        2D phase portrait showing trajectory.
        
        Args:
            result: Simulation result
            var1: First variable (x-axis)
            var2: Second variable (y-axis)
            nullclines: Whether to show nullclines
            title: Optional title
            save_path: Path to save
        
        Returns:
            Figure object
        """
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Trajectory
        x = result[var1]
        y = result[var2]
        
        # Color by time
        t = result['t']
        t_norm = (t - t[0]) / (t[-1] - t[0])
        
        for i in range(len(x) - 1):
            ax.plot(x[i:i+2], y[i:i+2], 
                   color=cm.viridis(t_norm[i]),
                   linewidth=2, alpha=0.7)
        
        # Start and end markers
        ax.plot(x[0], y[0], 'go', markersize=10, label='Start', zorder=5)
        ax.plot(x[-1], y[-1], 'rs', markersize=10, label='End', zorder=5)
        
        ax.set_xlabel(Visualizer.STATE_LABELS[var1])
        ax.set_ylabel(Visualizer.STATE_LABELS[var2])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if title:
            ax.set_title(title, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_bifurcation_diagram(continuation_result: Dict,
                                state_var: str = 'D',
                                param_label: str = 'Input Rate [μM/min]',
                                title: Optional[str] = None,
                                save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot bifurcation diagram.
        
        Args:
            continuation_result: Output from BifurcationAnalyzer.continuation_1d()
            state_var: State variable to plot (vertical axis)
            param_label: Label for parameter axis
            title: Optional title
            save_path: Save path
        
        Returns:
            Figure object
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        var_idx = {'N': 0, 'R_a': 1, 'R_d': 2, 'R_T': 3, 'D': 4}[state_var]
        
        param_vals = continuation_result['param_values']
        steady_states = continuation_result['steady_states'][:, var_idx]
        
        if 'stability' in continuation_result:
            stable = continuation_result['stability']
            
            # Plot stable branch
            ax.plot(param_vals[stable], steady_states[stable],
                   'b-', linewidth=2.5, label='Stable')
            
            # Plot unstable branch
            if np.any(~stable):
                ax.plot(param_vals[~stable], steady_states[~stable],
                       'r--', linewidth=2.5, label='Unstable')
        else:
            ax.plot(param_vals, steady_states, 'b-', linewidth=2.5)
        
        ax.set_xlabel(param_label)
        ax.set_ylabel(Visualizer.STATE_LABELS[state_var])
        ax.grid(True, alpha=0.3)
        
        if 'stability' in continuation_result:
            ax.legend()
        
        if title:
            ax.set_title(title, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_comparison(results: Dict[str, Dict],
                       variable: str = 'D',
                       title: Optional[str] = None,
                       save_path: Optional[str] = None) -> plt.Figure:
        """
        Compare multiple simulation results.
        
        Args:
            results: Dictionary {label: result} of simulation results
            variable: Variable to compare
            title: Optional title
            save_path: Save path
        
        Returns:
            Figure object
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
        
        for (label, result), color in zip(results.items(), colors):
            ax.plot(result['t'], result[variable],
                   linewidth=2, label=label, color=color)
        
        ax.set_xlabel('Time [min]')
        ax.set_ylabel(Visualizer.STATE_LABELS[variable])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if title:
            ax.set_title(title, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_parameter_sensitivity(param_name: str,
                                   param_values: np.ndarray,
                                   metric_values: np.ndarray,
                                   metric_label: str,
                                   critical_threshold: Optional[float] = None,
                                   title: Optional[str] = None,
                                   save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot parameter sensitivity analysis.
        
        Args:
            param_name: Parameter name
            param_values: Array of parameter values
            metric_values: Corresponding metric (e.g., max dopamine)
            metric_label: Label for metric axis
            critical_threshold: Optional threshold line
            title: Optional title
            save_path: Save path
        
        Returns:
            Figure object
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(param_values, metric_values, 'b-', linewidth=2.5)
        
        if critical_threshold is not None:
            ax.axhline(critical_threshold, color='r', linestyle='--',
                      linewidth=2, label=f'Threshold: {critical_threshold}')
            ax.legend()
        
        ax.set_xlabel(param_name)
        ax.set_ylabel(metric_label)
        ax.grid(True, alpha=0.3)
        
        if title:
            ax.set_title(title, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_multi_timescale(result: Dict,
                            title: Optional[str] = None,
                            save_path: Optional[str] = None) -> plt.Figure:
        """
        Specialized plot showing all three timescales.
        
        Args:
            result: Long simulation result showing all timescales
            title: Optional title
            save_path: Save path
        
        Returns:
            Figure object
        """
        fig = plt.figure(figsize=(14, 8))
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
        
        t = result['t']
        
        # Fast dynamics (top row, left)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(t, result['N'], color=Visualizer.STATE_COLORS['N'], linewidth=2)
        ax1.set_ylabel('Nicotine [μM]')
        ax1.set_title('FAST: Nicotine Kinetics', fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Fast dynamics (top row, right)
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(t, result['D'], color=Visualizer.STATE_COLORS['D'], linewidth=2)
        ax2.set_ylabel('Dopamine [AU]')
        ax2.set_title('FAST: Dopamine Response', fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Medium dynamics (middle row, spanning)
        ax3 = fig.add_subplot(gs[1, :])
        ax3.plot(t, result['R_a'], color=Visualizer.STATE_COLORS['R_a'], 
                linewidth=2, label='Active')
        ax3.plot(t, result['R_d'], color=Visualizer.STATE_COLORS['R_d'], 
                linewidth=2, linestyle='--', label='Desensitised')
        ax3.set_ylabel('Receptor Fraction')
        ax3.set_title('MEDIUM: Receptor Desensitisation', fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Slow dynamics (bottom row, spanning)
        ax4 = fig.add_subplot(gs[2, :])
        ax4.plot(t, result['R_T'], color=Visualizer.STATE_COLORS['R_T'], linewidth=2)
        ax4.set_ylabel('Total Receptors')
        ax4.set_xlabel('Time [min]')
        ax4.set_title('SLOW: Receptor Upregulation (Allostasis)', fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold', y=0.998)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig


def plot_vector_field_2d(model, var1_idx: int, var2_idx: int,
                         var1_range: Tuple[float, float],
                         var2_range: Tuple[float, float],
                         I_const: float,
                         n_points: int = 20,
                         trajectory: Optional[Dict] = None,
                         title: Optional[str] = None,
                         save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot 2D vector field with optional trajectory overlay.
    
    Args:
        model: NicotineRewardModel instance
        var1_idx: Index of first variable
        var2_idx: Index of second variable
        var1_range: Range for first variable
        var2_range: Range for second variable
        I_const: Constant input
        n_points: Grid resolution
        trajectory: Optional simulation result to overlay
        title: Optional title
        save_path: Save path
    
    Returns:
        Figure object
    """
    from stability import phase_portrait_2d
    
    # Compute vector field
    field = phase_portrait_2d(model, var1_idx, var2_idx,
                             var1_range, var2_range, I_const, n_points)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot streamlines
    speed = np.sqrt(field['dvar1']**2 + field['dvar2']**2)
    ax.streamplot(field['var1_grid'], field['var2_grid'],
                 field['dvar1'], field['dvar2'],
                 color=speed, cmap='viridis', density=1.5,
                 linewidth=1, arrowsize=1.5)
    
    # Overlay trajectory if provided
    if trajectory is not None:
        var_names = ['N', 'R_a', 'R_d', 'R_T', 'D']
        var1_name = var_names[var1_idx]
        var2_name = var_names[var2_idx]
        
        ax.plot(trajectory[var1_name], trajectory[var2_name],
               'r-', linewidth=2.5, alpha=0.8, label='Trajectory')
        ax.plot(trajectory[var1_name][0], trajectory[var2_name][0],
               'go', markersize=12, label='Start')
        ax.plot(trajectory[var1_name][-1], trajectory[var2_name][-1],
               'rs', markersize=12, label='End')
        ax.legend()
    
    ax.set_xlabel(Visualizer.STATE_LABELS[var_names[var1_idx]])
    ax.set_ylabel(Visualizer.STATE_LABELS[var_names[var2_idx]])
    
    if title:
        ax.set_title(title, fontweight='bold')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


if __name__ == "__main__":
    """Demonstration of visualization tools"""
    from model import NicotineRewardModel, bolus_input
    
    print("Generating demonstration plots...")
    
    # Create model and simulate
    model = NicotineRewardModel()
    I = bolus_input(amplitude=5.0, t_start=10, duration=5)
    result = model.simulate(t_span=(0, 300), input_func=I)
    
    # Time series
    vis = Visualizer()
    fig1 = vis.plot_timeseries(result, 
                               variables=['N', 'R_a', 'R_d', 'D'],
                               title='Single Cigarette Response')
    
    # Phase portrait
    fig2 = vis.plot_phase_portrait_2d(result, 'R_a', 'D',
                                      title='Receptor Activation vs Dopamine')
    
    plt.show()
    print("Visualization demonstration complete.")
