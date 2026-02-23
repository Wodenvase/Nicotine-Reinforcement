"""
Nicotine Reward System Model

Computational neuropharmacology model of nicotine feedback dynamics.
"""

__version__ = "1.0.0"
__author__ = "Systems Neuroscience Lab"

from .model import NicotineRewardModel, bolus_input, continuous_input, repeated_bolus
from .stability import StabilityAnalyzer
from .bifurcation import BifurcationAnalyzer
from .visualization import Visualizer

__all__ = [
    'NicotineRewardModel',
    'bolus_input',
    'continuous_input', 
    'repeated_bolus',
    'StabilityAnalyzer',
    'BifurcationAnalyzer',
    'Visualizer'
]
