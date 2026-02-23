# Feedback Loop Dynamics of Nicotine in the Brain Reward System

## Core Question

How do fast pharmacological effects of nicotine interact with slower neural and receptor-level feedback to produce reinforcement, tolerance, and withdrawal?

## Scientific Framing

This project models nicotine dependence as a **nonlinear feedback system with multiple interacting timescales**, focusing on stability analysis and emergent behaviour in the mesolimbic dopamine pathway.

## System Overview

### Biological Substrate

Nicotine acts primarily on **α4β2 nicotinic acetylcholine receptors (nAChRs)** in the VTA → NAc dopamine pathway.

### Key Feedback Loops

1. **Fast positive feedback (reinforcement)**
   - Nicotine → nAChR activation → dopamine release → behavioural reinforcement

2. **Slow negative feedback (tolerance)**
   - Prolonged nicotine → receptor desensitisation → reduced response

3. **Ultra-slow positive feedback (dependence)**
   - Chronic exposure → receptor upregulation → allostatic shift

## Mathematical Model

### State Variables

- `N(t)` — brain nicotine concentration (μM)
- `R_a(t)` — active receptors (fraction)
- `R_d(t)` — desensitised receptors (fraction)
- `R_T(t)` — total receptor pool (slow adaptation)
- `D(t)` — dopamine level (arbitrary units)

### Multi-Timescale ODE System

**Nicotine pharmacokinetics (fast):**
```
dN/dt = I(t) - k_N · N
```

**Receptor activation & desensitisation (medium):**
```
dR_a/dt = k_on · N · (R_T - R_a - R_d) - k_off · R_a - k_des · R_a
dR_d/dt = k_des · R_a - k_res · R_d
```

**Dopamine dynamics (fast):**
```
dD/dt = k_D · R_a - k_clear · D
```

**Receptor upregulation (ultra-slow):**
```
dR_T/dt = ε · (D - D_0)
```

Where ε ≪ 1 captures allostatic adaptation.

## Analysis Components

### ✓ Timescale Separation
- Fast: nicotine clearance, dopamine dynamics (~minutes)
- Medium: receptor desensitisation (~hours)
- Slow: receptor upregulation (~days-weeks)

### ✓ Stability Analysis
- Steady states under constant nicotine exposure
- Jacobian eigenvalue analysis
- Loss of stability → withdrawal/craving dynamics

### ✓ Bifurcation Analysis
- Parameter continuation in intake rate `I(t)`
- Transition from casual use → dependence
- Hysteresis effects

### ✓ Control Theory Perspective
- Cigarettes vs nicotine patches = different `I(t)` profiles
- Spike vs sustained delivery → different reinforcement

## What This Model Shows

### ✅ Emphasis
- Qualitative dynamical regimes
- Parameter sensitivity and robustness
- Control strategies for harm reduction
- Emergent timescale separation

### ❌ Limitations
- Not predictive of individual addiction risk
- Not for clinical outcome prediction
- Simplified to essential mechanisms

## Key References

### Receptor Pharmacology
1. **Dani JA, Bertrand D** (2007). Nicotinic acetylcholine receptors and nicotinic cholinergic mechanisms of the central nervous system. *Annu Rev Pharmacol Toxicol* 47:699-729.

2. **Picciotto MR, Zoli M** (2002). Nicotinic receptors in aging and dementia. *J Neurobiol* 53(4):641-655.

3. **Quick MW, Lester RAJ** (2002). Desensitization of neuronal nicotinic receptors. *J Neurobiol* 53(4):457-478.

### Dopamine & Reward
4. **Schultz W** (2002). Getting formal with dopamine and reward. *Neuron* 36(2):241-263.

5. **Di Chiara G** (2000). Role of dopamine in the behavioural actions of nicotine related to addiction. *Eur J Pharmacol* 393(1-3):295-314.

### Systems/Mathematical Modeling
6. **Koob GF, Le Moal M** (2008). Neurobiological mechanisms for opponent motivational processes in addiction. *Phil Trans R Soc B* 363:3113-3123.

7. **Gutkin BS, Dehaene S, Changeux JP** (2006). A neurocomputational hypothesis for nicotine addiction. *PNAS* 103(4):1106-1111.

8. **Graupner M, Gutkin B** (2009). Modeling nicotinic modulation of dopamine in the basal ganglia. *BMC Neuroscience* 10(Suppl 1):P115.

### Pharmacokinetics
9. **Benowitz NL, Hukkanen J, Jacob P III** (2009). Nicotine chemistry, metabolism, kinetics and biomarkers. *Handb Exp Pharmacol* 192:29-60.

### Allostatic Theory
10. **Koob GF, Le Moal M** (1997). Drug abuse: Hedonic homeostatic dysregulation. *Science* 278(5335):52-58.

## Project Structure

```
Nicotine/
├── README.md                    # This file
├── RESULTS.md                   # 🔬 Computational results & scientific inferences
├── Project.toml                 # Julia dependencies
├── src/                         # Python implementation
│   ├── model.py                # Core ODE system
│   ├── stability.py            # Stability analysis tools
│   ├── bifurcation.py          # Bifurcation analysis
│   └── visualization.py        # Plotting utilities
├── julia/                       # Julia implementation (faster)
│   ├── NicotineModel.jl        # Core model
│   ├── single_dose.jl          # Single cigarette analysis
│   └── bifurcation_analysis.jl # Full bifurcation study
├── examples/                    # Python examples
│   ├── single_dose.py          # Acute administration
│   ├── chronic_exposure.py     # Steady-state regime
│   ├── withdrawal.py           # Cessation dynamics
│   ├── patch_vs_cigarette.py  # Delivery mode comparison
│   └── bifurcation_analysis.py # Full bifurcation diagram
├── requirements.txt             # Python dependencies
└── tests/                       # Unit tests
    └── test_model.py
```

## Installation

### Julia (Recommended - 3-4× faster)
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run analyses
cd julia
julia --project=.. single_dose.jl
julia --project=.. bifurcation_analysis.jl
```

### Python
```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run examples
python examples/single_dose.py
python examples/bifurcation_analysis.py
```

## Quick Start

### Single Cigarette Response (Julia)
```julia
include("julia/NicotineModel.jl")

params = Parameters()
I = bolus_input(5.0, 10.0, 5.0)  # amplitude, start, duration
sol = simulate((0.0, 300.0), I, params)

# Analyze timescales
u_final = sol.u[end]
ts = timescale_separation(u_final, 0.0, params)
println("Timescale ratio: $(ts.ratio)×")  # ~100-1000×
```

### Bifurcation Analysis
```julia
# Vary intake rate
I_range = range(0.0, 1.5, length=60)

# Compute steady states
for I_val in I_range
    sol = simulate((0.0, 100000.0), continuous_input(I_val), params)
    # Extract dopamine, R_T, etc.
end
```

## Key Results (See RESULTS.md for details)

### 🔬 Scientific Findings

1. **Multi-timescale separation confirmed**: 100-1000× ratio between fast (minutes) and slow (weeks) dynamics
   - Fast: nicotine clearance (τ ~ 2-20 min)
   - Medium: receptor desensitisation (τ ~ 20-100 min)  
   - Slow: receptor upregulation (τ ~ days-weeks)

2. **Bifurcation analysis**: Smooth continuous transition from casual use → dependence
   - No threshold or "addiction switch"
   - Regime boundaries: Casual (I < 0.3), Regular (0.3-0.8), Dependent (I > 0.8 μM/min)
   - Supports allostatic theory (Koob & Le Moal, 2008)

3. **Withdrawal dynamics**: Dopamine falls 15% below baseline, nadir at days 2-3
   - Explains acute withdrawal symptoms
   - Extended recovery (weeks) explains relapse vulnerability

4. **Delivery kinetics matter**: Cigarettes produce 1.4× higher dopamine peaks, 9× greater variability vs patches
   - Reinforcement driven by phasic spikes, not static levels
   - Control theory: patches act as low-pass filters

5. **Stability**: System remains stable across all intake rates (no oscillations)
   - Craving/relapse not from inherent dynamics but perturbation response in high-gain state

### 📊 Computational Performance

| Analysis | Python | Julia | Speedup |
|----------|--------|-------|---------|
| Single dose (300 min) | 0.8 s | 0.3 s | 2.7× |
| Chronic (30 days) | 15 s | 4 s | 3.8× |
| Bifurcation (60 pts) | 180 s | 45 s | 4.0× |

**Conclusion**: Julia's DifferentialEquations.jl provides superior performance for stiff ODE systems.

### 🧠 Clinical Implications

1. **NRT design**: Slow-release formulations reduce reinforcement by input shaping
2. **Cessation timing**: Intensive support needed weeks 1-4 (not just first days)
3. **Gradual reduction**: Tapering allows R_T downregulation, minimizing withdrawal
4. **Pharmacotherapy**: Varenicline mimics patch kinetics at receptor level

### 🎯 Model Validation

✓ Nicotine t₁/₂: 11.7 min (brain) vs 120 min (plasma) - appropriate for CNS model  
✓ Receptor desensitisation: 50% in 60 min (matches Quick & Lester 2002)  
✓ Receptor upregulation: 15-30% over 4 weeks (matches Picciotto & Zoli 2002)  
✓ Dopamine response: 1.5-2× baseline (matches Di Chiara 2000)

## Mathematical Depth

This is **not just simulation**—we provide:

✅ **Timescale separation**: Eigenvalue analysis, timescale ratios  
✅ **Stability analysis**: Jacobian computation, linearization  
✅ **Bifurcation theory**: One-parameter continuation, regime classification  
✅ **Control perspective**: Input shaping analysis (cigarettes vs patches)  
✅ **Validation**: Parameters grounded in neuropharmacology literature

---

*Computational neuropharmacology model — Systems neuroscience approach to addiction dynamics*

**For detailed results, mechanisms, and inferences → See [RESULTS.md](RESULTS.md)**
