"""
    NicotineRewardModel

Multi-timescale ODE model for nicotine feedback dynamics in the brain reward system.

Julia implementation for high-performance computational neuropharmacology.

State variables:
- N: brain nicotine concentration [μM]
- Ra: active α4β2 nAChRs [fraction]
- Rd: desensitised receptors [fraction]  
- RT: total receptor pool
- D: dopamine level [AU]
"""

using DifferentialEquations
using LinearAlgebra
using Statistics

"""
    Parameters

Default physiologically-based parameter set.
"""
struct Parameters
    # Nicotine pharmacokinetics
    k_N::Float64        # clearance rate [1/min]
    
    # Receptor kinetics
    k_on::Float64       # binding rate [1/(μM·min)]
    k_off::Float64      # unbinding rate [1/min]
    k_des::Float64      # desensitisation rate [1/min]
    k_res::Float64      # resensitisation rate [1/min]
    
    # Dopamine dynamics
    k_D::Float64        # release rate [AU/min]
    k_clear::Float64    # clearance [1/min]
    
    # Allostatic adaptation
    epsilon::Float64    # upregulation rate [1/min]
    D_0::Float64        # baseline dopamine [AU]
    
    # Initial conditions
    RT_0::Float64       # initial receptor density
end

function Parameters()
    Parameters(
        0.06,      # k_N
        0.5,       # k_on
        0.1,       # k_off
        0.05,      # k_des
        0.01,      # k_res
        1.0,       # k_D
        0.5,       # k_clear
        0.0001,    # epsilon
        0.5,       # D_0
        1.0        # RT_0
    )
end

"""
    nicotine_dynamics!(du, u, p, t)

ODE system for nicotine reward dynamics.

# Arguments
- `du`: derivatives (output)
- `u`: state vector [N, Ra, Rd, RT, D]
- `p`: tuple (params, input_func)
- `t`: time
"""
function nicotine_dynamics!(du, u, p, t)
    params, input_func = p
    N, Ra, Rd, RT, D = u
    
    # Current input
    I = input_func(t)
    
    # Available receptors
    R_available = max(0.0, RT - Ra - Rd)
    
    # Derivatives
    du[1] = I - params.k_N * N                                              # dN/dt
    du[2] = params.k_on * N * R_available - params.k_off * Ra - params.k_des * Ra  # dRa/dt
    du[3] = params.k_des * Ra - params.k_res * Rd                          # dRd/dt
    du[4] = params.epsilon * (D - params.D_0)                               # dRT/dt
    du[5] = params.k_D * Ra - params.k_clear * D                           # dD/dt
end

"""
    baseline_state(params)

Return baseline equilibrium state (no nicotine).
"""
function baseline_state(params::Parameters)
    [0.0, 0.0, 0.0, params.RT_0, params.D_0]
end

"""
    simulate(tspan, input_func, params; u0=nothing, saveat=nothing)

Simulate the model over a time interval.

# Arguments
- `tspan`: time span (t_start, t_end)
- `input_func`: nicotine input function I(t)
- `params`: model parameters
- `u0`: initial state (uses baseline if nothing)
- `saveat`: time points for output

# Returns
- ODESolution object with time and state history
"""
function simulate(tspan, input_func, params::Parameters; 
                 u0=nothing, saveat=nothing)
    if isnothing(u0)
        u0 = baseline_state(params)
    end
    
    prob = ODEProblem(nicotine_dynamics!, u0, tspan, (params, input_func))
    
    if isnothing(saveat)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-10)
    else
        sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-10, saveat=saveat)
    end
    
    return sol
end

# === Input functions ===

"""
    bolus_input(amplitude, t_start, duration)

Single bolus dose (e.g., cigarette).
"""
function bolus_input(amplitude, t_start, duration)
    return t -> (t_start <= t < t_start + duration) ? amplitude : 0.0
end

"""
    continuous_input(rate)

Constant infusion (e.g., nicotine patch).
"""
function continuous_input(rate)
    return t -> rate
end

"""
    repeated_bolus(amplitude, interval, n_doses, duration=5.0)

Repeated dosing pattern.
"""
function repeated_bolus(amplitude, interval, n_doses, duration=5.0)
    return function(t)
        for i in 0:(n_doses-1)
            t_start = i * interval
            if t_start <= t < t_start + duration
                return amplitude
            end
        end
        return 0.0
    end
end

# === Analysis functions ===

"""
    compute_jacobian(u, I_const, params)

Compute Jacobian matrix at state u.
"""
function compute_jacobian(u, I_const, params::Parameters)
    N, Ra, Rd, RT, D = u
    R_available = max(0.0, RT - Ra - Rd)
    
    J = zeros(5, 5)
    
    # dN/dt derivatives
    J[1, 1] = -params.k_N
    
    # dRa/dt derivatives
    J[2, 1] = params.k_on * R_available
    J[2, 2] = -params.k_on * N - params.k_off - params.k_des
    J[2, 3] = -params.k_on * N
    J[2, 4] = params.k_on * N
    
    # dRd/dt derivatives
    J[3, 2] = params.k_des
    J[3, 3] = -params.k_res
    
    # dRT/dt derivatives
    J[4, 5] = params.epsilon
    
    # dD/dt derivatives
    J[5, 2] = params.k_D
    J[5, 5] = -params.k_clear
    
    return J
end

"""
    stability_analysis(u, I_const, params)

Perform eigenvalue analysis at state u.

Returns (eigenvalues, timescales, is_stable)
"""
function stability_analysis(u, I_const, params::Parameters)
    J = compute_jacobian(u, I_const, params)
    eigenvalues = eigvals(J)
    
    real_parts = real.(eigenvalues)
    is_stable = all(real_parts .< 0)
    
    timescales = 1.0 ./ abs.(real_parts)
    timescales[isinf.(timescales)] .= NaN
    
    return eigenvalues, timescales, is_stable
end

"""
    timescale_separation(u, I_const, params)

Quantify timescale separation.
"""
function timescale_separation(u, I_const, params::Parameters)
    _, timescales, _ = stability_analysis(u, I_const, params)
    
    valid_timescales = filter(!isnan, timescales)
    if isempty(valid_timescales)
        return (fast=NaN, slow=NaN, ratio=NaN)
    end
    
    sorted_ts = sort(valid_timescales)
    fast = sorted_ts[1]
    slow = sorted_ts[end]
    
    return (fast=fast, slow=slow, ratio=slow/fast)
end

println("✓ NicotineRewardModel loaded")
println("  Use: simulate(tspan, input_func, Parameters())")
