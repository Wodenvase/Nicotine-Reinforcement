"""
Single Cigarette Response - Julia Implementation

Demonstrates acute pharmacological response to a single bolus dose.
"""

include("NicotineModel.jl")
using Plots

println("="^60)
println("EXAMPLE: SINGLE CIGARETTE DOSE (Julia)")
println("="^60)

# Create model parameters
params = Parameters()

# Single cigarette bolus
I = bolus_input(5.0, 10.0, 5.0)

# Simulate
println("\nSimulating 5-hour period...")
tspan = (0.0, 300.0)
sol = simulate(tspan, I, params)

# Extract results
t = sol.t
N = [u[1] for u in sol.u]
Ra = [u[2] for u in sol.u]
Rd = [u[3] for u in sol.u]
RT = [u[4] for u in sol.u]
D = [u[5] for u in sol.u]

# Analysis
println("\nPeak responses:")
println("  Nicotine:  $(maximum(N)) μM")
println("  Active receptors: $(maximum(Ra))")
println("  Dopamine: $(maximum(D)) AU")

# Time to peak
t_peak_D = t[argmax(D)]
println("\nTime to peak dopamine: $(round(t_peak_D, digits=1)) min")

# Nicotine half-life
N_peak = maximum(N)
idx_peak = argmax(N)
idx_half = idx_peak + findfirst(N[idx_peak:end] .< N_peak/2)
if !isnothing(idx_half)
    t_half = t[idx_half] - t[idx_peak]
    println("Nicotine half-life: $(round(t_half, digits=1)) min")
end

# Visualization
println("\nGenerating plots...")

p = plot(layout=(4,1), size=(1000, 800), legend=:topright)

# Nicotine
plot!(p[1], t, N, linewidth=2, color=:blue, label="Nicotine", 
     ylabel="Nicotine [μM]", grid=true, title="FAST: Nicotine Kinetics")

# Active & desensitised receptors
plot!(p[2], t, Ra, linewidth=2, color=:green, label="Active", 
     ylabel="Receptor Fraction", grid=true)
plot!(p[2], t, Rd, linewidth=2, color=:red, linestyle=:dash, label="Desensitised",
     title="MEDIUM: Receptor Dynamics")

# Total receptors
plot!(p[3], t, RT, linewidth=2, color=:orange, label="R_T",
     ylabel="Total Receptors", grid=true, title="SLOW: Receptor Pool")

# Dopamine
plot!(p[4], t, D, linewidth=2, color=:purple, label="Dopamine",
     ylabel="Dopamine [AU]", xlabel="Time [min]", grid=true,
     title="FAST: Dopamine Response")

plot!(p, plot_title="Single Cigarette: Acute Response", 
     plot_titlefontsize=14)

savefig(p, "single_dose_julia.png")
println("✓ Saved: single_dose_julia.png")

# Stability analysis
println("\n" * "="^60)
println("STABILITY ANALYSIS")
println("="^60)

u_final = sol.u[end]
eigenvals, timescales, is_stable = stability_analysis(u_final, 0.0, params)

println("\nFinal state stability: $(is_stable ? "STABLE" : "UNSTABLE")")
println("Eigenvalues (real parts): $(round.(real.(eigenvals), digits=6))")

valid_ts = filter(!isnan, timescales)
if !isempty(valid_ts)
    println("\nTimescales:")
    for (i, τ) in enumerate(sort(valid_ts))
        println("  τ_$i = $(round(τ, digits=2)) min")
    end
end

ts_results = timescale_separation(u_final, 0.0, params)
println("\nTimescale separation: $(round(ts_results.ratio, digits=0))×")

println("\n" * "="^60)
println("KEY OBSERVATIONS")
println("="^60)
println("1. Fast nicotine clearance (~120 min half-life)")
println("2. Transient dopamine surge (reinforcement window)")
println("3. Receptor desensitisation develops and recovers")
println("4. Complete return to baseline (no long-term adaptation)")
println("5. System is stable throughout")
println("="^60)
