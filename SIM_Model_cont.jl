# Continuous Time SIM Model, Julia version
# Created by Marco Veronese Passarella
# 1 August 2025

# ───────────────────────────────────────────────
# RECALL PACKAGES
# ───────────────────────────────────────────────
using DifferentialEquations
using Plots

# ───────────────────────────────────────────────
# PARAMETERS & SETTINGS
# ───────────────────────────────────────────────
alpha1 = 0.6
alpha2 = 0.4
theta = 0.2
gov_before = 0.0  # Initial government expenditure
gov_after = 20.0  # Government expenditure after shock
t_switch = 15.0 # Time when shock occurs
tspan = (0.0, 65.0) # Simulation time span

# ───────────────────────────────────────────────
# MODEL DYNAMICS (Continuous-time ODE system)
# ───────────────────────────────────────────────
function sim_model!(du, u, p, t)
    y, h_h, h_s = u
    alpha1, alpha2, theta = p
    
    # Government expenditure (step function)
    gov = t >= t_switch ? gov_after : gov_before
    
    # Auxiliary variables
    yd = (1 - theta) * y
    c = alpha1 * yd + alpha2 * h_h
    
    # Differential equations
    du[1] = c + gov - y  # dy/dt: Income adjusts to demand
    du[2] = yd - c      # dh_h/dt: Change in household cash
    du[3] = gov - theta * y  # dh_s/dt: Government deficit
    
    # Note: In continuous time, we don't need iterations as the ODE solver handles
    # the simultaneous determination of variables through differential equations
end

# ───────────────────────────────────────────────
# INITIAL CONDITIONS
# ───────────────────────────────────────────────
u0 = [0.0, 0.0, 0.0]  # [y, h_h, h_s] initial values
p = (alpha1, alpha2, theta)

# ───────────────────────────────────────────────
# SOLVE THE ODE SYSTEM
# ───────────────────────────────────────────────
prob = ODEProblem(sim_model!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# ───────────────────────────────────────────────
# POST-PROCESSING FOR PLOTTING
# ───────────────────────────────────────────────
t = sol.t
y = sol[1,:]
h_h = sol[2,:]
h_s = sol[3,:]

# Calculate steady-state solution and other variables
y_star = zeros(length(t))
c = zeros(length(t))
yd = zeros(length(t))

for (i, ti) in enumerate(t)
    gov = ti >= t_switch ? gov_after : gov_before
    y_star[i] = gov / theta
    yd[i] = (1 - theta) * y[i]
    c[i] = alpha1 * yd[i] + alpha2 * h_h[i]
end

# ───────────────────────────────────────────────
# PLOT RESULTS
# ───────────────────────────────────────────────
plot(t, y_star,
    label="Steady state solution Y*",
    linewidth=2,
    linestyle=:solid,
    color=:blue,
    xlabel="Time",
    ylabel="",
    ylim=(0,130),
    title="Figure 3.1: Impact of Y and Y* of a permanent increase in G",
    titlefontsize=10,     
    legendfontsize=10,     
    guidefontsize=10,      
    tickfontsize=10)       

plot!(t, y,
    label="Income Y",
    linewidth=2,
    linestyle=:solid,
    color=:green)
