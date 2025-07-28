# SIM Model, Julia version
# Created by Marco Veronese Passarella 
# 28 July 2025

# ───────────────────────────────────────────────
# RECALL PACKAGES
# ───────────────────────────────────────────────
using Plots

# ───────────────────────────────────────────────
# PARAMETERS & SETTINGS
# ───────────────────────────────────────────────
nPeriods = 65
alpha1 = 0.6
alpha2 = 0.4
theta = 0.2

# ───────────────────────────────────────────────
# VARIABLES 
# ───────────────────────────────────────────────
y   = zeros(nPeriods)   # Income
y_star = zeros(nPeriods) # Steady-state income
c   = zeros(nPeriods)   # Consumption
g   = zeros(nPeriods)   # Government expenditures
t   = zeros(nPeriods)   # Taxes
yd  = zeros(nPeriods)   # Disposable income
h_h = zeros(nPeriods)   # Cash demand by households
h_s = zeros(nPeriods)   # Cash supply
n   = zeros(nPeriods)   # Employment
w   = ones(nPeriods)    # Wages (constant = 1)

# ───────────────────────────────────────────────
# MODEL DYNAMICS
# ───────────────────────────────────────────────
for i in 2:nPeriods # Time loop
    
    for iter in 1:20   # Iterations loop

        # Government expenditures shock after period 15
        if i >= 15
            g[i] = 20
        end

        # Consumption function
        c[i] = alpha1 * yd[i] + alpha2 * h_h[i-1]
        
        # Output determination
        y[i] = c[i] + g[i]

        # Steady-state solution
        y_star[i] = g[i] / theta

        # Employment
        n[i] = y[i] / w[i]

        # Taxes
        t[i] = theta * w[i] * n[i]

        # Disposable income (accounting identity)
        yd[i] = w[i] * n[i] - t[i]

        # Cash supply (government deficit)
        h_s[i] = h_s[i-1] + g[i] - t[i]

        # Cash held by households
        h_h[i] = h_h[i-1] + yd[i] - c[i]
    end
end

# ───────────────────────────────────────────────
# PLOT RESULTS
# ───────────────────────────────────────────────
plot(2:nPeriods, y_star[2:end],
    label="Steady state solution Y*",
    linewidth=2,
    linestyle=:solid,
    color=:blue,
    xlabel="Periods",
    ylabel="",
    ylim=(0,130),
    title="Figure 3.1: Impact of Y and Y* of a permanent increase in G",
    titlefontsize=10,     
    legendfontsize=10,     
    guidefontsize=10,      
    tickfontsize=10)       

plot!(2:nPeriods, y[2:end],
    label="Income Y",
    linewidth=2,
    linestyle=:solid,
    color=:green)