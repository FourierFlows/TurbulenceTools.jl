#using TurbulenceTools.TwoDTurbTools.StochasticForcingProblems
using TurbulenceTools.VerticallyCosineTools.StochasticWaveTurbProblems
using JLD2

L = 1600e3
f = 1e-4
N = 5e-3

# σ^2 = f^2 + (N*k/m)^2
σ = 2f
nkw = 16
k = nkw*2π/L
m = N*k/sqrt(σ^2-f^2)
ε = 0.1

tσ = 2π/σ
dt = 1e-16*tσ
tf = 100dt
nt = round(Int, tf/dt)
ns = nt

filename = joinpath("data", "equilq_n256_Ro1.0e-01_Rmu1.0e-04_ki16.jld2")
plotname = "test"
savename = "$plotname.jld2"

@printf("Running stochastic wave turb for %d steps...\n", nt)

prob, diags = startfromfile(filename; stepper="FilteredRK4", f=f, N=N, m=m,
  ε=ε, nkw=nkw, tf=tf, ns=nt, withplot=true, plotname=plotname)
