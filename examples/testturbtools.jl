include("turbtools.jl")

using TurbulenceTools

for ki = [8, 16]
  for μ = [0, 1e-3, 1e-2, 1e-1, 1]
    runsteadyforcingproblem(tf=100, dt=1e-3, withplot=true, ns=10, n=256,
      ν=1e-3, μ=1e-2, stepper="FilteredRK4")
    runstochasticforcingproblem(tf=100, dt=1e-3, withplot=true, ns=10, n=256,
      ν=1e-3, μ=1e-2, stepper="FilteredRK4")
  end
end
