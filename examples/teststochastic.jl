using TurbulenceTools

for ki = [8, 16]
  for μ = [0, 1e-2]

    runstochasticforcingproblem(tf=10, dt=5e-3, withplot=true, ns=10, n=256,
      ν=5e-3, μ=μ, stepper="FilteredRK4")

  end
end

