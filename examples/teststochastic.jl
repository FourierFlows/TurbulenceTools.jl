using TurbulenceTools

n = 128
ν = 5e-3
tf = 10/ν

for ki = [8, 16]
  for μ = [0, 1e-2]

    # 1/1e-2 = 100...
    # 1/5e-3 = 200...
    # 1/2e-3 = 500...

    @printf "Running with ki=%d and μ=%4.0e...\n" ki μ

    plotname = @sprintf("teststochastic_n%d_ki%02d_mu%4.0e", n, ki, μ)
    runstochasticforcingproblem(tf=tf, dt=1e-2, withplot=true, ns=10, n=n,
      ν=ν, μ=μ, stepper="FilteredRK4", plotname=plotname)

  end
end
