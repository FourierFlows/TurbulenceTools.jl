using TurbulenceTools

 n = 256
 ν = 1e-6
tf = 200

fi = 1.0
 μ = 1e0

for ki = [16, 32, 64]

  # 1/1e-2 = 100...
  # 1/5e-3 = 200...
  # 1/2e-3 = 500...

  @printf "Running with fi=%.2f, ki=%d, and μ=%4.0e...\n" fi ki μ

  plotname = @sprintf("forcing_fi%02d_n%d_ki%02d_mu%4.0e", Int(fi*10), n, ki, μ)

  runstochasticforcingproblem(fi=fi, ki=ki, tf=tf, dt=5e-3, withplot=true, 
    ns=10, n=n, ν=ν, μ=μ, stepper="FilteredRK4", plotname=plotname,
    output=true, filename=plotname)

end
