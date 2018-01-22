using PyPlot, FourierFlows, JLD2
import FourierFlows.VerticallyFourierBoussinesq
import FourierFlows.VerticallyFourierBoussinesq: totalenergy, mode0energy,
  mode1energy, mode0dissipation, mode1dissipation, mode0drag, mode1drag

@load "testturb.jld2" q

nx, ny = size(q) 
L =  2π   # Domain
nu0, nnu0 = 1e-5,   1  # Viscosity
nu1, nnu1 = 1e-3,   1  # Viscosity
mu0, nmu0 =  0.1,   0
mu1, nmu1 =  0.1,   0
  f, N, m =  1.0, 1.0, 8.0
   uw, kw = 0.02,  16 

# ε = U / σ L 
# σ = f sqrt[ 1 + (Nk/m)^2 ] = sqrt(2)-4
# tσ = 1/σ = 0.25 - 0.5
# tZ = Z^{-1}

 σ = f*sqrt(1 + (N*kw/m)^2)
tσ = 2π/σ
dt = tσ/200
ti = tσ/4
ni = round(Int, ti/dt)
ns = 10
nt = ns*ni

@printf "σ/f: %.2f, ε: %.3f\n" σ/f maximum(abs.(q))/σ

function makesquare!(axs)
  for ax in axs
    ax[:set_aspect](1, adjustable="box")
    ax[:set_aspect](1, adjustable="box")
  end
  nothing
end

function primp(diags)
  E, E0, E1, D0, D1, R0, R1 = diags
  ii = 1:E.count-1
  ii₊₁ = 2:E.count
  dEdt = (E[ii₊₁] - E[ii]) / dt
  E.time[ii], dEdt, E[ii], E0[ii], E1[ii], D0[ii], D1[ii], R0[ii], R1[ii]
end

function makeplot!(axs, prob, diags)
  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.Z)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, real.(prob.vars.u))
  makesquare!(axs[1:2])

  t, dEdt, E, E0, E1, D0, D1, R0, R1 = primp(diags)
  sca(axs[3]); cla()
  plot(t, E,    label=L"\mathcal{E}") 
  plot(t, E0,   label=L"E") 
  plot(t, E1,   label=L"e")
  legend()
  xlabel(L"t")
  ylabel("Energy")

  sca(axs[4]); cla()
  plot(t, dEdt, "k-", label=L"E_t")
  plot(t, -D0,   label=L"D")
  plot(t, -D1,   label=L"d")
  plot(t, -R0,   label=L"R")
  plot(t, -R1,   label=L"r")
  plot(t, -D0 - D1 - R0 - R1, "y:", label=L"-D-d-R-r")
  legend()
  xlabel(L"t")
  ylabel("Energy tendency")

  axs[1][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)
  axs[2][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)

  tight_layout()
  nothing
end

prob = VerticallyFourierBoussinesq.Problem(f=f, N=N, m=m,
  nx=nx, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, 
  mu1=mu1, nmu1=nmu1, dt=dt, stepper="RK4")

VerticallyFourierBoussinesq.set_Z!(prob, q)
VerticallyFourierBoussinesq.set_planewave!(prob, uw, kw)

E = Diagnostic(totalenergy, prob, nsteps=nt)
E0 = Diagnostic(mode0energy, prob, nsteps=nt)
E1 = Diagnostic(mode1energy, prob, nsteps=nt)
D0 = Diagnostic(mode0dissipation, prob, nsteps=nt)
D1 = Diagnostic(mode1dissipation, prob, nsteps=nt)
R0 = Diagnostic(mode0drag, prob, nsteps=nt)
R1 = Diagnostic(mode1drag, prob, nsteps=nt)
diags = [E, E0, E1, D0, D1, R0, R1]

fig, axs = subplots(ncols=2, nrows=2, figsize=(8,8))
for i = 1:ns
  @time stepforward!(prob, diags, ni)
  makeplot!(axs, prob, diags)
  pause(0.1)
end
