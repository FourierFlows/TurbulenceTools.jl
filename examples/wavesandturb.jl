using PyPlot, FourierFlows, JLD2, TurbulenceTools,
      FourierFlows.VerticallyCosineBoussinesq

@load "vorticity.jld2" q

L = 1600e3
n = 256

# Initial turbulence parameters:
Ro = 0.1      # Rossby number
Rμ = 1e-4     # Ratio of drag and inertial time-scale
Rν = 2e1      # Ratio of viscous and inertial time-scale
 f = 1e-4     # Planetary vorticity

tq = 1/(f*Ro) # Eddy turnover time-scale
tμ = 1/(f*Rμ) # Bottom drag time-scale
tν = 1/(f*Rν) # Viscous time-scale

 μ = 1/tμ
kν = 2π/L * n/3 
 ν = 1/(kν*tν)
ki = 16*2π/L
fi = f*Ro/ki * sqrt(μ) # P = fi²

nu0 = ν
nu1 = ν
mu0 = μ
mu1 = μ

nnu0 = 2
nnu1 = 2
nmu0 = 0
nmu1 = 0


f = 1e-4
N = 5e-3

ε = 0.05
α = 1
σ = f*sqrt(α+1)
nkw = 16
k = nkw*2π/L
m = N*k/sqrt(σ^2-f^2)


tq0 = 1/maximum(abs.(q))

# ε = U / σ L = q / σ = u k / σ. L \approx U / q => u = ε σ U / q?
# σ = f sqrt[ 1 + (Nk/m)^2 ] = sqrt(2)-4
# tσ = 1/σ = 0.25 - 0.5
# tZ = Z^{-1}

 σ = sqrt(f^2 + (N*k/m)^2)
tσ = 2π/σ
dt = tσ/20
ni = round(Int, 10tσ/dt)
#dt = tq0/5 #tσ/100
#ni = 100 #round(Int, tσ/dt)
ns = 400
nt = ns*ni

# Equilibration e1 ~ u*aᵤ/μ
# u² ~ e1 ~ u*aᵤ/μ => aᵤ ~ u*μ
q₀ = f*Ro
U₀ = 0.1
u₀ = ε*σ*U₀/q₀
iᵤ = round(Int, k*L/2π) + 1
aᵤ = u₀*μ * n^2/2
kii = ki*L/2π
amplitude = fi*ki/sqrt(dt) * n^2/4
function calcF!(F, sol, t, s, v, p, g)

  F[iᵤ, 1, 2] = aᵤ*exp(im*σ*t)

  if t == s.t # not a substep
    @views F[:, :, 1] .= 0.0
    # Vorticity forcing
    θk = 2π*rand() 
    phase = 2π*im*rand()
    i₁ = round(Int, abs(kii*cos(θk))) + 1
    j₁ = round(Int, abs(kii*sin(θk))) + 1  # j₁ >= 1
    j₂ = g.nl + 2 - j₁                    # e.g. j₁ = 1 => j₂ = nl+1
    if j₁ != 1  # apply forcing to l = (+/-)l★ mode
      F[i₁, j₁, 1] = amplitude*exp(phase)
      F[i₁, j₂, 1] = amplitude*exp(phase)
    else        # apply forcing to l=0 mode
      F[i₁, 1, 1] = 2amplitude*exp(phase)
    end
  end

  nothing
end

function makeplot!(axs, prob, diags)
  close("all")
  fig, axs = subplots(ncols=3, figsize=(12,4))
  plotwaveturb!(axs[1:2], prob)
  plotenergy!(axs[3], diags)
  tight_layout()
  nothing
end

function plotwaveturb!(axs, prob)
  sp = sqrt.(prob.vars.u.^2 + prob.vars.v.^2)
  εu = sp * k / σ
  RoZ = prob.vars.Z/f
  maxRo = maximum(RoZ)

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, RoZ, cmap="RdBu_r",
    vmin=-maxRo/2, vmax=maxRo/2)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, εu)

  makesquare!(axs)
  ticksoff!(axs)

  for ax in axs; sca(ax); colorbar(orientation="horizontal"); end

  nothing
end

function plotenergy!(ax, diags)

  E, E0, E1 = diags
  ii = 1:E.count
  t = E.time[ii]

  sca(ax)
  plot(t, E[ii]/E[1], "k-",  label="Total")
  plot(t, E0[ii]/E[1], "--", label="Barotropic")
  plot(t, E1[ii]/E[1], "--", label="Baroclinic")
  xlabel(L"t")
  ylabel("Energy")
  legend()

  nothing
end
  

prob = VerticallyCosineBoussinesq.Problem(f=f, N=N, m=m,
  nx=n, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, 
  mu1=mu1, nmu1=nmu1, dt=dt, stepper="FilteredRK4", calcF=calcF!)

VerticallyCosineBoussinesq.set_Z!(prob, q)
VerticallyCosineBoussinesq.set_planewave!(prob, u₀, kw)

E = Diagnostic(totalenergy, prob, nsteps=nt)
E0 = Diagnostic(mode0energy, prob, nsteps=nt)
E1 = Diagnostic(mode1energy, prob, nsteps=nt)
diags = [E, E0, E1]

@printf "1/m: %.1f m, σ/f: %.3f, max Ro: %.3f\n" 1/m σ/f maximum(abs.(q))/f

for i = 1:ns
  @time stepforward!(prob, diags, ni)
  makeplot!(axs, prob, diags)
  pause(0.1)
end
