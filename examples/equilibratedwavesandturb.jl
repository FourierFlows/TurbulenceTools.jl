using PyPlot, FourierFlows, JLD2, TurbulenceTools,
      FourierFlows.VerticallyCosineBoussinesq

#@load "data/vorticity2_ki8.jld2" q

filename = "equilq_n256_Ro2.0e-01_Rmu1.0e-04_Rnu_1.0e+00_ki16.jld2"
filepath = joinpath("data", filename)

file = jldopen(filepath)
step = keys(file["timeseries/sol"])[end]
qh = file["timeseries/sol/$step"]
nx = 2*(size(qh)[1]-1)
q = irfft(qh, nx)
close(file)

srand(100) # Set random number seed for reproducibility

include("sampleparams.jl")

nu0 = ν
nu1 = ν
mu0 = μ
mu1 = μ

nnu0 = 1
nnu1 = 1
nmu0 = 0
nmu1 = 0

N = 5e-3

ε = 0.2 
α = 3.0
σ = f*sqrt(α+1)
nkw = 16 
k = nkw*2π/L
m = N*k/sqrt(σ^2-f^2)

filename = @sprintf("waveturb_k%02d_ep%03d_alpha%03d.jld2", nkw, 100ε, 100α)

tq0 = 1/maximum(abs.(q))

# ε = U / σ L = q / σ = u k / σ. L \approx U / q => u = ε σ U / q?
 σ = sqrt(f^2 + (N*k/m)^2)
tσ = 2π/σ
dt = tσ/20
ni = round(Int, 10tσ/dt)
ns = 1000
nt = ns*ni

# Equilibration μ*e1 ~ Γ*u*aᵤ, where Γ is fraction of wave energy in k₁ mode
# u² ~ e1 ~ aᵤ/u*μ => aᵤ ~ u*μ/Γ
q₀ = f*Ro
U₀ = 0.1
u₀ = ε*σ*U₀/q₀
iᵤ = round(Int, k*L/2π) + 1
 Γ = 0.01

kii = ki*L/2π
amplitude = fi*ki/sqrt(dt) * n^2/4

aᵤ = u₀ * μ/Γ * n^2/2
aᵥ = aᵤ * (-im*f/σ)
aᵣ = aᵤ * (σ^2 - f^2)/σ

function calcF!(F, sol, t, s, v, p, g)

  F[iᵤ, 1, 2] = aᵤ*exp(-im*σ*t) # u
  F[iᵤ, 1, 3] = aᵥ*exp(-im*σ*t) # v
  F[iᵤ, 1, 4] = aᵣ*exp(-im*σ*t) # p

  if t == s.t # not a substep
    @views @. F[:, :, 1] = 0
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

function makeplot!(prob, diags)
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
  RoQ = mode0apv(prob)/f
  maxRo = maximum(RoQ)

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, RoQ, cmap="RdBu_r",
    vmin=-3maxRo/4, vmax=3maxRo/4)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, εu)

  makesquare!(axs)
  ticksoff!(axs)

  for ax in axs
    sca(ax)
    divider = 
    colorbar(orientation="horizontal")
  end

  nothing
end

function plotenergy!(ax, diags)
  E, E0, E1 = diags
  ii = 1:E.count
  t = E.time[ii]

  sca(ax)
  #plot(t, E[ii]/E[1], "k-",  label="Total")
  plot(t, E0[ii]/E0[1], label="\$ \\Delta \$Barotropic")
  plot(t, E1[ii]/E1[1], label="\$ \\Delta \$Baroclinic")
  xlabel(L"t")
  ylabel("Energy")
  legend()
  nothing
end

prob = VerticallyCosineBoussinesq.Problem(f=f, N=N, m=m,
  nx=n, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, 
  mu1=mu1, nmu1=nmu1, dt=dt, stepper="FilteredRK4", calcF=calcF!)

VerticallyCosineBoussinesq.set_Z!(prob, q)
#VerticallyCosineBoussinesq.set_planewave!(prob, u₀, nkw)

E = Diagnostic(totalenergy, prob, nsteps=nt)
E0 = Diagnostic(mode0energy, prob, nsteps=nt)
E1 = Diagnostic(mode1energy, prob, nsteps=nt)
diags = [E, E0, E1]

getsol(prob) = deepcopy(prob.state.sol)
output = Output(prob, filename, (:sol, getsol))

@printf "1/m: %.1f m, σ/f: %.3f, max Ro: %.3f\n" 1/m σ/f maximum(abs.(q))/f

for i = 1:ns
  @time stepforward!(prob, diags, ni)
  saveoutput(output)
  makeplot!(prob, diags)
  plotname = @sprintf("%s_%09d.png", filename[1:end-5], prob.step)
  plotpath = joinpath("plots", plotname)
  savefig(plotpath, dpi=240)
end
