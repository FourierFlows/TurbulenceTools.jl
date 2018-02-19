module VerticallyCosineTools

using TurbulenceTools, FourierFlows, FourierFlows.VerticallyCosineBoussinesq, 
      PyPlot, JLD2

import FourierFlows.VerticallyCosineBoussinesq
import FourierFlows: jacobian, jacobianh

#import TurbulenceTools.TwoDTurb

export startfromfile, runproblem, makeproblem, makeplot!,
       loadgridparams, loadtimestep, loadparams, loadforcingparams, loadlastsolution

calcN_forced! = VerticallyCosineBoussinesq.calcN_forced!
calcN! = VerticallyCosineBoussinesq.calcN!

function usigvsig(prob, σ; forced=false)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  @views @. v.uh = s.sol[:, :, 2]
  @views @. v.vh = s.sol[:, :, 3]

  # RHS[:, :, 4] = ph_t
  RHS = zeros(eltype(s.sol), size(s.sol))
  if forced
    calcN_forced!(RHS, s.sol, prob.t, s, v, p, g)
  else
    calcN!(RHS, s.sol, prob.t, s, v, p, g)
  end
  @. RHS += prob.eqn.LC * s.sol

  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)

  @views ut = irfft(RHS[:, :, 2], g.nx)
  @views vt = irfft(RHS[:, :, 3], g.nx)

  usig = @views exp(im*σ*prob.t)/(2*p.f) * (v.u .+ im/σ*ut)
  vsig = @views exp(im*σ*prob.t)/(2*p.f) * (v.v .+ im/σ*vt)

  usig, vsig
end

function waveapv(usig, vsig, f, sig, g)
  
  j1 = real.(
    im.*jacobian(conj.(usig), usig, g) .+ im.*jacobian(conj.(vsig), vsig, g))
  j2 = real.(
    jacobian(conj.(vsig), usig, g) .+ jacobian(vsig, conj.(usig), g))

  u2h = rfft(abs2.(usig))
  v2h = rfft(abs2.(vsig))

  uv = @. real(usig*conj(vsig) + conj(usig)*vsig)
  uvh = rfft(uv) 

  plaph = @. -(g.kr^2*u2h + g.l^2*v2h + g.kr*g.l*uvh)
  plap = irfft(plaph, g.nx)

  @. j1/sig + f/(2*sig^2)*(j2 + plap)
end

function waveapv(prob, sig)
  usig, vsig = usigvsig(prob, sig)
  waveapv(usig, vsig, prob.params.f, sig, prob.grid)
end

function waveinducedspeed(prob, sig)
  g = prob.grid
  qw = waveapv(prob, sig)
  qwh = rfft(qw)
  uwh = @.  im*g.kr * g.invKKrsq * qwh
  vwh = @. -im*g.l  * g.invKKrsq * qwh
  uw = irfft(uwh, g.nx)
  vw = irfft(vwh, g.nx)
  @. sqrt(uw^2+vw^2)
end

"""
    loadgridparams(filename)

Returns nx, Lx, ny, Ly from the FourierFlows output stored in filename.
"""
function loadgridparams(filename)
  file = jldopen(filename)
  nx = file["grid/nx"]
  ny = file["grid/ny"]
  Lx = file["grid/Lx"]
  Ly = file["grid/Ly"]
  close(file)
  nx, Lx, ny, Ly
end

"""
    loadtimestep(filename)

Returns dt from the FourierFlows output stored in filename.
"""
function loadtimestep(filename)
  file = jldopen(filename)
  dt = file["timestepper/dt"]
  dt
end 

"""
    loadparams(filename)

Returns ν, nν, μ, nμ from the FourierFlows output stored in filename.
"""
function loadparams(filename)
  file = jldopen(filename)
   ν = file["params/ν"]
  nν = file["params/nν"]
   μ = file["params/μ"]
  nμ = file["params/nμ"]
  close(file)
  ν, nν, μ, nμ
end

"""
    loadforcingparams(filename)

Returns fi, ki from the FourierFlows output stored in filename.
"""
function loadforcingparams(filename)
  file = jldopen(filename)
  fi = file["forcingparams/fi"]
  ki = file["forcingparams/ki"]
  fi, ki
end

"""
    loadlastsolution(filename)

Returns the value of :sol with the highest timestep in the timeseris stored 
in the FourierFlows output file filename.
"""
function loadlastsolution(filename)
  file = jldopen(filename)
  laststep = parse(keys(file["timeseries/sol"])[end])
  sol = file["timeseries/sol/$laststep"]
  t = file["timeseries/t/$laststep"]
  laststep, t, sol
end

#loadgridparams(filename) = TwoDTurb.loadgridparams(filename)
#loadtimestep(filename) = TwoDTurb.loadtimestep(filename)
#loadparams(filename) = TwoDTurb.loadparams(filename)
#loadforcingparams(filename) = TwoDTurb.loadforcingparams(filename)
#loadlastsolution(filename) = TwoDTurb.loadlastsolution(filename)

"""
    cfl(prob)

Returns the CFL number defined by CFL = max([max(U)*dx/dt max(V)*dy/dt]).
"""
function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.U)/prob.grid.dx, maximum(prob.vars.V)/prob.grid.dy,
     maximum(prob.vars.u)/prob.grid.dx, maximum(prob.vars.v)/prob.grid.dy  ])
end

"""
    prob, diags = startfromfile(filename; kwargs...)

Start a stochastically-forced wave-turbulence interaction problem from the
vorticity field stored in the file with path filename.
"""
function startfromfile(filename; stepper="RK4", f=1.0, N=1.0, m=1.0,
                       nu1=nothing, nnu1=nothing, mu1=nothing, nmu1=nothing, 
                       ε=0.1, nkw=16, tf=1, ns=1, withplot=false, dt=nothing,
                       plotname=nothing, message=nothing)

  # Extract two-dimensional turbulence parameters
  file = jldopen(filename, "r")
    fi = file["forcingparams/fi"]
    ki = file["forcingparams/ki"]
   nu0 = file["params/ν"]
  nnu0 = file["params/nν"]
   mu0 = file["params/μ"]
  nmu0 = file["params/nμ"]
    nx = file["grid/nx"]
    Lx = file["grid/Lx"]
  
  if dt == nothing
    dt = file["timestepper/dt"]
  end

  laststep = keys(file["timeseries/sol"])[end]
  qh = file["timeseries/sol/$laststep"]
  close(file)

  if  nu1 == nothing;  nu1=nu0;  end
  if  mu1 == nothing;  mu1=mu0;  end
  if nnu1 == nothing; nnu1=nnu0; end
  if nmu1 == nothing; nmu1=nmu0; end
  
  g = TwoDGrid(nx, Lx)
  Uh =  im * g.l  .* g.invKKrsq .* q0h
  Vh = -im * g.kr .* g.invKKrsq .* q0h
  U = irfft(Uh, nx)
  V = irfft(Vh, nx)
  q = irfft(qh, nx)
  
  U₀ = 0.1
  u₀ = ε*σ*U₀/q₀
  iᵤ = round(Int, k*L/2π) + 1
   Γ = 0.01



  prob, diags, nt = makeproblem(q0, u0; n=nx, L=Lx, nu0=nu0, nnu0=nnu0, 
    nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1, dt=dt,
    f=f, N=N, m=m, ε=ε, nkw=nkw, fi=fi, ki=ki, tf=tf, stepper=stepper,
    q0=q0)

  fileprefix = filename[1:end-5]
  newfilename = "$fileprefix-waveturb.jld2"
  rm(newfilename, force=true)
  output = getbasicoutput(prob; filename=newfilename, filedir=".")

  runproblem(prob, diags, nt; ns=ns, withplot=withplot, output=output,
    plotname=plotname, message=message)

  prob, diags
end

"""
    runproblem(prob, diags, nt; kwargs...)

Run a stochastically forced wave-turbulence problem. 
"""
function runproblem(prob, diags, nt; ns=1, withplot=false, output=nothing,
                    plotname=nothing, message=nothing)

  if withplot; 
    fig, axs = subplots(ncols=3, figsize=(12, 4))
    makeplot!(axs, prob, diags)
    test = readline()
  end

  nint = round(Int, nt/ns)
  for i = 1:ns
    tic()
    stepforward!(prob, diags, nint)
    tc = toq()
    updatevars!(prob)  
    
    @printf(
      "step: %04d, t: %.2e, cfl: %.3f, tc: %.2f s\n",
      prob.step, prob.t, cfl(prob), tc)

    if message != nothing; println(message(prob)); end

    if withplot     
      makeplot!(axs, prob, diags)
      if plotname != nothing
        plotdir = "plots"
        fullplotname = joinpath(plotdir,
          @sprintf("%s_%d.png", plotname, prob.step))
        if !isdir(plotdir); mkdir(plotdir); end
        savefig(fullplotname, dpi=240)
      end
    end

    if output != nothing
      println(prob.step, " ", output.filename)
      saveoutput(output)
    end
  end

  updatevars!(prob)
  nothing
end

    
"""
    makeproblem(; kwargs...)

Returns a Problem, vector of Diagnostics, and number of timesteps for an
initiated stochastically-forced wave-turbulence interaction problem.
"""
function makeproblem(q0, u0; n=128, L=2π, nu0=1e-6, nnu0=1, dt=1.0,
  nu1=1e-6, nnu1=1, mu0=1e-6, nmu0=1, mu1=1e-6, nmu1=1, f=1.0, N=1.0, m=4.0, 
  ε=0.1, nkw=16, fi=1.0, ki=8, tf=1, stepper="RK4")

  kii = ki*L/2π
  amplitude = fi*ki/sqrt(dt) * n^2/2

  aᵤ = u₀ * μ/Γ * n^2/2
  aᵥ = aᵤ * (-im*f/σ)
  aᵣ = aᵤ * (σ^2 - f^2)/σ

  q₀ = f*Ro
  U₀ = 0.1
  u₀ = ε*σ*U₀/q₀
  iᵤ = round(Int, k*L/2π) + 1
   Γ = 0.01

  function calcF!(F, sol, t, s, v, p, g)

    if t == s.t # not a substep
      F .= 0.0
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

  nt = round(Int, tf/dt)
  prob = VerticallyCosineBoussinesq.Problem(f=f, N=N, m=m, nx=n, Lx=L, nu0=nu0, 
    nnu0=nnu0, nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1, 
    dt=dt, stepper=stepper, calcF=calcF!)

  # wave nonlinearity is measured by ε = uw*ki/σ
  kw = nkw*2π/L
  σ = f*sqrt(1 + (N*kw/m)^2) # wave frequency
  uw = ε*σ/ki # wave velocity defined in terms of wave nonlinearity
    
  # For now, initial wave condition is a plane wave.
  #set_planewave!(prob, uw, kw)
  set_Z!(prob, q0)

  diags = getdiags(prob, nt)

  prob, diags, nt
end


"""
    getdiags(prob, nt)

Returns a vector of diagnostics for a stochastically-forced wave-turbulence
interaction problem in the rigid lid Boussinesq system truncated to two 
vertical modes with constant stratification and rotation.
"""
function getdiags(prob, nt)
  E = Diagnostic(totalenergy, prob, nsteps=nt)
  E0 = Diagnostic(mode0energy, prob, nsteps=nt)
  E1 = Diagnostic(mode1energy, prob, nsteps=nt)
  #D0 = Diagnostic(mode0dissipation, prob, nsteps=nt)
  #D1 = Diagnostic(mode1dissipation, prob, nsteps=nt)
  #R0 = Diagnostic(mode0drag, prob, nsteps=nt)
  #R1 = Diagnostic(mode1drag, prob, nsteps=nt)
  [E, E0, E1]
end

"""
    makeplot!(axs, prob, diags)

Makes a plot of stochastically-forced waves and turbulence.
"""
function makeplot!(axs, prob, diags)
  E, E0, E1 = diags
  t = E.time

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.Z)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, real.(prob.vars.u))
  makesquare!(axs[1:2])

  #=
  sca(axs[3]); cla()
  plot(t, E,  label=L"\mathcal{E}") 
  plot(t, E0, label=L"E") 
  plot(t, E1, label=L"e")
  legend()
  xlabel(L"t")
  ylabel("Energy")
  =#

  axs[1][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)
  axs[2][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)

  tight_layout()
  nothing
end

end # module
