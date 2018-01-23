module StochasticWaveTurbProblems

export startfromfile, makestochasticforcingproblem, makeplot!


"""
    startfromfile(filename; kwargs...)

Start a stochastically-forced wave-turbulence interaction problem from the
vorticity field stored in the file with path filename.
"""
function startfromfile(filename; stepper="RK4", f=1.0, N=1.0, m=1.0,
                       nu1=nothing, nnu1=nothing, mu1=nothing, nmu1=nothing, 
                       ε=0.1, nkw=16, tf=1, ns=1, withplot=false, 
                       plotname=nothing, message=nothing)

  # Extract two-dimensional turbulence parameters
  jldopen(filename, "r") do file
      fi = file["forcingparams/fi"]
      ki = file["forcingparams/ki"]
     nu0 = file["params/ν"]
    nnu0 = file["params/nν"]
     mu0 = file["params/μ"]
    nmu0 = file["params/nμ"]
      nx = file["grid/nx"]
      Lx = file["grid/Lx"]

    laststep = keys(file["timeseries/sol"])[end]
    qh = file["timeseries/sol/$laststep"]
  end

  if  nu1 == nothing;  nu1=nu0;  end
  if  mu1 == nothing;  mu1=mu0;  end
  if nnu1 == nothing; nnu1=nnu0; end
  if nmu1 == nothing; nmu1=nmu0; end
  
  q0 = irfft(q0, nx)

  prob, diags, nt = makeproblem(; n=nx, L=Lx, nu0=nu0, nnu0=nnu0, 
    nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1
    f=f, N=N, m=m, ε=ε, nkw=nkw, fi=fi, ki=ki, tf=tf, stepper=stepper,
    q0=q0)

  fileprefix = filename[1:end-5]
  newfilename = "$fileprefix-waveturb.jld2"
  output = getbasicouput(prob; filename=newfilename)

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

  if withplot; fig, axs = subplots(ncols=3, figsize=(12, 4)); end

  nint = round(Int, nt/ns)
  for i = 1:ns
    tic()
    stepforward!(prob, diags, nint)
    tc = toq()
    updatevars!(prob)  
    
    @printf(
      "step: %04d, t: %.2e, cfl: %.3f, tc: %.2f s, <res>: %.3e, <I>: %.2f\n", 
      prob.step, prob.t, cfl(prob), tc, resnorm, avgI)

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

    if output != nothing; saveoutput(output); end
  end

  updatevars!(prob)
  nothing
end

    
"""
    makeproblem(; kwargs...)

Returns a Problem, vector of Diagnostics, and number of timesteps for an
initiated stochastically-forced wave-turbulence interaction problem.
"""
function makeproblem(; n=128, L=2π, nu0=1e-6, nnu0=1,
  nu1=1e-6, nnu1=1, mu0=1e-6, nmu0=1, mu1=1e-6, nmu1=1, f=1.0, N=1.0, m=4.0, 
  ε=0.1, nkw=16, fi=1.0, ki=8, tf=1, stepper="RK4", q0=nothing)

  kii = ki*L/2π
  amplitude = fi*ki/sqrt(dt) * n^2/4
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
  prob = Problem(f=f, N=N, m=m, nx=n, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, 
    nnu1=nnu1, mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1, dt=dt, 
    stepper=stepper, calcF=calcF!)

  # wave nonlinearity is measured by ε = uw*ki/σ
  kw = nkw*2π/L
  σ = f*sqrt(1 + (N*kw/m)^2) # wave frequency
  uw = ε*σ/ki # wave velocity defined in terms of wave nonlinearity
    
  set_planewave!(prob, uw, kw)

  if q0 != nothing; set_Z!(prob, q0); end

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
  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.Z)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, real.(prob.vars.u))
  makesquare!(axs[1:2])

  sca(axs[3]); cla()
  plot(t, E,    label=L"\mathcal{E}") 
  plot(t, E0,   label=L"E") 
  plot(t, E1,   label=L"e")
  legend()
  xlabel(L"t")
  ylabel("Energy")

  axs[1][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)
  axs[2][:tick_params]( 
    bottom=false, left=false, labelbottom=false, labelleft=false)

  tight_layout()
  nothing
end


end # module
