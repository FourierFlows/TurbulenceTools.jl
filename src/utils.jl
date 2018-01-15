import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

export cfl, getresidual, getdiags, getbasicoutput, savediags

function cfl(prob)
  prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])
end


function getresidual(prob, E, I, D, R; i₀=1)
  dEdt₀ = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count
  dEdt₁ = I[ii] - D[ii] - R[ii] # dEdt = I - D - R
  dEdt₀ - dEdt₁
end

function getresidual(prob, diags; kwargs...)
  E, Z, D, I, R = diags
  getresidual(prob, E, I, D, R; kwargs...)
end


function getdiags(prob, nt)
  E = Diagnostic(energy,      prob, nsteps=nt)
  Z = Diagnostic(enstrophy,   prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  I = Diagnostic(injection,   prob, nsteps=nt)
  diags = [E, Z, D, I, R]

  diags
end


function getbasicoutput(prob; filename="default")
  getsol(prob) = deepcopy(prob.state.sol)
  Output(prob, filename, (:sol, getsol))
end


function savediags(out, diags)
  E, Z, D, I, R = diags
  savediagnostic(E, "energy", out.filename)
  savediagnostic(Z, "enstrophy", out.filename)
  savediagnostic(D, "dissipation", out.filename)
  savediagnostic(I, "injection", out.filename)
  savediagnostic(R, "drag", out.filename)
  nothing
end
