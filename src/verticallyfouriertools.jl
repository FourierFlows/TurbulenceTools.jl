module VerticallyFourierTools

using JLD2
using TurbulenceTools, FourierFlows, FourierFlows.VerticallyFourierBoussinesq 

import FourierFlows.VerticallyFourierBoussinesq: mode0apv
import TurbulenceTools.TwoDTurbTools
import FourierFlows: jacobian

export loadgridparams, loadtimestep, loadparams, loadsolution, waveapv, waveinducedspeed, pvinducedspeed,
       loadlastsolution, loadfirstsolution

loadgridparams = TwoDTurbTools.loadgridparams
loadtimestep = TwoDTurbTools.loadtimestep

"""
    loadparams(filename)

Returns ν, nν, μ, nμ from the FourierFlows output stored in filename.
"""
function loadparams(filename)
  file = jldopen(filename)
   nu0 = file["params/nu0"]
  nnu0 = file["params/nnu0"]
   mu0 = file["params/mu0"]
  nmu0 = file["params/nmu0"]
   nu1 = file["params/nu1"]
  nnu1 = file["params/nnu1"]
   mu1 = file["params/mu1"]
  nmu1 = file["params/nmu1"]
     f = file["params/f"]
     N = file["params/N"]
     m = file["params/m"]
  close(file)
  nu0, nnu0, mu0, nmu0, nu1, nnu1, mu1, nmu1, f, N, m
end

"""
    loadforcingparams(filename)

Returns fi, ki from the FourierFlows output stored in filename.
"""
function loadforcingparams(filename)
  file = jldopen(filename)
  fi = file["forcingparams/fi"]
  ki = file["forcingparams/ki"]
  close(file)
  fi, ki
end

"""
    loadsolution(filename, step)

Returns `(solc, solr)` from the specified `step` in `filename`.
"""
function loadsolution(filename, step)
  file = jldopen(filename)
  solc = file["timeseries/solc/$step"]
  solr = file["timeseries/solr/$step"]
     t = file["timeseries/t/$step"]
  close(file)
  t, solc, solr
end

function getsteps(filename)
  file = jldopen(filename)
  steps = parse.(keys(file["timeseries/t"]))
  close(file)
  steps
end

 loadlastsolution(filename) = loadsolution(filename, getsteps(filename)[end])
loadfirstsolution(filename) = loadsolution(filename, getsteps(filename)[1])

function waveapv(usig, vsig, f, sig, g)
  j1 = @. real(im*$jacobian(conj.(usig), usig, g) + im*$jacobian(conj.(vsig), vsig, g))
  j2 = @. real(   $jacobian(conj.(vsig), usig, g) +    $jacobian(vsig, conj.(usig), g))

  u2h = rfft(abs2.(usig))
  v2h = rfft(abs2.(vsig))

  uv = @. real(usig*conj(vsig) + conj(usig)*vsig)
  uvh = rfft(uv) 

  plaph = @. -(g.kr^2*u2h + g.l^2*v2h + g.kr*g.l*uvh)
  plap = irfft(plaph, g.nx)

  @. 2j1/sig + f/sig^2*(j2 + plap)
end

function waveapv(prob, sig)
  usig = @. exp(-im*sig*prob.state.t) * prob.vars.u
  vsig = @. exp(-im*sig*prob.state.t) * prob.vars.v
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

function pvinducedspeed(prob)
  g = prob.grid
  Q = mode0apv(prob)
  Qh = rfft(Q)
  Uwh = @.  im*g.kr * g.invKKrsq * Qh
  Vwh = @. -im*g.l  * g.invKKrsq * Qh
  Uw = irfft(Uwh, g.nx)
  Vw = irfft(Vwh, g.nx)
  @. sqrt(Uw^2 + Vw^2)
end

end # module
