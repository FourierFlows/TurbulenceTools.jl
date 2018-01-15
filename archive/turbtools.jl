# Tools for a simple turbulence code.
__precompile__()

import Base: eltype

struct Grid{T}
  nx::Int
  Lx::T
  dx::T
  ny::Int
  Ly::T
  dy::T
  nk::Int
  nl::Int
  x::Array{T,2}
  y::Array{T,2}
  k::Array{T,2}
  l::Array{T,2}
  X::Array{T,2}
  Y::Array{T,2}
  Ksq::Array{T,2}
  invKsq::Array{T,2}
  rfftplan::Base.DFT.FFTW.rFFTWPlan{T, -1, false, 2}
  irfftplan::Base.DFT.ScaledPlan{Complex{T},
    Base.DFT.FFTW.rFFTWPlan{Complex{T}, 1, false, 2}, T}
end

function Grid(n, L; effort=FFTW.MEASURE, nthreads=Sys.CPU_CORES)
  T = typeof(L)
  Δ = T(L/n)

  x = reshape(linspace(T(-L/2), T(L/2-Δ), n), (n, 1))
  y = reshape(x, (1, n))

  X = [ x[i] for i=1:n, j=1:n ]
  Y = [ y[j] for i=1:n, j=1:n ]

  nk = Int(n/2+1)
  i₁ = 0:Int(n/2)
  i₂ = Int(-n/2+1):-1

  k = reshape(T(2π/L)*i₁, (nk, 1))
  l = reshape(T(2π/L)*cat(1, i₁, i₂), (1, n))

  Ksq = k.^2 .+ l.^2
  invKsq = 1./Ksq
  invKsq[1, 1] = 0

  # FFT plans
  FFTW.set_num_threads(nthreads)
  rfftplan  = plan_rfft(Array{T,2}(n, n); flags=effort)
  irfftplan = plan_irfft(Array{Complex{T},2}(nk, n), n; flags=effort)

  Grid{T}(n, L, Δ, n, L, Δ, nk, n, x, y, k, l, X, Y, Ksq, invKsq,
          rfftplan, irfftplan)
end

eltype{T}(::Type{Grid{T}}) = T
eltype(g::Grid) = eltype(typeof(g))


macro physvars(g, vars...)
  expr = Expr(:block)
  append!(expr.args, 
    [:( $(esc(var)) = zeros(
      eltype($(esc(g))), $(esc(g)).nx, $(esc(g)).ny);) for var in vars]
  )
  expr
end

macro transvars(g, vars...)
  expr = Expr(:block)
  append!(expr.args, 
    [:( $(esc(var)) = zeros(
      Complex{eltype($(esc(g)))}, $(esc(g)).nk, $(esc(g)).nl);) 
      for var in vars]
  )
  expr
end


struct Vars{T}
  u::Array{T,2}
  v::Array{T,2}
  q::Array{T,2}
  uq::Array{T,2}
  vq::Array{T,2}
  uh::Array{Complex{T},2}
  vh::Array{Complex{T},2}
  qh::Array{Complex{T},2}
  qth::Array{Complex{T},2}
  psih::Array{Complex{T},2}
  uqh::Array{Complex{T},2}
  vqh::Array{Complex{T},2}
end

function Vars(g)
  T = eltype(g)
  @physvars g u v q uq vq
  @transvars g uh vh qh qth psih uqh vqh
  Vars(u, v, q, uq, vq, uh, vh, qh, qth, psih, uqh, vqh)
end


struct Params{T}
  ν::T
  nν::Int
  μ::T
  dt::T
end

struct StochasticParams{T}
  ν::T
  nν::Int
  μ::T
  F::T
  nkF::Int
  dt::T
  sqrtdt::T
end

function StochasticParams(ν, nν, μ, F, nkF, dt)
  sqrtdt = sqrt(dt)
  StochasticParams{typeof(ν)}(ν, nν, μ, F, nkF, dt, sqrtdt)
end

abstract type AbstractProblem end

mutable struct Problem{T} <: AbstractProblem
  grid::Grid{T}
  vars::Vars{T}
  params::Params{T}
  t::T
  step::Int
end  

mutable struct StochasticProblem{T} <: AbstractProblem
  grid::Grid{T}
  vars::Vars{T}
  params::StochasticParams{T}
  t::T
  step::Int
end  

function Problem(n, L, ν, nν, μ, dt; kwargs...) 
  g = Grid(n, L; kwargs...)
  v = Vars(g)
  p = Params(ν, nν, μ, dt)
  Problem{typeof(L)}(g, v, p, 0.0, 0)
end

function StochasticProblem(n, L, ν, nν, μ, dt; kwargs...)
  g = Grid(n, L; kwargs...)
  v = Vars(g)
  p = StochasticParams(ν, nν, μ, F, nkF, dt)
  StochasticProblem{typeof(L)}(g, v, p, 0.0, 0)
end

function calcNL!(v, p, g)
  # Calculate right hand side of vorticity equation.
  @. v.qth = v.qh 
  A_mul_B!(v.q, g.irfftplan, v.qth)

  @. v.uh =  im*g.l*g.invKsq*v.qh
  @. v.vh = -im*g.k*g.invKsq*v.qh

  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)

  @. v.uq = v.u*v.q
  @. v.vq = v.v*v.q

  A_mul_B!(v.uqh, g.rfftplan, v.uq)
  A_mul_B!(v.vqh, g.rfftplan, v.vq)
end

function stepforward!(v, p::Params, g)
  calcNL!(v, p, g)
  
  @. v.qh += p.dt*(
    -im*g.k*v.uqh - im*g.l*v.vqh - (p.ν*g.Ksq^(0.5*p.nν) + p.μ)*v.qh
  )
end












function stepforward!(v, p::StochasticParams, g)
  calcNL!(v, p, g)
  
  # Step forward
  @. v.qh += p.dt*(
    -im*g.k*v.uqh - im*g.l*v.vqh - (p.ν*g.Ksq^(0.5*p.nν) + p.μ)*v.qh
    #- p.μ*v.qh # linear drag only
  )

  # Add stochastic forcing q_t = f*cos(kf*x + lf*y + θ)
  θ = 2π*rand() 
  ξ = 2π*rand() 
  i₁ = round(Int, abs(p.nkF*cos(θ))) + 1
  j₁ = round(Int, abs(p.nkF*sin(θ))) + 1 # j₁ >= 1
  j₂ = g.nl - j₁ + 2 # e.g. j₁ = 1 => j₂ = nl+1

  norm = g.nx*g.ny/4
  if j₁ != 1
    v.qh[i₁, j₁] += norm*sqrt(p.dt)*p.F * exp(im*ξ)
    v.qh[i₁, j₂] += norm*sqrt(p.dt)*p.F * exp(im*ξ)
  else
    v.qh[i₁, j₁] += 2*norm*sqrt(p.dt)*p.F * exp(im*ξ)
  end

  nothing 
end











function stepforward!(prob::AbstractProblem, nsteps)
  for i = 1:nsteps
    stepforward!(prob.vars, prob.params, prob.grid)
    prob.t += prob.params.dt
    prob.step += 1
  end
end
