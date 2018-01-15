include("./turbtools.jl")

using PyPlot

  n = 128
  L = 2π
  ν = 1e-8       # Laplacian viscosity
 nν = 4
  μ = 0.0        # Large-scale drag
 dt = 1e-2       # Time step
  F = 1e-2       # Stochastic forcing
nkF = 16
 nt = 100000     # Number of time steps

prob = StochasticProblem(n, L, ν, nν, μ, dt)

#q₀ = rand(n, n)
q₀ = zeros(n, n)
qh₀ = rfft(q₀)
@. prob.vars.qh = qh₀

# Step forward
plots = 20
fig = figure(); tic()
for i = 1:plots
  stepforward!(prob, round(Int, nt/plots))

  cfl = maximum(prob.vars.u)*prob.grid.dx/prob.params.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  clf(); imshow(prob.vars.q); pause(0.01)
end


# Final plot
close("all")
fig = figure()
q = irfft(prob.vars.qh, n)
X, Y = prob.grid.X, prob.grid.Y

pcolormesh(X, Y, q)

colorbar()
xlabel(L"x")
ylabel(L"y")
title("\$q(x, y, t = $(prob.t))\$")

show()
