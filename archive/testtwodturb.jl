include("./turbtools.jl")

using PyPlot

 T = Float64
 n = 512
 L = T(2π)
 ν = T(8e-5)       # Laplacian viscosity
 μ = T(0)          # Laplacian viscosity
nν = 2
dt = T(1e-2)       # Time step
nsteps = 5000      # Number of time steps

prob = Problem(n, L, ν, nν, μ, dt)

q₀ = rand(n, n)
qh₀ = rfft(q₀)
@. prob.vars.qh = qh₀

# Step forward
plots = 5
fig = figure(); tic()
for i = 1:plots
  stepforward!(prob, round(Int, nsteps/plots))

  cfl = maximum(prob.vars.u)*prob.grid.dx/prob.params.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  #clf(); imshow(prob.vars.q); pause(0.01)
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
