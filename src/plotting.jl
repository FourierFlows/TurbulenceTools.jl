using PyPlot

export makeplot

function makeplot(prob, diags; stochasticforcing=false)

  TwoDTurb.updatevars!(prob)  
  E, Z, D, I, R = diags

  close("all")
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count

  # dEdt = I - D - R?
  dEdt₁ = I[ii] - D[ii] - R[ii]
  residual = dEdt - dEdt₁

  plot(E.time[ii], I[ii], label="injection (\$I\$)")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  plot(E.time[ii], residual, "c-", label="residual")

  if !stochasticforcing
    plot(E.time[ii], dEdt₁, label=L"I-D-R")
    plot(E.time[ii], dEdt, "k:", label=L"E_t")
  end

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()
  pause(0.1)

  nothing
end
