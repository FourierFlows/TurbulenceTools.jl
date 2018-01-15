using PyPlot

n = 512
L = 2π
Δ = L/n

x = reshape(linspace(-L/2, L/2-Δ, n), (n, 1))
y = reshape(x, (1, n))

X = [ x[i] for i=1:n, j=1:n ]
Y = [ y[j] for i=1:n, j=1:n ]

nk = Int(n/2+1)
i₁ = 0:Int(n/2)
i₂ = Int(-n/2+1):-1

k = reshape(2π/L*i₁, (nk, 1))
l = reshape(2π/L*cat(1, i₁, i₂), (1, n))

K = [ k[i] for i=1:nk, j=1:n ]
L = [ l[i] for i=1:nk, j=1:n ]

nk = 4
θ = 2π*rand() 
ξ = 2π*rand() 
i₁ = round(Int, abs(nk*cos(θ))) + 1
j₁ = round(Int, abs(nk*sin(θ))) + 1 # j₁ >= 1
j₂ = n - j₁ + 2 # j₁ = 1 => lf = 0; etc.


fh = zeros(K)

#fh[i₁, j₁] = 2*(0.5*n)^2
if j₁ != 1
  fh[i₁, j₁] = (0.5*n)^2 * exp(im*ξ)
  fh[i₁, j₂] = (0.5*n)^2 * exp(im*ξ)
else
  fh[i₁, j₁] = 2*(0.5*n)^2 * exp(im*ξ)
end

f = irfft(fh, n)

close("all")
imshow(f)
colorbar()
