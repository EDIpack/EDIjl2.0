push!(LOAD_PATH, joinpath(@__DIR__, "../EDIjl2/src"))

using DelimitedFiles
using EDIjl2

function dens_bethe(x, d)
    root = sqrt(Complex(1 - (x / d)^2))
    dens_bethe = (2 / (pi * d)) * root
    return real(dens_bethe)
end

# Script
wmixing = 0.3
wband = 1.0

libpath = find_EDIpack()
EDIsolver = InitEDIjl(libpath)
read_input(EDIsolver,"inputED.conf")

#####################################################

Eband = range(-wband, wband, length=1000)
de = step(Eband)
Dband = dens_bethe.(Eband, wband)
wm = (pi / EDIsolver.beta) .* (2 .* (0:(EDIsolver.Lmats - 1)) .+ 1)


hloc = zeros(ComplexF64, 1, 1, 1, 1)
Gmats = zeros(ComplexF64, 1, 1, 1, 1,EDIsolver.Lmats)
Delta = zeros(ComplexF64, 1, 1, 1, 1,EDIsolver.Lmats)

set_hloc(EDIsolver, hloc)
global bath = init_solver(EDIsolver)

for iloop in 0:100

  converged = false  
  
  bath_old = copy(bath)
  
  solve(EDIsolver, bath)
  gimp = get_gimp(EDIsolver; axis="m")
  smats = get_sigma(EDIsolver; axis="m")

  zeta = wm * im .- smats[1, 1, 1, 1, :]

  Gmats[1,1,1,1,:] .= sum(Dband' ./ (zeta .- Eband'), dims=2) .* de
  Delta[1,1,1,1,:] .= 0.25 * wband * Gmats[1,1,1,1,:]

  writedlm("Delta_iw.dat", [wm imag.(Delta[1,1,1,1,:]) real.(Delta[1,1,1,1,:])])

  bath_old = copy(bath)
  global bath = chi2_fitgf(EDIsolver,Delta,bath_old)

  result, converged = check_convergence(Gmats)
  global bath = 0.3.*bath + (1.0-0.3).*bath_old
  
  
  if converged
    break
  end
  
end
