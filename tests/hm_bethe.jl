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

#global variables
ed = EDIjl2.global_env

EDIjl2.read_input("inputED.conf")

#####################################################

Eband = range(-wband, wband, length=1000)
de = step(Eband)
Dband = dens_bethe.(Eband, wband)
wm = (pi / ed.beta) .* (2 .* (0:(ed.Lmats - 1)) .+ 1)


hloc = zeros(ComplexF64, 1, 1, 1, 1)
Gmats = zeros(ComplexF64, 1, 1, 1, 1,ed.Lmats)
Delta = zeros(ComplexF64, 1, 1, 1, 1,ed.Lmats)

EDIjl2.set_hloc(hloc)
global bath = EDIjl2.init_solver()

for iloop in 0:100

  converged = false  
  
  bath_old = copy(bath)
  
  EDIjl2.solve(bath)
  
  dens = EDIjl2.get_dens()
  mag = EDIjl2.get_mag()
  docc = EDIjl2.get_docc()
  phisc = EDIjl2.get_phi()
  println("Density = ", dens)
  println("Magnetization = ", mag)
  println("Double occupation = ", docc)
  println("Superconductive phi = ", phisc)
  
  gimp = EDIjl2.get_gimp(axis="m")
  smats = EDIjl2.get_sigma(axis="m")

  zeta = wm * im .- smats[1, 1, 1, 1, :]

  Gmats[1,1,1,1,:] .= sum(Dband' ./ (zeta .- Eband'), dims=2) .* de
  Delta[1,1,1,1,:] .= 0.25 * wband * Gmats[1,1,1,1,:]

  writedlm("Delta_iw.dat", [wm imag.(Delta[1,1,1,1,:]) real.(Delta[1,1,1,1,:])])

  bath_old = copy(bath)
  global bath = EDIjl2.chi2_fitgf(Delta,bath_old)

  result, converged = EDIjl2.check_convergence(Gmats)
  global bath = 0.3.*bath + (1.0-0.3).*bath_old
  
  
  if converged
    break
  end
  
end
