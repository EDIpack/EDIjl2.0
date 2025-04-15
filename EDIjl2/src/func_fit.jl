using Libdl
using Base.Threads

function chi2_fitgf(args...; ispin=1, iorb=nothing, fmpi=true, link::Link=global_env)
    """
    Fits the Weiss field or Hybridization function with a discrete set of levels.
    """
    
    fit_single_normal_n3 = Libdl.dlsym(link.library, "chi2_fitgf_single_normal_n3")
    fit_single_normal_n5 = Libdl.dlsym(link.library, "chi2_fitgf_single_normal_n5")
    fit_single_superc_n5 = Libdl.dlsym(link.library, "chi2_fitgf_single_superc_n3")
    fit_single_superc_n5 = Libdl.dlsym(link.library, "chi2_fitgf_single_superc_n5")
    
    fit_lattice_normal_n3 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_normal_n3")
    fit_lattice_normal_n4 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_normal_n4")
    fit_lattice_normal_n6 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_normal_n6")
    fit_lattice_superc_n3 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_superc_n3")
    fit_lattice_superc_n4 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_superc_n4")
    fit_lattice_superc_n6 = Libdl.dlsym(link.library, "chi2_fitgf_lattice_superc_n6")

    
    
    iorb = isnothing(iorb) ? 0 : iorb

    if length(args) == 2  # normal case
        g = args[1]
        bath = args[2]
        bath_copy = copy(bath)
        dim_g = collect(Int, size(g))
        dim_bath = collect(Int, size(bath_copy))

        if length(dim_bath) == 1  # single impurity
            if length(dim_g) == 3
                ccall(fit_single_normal_n3, Cvoid,
                      (Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                      g, dim_g, bath_copy, dim_bath, ispin, iorb, fmpi)
            elseif length(dim_g) == 5
                ccall(fit_single_normal_n5, Cvoid,
                      (Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                      g, dim_g, bath_copy, dim_bath, ispin, iorb, fmpi)
            else
                error("chi_fitgf_normal: Expected dim(g) = 3 or 5")
            end
        elseif length(dim_bath) == 2               #lattice   
            if link.has_ineq
              if length(dim_g) == 3
                  ccall(fit_single_normal_n3, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, bath_copy, dim_bath, ispin, iorb, fmpi)
              elseif length(dim_g) == 4
                  ccall(fit_single_normal_n4, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, bath_copy, dim_bath, ispin, iorb, fmpi)
              elseif length(dim_g) == 6
                  ccall(fit_single_normal_n6, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, bath_copy, dim_bath, ispin, iorb, fmpi)
              else
                  error("chi_fitgf_normal: Expected dim(g) = 3, 4 or 6")
              end
            else
              error("Can't use r-DMFT routines without installing edipack2ineq")
            end
        else
            error("chi_fitgf_normal: Expected dim(bath) = 1 or 2")
        end
    elseif length(args) == 3 #superconductive case
        g = args[1]
        f = args[2]
        bath = args[3]
        bath_copy = copy(bath)
        dim_g = collect(Int, size(g))
        dim_f = collect(Int, size(f))
        dim_bath = collect(Int, size(bath_copy))

        if length(dim_bath) == 1  # single impurity
            if length(dim_g) == 3
                ccall(fit_single_normal_n3, Cvoid,
                      (Ptr{ComplexF64}, Ptr{Cint}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                      g, dim_g, f, dim_f, bath_copy, dim_bath, ispin, iorb, fmpi)
            elseif length(dim_g) == 5
                ccall(fit_single_normal_n5, Cvoid,
                      (Ptr{ComplexF64}, Ptr{Cint}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                      g, dim_g, f, dim_f, bath_copy, dim_bath, ispin, iorb, fmpi)
            else
                error("chi_fitgf_superc: Expected dim(g) = 3 or 5")
            end
        elseif length(dim_bath) == 2              #lattice   
            if link.has_ineq
              if length(dim_g) == 3
                  ccall(fit_single_normal_n3, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, f, dim_f, bath_copy, dim_bath, ispin, iorb, fmpi)
              elseif length(dim_g) == 4
                  ccall(fit_single_normal_n4, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, f, dim_f, bath_copy, dim_bath, ispin, iorb, fmpi)
              elseif length(dim_g) == 6
                  ccall(fit_single_normal_n6, Cvoid,
                        (Ptr{ComplexF64}, Ptr{Cint}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Bool),
                        g, dim_g, f, dim_f, bath_copy, dim_bath, ispin, iorb, fmpi)
              else
                  error("chi_fitgf_superc: Expected dim(g) = 3, 4 or 6")
              end
            else
              error("Can't use r-DMFT routines without installing edipack2ineq")
            end
        else
            error("chi_fitgf_superc: Expected dim(bath) = 1 or 2")
        end
    else
        error("chi_fitgf: Expected g,(f),bath as arguments")
    end
    return bath_copy
end

