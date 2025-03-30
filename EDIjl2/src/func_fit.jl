using Libdl
using Base.Threads

function chi2_fitgf(link::Link, args...; ispin=1, iorb=nothing, fmpi=true)
    """
    Fits the Weiss field or Hybridization function with a discrete set of levels.
    """
    
    fit_single_normal_n3 = Libdl.dlsym(link.library, "chi2_fitgf_single_normal_n3")
    fit_single_normal_n5 =Libdl.dlsym(link.library, "chi2_fitgf_single_normal_n5")
    
    
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
        else
            error("chi_fitgf_normal: Expected dim(bath) = 1 or 2")
        end
    else
        error("chi_fitgf: Expected g,bath as arguments")
    end
    return bath_copy
end

