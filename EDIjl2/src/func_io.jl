using Libdl
using Base.Threads


function get_gimp(; ilat::Union{Int, Nothing}=nothing, ishape::Union{Int, Nothing}=nothing, axis::String="m", typ::String="n",  zeta::Union{Array{ComplexF64}, Nothing}=nothing, link::Link=global_env)

    # Get function pointers
    ed_get_gimp_site_n3 = Libdl.dlsym(link.library, "get_gimp_site_n3")
    ed_get_gimp_site_n5 = Libdl.dlsym(link.library, "get_gimp_site_n5")
    
    if link.has_ineq
        ed_get_gimp_lattice_n3 = Libdl.dlsym(link.library, "get_gimp_lattice_n3")
        ed_get_gimp_lattice_n4 = Libdl.dlsym(link.library, "get_gimp_lattice_n4")
        ed_get_gimp_lattice_n6 = Libdl.dlsym(link.library, "get_gimp_lattice_n6")
    end

    nspin_aux = link.Nspin
    norb_aux = link.Norb

    if zeta !== nothing
        if isa(zeta, Number) 
            zeta = [zeta]
        end
        nfreq = Cint(size(zeta, 1))
        zflag = 1
        if any(abs(real(zeta)) .> 1e-10)
            axis = "r"
        end
    else
        zeta = [0.0]
        if axis == "m"
            nfreq = link.Lmats
        else
            nfreq = link.Lreal
        end
        zflag = 0  
    end

    if ishape == nothing
        ishape = link.dim_hloc + 1
    end

    # axis  
    if axis == "m"
        axisint = 0
    elseif axis == "r"
        axisint = 1
    else
        error("get_gimp: axis can only be 'm' or 'r'")
    end

    # typ
    if typ == "n"
        typint = 0
    elseif typ == "a"
        typint = 1
    else
        error("get_gimp: typ can only be 'n' or 'a'")
    end
    
    if link.Nineq == 0
        if ilat !== nothing
            error("ilat is not defined in single-impurity DMFT")
        end
        if ishape == 3
            gimp = zeros(Complex{Float64}, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
            ccall(ed_get_gimp_site_n3, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Ptr{ComplexF64},Cint, Cint), gimp, axisint, typint, zeta, nfreq, zflag)
        elseif ishape == 5
            gimp = zeros(Complex{Float64}, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
            ccall(ed_get_gimp_site_n5, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Ptr{ComplexF64},Cint, Cint), gimp, axisint, typint, zeta, nfreq, zflag)
        else
            error("Shape(array) != 3,5 in get_gimp_site")
        end
        return gimp
    else
        if link.has_ineq
            if ishape == 3
                gimp = zeros(Complex{Float64}, link.Nineq * nspin_aux * norb_aux, link.Nineq * nspin_aux * norb_aux, nfreq)
                ccall(ed_get_gimp_lattice_n3, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), gimp, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            elseif ishape == 4
                gimp = zeros(Complex{Float64}, link.Nineq, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
                ccall(ed_get_gimp_lattice_n4, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), gimp, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            elseif ishape == 6
                gimp = zeros(Complex{Float64}, link.Nineq, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
                ccall(ed_get_gimp_lattice_n6, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), gimp, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            else
                error("Shape(array) != 3,4,6 in get_gimp_lattice")
            end
            if ilat !== nothing && ishape != 3
                return gimp[ilat]
            else
                return gimp
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end

using Libdl
using Base.Threads


function get_sigma(; ilat::Union{Int, Nothing}=nothing, ishape::Union{Int, Nothing}=nothing, axis::String="m", typ::String="n",  zeta::Union{Array{ComplexF64}, Nothing}=nothing, link::Link=global_env)

    # Get function pointers
    ed_get_sigma_site_n3 = Libdl.dlsym(link.library, "get_sigma_site_n3")
    ed_get_sigma_site_n5 = Libdl.dlsym(link.library, "get_sigma_site_n5")
    
    if link.has_ineq
        ed_get_sigma_lattice_n3 = Libdl.dlsym(link.library, "get_sigma_lattice_n3")
        ed_get_sigma_lattice_n4 = Libdl.dlsym(link.library, "get_sigma_lattice_n4")
        ed_get_sigma_lattice_n6 = Libdl.dlsym(link.library, "get_sigma_lattice_n6")
    end

    nspin_aux = link.Nspin
    norb_aux = link.Norb

    if zeta !== nothing
        if isa(zeta, Number) 
            zeta = [zeta]
        end
        nfreq = Cint(size(zeta, 1))
        zflag = 1
        if any(abs(real(zeta)) .> 1e-10)
            axis = "r"
        end
    else
        zeta = [0.0]
        if axis == "m"
            nfreq = link.Lmats
        else
            nfreq = link.Lreal
        end
        zflag = 0  
    end

    if ishape == nothing
        ishape = link.dim_hloc + 1
    end

    # axis  
    if axis == "m"
        axisint = 0
    elseif axis == "r"
        axisint = 1
    else
        error("get_sigma: axis can only be 'm' or 'r'")
    end

    # typ
    if typ == "n"
        typint = 0
    elseif typ == "a"
        typint = 1
    else
        error("get_sigma: typ can only be 'n' or 'a'")
    end
    
    if link.Nineq == 0
        if ilat !== nothing
            error("ilat is not defined in single-impurity DMFT")
        end
        if ishape == 3
            sigma = zeros(Complex{Float64}, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
            ccall(ed_get_sigma_site_n3, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Ptr{ComplexF64},Cint, Cint), sigma, axisint, typint, zeta, nfreq, zflag)
        elseif ishape == 5
            sigma = zeros(Complex{Float64}, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
            ccall(ed_get_sigma_site_n5, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Ptr{ComplexF64},Cint, Cint), sigma, axisint, typint, zeta, nfreq, zflag)
        else
            error("Shape(array) != 3,5 in get_sigma_site")
        end
        return sigma
    else
        if link.has_ineq
            if ishape == 3
                sigma = zeros(Complex{Float64}, link.Nineq * nspin_aux * norb_aux, link.Nineq * nspin_aux * norb_aux, nfreq)
                ccall(ed_get_sigma_lattice_n3, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), sigma, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            elseif ishape == 4
                sigma = zeros(Complex{Float64}, link.Nineq, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
                ccall(ed_get_sigma_lattice_n4, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), sigma, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            elseif ishape == 6
                sigma = zeros(Complex{Float64}, link.Nineq, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
                ccall(ed_get_sigma_lattice_n6, Cvoid, (Ptr{ComplexF64}, Cint, Cint, Cint, Ptr{ComplexF64},Cint, Cint), sigma, link.Nineq, axisint, typint, zeta, nfreq, zflag)
            else
                error("Shape(array) != 3,4,6 in get_sigma_lattice")
            end
            if ilat !== nothing && ishape != 3
                return sigma[ilat]
            else
                return sigma
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end




function get_dens(; ilat::Union{Nothing, Int}=nothing, iorb::Union{Nothing, Int}=nothing, link::Link=global_env)
    """
    This function returns the value of the charge density.
    """

    aux_norb = link.Norb

    # Load the function ed_get_dens_n1
    ed_get_dens_n1_wrap = Libdl.dlsym(link.library, "ed_get_dens_n1")

    # If has_ineq, load ed_get_dens_n2
    ed_get_dens_n2_wrap = nothing
    if link.has_ineq
        ed_get_dens_n2_wrap = Libdl.dlsym(link.library, "ed_get_dens_n2")
        if isnothing(ed_get_dens_n2_wrap)
            error("ed_get_dens_n2 function not found in the shared library")
        end
    end

    if link.Nineq == 0
        densvec = zeros(Float64, aux_norb)
        ccall(ed_get_dens_n1_wrap, Cvoid, (Ptr{Float64},), densvec)

        if ilat !== nothing
            error("ilat cannot be specified for single-impurity DMFT")
        elseif iorb !== nothing
            return densvec[iorb]  # Julia is 1-based
        else
            return densvec
        end
    else
        if link.has_ineq
            densvec = zeros(Float64, link.Nineq, aux_norb)
            ccall(ed_get_dens_n2_wrap, Cvoid, (Ptr{Float64}, Cint), densvec, link.Nineq)

            if ilat !== nothing && iorb !== nothing
                return densvec[ilat, iorb]
            elseif ilat === nothing && iorb !== nothing
                return densvec[:, iorb]
            elseif ilat !== nothing && iorb === nothing
                return densvec[ilat, :]
            else
                return densvec
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end




function get_mag(; icomp::Union{Nothing, String}=nothing, ilat::Union{Nothing, Int}=nothing, iorb::Union{Nothing, Int}=nothing, link::Link=global_env)
    """
    This function returns the value of the magnetization.
    """

    # Map string component to integer
    if icomp !== nothing
        icomp = lowercase(icomp)
        if icomp == "x"
            icomp = 0
        elseif icomp == "y"
            icomp = 1
        elseif icomp == "z"
            icomp = 2
        else
            error("icomp must be \"x\", \"y\", or \"z\"")
        end
    end

    # Read Norb from the shared library
    aux_norb = link.Norb

    # Load ed_get_mag_n2
    ed_get_mag_n2_wrap = Libdl.dlsym(link.library, "ed_get_mag_n2")
    if isnothing(ed_get_mag_n2_wrap)
        error("ed_get_mag_n2 function not found in the shared library")
    end

    # If has_ineq, load ed_get_mag_n3
    ed_get_mag_n3_wrap = nothing
    if link.has_ineq
        ed_get_mag_n3_wrap = Libdl.dlsym(link.library, "ed_get_mag_n3")
        if isnothing(ed_get_mag_n3_wrap)
            error("ed_get_mag_n3 function not found in the shared library")
        end
    end

    if link.Nineq == 0
        magvec = zeros(Float64, 3, aux_norb)
        ccall(ed_get_mag_n2_wrap, Cvoid, (Ptr{Float64},), magvec)

        if ilat !== nothing
            error("ilat cannot be specified for single-impurity DMFT")
        elseif iorb !== nothing && icomp !== nothing
            return magvec[icomp, iorb]
        elseif iorb !== nothing && icomp === nothing
            return magvec[:, iorb]
        elseif iorb === nothing && icomp !== nothing
            return magvec[icomp, :]
        else
            return magvec
        end
    else
        if link.has_ineq
            magvec = zeros(Float64, link.Nineq, 3, aux_norb)
            ccall(ed_get_mag_n3_wrap, Cvoid, (Ptr{Float64}, Cint), magvec, link.Nineq)

            if ilat !== nothing
                if iorb !== nothing && icomp !== nothing
                    return magvec[ilat, icomp, iorb]
                elseif iorb === nothing && icomp !== nothing
                    return magvec[ilat, icomp, :]
                elseif iorb !== nothing && icomp === nothing
                    return magvec[ilat, :, iorb]
                elseif iorb === nothing && icomp === nothing
                    return magvec[ilat, :, :]
                end
            else
                if iorb !== nothing && icomp !== nothing
                    return magvec[:, icomp, iorb]
                elseif iorb === nothing && icomp !== nothing
                    return magvec[:, icomp, :]
                elseif iorb !== nothing && icomp === nothing
                    return magvec[:, :, iorb]
                elseif iorb === nothing && icomp === nothing
                    return magvec
                end
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end



function get_docc(; ilat::Union{Nothing, Int}=nothing, iorb::Union{Nothing, Int}=nothing, link::Link=global_env)
    """
    This function returns the value of the double occupation.
    """

    # Read Norb from the shared library
    aux_norb = link.Norb

    # Load ed_get_docc_n1
    ed_get_docc_n1_wrap = Libdl.dlsym(link.library, "ed_get_docc_n1")
    if isnothing(ed_get_docc_n1_wrap)
        error("ed_get_docc_n1 function not found in the shared library")
    end

    # If has_ineq, load ed_get_docc_n2
    ed_get_docc_n2_wrap = nothing
    if link.has_ineq
        ed_get_docc_n2_wrap = Libdl.dlsym(link.library, "ed_get_docc_n2")
        if isnothing(ed_get_docc_n2_wrap)
            error("ed_get_docc_n2 function not found in the shared library")
        end
    end

    if link.Nineq == 0
        doccvec = zeros(Float64, aux_norb)
        ccall(ed_get_docc_n1_wrap, Cvoid, (Ptr{Float64},), doccvec)

        if ilat !== nothing
            error("ilat cannot be specified for single-impurity DMFT")
        elseif iorb !== nothing
            return doccvec[iorb]
        else
            return doccvec
        end
    else
        if link.has_ineq
            doccvec = zeros(Float64, link.Nineq, aux_norb)
            ccall(ed_get_docc_n2_wrap, Cvoid, (Ptr{Float64}, Cint), doccvec, link.Nineq)

            if ilat !== nothing && iorb !== nothing
                return doccvec[ilat, iorb]
            elseif ilat === nothing && iorb !== nothing
                return doccvec[:, iorb]
            elseif ilat !== nothing && iorb === nothing
                return doccvec[ilat, :]
            else
                return doccvec
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end



function get_phi(; ilat::Union{Nothing, Int}=nothing, iorb::Union{Nothing, Int}=nothing, jorb::Union{Nothing, Int}=nothing, link::Link=global_env)
    """
    This function returns the value of the superconductive order parameter: \\phi = \\langle c_{\\uparrow} c_{\\downarrow} \\rangle.
    """

    # Read Norb from the shared library
    aux_norb = link.Norb

    # Load ed_get_phisc_n2
    ed_get_phisc_n2_wrap = Libdl.dlsym(link.library, "ed_get_phisc_n2")
    if isnothing(ed_get_phisc_n2_wrap)
        error("ed_get_phisc_n2 function not found in the shared library")
    end

    # If has_ineq, load ed_get_phisc_n3
    ed_get_phisc_n3_wrap = nothing
    if link.has_ineq
        ed_get_phisc_n3_wrap = Libdl.dlsym(link.library, "ed_get_phisc_n3")
        if isnothing(ed_get_phisc_n3_wrap)
            error("ed_get_phisc_n3 function not found in the shared library")
        end
    end

    if link.Nineq == 0
        phivec = zeros(Float64, aux_norb, aux_norb)
        ccall(ed_get_phisc_n2_wrap, Cvoid, (Ptr{Float64},), phivec)

        if ilat !== nothing
            error("ilat cannot be specified for single-impurity DMFT")
        elseif iorb !== nothing && jorb !== nothing
            return phivec[iorb, jorb]
        elseif iorb !== nothing && jorb === nothing
            return phivec[iorb, :]
        elseif jorb !== nothing && iorb === nothing
            return phivec[:, jorb]
        else
            return phivec
        end
    else
        if link.has_ineq
            phivec = zeros(Float64, link.Nineq, aux_norb, aux_norb)
            ccall(ed_get_phisc_n3_wrap, Cvoid, (Ptr{Float64}, Cint), phivec, link.Nineq)

            if ilat !== nothing
                if iorb !== nothing && jorb !== nothing
                    return phivec[ilat, iorb, jorb]
                elseif iorb !== nothing && jorb === nothing
                    return phivec[ilat, iorb, :]
                elseif jorb !== nothing && iorb === nothing
                    return phivec[ilat, :, jorb]
                else
                    return phivec[ilat, :, :]
                end
            else
                if iorb !== nothing && jorb !== nothing
                    return phivec[:, iorb, jorb]
                elseif iorb !== nothing && jorb === nothing
                    return phivec[:, iorb, :]
                elseif jorb !== nothing && iorb === nothing
                    return phivec[:, :, jorb]
                else
                    return phivec
                end
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end




function get_eimp(; ilat::Union{Nothing, Int}=nothing, ikind::Union{Nothing, Int}=nothing, link::Link=global_env)
    """
    This function returns the value of the local energy components.
    """

    # Load ed_get_eimp_n1
    ed_get_eimp_n1_wrap = Libdl.dlsym(link.library, "ed_get_eimp_n1")
    if isnothing(ed_get_eimp_n1_wrap)
        error("ed_get_eimp_n1 function not found in the shared library")
    end

    # If has_ineq, load ed_get_eimp_n2
    ed_get_eimp_n2_wrap = nothing
    if link.has_ineq
        ed_get_eimp_n2_wrap = Libdl.dlsym(link.library, "ed_get_eimp_n2")
        if isnothing(ed_get_eimp_n2_wrap)
            error("ed_get_eimp_n2 function not found in the shared library")
        end
    end

    if link.Nineq == 0
        eimp_vec = zeros(Float64, 4)
        ccall(ed_get_eimp_n1_wrap, Cvoid, (Ptr{Float64},), eimp_vec)

        if ilat !== nothing
            error("ilat cannot be specified for single-impurity DMFT")
        elseif ikind !== nothing
            return eimp_vec[ikind]
        else
            return eimp_vec
        end
    else
        if link.has_ineq
            eimp_vec = zeros(Float64, link.Nineq, 4)
            ccall(ed_get_eimp_n2_wrap, Cvoid, (Ptr{Float64}, Cint), eimp_vec, link.Nineq)

            if ilat !== nothing && ikind !== nothing
                return eimp_vec[ilat, ikind]
            elseif ilat === nothing && ikind !== nothing
                return eimp_vec[:, ikind]
            elseif ilat !== nothing && ikind === nothing
                return eimp_vec[ilat, :]
            else
                return eimp_vec
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end


function get_g0and(zeta::AbstractArray{Complex{T}}, bath::AbstractArray{T}; ishape::Union{Nothing, Int}=nothing, typ::String="n", link::Link=global_env) where T
    """
    This function calculates the value of the Anderson Impurity Model's noninteracting Green's function.
    """

    # Load functions dynamically
    ed_get_g0and_n3_wrap = Libdl.dlsym(link.library, "get_g0and_n3")
    if isnothing(ed_get_g0and_n3_wrap)
        error("get_g0and_n3 function not found in the shared library")
    end

    ed_get_g0and_n5_wrap = Libdl.dlsym(link.library, "get_g0and_n5")
    if isnothing(ed_get_g0and_n5_wrap)
        error("get_g0and_n5 function not found in the shared library")
    end

    # Get the auxiliary values for norb and nspin
    norb_aux = link.Norb
    nspin_aux = link.Nspin

    # Ensure zeta is of complex type
    zeta = Complex{T}[complex(v) for v in zeta]

    nfreq = size(zeta, 1)
    dimbath = size(bath, 1)

    # Determine the axis type (real or imaginary)
    axis = ""
    if any(abs(real(zeta)) .> 1e-10)
        axis = "r"
    elseif any(abs(imag(zeta)) .> 1e-10)
        axis = "m"
    else
        throw(ValueError("get_g0and: frequencies can only be purely real or purely imaginary"))
    end

    if isnothing(ishape)
        ishape = link.dim_hloc + 1
    end

    # Case for 3-dimensional shape
    if ishape == 3
        G0and = zeros(Complex{T}, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
        DimG0and = Int64[nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq]

        ccall(ed_get_g0and_n3_wrap, Cvoid, 
            (Ptr{Complex{T}}, Cint, Ptr{T}, Cint, Ptr{Complex{T}}, Ptr{Int64}, Cstring, Cstring),
            zeta, nfreq, bath, dimbath, G0and, DimG0and, axis, typ
        )
    # Case for 5-dimensional shape
    elseif ishape == 5
        G0and = zeros(Complex{T}, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
        DimG0and = Int64[nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq]

        ccall(ed_get_g0and_n5_wrap, Cvoid, 
            (Ptr{Complex{T}}, Cint, Ptr{T}, Cint, Ptr{Complex{T}}, Ptr{Int64}, Cstring, Cstring),
            zeta, nfreq, bath, dimbath, G0and, DimG0and, axis, typ
        )
    else
        throw(ValueError("Shape(array) != 3,5 in get_g0and"))
    end

    return G0and
end




function get_delta(zeta::AbstractArray{Complex{T}}, bath::AbstractArray{T}; ishape::Union{Nothing, Int}=nothing, typ::String="n", link::Link=global_env) where T
    """
    This function calculates the value of the Anderson Impurity Model's hybridization function.
    """

    # Load functions dynamically
    ed_get_delta_n3_wrap = Libdl.dlsym(link.library, "get_delta_n3")
    if isnothing(ed_get_delta_n3_wrap)
        error("get_delta_n3 function not found in the shared library")
    end

    ed_get_delta_n5_wrap = Libdl.dlsym(link.library, "get_delta_n5")
    if isnothing(ed_get_delta_n5_wrap)
        error("get_delta_n5 function not found in the shared library")
    end

    # Get the auxiliary values for norb and nspin
    norb_aux = link.Norb
    nspin_aux = link.Nspin

    # Ensure zeta is of complex type
    zeta = Complex{T}[complex(v) for v in zeta]

    nfreq = size(zeta, 1)
    dimbath = size(bath, 1)

    # Determine the axis type (real or imaginary)
    axis = ""
    if any(abs(real(zeta)) .> 1e-10)
        axis = "r"
    elseif any(abs(imag(zeta)) .> 1e-10)
        axis = "m"
    else
        throw(ValueError("get_delta: frequencies can only be purely real or purely imaginary"))
    end

    if isnothing(ishape)
        ishape = link.dim_hloc + 1
    end

    # Case for 3-dimensional shape
    if ishape == 3
        Delta = zeros(Complex{T}, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq)
        DimDelta = Int64[nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq]

        ccall(ed_get_delta_n3_wrap, Cvoid, 
            (Ptr{Complex{T}}, Cint, Ptr{T}, Cint, Ptr{Complex{T}}, Ptr{Int64}, Cstring, Cstring),
            zeta, nfreq, bath, dimbath, Delta, DimDelta, axis, typ
        )
    # Case for 5-dimensional shape
    elseif ishape == 5
        Delta = zeros(Complex{T}, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq)
        DimDelta = Int64[nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq]

        ccall(ed_get_delta_n5_wrap, Cvoid, 
            (Ptr{Complex{T}}, Cint, Ptr{T}, Cint, Ptr{Complex{T}}, Ptr{Int64}, Cstring, Cstring),
            zeta, nfreq, bath, dimbath, Delta, DimDelta, axis, typ
        )
    else
        throw(ValueError("Shape(array) != 3,5 in get_delta"))
    end

    return Delta
end


function get_chi(chan::String="spin", zeta::AbstractArray{Complex{T}}=nothing, axis::String=nothing, ilat::Union{Nothing, Int}=nothing, link::Link=global_env) where T
    """
    This function calculates the value of the Anderson Impurity Model's response function χ.
    """

    # Load functions dynamically
    ed_get_spinchi_wrap = Libdl.dlsym(link.library, "ed_get_spinchi")
    if isnothing(ed_get_spinchi_wrap)
        error("ed_get_spinchi function not found in the shared library")
    end

    ed_get_denschi_wrap = Libdl.dlsym(link.library, "ed_get_denschi")
    if isnothing(ed_get_denschi_wrap)
        error("ed_get_denschi function not found in the shared library")
    end

    ed_get_pairchi_wrap = Libdl.dlsym(link.library, "ed_get_pairchi")
    if isnothing(ed_get_pairchi_wrap)
        error("ed_get_pairchi function not found in the shared library")
    end

    ed_get_exctchi_wrap = Libdl.dlsym(link.library, "ed_get_exctchi")
    if isnothing(ed_get_exctchi_wrap)
        error("ed_get_exctchi function not found in the shared library")
    end

    # Get the auxiliary values for norb and nspin
    aux_norb = link.Norb
    aux_Lmats = link.Lmats
    aux_Lreal = link.Lreal
    aux_Ltau = link.Ltau
    edmode = link.get_ed_mode()

    if edmode != 1
        throw(ValueError("Susceptibility calculation not supported for ed_mode not = normal"))
    end

    # Initialize flags and variables
    zetaflag = 1
    Nsites = if link.Nineq == 0 1 else link.Nineq end
    latticeflag = if link.Nineq > 0 1 else 0 end

    if axis == nothing
        throw(ValueError("Axis is required"))
    end

    if zeta == nothing
        if axis == "m"
            zeta = Complex{T}[0.0]
            zetaflag = 0
            axisflag = 0
            nfreq = aux_Lmats
        elseif axis == "r"
            zeta = Complex{T}[0.0]
            zetaflag = 0
            axisflag = 1
            nfreq = aux_Lreal
        elseif axis == "t"
            zeta = Complex{T}[0.0]
            zetaflag = 0
            axisflag = 2
            nfreq = aux_Ltau
        end
    else
        axisflag = axis == "m" ? 0 : axis == "r" ? 1 : axis == "t" ? 2 : throw(ValueError("axis can only be m, r, or t"))
        nfreq = size(zeta, 1)
    end

    chi = Complex{T}[]

    if chan in ["spin", "s"]
        chi = zeros(Complex{T}, Nsites, aux_norb, aux_norb, nfreq)
        ccall(ed_get_spinchi_wrap, Cvoid,
            (Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Cint, Cint, Cint, Cint),
            chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag
        )
    elseif chan in ["dens", "d"]
        chi = zeros(Complex{T}, Nsites, aux_norb, aux_norb, nfreq)
        ccall(ed_get_denschi_wrap, Cvoid,
            (Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Cint, Cint, Cint, Cint),
            chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag
        )
    elseif chan in ["pair", "p"]
        chi = zeros(Complex{T}, Nsites, aux_norb, aux_norb, nfreq)
        ccall(ed_get_pairchi_wrap, Cvoid,
            (Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Cint, Cint, Cint, Cint),
            chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag
        )
    elseif chan in ["exct", "e"]
        chi = zeros(Complex{T}, Nsites, 3, aux_norb, aux_norb, nfreq)
        ccall(ed_get_exctchi_wrap, Cvoid,
            (Ptr{Complex{T}}, Ptr{Complex{T}}, Cint, Cint, Cint, Cint, Cint),
            chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag
        )
    else
        throw(ValueError("Invalid channel"))
    end

    # Return based on real-space DMFT flag
    if link.Nineq == 0
        return chi[1]
    elseif isnothing(ilat)
        return chi
    else
        return chi[ilat]
    end
end
