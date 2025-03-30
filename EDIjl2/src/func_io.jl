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
