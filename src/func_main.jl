using Libdl



function init_solver(link::Link; bath::Union{Nothing, Array{Float64, 2}, Array{Float64, 1}}=nothing, Nb::Union{Nothing, Int}=nothing, Nlat::Union{Nothing, Int}=nothing)
    
    init_solver_site = Libdl.dlsym(link.library, "init_solver_site")
    init_solver_ineq = Libdl.dlsym(link.library, "init_solver_ineq")

    if isnothing(bath)
        if isnothing(Nb) && isnothing(Nlat)
            Nb = ccall(dlsym(link.library, :get_bath_dimension), Int, ())
            bath = zeros(Float64, Nb)
        elseif isnothing(Nb) && !isnothing(Nlat)
            if link.has_ineq
                Nb = ccall(dlsym(link.library, :get_bath_dimension), Int, ())
                bath = zeros(Float64, Nlat, Nb)
            else
                error("Can't use r-DMFT routines without installing edipack2ineq")
            end
        elseif !isnothing(Nb) && isnothing(Nlat)
            bath = zeros(Float64, Nb)
        elseif !isnothing(Nb) && !isnothing(Nlat)
            if link.has_ineq
                bath = zeros(Float64, Nlat, Nb)
            else
                error("Can't use r-DMFT routines without installing edipack2ineq")
            end
        end
    else
        if !isnothing(Nb) || !isnothing(Nlat)
            println("INIT_SOLVER WARNING: Bath vector provided, Nb and/or Nlat are discarded")
        end
    end
    
    dim_bath = collect(Cint, size(bath))
    
    if length(dim_bath) < 2
        ccall(init_solver_site, Cvoid, (Ptr{Float64}, Ptr{Cint}), bath, dim_bath)
        link.Nineq = Cint(0)
    else
        if link.has_ineq
            ccall(init_solver_ineq, Cvoid, (Ptr{Float64}, Ptr{Cint}), bath, dim_bath)
            link.Nineq = Cint(size(bath, 1))
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end

    return bath
end


function solve(link::Link, bath::Array{Float64}, flag_gf::Bool=true, flag_mpi::Bool=true, mpi_lanc::Bool=false)
    solve_site = Libdl.dlsym(link.library, "solve_site")
    solve_ineq = Libdl.dlsym(link.library, "solve_ineq")

    dim_bath = collect(Cint, size(bath))
    
    if length(dim_bath) < 2
        ccall(solve_site, Nothing,
              (Ptr{Float64}, Ptr{Int}, Cint, Cint),
              bath, dim_bath, flag_gf, flag_mpi)
    else
        if has_ineq
            ccall(solve_ineq, Nothing,
                  (Ptr{Float64}, Ptr{Int}, Cint, Cint),
                  bath, dim_bath, flag_gf, mpi_lanc)
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    end
end

function finalize_solver(link::Link)
    """
    Cleans up the ED environment, deallocates relevant arrays, 
    and allows for another call to `init_solver`.
    """
    
    finalize_solver = Libdl.dlsym(link.library, "finalize_solver")
    
    if link.Nineq === nothing
        println("ED environment is not initialized yet")
        return
    end

    ccall(finalize_solver, Cvoid, (Cint,), link.Nineq)
    
    link.Nineq = nothing
    link.dim_hloc = 0
    link.Nsym = nothing

    if hasproperty(link, :oldfunc)
        delete!(link, :oldfunc)
    end
    if hasproperty(link, :gooditer)
        delete!(link, :gooditer)
    end
    if hasproperty(link, :whichiter)
        delete!(link, :whichiter)
    end

    println("ED environment finalized")
end
