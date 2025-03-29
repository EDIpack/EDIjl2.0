using Libdl
using Base.Threads


function set_hloc(link::Link, hloc::Array{ComplexF64}, Nlat::Union{Int, Nothing}=nothing)
    # Get function pointers
    ed_set_Hloc_single_N2 = Libdl.dlsym(link.library, "ed_set_Hloc_single_N2")
    ed_set_Hloc_single_N4 = Libdl.dlsym(link.library, "ed_set_Hloc_single_N4")
    
    if link.has_ineq
        ed_set_Hloc_lattice_N2 = Libdl.dlsym(link.library, "ed_set_Hloc_lattice_N2")
        ed_set_Hloc_lattice_N3 = Libdl.dlsym(link.library, "ed_set_Hloc_lattice_N3")
        ed_set_Hloc_lattice_N5 = Libdl.dlsym(link.library, "ed_set_Hloc_lattice_N5")
    end

    # Ensure the array is Fortran-ordered
    hloc = copy(hloc)  # Julia uses column-major order, so this should be fine
    dim_hloc = collect(Int, size(hloc))

    if Nlat !== nothing
        if link.has_ineq
            if length(dim_hloc) == 2
                ccall(ed_set_Hloc_lattice_N2, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), hloc, dim_hloc, Nlat)
            elseif length(dim_hloc) == 3
                ccall(ed_set_Hloc_lattice_N3, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), hloc, dim_hloc, Nlat)
            elseif length(dim_hloc) == 5
                ccall(ed_set_Hloc_lattice_N5, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), hloc, dim_hloc, Nlat)
            else
                error("ed_set_Hloc_lattice: dimension must be 2, 3, or 5")
            end
        else
            error("Can't use r-DMFT routines without installing edipack2ineq")
        end
    else
        if length(dim_hloc) == 2
            ccall(ed_set_Hloc_single_N2, Cvoid, (Ptr{ComplexF64}, Ptr{Int}), hloc, dim_hloc)
        elseif length(dim_hloc) == 4
            ccall(ed_set_Hloc_single_N4, Cvoid, (Ptr{ComplexF64}, Ptr{Int}), hloc, dim_hloc)
        else
            error("ed_set_Hloc_site: dimension must be 2 or 4")
        end
    end
end

function search_variable(link::Link, var::Float64, ntmp::Float64, converged::Bool)
    search_variable_wrap = Libdl.dlsym(link.library, :search_variable)

    var_arr = [var]
    ntmp_arr = [ntmp]
    conv_arr = [Int(converged)]  # Convert Bool to Int for C compatibility

    ccall(search_variable_wrap, Cvoid, 
          (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}), 
          var_arr, ntmp_arr, conv_arr)

    return var_arr[1], conv_arr[1] != 0  # Convert Int back to Bool
end


function check_convergence(link::Link, func::Array{ComplexF64}; threshold::Union{Nothing, Float64}=nothing, N1::Union{Nothing, Int}=nothing, N2::Union{Nothing, Int}=nothing)
    func = Array{ComplexF64}(func)  # Ensure it's an array
    err = 1.0
    conv_bool = false
    outfile = "error.err"

    # Retrieve default values if needed
    threshold = isnothing(threshold) ? unsafe_load(dlsym(link.library, :dmft_error)) : threshold
    N1 = isnothing(N1) ? unsafe_load(dlsym(link.library, :Nsuccess)) : N1
    N2 = isnothing(N2) ? unsafe_load(dlsym(link.library, :Nloop)) : N2

    # Initialize memory if first iteration
    if !haskey(link, :oldfunc)
        link.oldfunc = zeros(ComplexF64, size(func))
        link.whichiter = 0
        link.gooditer = 0
    end

    # Compute denominator and numerator in parallel
    denominator = Vector{Float64}(undef, size(func, 1))
    numerator = Vector{Float64}(undef, size(func, 1))

    @threads for i in 1:size(func, 1)
        denominator[i] = sum(abs, func[i, :, :, :])
        numerator[i] = sum(abs, func[i, :, :, :] .- link.oldfunc[i, :, :, :])
    end

    valid = denominator .!= 0
    errvec = real.(numerator[valid] ./ denominator[valid])

    # First iteration handling
    if link.whichiter == 0
        errvec .= 1.0
    end

    errvec = errvec[.!isnan.(errvec)]
    errmax = maximum(errvec)
    errmin = minimum(errvec)
    err = mean(errvec)
    link.oldfunc .= func

    if err < threshold
        link.gooditer += 1
    else
        link.gooditer = 0
    end

    link.whichiter += 1
    conv_bool = ((err < threshold) && (link.gooditer > N1) && (link.whichiter < N2)) || (link.whichiter >= N2)

    # Append errors to files using multithreading
    @threads for (filename, value) in [(outfile, err), (outfile * ".max", errmax), (outfile * ".min", errmin)]
        open(filename, "a") do file
            println(file, "$(link.whichiter) $(value)")
        end
    end

    if length(errvec) > 1
        open(outfile * ".distribution", "a") do file
            println(file, "$(link.whichiter) ", join(errvec, " "))
        end
    end

    # Print convergence message
    colorprefix = if conv_bool
        "\033[1;32m"
    elseif (err < threshold) && (link.gooditer <= N1)
        "\033[1;33m"
    else
        "\033[1;31m"
    end

    println(colorprefix * "Error: $(err) \033[0m")

    if link.whichiter >= N2
        println("Not converged after $N2 iterations.")
        open("ERROR.README", "a") do file
            println(file, "Not converged after $N2 iterations.")
        end
    end

    return err, conv_bool
end
