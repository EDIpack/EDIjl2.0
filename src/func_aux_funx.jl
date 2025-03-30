using Libdl
using Base.Threads
using Printf

# ANSI color codes
const BOLD = "\033[1m"
const GREEN = "\033[92m"
const YELLOW = "\033[93m"
const RED = "\033[91m"
const COLOREND = "\033[0m"


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
    link.dim_hloc = Cint(length(dim_hloc))

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


using Printf

# Define ANSI color codes
const BOLD = "\033[1m"
const GREEN = "\033[92m"
const YELLOW = "\033[93m"
const RED = "\033[91m"
const COLOREND = "\033[0m"


function check_convergence(arr1::AbstractArray; threshold::Float64=0.00001, N1::Int=2, N2::Int=100)
    global arr_old, gooditer, whichiter

    # If arr_old does not exist, initialize it
    if !@isdefined(arr_old)
        arr_old = copy(arr1)
        gooditer = 0
        whichiter = 0
        err = 1.0
        errmin = 1.0
        errmax = 1.0
        conv_bool = false
    else
        # Ensure both arrays have the same shape
        if size(arr1) != size(arr_old)
            throw(ArgumentError("Arrays must have the same size"))
        end

        # Compute the absolute difference between arr1 and arr_old
        denominator = sum(abs.(arr1), dims=ndims(arr1))
        diff = abs.(arr1 .- arr_old)

        # Sum along the last dimension
        reduced_diff = sum(diff, dims=ndims(arr1))
        
        # Compute the average error using Base functions
        err = sum(reduced_diff ./ denominator) / length(reduced_diff)  
        errmin = minimum(reduced_diff ./ denominator)
        errmax = maximum(reduced_diff ./ denominator)

        # Check convergence
        if err < threshold
            gooditer += 1
        else
            gooditer = 0
        end
        whichiter += 1
    end

    # Determine color based on error
    conv_bool = err < threshold && (gooditer >= N1)
    
    colorprefix = if conv_bool 
        BOLD * GREEN
    elseif (err < threshold) && (gooditer < N1)
        BOLD * YELLOW
    else
        BOLD * RED
    end

    # Print convergence message
    println(colorprefix * "Max Error = " * COLOREND * @sprintf("%.6e", errmax))
    println(colorprefix * "Avg Error = " * COLOREND * @sprintf("%.6e", err))
    println(colorprefix * "Min Error = " * COLOREND * @sprintf("%.6e", errmin))

    if whichiter >= N2
        println(colorprefix * "Not converged after $N2 iterations." * COLOREND)
        open("ERROR.README", "a") do file
            write(file, "Not converged after $N2 iterations.\n")
        end
    end

    println()  # Add a blank line for readability

    # Update arr_old for the next iteration
    arr_old .= arr1  # Update in-place to avoid reallocation

    return err, conv_bool
end

