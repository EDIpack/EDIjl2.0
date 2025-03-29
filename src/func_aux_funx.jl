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


function check_convergence(arr1::AbstractArray, arr2::AbstractArray; threshold::Float64=0.00001)
    # Ensure both arrays have the same shape
    if size(arr1) != size(arr2)
        throw(ArgumentError("Arrays must have the same size"))
    end

    # Get the number of dimensions (ndim) of the arrays
    dims = ndims(arr1)

    # Check if both arrays have the same number of dimensions
    if dims != ndims(arr2)
        throw(ArgumentError("Arrays must have the same number of dimensions"))
    end

    # Compute the absolute difference between arr1 and arr2
    denominator = sum(abs.(arr1), dims=dims)
    diff = abs.(arr1 .- arr2)

    # Sum along the last dimension to reduce to a (d-1)-dimensional array
    reduced_diff = sum(diff, dims=dims)  # This reduces along the last dimension
    result = maximum(reduced_diff ./ denominator)

    println( "Error = ",result)

    # Check if the result is less than the threshold
    is_below_threshold = all(result .< threshold)

    return result, is_below_threshold
end
