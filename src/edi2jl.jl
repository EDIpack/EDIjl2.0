using Libdl
using DelimitedFiles

mutable struct Link
    library::Ptr{Cvoid}
    has_ineq::Union{Bool, Nothing}
    Nineq::Union{Cint, Nothing}
    dim_hloc::Cint
end

const var_types = Dict(
    "Uloc" => [Cdouble,5],    # Array of float
    "beta" => [Cdouble,1],    # Scalar float
    "Norb" => [Cint,1],       # Scalar int
    "Nspin" => [Cint,1],      # Scalar int
    "Lmats" => [Cint,1],      # Scalar int
    "Lreal" => [Cint,1],      # Scalar int
    "dmft_error" => [Cint,1]  
)


function InitLink(libpath::String)
    lib = try
        Libdl.dlopen(libpath)
    catch e
        println("Cannot init Link class: invalid library: ", e)
        return Link(lib, nothing, nothing ,0)
    end

    has_ineq_ptr = try
        Libdl.dlsym(lib, "has_ineq")
    catch
        nothing
    end

    has_ineq = isnothing(has_ineq_ptr) ? nothing : Bool(unsafe_load(Ptr{Cint}(has_ineq_ptr)))
    return Link(lib, has_ineq, nothing ,0)
end

# Function to get a variable's pointer from the shared library
function get_variable_ptr(obj::Link, varname::String)
    try
        ptr = Libdl.dlsym(obj.library, varname)
        if ptr == nothing
            error("Symbol $varname not found in library!")
        end
        return ptr
    catch
        error("Variable $varname not found in library")
    end
end

# General get_variable function to handle any type based on the dictionary
function get_variable(obj::Link, name::Symbol)
    varname = String(name)
    
    # Look up the type from the dictionary
    if haskey(var_types, varname)
        ctype = var_types[varname][1]
        dimension = var_types[varname][2]
    else
        error("Unknown variable type for $varname")
    end

    ptr = get_variable_ptr(obj, varname)
    if ptr == nothing
        error("Failed to retrieve pointer for $varname")
    end

    # Retrieve the variable value based on its type
    if dimension == 1
      return unsafe_load(Ptr{ctype}(ptr))
    else
      return unsafe_wrap(Array, Base.unsafe_convert(Ptr{ctype}, ptr), dimension)
    end
end

# General set_variable function to handle any type based on the dictionary
function set_variable(obj::Link, name::Symbol, val)
    varname = String(name)
    
    # Look up the type from the dictionary
    if haskey(var_types, varname)
        ctype = var_types[varname][1]
        dimension = var_types[varname][2]
    else
        error("Unknown variable type for $varname")
    end

    ptr = get_variable_ptr(obj, varname)
    if ptr == nothing
        error("Failed to retrieve pointer for $varname")
    end

    # Set the variable value based on its type
    if dimension == 1
      unsafe_store!(Ptr{ctype}(ptr), ctype(val))
    else  # Array
      for (i, comp) in enumerate(val)
          if i <= dimension
              unsafe_store!(Ptr{ctype}(ptr) + (i - 1)*sizeof(ctype), ctype(comp))
          end        
      end
    end
end


# Accessing properties based on the dictionary
function Base.getproperty(obj::Link, name::Symbol)
    if name in fieldnames(Link)
        return getfield(obj, name)
    else
        # Use the dictionary to handle types dynamically
        return get_variable(obj, name)
    end
end

function Base.setproperty!(obj::Link, name::Symbol, val)
    if name in fieldnames(Link)
        setfield!(obj, name, val)
    else
        # Use the dictionary to handle types dynamically
        set_variable(obj, name, val)
    end
end

function get_bath_type(link::Link)::Int
    get_bath_type_ptr = try
        Libdl.dlsym(link.library, "get_bath_type")
    catch
        error("Function get_bath_type not found")
    end
    return ccall(get_bath_type_ptr, Cint, ())
end

function get_ed_mode(link::Link)::Int
    get_ed_mode_ptr = Libdl.dlsym(link.library, "get_ed_mode", false)
    isnothing(get_ed_mode_ptr) && error("Function get_ed_mode not found")
    return ccall(get_ed_mode_ptr, Cint, ())
end

# Load shared library
function find_library()
    libext = Sys.isapple() ? ".dylib" : ".so"
    libname = "libedipack2_cbinding" * libext
    search_paths = [get(ENV, "EDIPACK_PATH", ""), get(ENV, "LD_LIBRARY_PATH", ""), get(ENV, "DYLD_LIBRARY_PATH", "")]
    for path in split(join(search_paths, ":"), ":")
        libpath = joinpath(path, libname)
        if isfile(libpath)
            return libpath
        end
    end
    error("Library loading failed. Check installation of edipack2")
end

# Initialize global environment
libpath = find_library()
global_env = InitLink(libpath)

include("func_read_input.jl")
include("func_aux_funx.jl")
include("func_main.jl")
include("func_bath.jl")
include("func_io.jl")
include("func_fit.jl")

function dens_bethe(x, d)
    root = sqrt(Complex(1 - (x / d)^2))
    dens_bethe = (2 / (pi * d)) * root
    return real(dens_bethe)
end

# Script
wmixing = 0.3
wband = 1.0
read_input(global_env, "inputED.conf")


#####################################################

Eband = range(-wband, wband, length=1000)
de = step(Eband)
Dband = dens_bethe.(Eband, wband)
wm = (pi / global_env.beta) .* (2 .* (0:(global_env.Lmats - 1)) .+ 1)

println(global_env.Nspin)


hloc = zeros(ComplexF64, 1, 1, 1, 1)
Gmats = zeros(ComplexF64, 1, 1, 1, 1,global_env.Lmats)
Delta = zeros(ComplexF64, 1, 1, 1, 1,global_env.Lmats)

set_hloc(global_env, hloc)
bath = init_solver(global_env)

for iloop in 0:100

  converged = false  
  
  bath_old = copy(bath)
  
  solve(global_env, bath)
  gimp = get_gimp(global_env; axis="m")
  smats = get_sigma(global_env; axis="m")

  zeta = wm * im .- smats[1, 1, 1, 1, :]

  Gmats[1,1,1,1,:] .= sum(Dband' ./ (zeta .- Eband'), dims=2) .* de
  Delta[1,1,1,1,:] .= 0.25 * wband * Gmats[1,1,1,1,:]

  writedlm("Delta_iw.dat", [wm imag.(Delta[1,1,1,1,:]) real.(Delta[1,1,1,1,:])])

  bath_old = copy(bath)
  global bath = chi2_fitgf(global_env,Delta,bath_old)

  result, converged = check_convergence(Gmats)
  global bath = 0.3.*bath + (1.0-0.3).*bath_old
  
  
  if converged
    break
  end
  
end

