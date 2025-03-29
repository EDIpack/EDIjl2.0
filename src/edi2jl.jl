using Libdl

mutable struct Link
    library::Ptr{Cvoid}
    has_ineq::Union{Bool, Nothing}
    Nineq::Union{Cint, Nothing}
    dim_hloc::Cint
end

function Link(libpath::String)
    lib = try
        Libdl.dlopen(libpath)
    catch e
        println("Cannot init Link class: invalid library: ", e)
        return Link(C_NULL, nothing)
    end

    has_ineq_ptr = try
        Libdl.dlsym(lib, "has_ineq")
    catch
        nothing
    end

    has_ineq = isnothing(has_ineq_ptr) ? nothing : Bool(unsafe_load(Ptr{Cint}(has_ineq_ptr)))
    return Link(lib, has_ineq, nothing ,0)
end

function get_variable_ptr(obj::Link, varname::String)
    try
        return Libdl.dlsym(obj.library, varname)
    catch
        error("Variable $varname not found in library")
    end
end

function get_variable(obj::Link, name::Symbol, ctype)
    ptr = get_variable_ptr(obj, String(name))
    return unsafe_load(Ptr{ctype}(ptr))
end

function set_variable(obj::Link, name::Symbol, ctype, val)
    ptr = get_variable_ptr(obj, String(name))
    unsafe_store!(Ptr{ctype}(ptr), val)
end

function Base.getproperty(obj::Link, name::Symbol)
    if name in fieldnames(Link)
        return getfield(obj, name)
    else
        return get_variable(obj, name, Cint)  # Change `Cint` to the correct type as needed
    end
end

function Base.setproperty!(obj::Link, name::Symbol, val)
    if name in fieldnames(Link)
        setfield!(obj, name, val)
    else
        set_variable(obj, name, Cint, val)  # Change `Cint` to the correct type as needed
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
global_env = Link(libpath)

include("func_read_input.jl")
include("func_aux_funx.jl")
include("func_main.jl")
include("func_bath.jl")
include("func_io.jl")



println(global_env.Nspin)

read_input(global_env, "inputED.conf")
hloc = rand(ComplexF64, 1, 1, 1, 1)
println(global_env.Nspin)
global_env.Nspin=1
global_env.Norb=1
println(global_env.Nspin)

set_hloc(global_env, hloc)
bath = init_solver(global_env)
solve(global_env, bath)
g = get_gimp(global_env; axis="m")
smats = get_sigma(global_env; axis="m")


