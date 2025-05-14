using Libdl

function reset_umatrix(link::Link)
    """
    This function resets to 0 all the interaction coefficients.
    """

    # Load reset_umatrix function dynamically
    reset_umatrix_wrap = Libdl.dlsym(link.library, "reset_umatrix")
    if isnothing(reset_umatrix_wrap)
        error("reset_umatrix function not found in the shared library")
    end

    # Call the reset function
    ccall(reset_umatrix_wrap, Cvoid, ())
end


function add_twobody_operator(oi::Int, si::String, oj::Int, sj::String, ok::Int, sk::String, ol::Int, sl::String, Uijkl::Float64, link::Link)
    """
    This function lets the user add an interaction term on-the-fly.
    The input parameters are the spin and orbital indices of the second
    quantized operators and the interaction coefficient.
    """

    # Load add_twobody_operator function dynamically
    add_twobody_operator_wrap = Libdl.dlsym(link.library, "add_twobody_operator")
    if isnothing(add_twobody_operator_wrap)
        error("add_twobody_operator function not found in the shared library")
    end

    # Map spin values from strings to integers
    mapping = Dict("u" => 1, "d" => 2)
    spinvector = [mapping[si], mapping[sj], mapping[sk], mapping[sl]]

    # Adjust orbital indices by adding 1
    orbvector = [oi, oj, ok, ol]

    # Call the function
    ccall(add_twobody_operator_wrap, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cdouble),
        orbvector[1], spinvector[1], orbvector[2], spinvector[2],
        orbvector[3], spinvector[3], orbvector[4], spinvector[4], Uijkl
    )
end
