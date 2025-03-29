using Libdl



function get_bath_dimension(link::Link)
    """
    This function returns the correct dimension for the bath to be allocated \
    (for each impurity) given the parameters of the system.
    """
    # Call the C function using ccall
    return ccall(dlsym(link.library, :get_bath_dimension), Cint, ())
end


