using Libdl



function read_input(link::Link, input_string::String)
    """
    This function reads from the input file of EDIpack2. If the file does not
    exist, a template file is generated with default parameters.
    """
    
   # Load the function `read_input` from the library
    read_input_wrap = Libdl.dlsym(link.library, "read_input")
    
    if isnothing(read_input_wrap)
        error("read_input function not found in the shared library")
    end
    
    # Convert the input string to a C string (null-terminated)
    #c_string = Cstring(input_string)
    
    # Call the function
    ccall(read_input_wrap, Cvoid, (Cstring,), input_string)
end

