using Libdl



function get_bath_dimension(;link::Link=global_env)
    """
    This function returns the correct dimension for the bath to be allocated \
    (for each impurity) given the parameters of the system.
    """

    if get_bath_type(link) > 2  # replica/general
        if self.Nsym === nothing
            error("get_bath_dimension: no replica/general matrix is initialized")
        else
            bbathdim = ccall(dlsym(link.library, :get_bath_dimension_symmetries), Cint, (Ptr{Int},), self.Nsym)
        end
    else
        bathdim = ccall(dlsym(link.library, :get_bath_dimension_direct), Cint, ())
    end
    return bathdim
   
end






function set_hreplica(hvec::Array{ComplexF64}; lambdavec::Array{Float64}, link::Link=global_env)
    # Get function pointers
    init_hreplica_symmetries_d5 = Libdl.dlsym(link.library, "init_Hreplica_symmetries_d5")
    init_hreplica_symmetries_d3 = Libdl.dlsym(link.library, "init_Hreplica_symmetries_d3")
    
    if link.has_ineq
      init_hreplica_symmetries_lattice_d5 = Libdl.dlsym(link.library, "init_Hreplica_symmetries_lattice_d5")
      init_hreplica_symmetries_lattice_d3 = Libdl.dlsym(link.library, "init_Hreplica_symmetries_lattice_d3")
    end

    # Ensure the array is Fortran-ordered
    aux_norb =  link.norb
    aux_nspin = link.nspin
    dim_hvec = size(hvec)
    dim_lambdavec = size(lambdavec)
    self.Nsym = dim_lambdavec[1]

    if len(dim_hvec) == 3
        if len(dim_lambdavec) == 2
          ccall(init_hreplica_symmetries_d3, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
        elseif length(dim_lambdavec) == 3
          if link.has_ineq
            ccall(init_hreplica_symmetries_lattice_d3, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
          else
            error("Can't use r-DMFT routines without installing edipack2ineq")  
          end
        else
          error("Size(lambdavec) != 2 or 3  in set_Hreplica")
        end
    else
        if len(dim_lambdavec) == 2
          ccall(init_hreplica_symmetries_d5, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
        elseif length(dim_lambdavec) == 3
          if link.has_ineq
            ccall(init_hreplica_symmetries_lattice_d5, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
          else
            error("Can't use r-DMFT routines without installing edipack2ineq")  
          end
        else
          error("Size(lambdavec) != 2 or 3  in set_Hreplica")
        end
    end
end



function set_hgeneral(hvec::Array{ComplexF64}; lambdavec::Array{Float64}, link::Link=global_env)
    # Get function pointers
    init_hgeneral_symmetries_d5 = Libdl.dlsym(link.library, "init_Hgeneral_symmetries_d5")
    init_hgeneral_symmetries_d3 = Libdl.dlsym(link.library, "init_Hgeneral_symmetries_d3")
    
    if link.has_ineq
      init_hgeneral_symmetries_lattice_d5 = Libdl.dlsym(link.library, "init_Hgeneral_symmetries_lattice_d5")
      init_hgeneral_symmetries_lattice_d3 = Libdl.dlsym(link.library, "init_Hgeneral_symmetries_lattice_d3")
    end

    # Ensure the array is Fortran-ordered
    aux_norb =  link.norb
    aux_nspin = link.nspin
    dim_hvec = size(hvec)
    dim_lambdavec = size(lambdavec)
    self.Nsym = dim_lambdavec[1]

    if len(dim_hvec) == 3
        if len(dim_lambdavec) == 2
          ccall(init_hgeneral_symmetries_d3, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
        elseif length(dim_lambdavec) == 3
          if link.has_ineq
            ccall(init_hgeneral_symmetries_lattice_d3, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
          else
            error("Can't use r-DMFT routines without installing edipack2ineq")  
          end
        else
          error("Size(lambdavec) != 2 or 3  in set_Hgeneral")
        end
    else
        if len(dim_lambdavec) == 2
          ccall(init_hgeneral_symmetries_d5, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
        elseif length(dim_lambdavec) == 3
          if link.has_ineq
            ccall(init_hgeneral_symmetries_lattice_d5, Cvoid, (Ptr{ComplexF64}, Ptr{Int},Ptr{ComplexF64}, Ptr{Int}), hvec, dim_hvec, lambdavec, dim_lambdavec)
          else
            error("Can't use r-DMFT routines without installing edipack2ineq")  
          end
        else
          error("Size(lambdavec) != 2 or 3  in set_Hgeneral")
        end
    end
end


function break_symmetry_bath(bath::Array{Float64}, field::Float64, sign::Union{Float64, Array{Float64}}, save::Bool=true, link::Link=global_env)

  # Get function pointers
  break_symmetry_bath_site = Libdl.dlsym(link.library, "break_symmetry_bath_site")

  if link.has_ineq
    break_symmetry_bath_ineq = Libdl.dlsym(link.library, "break_symmetry_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(break_symmetry_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Float64, Float64, Cint), bath, bath_shape, field, sign, save_int)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          if sign isa Float64
            sign = sign * ones(bath_shape[1])
          end
          ccall(break_symmetry_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Float64, Ptr{Float64}, Cint), bath, bath_shape, field, sign, save_int)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end






function spin_symmetrize_bath(bath::Array{Float64}, save::Bool=true, link::Link=global_env)

  # Get function pointers
  spin_symmetrize_bath_site = Libdl.dlsym(link.library, "spin_symmetrize_bath_site")

  if link.has_ineq
    spin_symmetrize_bath_ineq = Libdl.dlsym(link.library, "spin_symmetrize_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(spin_symmetrize_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), bath, bath_shape, save_int)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          ccall(spin_symmetrize_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), bath, bath_shape, save_int)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end




function orb_symmetrize_bath(bath::Array{Float64}, orb1::Cint, orb2::Cint, save::Bool=true, link::Link=global_env)

  # Get function pointers
  orb_symmetrize_bath_site = Libdl.dlsym(link.library, "orb_symmetrize_bath_site")

  if link.has_ineq
    orb_symmetrize_bath_ineq = Libdl.dlsym(link.library, "orb_symmetrize_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(orb_symmetrize_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint, Cint, Cint), bath, bath_shape, orb1, orb2, save_int)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          ccall(orb_symmetrize_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint, Cint, Cint), bath, bath_shape, orb1, orb2, save_int)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end






function orb_equality_bath(bath::Array{Float64}, orb1::Cint, orb2::Cint, save::Bool=true, link::Link=global_env)

  # Get function pointers
  orb_equality_bath_site = Libdl.dlsym(link.library, "orb_equality_bath_site")

  if link.has_ineq
    orb_equality_bath_ineq = Libdl.dlsym(link.library, "orb_equality_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(orb_equality_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint, Cint), bath, bath_shape, indx, save_int)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          ccall(orb_equality_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint, Cint), bath, bath_shape, indx, save_int)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end



function ph_symmetrize_bath(bath::Array{Float64}, save::Bool=true, link::Link=global_env)

  # Get function pointers
  ph_symmetrize_bath_site = Libdl.dlsym(link.library, "ph_symmetrize_bath_site")

  if link.has_ineq
    ph_symmetrize_bath_ineq = Libdl.dlsym(link.library, "orb_equality_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(ph_symmetrize_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), bath, bath_shape, save_int)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          ccall(ph_symmetrize_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}, Cint), bath, bath_shape, save_int)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end


function save_array_as_bath(bath::Array{Float64}, link::Link=global_env)

  # Get function pointers
  save_array_as_bath_site = Libdl.dlsym(link.library, "save_array_as_bath_site")

  if link.has_ineq
    save_array_as_bath_ineq = Libdl.dlsym(link.library, "save_array_as_bath_ineq")
  end
  
  if save
      save_int = 1
  else
      save_int = 0
  end

  bath_shape = size(bath)

  if len(bath_shape) == 1
      if len(dim_lambdavec) == 2
        ccall(save_array_as_bath_site, Cvoid, (Ptr{ComplexF64}, Ptr{Int}), bath, bath_shape)
      elseif length(dim_lambdavec) == 3
        if link.has_ineq
          ccall(save_array_as_bath_ineq, Cvoid, (Ptr{ComplexF64}, Ptr{Int}), bath, bath_shape)
        else
          error("Can't use r-DMFT routines without installing edipack2ineq")  
        end
      else
        error("Size(bath) != 1 or 2  in set_Hgeneral")
      end
  end

end
