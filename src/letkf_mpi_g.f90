module letkf_mpi_g
  use mpi
  use global
  
  implicit none
  
  integer, save :: pe_root, pe_rank, pe_size
  integer, save :: mpi_comm_letkf

  logical, save :: pe_isroot

  integer, save, allocatable :: ens_list(:)
  !! A list of ensemble numbers (from 0 to MEM) that this
  !! process is responsible for the I/O
  integer, save :: ens_count

  integer, save, allocatable :: ens_map(:)
  !! for each index MEM, the number of the PE
  !! respsonsible for its I/O
  
  
  integer, save, allocatable :: ij_list(:)
  !! A list of grid point locations that this process
  !! is responsible for doing core LETKF algorithm on
  integer, save :: ij_count


  real, allocatable, save :: load_weights(:)
  
  integer, allocatable, save :: scatterv_count(:)
  integer, allocatable, save :: scatterv_displ(:)


  logical, save :: autotune
  character(len=:), allocatable, save :: autotune_outfile
end module letkf_mpi_g
