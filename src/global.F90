module global

  use cmfd_header, only: cmfd_obj
  use timing,      only: Timer

  implicit none
  save

  ! Main object
  type(cmfd_obj) :: cmfd

  ! Timing objects
  type(Timer) :: time_total  ! timer for whole calculation
  type(Timer) :: time_mat    ! timer for mat building
  type(Timer) :: time_power  ! timer for power iteration

  ! petsc error code
  integer :: ierr

  ! mpi parametesr
  logical :: master = .false. ! am i master
  integer :: rank             ! rank of processor
  integer :: n_procs          ! number of processors

  ! solver type
  character(len=25) :: solver_type

end module global
