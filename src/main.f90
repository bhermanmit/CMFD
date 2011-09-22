program main

  use global
  use cmfd_execute,  only: compute_diffcoef,cmfd_solver 
  use cmfd_utils,    only: read_input

  implicit none

  call read_input()
  call compute_diffcoef()
  call cmfd_solver()

end program main
