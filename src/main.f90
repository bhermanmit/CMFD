program main

  use global
  use cmfd_execute,  only: compute_diffcoef 
  use cmfd_utils,    only: read_input

  implicit none

  call read_input()
  call compute_diffcoef

end program main
