program main

  use global
  use cmfd_execute,  only: compute_diffcoef,cmfd_solver 
  use cmfd_utils,    only: read_input

  implicit none

  print *,"Reading input..."
  call read_input()
  print *,"Computing coupling coefficient"
  call compute_diffcoef()
  print *,"Solving linear system"
  call cmfd_solver()

end program main
