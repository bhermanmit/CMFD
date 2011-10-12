program main

  use global
  use cmfd_execute,  only: compute_diffcoef,cmfd_solver 
  use cmfd_utils,    only: read_input
  use timing,        only: timer_start, timer_stop

  implicit none

  call timer_start(time_total)
  print *,"Reading input..."
  call read_input()
  print *,"Computing coupling coefficient"
  call compute_diffcoef()
  print *,"Solving linear system"
  call cmfd_solver()
  call timer_stop(time_total)
  print *,"Total Execution time (s):",time_total%elapsed

end program main
