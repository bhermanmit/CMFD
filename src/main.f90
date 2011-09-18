program main

  use global
  use cmfd_execute,  only: allocate_cmfd,read_input 

  implicit none

  call read_input()
  call allocate_cmfd()

end program main
