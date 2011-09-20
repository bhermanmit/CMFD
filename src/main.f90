program main

  use global
  use cmfd_execute,  only: allocate_cmfd
  use cmfd_utils,    only: read_input

  implicit none

  call read_input()
  call allocate_cmfd()

end program main
