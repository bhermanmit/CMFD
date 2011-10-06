module global

  use cmfd_header, only: cmfd_obj
  use timing,      only: Timer

  implicit none
  save

  ! Main object
  type(cmfd_obj) :: cmfd

  ! Timing object
  type(Timer) :: time_mat    ! timer for mat building

end module global
