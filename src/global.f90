module global

  use cmfd_header, only: cmfd_obj

  implicit none
  save

  ! Main object
  type(cmfd_obj) :: cmfd

end module global
