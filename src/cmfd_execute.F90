!==============================================================================!
! MODULE: cmfd_execute
!
!> @author Bryan Herman
!>
!> @brief Routine for running the eigenvalue solve
!==============================================================================!

module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use global
  use cmfd_power_solver, only: cmfd_power_execute
  use cmfd_slepc_solver, only: cmfd_slepc_execute
  use cmfd_snes_solver,  only: cmfd_snes_execute
 
  implicit none

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

    ! intialize the cmfd object 
    call set_up_cmfd()

    ! execute solver
    select case (trim(solver_type))

      case('power')
        call cmfd_power_execute()
      case('slepc')
        call cmfd_slepc_execute()
      case('snes')
        call cmfd_snes_execute()
      case DEFAULT
        call cmfd_power_execute()

    end select

    ! write vtk file
    !if(.not. cmfd_only) call write_cmfd_vtk()

  end subroutine execute_cmfd

end module cmfd_execute
