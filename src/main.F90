program main

  use global
  use cmfd_execute,  only: execute_cmfd 
  use timing,        only: timer_start, timer_stop

  implicit none

#include <finclude/petsc.h90>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>

  ! initialize
  call initialize()

  ! begin total timer
  call timer_start(time_total)

  ! execute cmfd
  call execute_cmfd()

  ! stop timer
  call timer_stop(time_total)

  ! finalize run
  call finalize()

contains

!==============================================================================
! INITIALIZE
!==============================================================================

  subroutine initialize()

    integer :: i                ! loop counter
    integer :: argc             ! number of command line arguments
    character(len=25) :: argv   ! a command line argument

    ! initialize slepc/petsc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

    ! get mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD,n_procs,ierr)

    ! find master
    if (rank == 0) master = .true.

    ! print out number of procs
    if (master) print *,'Number of Processors:',n_procs

    ! get number of command line arguments
    argc = COMMAND_ARGUMENT_COUNT()

    ! begin loop around command line args
    do i = 1,argc

      ! get that argument
      call GET_COMMAND_ARGUMENT(i,argv)

      ! begin case structure
      select case(trim(argv))

        ! solver
        case('--solver_type')

          ! get next argument
          call GET_COMMAND_ARGUMENT(i+1,argv)

          ! set global var
          solver_type = trim(argv)

        ! do nothing here
        case DEFAULT

      end select

    end do

  end subroutine initialize

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize slepc/petsc
    call SlepcFinalize(ierr)

    ! print output
    if (master) then
      print *,"Total Execution time (s):",time_total%elapsed
      print *,"K-effective:",cmfd%keff
    end if

  end subroutine finalize

end program main
