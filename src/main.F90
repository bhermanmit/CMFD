program main

  use global
  use cmfd_execute,  only: execute_cmfd 
  use timing,        only: timer_start, timer_stop

  implicit none

#include <finclude/petsc.h90>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>

  ! begin total timer
  call timer_start(time_total)

  ! execute cmfd
  call execute_cmfd()

  ! stop timer
  call timer_stop(time_total)

contains

!==============================================================================
! INITIALIZE
!==============================================================================

  subroutine initialize()

    ! initialize slepc/petsc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

    ! get mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD,n_procs,ierr)

    ! find master
    if (rank == 0) master = .true.

    ! print out number of procs
    if (master) print *,'Number of Processors:',n_procs

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
