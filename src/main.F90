!==============================================================================!
!> @mainpage CMFD: Coarse Mesh Finite Difference Diffusion
!>
!> @section Overview
!>
!> This program solve the neutron diffusion equation in three-dimensions and
!> arbitrary number of energy groups. CMFD utilizing a number of external
!> packages which must be downloaded and configured before compiling. A list
!> of them is as follows:
!>  - PETSc
!>  - MPICH2
!>  - Fortran BLAS/LAPACK
!>  - HYPRE
!>  - SLEPc
!>
!> The packages PETSc/MPICH2/BLAS/LAPACK/HYPRE can be downloaded from 
!>  http://www.mcs.anl.gov/petsc/
!>
!> The package SLEPc can be downloaded from http://www.grycap.upv.es/slepc/
!>
!> Successful configuration has been completed with Intel v.12.0 compiler.
!> GNU has not been working when configuring HYPRE
!>
!> @section Building PETSc
!>
!> PETSc has been successfully configured and built with the following command:
!>
!> @verbatim
!>   ./configure --with-cc=icc --with-fc=ifort --download-mpich=1 
!> --download-f-blas-lapack=1 --download-hypre=1 --with-debugging=<0 or 1>
!> @endverbatim
!>
!> @section Compiling
!> 
!> Before compiling, the following environmental variables need to be set:
!>  - PETSC_DIR => path to PETSc build directory
!>  - PETSC_ARCH => name of PETSc build to use (located in petsc build dir)
!>  - SLEPC_DIR => path to SLEPc build directory
!>
!> Compiling is as easy as running the Makefile with:
!>
!> @verbatim
!>   make cmfd
!> @endverbatim
!>
!> @section Running
!>
!> To run CMFD code, execute the following:
!>
!> @verbatim
!>   $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -np <nprocs> ~/path-to-cmfd/cmfd [options]
!> @endverbatim
!>
!> <nprocs> is the number of processors to run with MPI
!>
!> Both SLEPc/PETSc and CMFD command line options can be set for [options].
!> Please refer to SLEPC/PETSC manuals for package specific command line options or
!> just type -help.  CMFD specific command line options are as follows:
!>   - --solver_type = <solver> where <solver> is <b>power</b>, <b>slepc</b>, <b>snes</b>
!>     - <b>power</b> uses the manual power iteration routine
!>     - <b>slepc</b> uses slepc to perform the eigenvalue calculation
!>     - <b>snes</b> uses JFNK nonlinear method to solve eigenvalue calculation
!>
!==============================================================================!

program main

  use constants
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

    integer :: i                    ! loop counter
    integer :: argc                 ! number of command line arguments
    character(len=25) :: argv       ! a command line argument
    character(len=10) :: today_date
    character(len=8)  :: today_time

     !========================================
     ! INITIALIZE MPI

     ! initialize slepc/petsc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

    ! get mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD,n_procs,ierr)

    ! find master
    if (rank == 0) master = .true.

   
    !=========================================
    ! WRITE HEADING INFO

    if (master) then

      ! write header
      write(*, FMT='(/11(A/))') &
    & '    ,o888888o.           ,8.       ,8.          8 8888888888   8 888888888o.      ', &
    & '   8888     `88.        ,888.     ,888.         8 8888         8 8888    `^888.   ', &
    & ',8 8888       `8.      .`8888.   .`8888.        8 8888         8 8888        `88. ', &
    & '88 8888               ,8.`8888. ,8.`8888.       8 8888         8 8888         `88 ', &
    & '88 8888              ,8"8.`8888,8^8.`8888.      8 888888888888 8 8888          88 ', &
    & '88 8888             ,8" `8.`8888" `8.`8888.     8 8888         8 8888          88 ', &
    & '88 8888            ,8"   `8.`88"   `8.`8888.    8 8888         8 8888         ,88 ', &
    & '`8 8888       .8" ,8"     `8.`"     `8.`8888.   8 8888         8 8888        ,88" ', &
    & '   8888     ,88" ,8"       `8        `8.`8888.  8 8888         8 8888    ,o88P"   ', &
    & '    `8888888P"  ,8"         `         `8.`8888. 8 8888         8 888888888P"      ', &
    & '__________________________________________________________________________________'

      ! Write version information
      write(*, FMT=*) &
           '     Developed At:  Massachusetts Institute of Technology'
      write(*, FMT='(6X,"Version:",7X,I1,".",I1,".",I1)') &
           VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE

      ! Write the date and time
      call get_today(today_date, today_time)
      write(*, FMT='(6X,"Date/Time:",5X,A,1X,A)') &
           trim(today_date), trim(today_time)

      ! Write information on number of processors
      write(*, FMT='(1X,A,I0)') '     MPI Processes: ',n_procs

    end if

    !=======================================
    ! PROCESS COMMAND LINE ARGS

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
      write(*, FMT='(/,/,A)') 'RESULTS'
      write(*, FMT='(A)')     '**********************'
      write(*, FMT='(/,"Total Execution time (s): ",F0.4)') time_total%elapsed
      write(*, FMT='("K-effective: ",F0.8,/)') cmfd%keff
    end if

  end subroutine finalize

!===============================================================================
! GET_TODAY determines the date and time at which the program began execution
! and returns it in a readable format
!===============================================================================

  subroutine get_today(today_date, today_time)

    character(10), intent(out) :: today_date
    character(8),  intent(out) :: today_time
    
    integer       :: val(8)
    character(8)  :: date_
    character(10) :: time_
    character(5)  :: zone
    
    call date_and_time(date_, time_, zone, val)
    ! val(1) = year (YYYY)
    ! val(2) = month (MM)
    ! val(3) = day (DD)
    ! val(4) = timezone
    ! val(5) = hours (HH)
    ! val(6) = minutes (MM)
    ! val(7) = seconds (SS)
    ! val(8) = milliseconds

    if (val(2) < 10) then
       if (val(3) < 10) then
          today_date = date_(6:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(6:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    else
       if (val(3) < 10) then
          today_date = date_(5:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(5:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    end if
    today_time = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end subroutine get_today

end program main
