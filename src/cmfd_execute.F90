module cmfd_execute

  use global

  implicit none

contains

!===============================================================================
! ALLOCATE_CMFD allocates all of the space for the cmfd object based on tallies.
!===============================================================================

  subroutine allocate_cmfd()

#include <finclude/petsc.h90>

    Mat         :: T      ! test matrix
    integer     :: i
    integer     :: j 
    integer     :: k
    integer     :: g
    integer     :: ierr
    integer     :: nx
    integer     :: ny
    PetscViewer viewer

    ! dimensions
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)

    ! initialize PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    ! set up loss matrix
    call MatCreate(PETSC_COMM_WORLD,T,ierr)
    call MatSetType(T,MATAIJ,ierr)
    call MatSetSizes(T,PETSC_DECIDE,PETSC_DECIDE,nx,ny,ierr)
    call MatSetFromOptions(T,ierr)

    ! begin matrix assembly
    call MatAssemblyBegin(T,MAT_FLUSH_ASSEMBLY,ierr)

    ! Begin loops to set values
    do i =1,nx
      do j=1,ny

        call MatSetValue(T,i-1,j-1,cmfd%totalxs(i,j,1,1),INSERT_VALUES,ierr)

      end do
    end do

    call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'matrixmat',FILE_MODE_WRITE,viewer,ierr)

    call MatView(T,viewer,ierr)

    call PetscViewerDestroy(viewer,ierr)

    call MatDestroy(T,ierr)

    call PetscFinalize(ierr)

  end subroutine allocate_cmfd

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections.
!===============================================================================

  subroutine compute_xs()

  end subroutine compute_xs

!===============================================================================
! COMPUTE_DIFFCOEF computes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_diffcoef()

  end subroutine compute_diffcoef

!===============================================================================
! COMPUTE_DHAT computes the nonlinear coupling coefficient
!===============================================================================

  subroutine compute_dhat()

  end subroutine compute_dhat

!===============================================================================
! ACCUMULATE accumulates cycle to cycle diffusion parameters
!===============================================================================

  subroutine accumulate()

  end subroutine accumulate

!===============================================================================
! CMFD_SOLVER in the main power iteration routine for the cmfd calculation
!===============================================================================

  subroutine cmfd_solver()

#include <finclude/petsc.h90>

    Mat         :: M      ! loss matrix
    Mat         :: F      ! production matrix
    Vec         :: phi_n  ! new flux eigenvector
    Vec         :: phi_o  ! old flux eigenvector
    Vec         :: S_n    ! new source vector
    Vec         :: S_o    ! old source vector
    PetscScalar :: k_n    ! new k-eigenvalue
    PetscScalar :: k_o    ! old k-eigenvlaue
    PetscScalar :: num    ! numerator for eigenvalue update
    PetscScalar :: den    ! denominator for eigenvalue update
    PetscInt    :: ierr   ! error flag
    KSP         :: krylov ! krylov solver
    PC          :: prec   ! preconditioner for krylov

    integer :: i       ! iteration counter
    logical :: iconv   ! is problem converged

    ! reset convergence flag
    iconv = .FALSE.

    ! initialize PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    ! initialize matrices and vectors
    call init_data(M,F,phi_n,phi_o,S_n,S_o,k_n,k_o,krylov,prec)

    ! set up M loss matrix
    call loss_matrix(M)

    ! set up F production matrix
    call prod_matrix(F)

    ! begin power iteration
    do i = 1,10000

      ! compute source vector
      call MatMult(F,phi_o,S_o)

      ! normalize source vector
      call VecScale(S_o,1.0/k_o)

      ! compute new flux vector
      call KSPSetOperators(krylov, M, M, SAME_NONZERO_PATTERN, ierr)
      call KSPSolve(krylov,S_o,phi_n,ierr)

      ! compute new source vector
      call MatMult(F,phi_n,S_n)

      ! compute new k-eigenvalue
      call VecSum(S_n,num)
      call VecSum(S_o,den)
      k_n = num/den

      ! check convergence
      call convergence(phi_n,phi_o,S_n,S_o,k_o,k_n,iconv)

      ! to break or not to break
      if (iconv) exit

      ! record old values
      call VecCopy(phi_n,phi_o)
      k_o = k_n

    end do

    ! finalize PETSc
    call PetscFinalize(ierr)

  end subroutine cmfd_solver

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data(M,F,phi_n,phi_o,S_n,S_o,k_n,k_o,krylov,prec)

#include <finclude/petsc.h90>

    ! arguments
    Mat         :: M          ! loss matrix
    Mat         :: F          ! production matrix
    Vec         :: phi_n      ! new flux eigenvector
    Vec         :: phi_o      ! old flux eigenvector
    Vec         :: S_n        ! new source vector
    Vec         :: S_o        ! old seource vector
    real(8)     :: k_n        ! new k-eigenvalue
    real(8)     :: k_o        ! old k-eigenvalue
    KSP         :: krylov     ! krylov solver
    PC          :: prec       ! preconditioner for krylov
 
    ! local variables
    integer             :: n           ! dimensions of matrix
    integer             :: nz_M        ! number of nonzeros in loss matrix
    integer             :: nz_F        ! number of nonzeros in prod matrix
    integer             :: ierr        ! error flag
    integer             :: nx          ! maximum number of x cells
    integer             :: ny          ! maximum number of y cells
    integer             :: nz          ! maximum number of z cells
    integer             :: ng          ! maximum number of groups
    integer             :: n_corner    ! number of corner cells
    integer             :: n_edge      ! number of edge cells
    integer             :: n_side      ! number of side cells
    integer             :: n_int       ! number of interior cells 
    real(8)             :: guess=1.0   ! initial guess
    real(8)          :: ktol=1.e-7  ! krylov tolerance


    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! calculate dimensions of matrix
    n = nx*ny*nz*ng

    ! calculate # of types of cells
    n_corner = 8
    n_edge = 4*(nx + ny + nz) - 24
    n_side = 2*(nx*ny + ny*nz + nz*nx) - 8*(nx + ny + nz) + 24
    n_int = nx*ny*nz - n_side - n_edge - n_corner

    ! calculate number of non-zeros in each matrix
    nz_M = ng*(nx*ny*nz + n_int*6 + n_side*5 + n_edge*4 + n_corner*3) +        &
   &       (ng**2 -ng)*nx*ny*nz
    nz_F = ng**2*nx*ny*nz

    ! set up loss matrix
    call MatCreate(PETSC_COMM_WORLD,M,ierr)
    call MatSetType(M,MATAIJ,ierr)
    call MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,nz_M,ierr)
    call MatSetFromOptions(M,ierr)

    ! set up production matrix
    call MatCreate(PETSC_COMM_WORLD,F,ierr)
    call MatSetType(F,MATAIJ,ierr)
    call MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,nz_F,ierr)
    call MatSetFromOptions(F,ierr)

    ! set up flux vectors
    call VecCreate(PETSC_COMM_WORLD,phi_n,ierr)
    call VecSetSizes(phi_n,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi_n,ierr)
    call VecCreate(PETSC_COMM_WORLD,phi_o,ierr)
    call VecSetSizes(phi_o,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi_o,ierr)

    ! set up source vectors
    call VecCreate(PETSC_COMM_WORLD,S_n,ierr)
    call VecSetSizes(S_n,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(S_n,ierr)
    call VecCreate(PETSC_COMM_WORLD,S_o,ierr)
    call VecSetSizes(S_o,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(S_o,ierr)

    ! set initial guess
    call VecSet(phi_n,guess)
    call VecSet(phi_o,guess)
    k_n = guess
    k_o = guess  

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_WORLD,krylov)
    call KSPSetTolerances(krylov,ktol,PETSC_DEFAULT_DOUBLE_PRECISION,          &
   &                      PETSC_DEFAULT_DOUBLE_PRECISION,                      &
   &                      PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetType(krylov,KSPGMRES,ierr)
    call KSPSetFromOptions(krylov,ierr)
    call KSPGetPC(krylov,Prec,ierr)
    call PCSetType(Prec,PCILU,ierr)

  end subroutine init_data 

!===============================================================================
! LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine loss_matrix(M)

    ! arguments
    Mat    :: M   ! loss matrix

    ! local variables

  end subroutine loss_matrix

!===============================================================================
! PROD_MATRIX creates the matrix representing production of neutrons
!===============================================================================

  subroutine prod_matrix(F)

    ! arguments
    Mat     :: F  ! production matrix

    ! local variables

  end subroutine prod_matrix

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine convergence(phi_n,phi_o,S_n,S_o,k_o,k_n,iconv)

    ! arguments
    Vec         :: phi_n  ! new flux eigenvector
    Vec         :: phi_o  ! old flux eigenvector
    Vec         :: S_n    ! new source vector
    Vec         :: S_o    ! old source vector
    PetscScalar :: k_n    ! new k-eigenvalue
    PetscScalar :: k_o    ! old k-eigenvalue
    logical     :: iconv  ! is the problem converged

    ! local variables
 
  end subroutine convergence

!===============================================================================
! SOURCE_PDF calculates the probability distribution of the cmfd fission source
!===============================================================================

  subroutine source_pdf()

  end subroutine source_pdf

!===============================================================================
! COUNT_SOURCE determines the number of source sites in each mesh box
!===============================================================================

  subroutine count_source()

  end subroutine count_source

!===============================================================================
! WEIGHT_FACTORS calculates the weight adjustment factors for next MC cycle
!===============================================================================

  subroutine weight_factors()

  end subroutine weight_factors

!===============================================================================
! ADJUST_WEIGHT adjusts the initial weight of the particle
!===============================================================================

  subroutine adjust_weight()

! should this do all source particles at once or called each time before a
! neutron is born in the MC cycle?  MC21 is the latter

  end subroutine adjust_weight

end module cmfd_execute
