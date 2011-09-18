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
    PetscScalar :: sd

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

    print *, 'HERE'
    ! Begin loops to set values
    do i =1,nx
      do j=1,ny

        print *, cmfd%totalxs(i,j,1,1)
        sd = cmfd%totalxs(i,j,1,1)
        call MatSetValue(T,i-1,j-1,sd,INSERT_VALUES,ierr)

      end do
    end do

    call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)

    call MatView(T,PETSC_VIEWER_STDOUT_SELF,ierr)

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
    PetscScalar :: k_n        ! new k-eigenvalue
    PetscScalar :: k_o        ! old k-eigenvalue
    KSP         :: krylov     ! krylov solver
    PC          :: prec       ! preconditioner for krylov
 
    ! local variables
    PetscInt             :: n           ! dimensions of matrix
    PetscInt             :: ierr        ! error flag
    PetscScalar          :: guess=1.0   ! initial guess
    PetscScalar          :: ktol=1.e-7  ! krylov tolerance


    ! get indices for allocation
!    n = cmfd%indices(1)*cmfd%indices(2)*cmfd%indices(3)*cmfd%indices(4)
    n = 10
    ! set up loss matrix
    call MatCreate(PETSC_COMM_WORLD,M,ierr)
    call MatSetType(M,MATAIJ,ierr)
    call MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
    call MatSetFromOptions(M,ierr)

    ! set up production matrix
    call MatCreate(PETSC_COMM_WORLD,F,ierr)
    call MatSetType(F,MATAIJ,ierr)
    call MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
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

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_input()

    use xml_data_cmfd_t

    ! local variables
    integer :: nx   ! number of volumes in x-direction
    integer :: ny   ! number of volumes in y-direction
    integer :: nz   ! number of volumes in z-direction
    integer :: ng   ! number of energy groups
    integer :: i    ! x iteration counter
    integer :: j    ! y iteration counter
    integer :: k    ! z iteration counter
    integer :: g    ! g iteration counter
    integer :: nvec ! index for location in xml vector

    ! read xml input file
    call read_xml_file_cmfd_t('cmfd.xml')

    ! get mesh and group indices
    nx = geometry%nx
    ny = geometry%ny
    nz = geometry%nz
    ng = geometry%ng

    ! record indices in object
    cmfd%indices(1) = nx
    cmfd%indices(2) = ny
    cmfd%indices(3) = nz
    cmfd%indices(4) = ng 

    ! allocate totxs vector
    allocate(cmfd%totalxs(nx,ny,nz,ng))

    print *, mat(1)%totxs

    ! read in total cross section vector
    GROUP: do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! get vector index
            nvec = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(g-1)

            ! move input to cmfd data object
            cmfd%totalxs(i,j,k,g) = mat(1)%totxs(nvec)

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

    write(*,*) "Hello!"
    print *, cmfd%totalxs

  end subroutine read_input

end module cmfd_execute
