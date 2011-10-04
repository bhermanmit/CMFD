module cmfd_execute

  use global

  implicit none

contains

!===============================================================================
! ALLOCATE_CMFD allocates all of the space for the cmfd object based on tallies.
!===============================================================================

  subroutine allocate_cmfd()

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

    ! local variables
    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: nxyz(3,2)          ! single vector containing boundary locations
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: xyz_idx            ! index for determining if x,y or z leakage
    integer :: dir_idx            ! index for determining - or + face of cell
    integer :: shift_idx          ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)        ! spatial indices of neighbour
    integer :: bound(6)           ! vector containing indices for boudary check
    real(8) :: albedo(6)          ! albedo vector with global boundaries
    real(8) :: cell_totxs         ! total cross section of current ijk cell
    real(8) :: cell_dc            ! diffusion coef of current cell
    real(8) :: cell_hxyz(3)       ! cell dimensions of current ijk cell
    real(8) :: neig_totxs         ! total xs of neighbor cell
    real(8) :: neig_dc            ! diffusion coefficient of neighbor cell
    real(8) :: neig_hxyz(3)       ! cell dimensions of neighbor cell
    real(8) :: dtilda             ! finite difference coupling parameter 

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/) 

    ! get boundary condition information
    albedo = cmfd%albedo

    ! geting loop over group and spatial indices
    GROUP:  do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! get cell data
            cell_dc = cmfd%diffcof(i,j,k,g)
            cell_hxyz = cmfd%hxyz(i,j,k,:)


            ! setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dtilda
                dtilda = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) +  &
               &         (1-albedo(l))*cell_hxyz(xyz_idx))

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor cell data
                neig_dc = cmfd%diffcof(neig_idx(1),neig_idx(2),neig_idx(3),g)
                neig_hxyz = cmfd%hxyz(neig_idx(1),neig_idx(2),neig_idx(3),:)
  
                ! compute dtilda
                dtilda = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc +     &
               &          cell_hxyz(xyz_idx)*neig_dc)

              end if
  
              ! record dtilda in cmfd object
              cmfd%dtilda(i,j,k,g,l) = dtilda

            end do LEAK

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

  end subroutine compute_diffcoef

!===============================================================================
! COMPUTE_DHAT computes the nonlinear coupling coefficient
!===============================================================================

  subroutine compute_dhat()

   ! local variables
    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: nxyz(3,2)          ! single vector containing boundary locations
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: xyz_idx            ! index for determining if x,y or z leakage
    integer :: dir_idx            ! index for determining - or + face of cell
    integer :: shift_idx          ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)        ! spatial indices of neighbour
    integer :: bound(6)           ! vector containing indices for boudary check
    real(8) :: cell_dtilda(6)     ! cell dtilda for each face
    real(8) :: cell_flux          ! flux in current cell
    real(8) :: current(3,2)       ! cell current at each face
    real(8) :: neig_flux          ! flux in neighbor cell
    real(8) :: dhat               ! dhat equivalence parameter

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! geting loop over group and spatial indices
    GROUP:  do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! get cell data
            cell_dtilda = cmfd%dtilda(i,j,k,g,:)
            cell_flux = cmfd%flux(i,j,k,g)
            current(1,:) = cmfd%currentX(i-1:i,j,k,g)
            current(2,:) = cmfd%currentY(i,j-1:j,k,g)
            current(3,:) = cmfd%currentZ(i,j,k-1:k,g)


            ! setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dhat
                dhat = (current(xyz_idx,dir_idx) - shift_idx*cell_dtilda(l)*   &
               &        cell_flux)/cell_flux

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor cell data
                neig_flux = cmfd%flux(neig_idx(1),neig_idx(2),neig_idx(3),g)

                ! compute dhat 
                dhat = (current(xyz_idx,dir_idx) + shift_idx*cell_dtilda(l)*   &
               &       (neig_flux - cell_flux))/(neig_flux + cell_flux)

              end if

              ! record dtilda in cmfd object
              cmfd%dhat(i,j,k,g,l) = dhat

            end do LEAK

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

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

    Mat         :: M       ! loss matrix
    Mat         :: F       ! production matrix
    Vec         :: phi_n   ! new flux eigenvector
    Vec         :: phi_o   ! old flux eigenvector
    Vec         :: S_n     ! new source vector
    Vec         :: S_o     ! old source vector
    real(8)     :: k_n     ! new k-eigenvalue
    real(8)     :: k_o     ! old k-eigenvlaue
    real(8)     :: num     ! numerator for eigenvalue update
    real(8)     :: den     ! denominator for eigenvalue update
    real(8)     :: one=1.0 ! one
    integer     :: ierr    ! error flag
    KSP         :: krylov  ! krylov solver
    PC          :: prec    ! preconditioner for krylov
    PetscViewer :: viewer  ! viewer for answer

    integer :: i       ! iteration counter
    logical :: iconv   ! is problem converged

    ! reset convergence flag
    iconv = .FALSE.

    ! initialize PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    ! initialize matrices and vectors
    call init_data(M,F,phi_n,phi_o,S_n,S_o,k_n,k_o,krylov,prec)

    print *,"Setting up matrices"
    ! set up M loss matrix
    call loss_matrix(M)

    ! set up F production matrix
    call prod_matrix(F)

    print *,"Beginning power iteration"
    ! begin power iteration
    do i = 1,10000

      ! compute source vector
      call MatMult(F,phi_o,S_o,ierr)

      ! normalize source vector
      call VecScale(S_o,one/k_o,ierr)

      print *,"Performing Krylov solve"
      ! compute new flux vector
      call KSPSetOperators(krylov, M, M, SAME_NONZERO_PATTERN, ierr)
      call KSPSolve(krylov,S_o,phi_n,ierr)

      ! compute new source vector
      call MatMult(F,phi_n,S_n,ierr)

      ! compute new k-eigenvalue
      call VecSum(S_n,num,ierr)
      call VecSum(S_o,den,ierr)
      k_n = num/den

      ! renormalize the old source
      call VecScale(S_o,k_o,ierr)

      ! check convergence
      call convergence(phi_n,phi_o,S_n,S_o,k_o,k_n,iconv)

      ! to break or not to break
      if (iconv) exit

      ! record old values
      call VecCopy(phi_n,phi_o,ierr)
      k_o = k_n
      print *,"keff:",k_n
    end do

    ! compute source pdf and record in cmfd object
    call source_pdf(S_n)

    ! output answers
    print *,'keff:',k_o
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'fluxvec.bin',FILE_MODE_WRITE, &
                               viewer,ierr)
    call VecView(phi_n,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

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
    Vec         :: S_o        ! old source vector
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
    call VecSet(phi_n,guess,ierr)
    call VecSet(phi_o,guess,ierr)
    k_n = guess
    k_o = guess  

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_WORLD,krylov,ierr)
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

    use cmfd_utils,  only: get_matrix_idx

#include <finclude/petsc.h90>

    ! arguments
    Mat    :: M   ! loss matrix

    ! local variables
    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: nxyz(3,2)          ! single vector containing boundary locations
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: h                  ! energy group when doing scattering
    integer :: cell_mat_idx       ! matrix index of current cell
    integer :: neig_mat_idx       ! matrix index of neighbor cell
    integer :: scatt_mat_idx      ! matrix index for h-->g scattering terms
    integer :: bound(6)           ! vector for comparing when looking for boundaries
    integer :: xyz_idx            ! index for determining if x,y or z leakage
    integer :: dir_idx            ! index for determining - or + face of cell
    integer :: neig_idx(3)        ! spatial indices of neighbour
    integer :: shift_idx          ! parameter to shift index by +1 or -1
    integer :: ierr               ! Persc error code
    real(8) :: totxs              ! total macro cross section
    real(8) :: scattxsgg          ! scattering macro cross section g-->g
    real(8) :: scattxshg          ! scattering macro cross section h-->g
    real(8) :: dtilda(6)          ! finite difference coupling parameter
    real(8) :: dhat(6)            ! nonlinear coupling parameter
    real(8) :: hxyz(3)            ! cell lengths in each direction
    real(8) :: jn                 ! direction dependent leakage coeff to neig
    real(8) :: jo(6)              ! leakage coeff in front of cell flux
    real(8) :: jnet               ! net leakage from jo
    real(8) :: val                ! temporary variable before saving to matrix 
    PetscViewer :: viewer         ! viewer to write out matrix to binary file

    ! initialize matrix for building
    call MatAssemblyBegin(M,MAT_FLUSH_ASSEMBLY,ierr)

    ! get maximum indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! begin iteration loops
    GROUP:  do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            print *, "Setting up (i,j,k,g) cell:",i,j,k,g

            ! get matrix index of cell
            cell_mat_idx = get_matrix_idx(i,j,k,g,nx,ny,nz)

            ! retrieve cell data
            totxs = cmfd%totalxs(i,j,k,g)
            scattxsgg = cmfd%scattxs(i,j,k,g,g)
            dtilda = cmfd%dtilda(i,j,k,g,:)
            dhat = cmfd%dhat(i,j,k,g,:)
            hxyz = cmfd%hxyz(i,j,k,:)

            ! create boundary vector 
            bound = (/i,i,j,j,k,k/)

            ! begin loop over leakages
            ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z 
            LEAK: do l = 1,6

              ! define (x,y,z) and (-,+) indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2

              ! calculate spatial indices of neighbor
              neig_idx = (/i,j,k/)                ! begin with i,j,k
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
              neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx) 

              ! check for global boundary
              if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

                ! compute leakage coefficient for neighbor
                jn = -dtilda(l) + shift_idx*dhat(l)

                ! get neighbor matrix index
                neig_mat_idx = get_matrix_idx(neig_idx(1),neig_idx(2),         &
               &                              neig_idx(3),g,nx,ny,nz) 

                ! compute value to bank
                val = jn/hxyz(xyz_idx)

                ! record value in matrix
                call MatSetValue(M,cell_mat_idx-1,neig_mat_idx-1,val,          &
               &                 INSERT_VALUES,ierr)

              end if

              ! compute leakage coefficient for target
              jo(l) = shift_idx*dtilda(l) + dhat(l) 

            end do LEAK

            ! calate net leakage coefficient for target
            jnet = (jo(2) - jo(1))/hxyz(1) + (jo(4) - jo(3))/hxyz(2) +         &
           &       (jo(6) - jo(5))/hxyz(3) 

            ! calculate loss of neutrons
            val = jnet + totxs - scattxsgg
            print *,totxs-scattxsgg

            ! record diagonal term
            call MatSetValue(M,cell_mat_idx-1,cell_mat_idx-1,val,INSERT_VALUES,&
           &                 ierr)

            ! begin loop over off diagonal in-scattering
            SCATTR: do h = 1,ng

              ! cycle though if h=g
              if (h == g) then
                cycle
              end if

              ! get matrix index of in-scatter
              scatt_mat_idx = get_matrix_idx(i,j,k,h,nx,ny,nz)

              ! get scattering macro xs
              scattxshg = cmfd%scattxs(i,j,k,h,g)

              ! record value in matrix (negate it)
              val = -scattxshg
              call MatSetValue(M,cell_mat_idx-1,scatt_mat_idx-1,val,            &
             &                 INSERT_VALUES,ierr)

            end do SCATTR

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

    ! finalize matrix assembly
    call MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr)

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'lossmat.bin',FILE_MODE_WRITE, &
                               viewer,ierr)
    call MatView(M,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine loss_matrix

!===============================================================================
! PROD_MATRIX creates the matrix representing production of neutrons
!===============================================================================

  subroutine prod_matrix(F)

    use cmfd_utils,  only: get_matrix_idx

#include <finclude/petsc.h90>

    ! arguments
    Mat     :: F  ! production matrix

    ! local variables
    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: h                  ! energy group when doing scattering
    integer :: gmat_idx           ! index in matrix for energy group g
    integer :: hmat_idx           ! index in matrix for energy group h
    integer :: ierr               ! Petsc error code
    real(8) :: nfissxs            ! nufission cross section h-->g
    real(8) :: val                ! temporary variable for nfissxs
    PetscViewer :: viewer         ! viewer to print out matrix

    ! initialize matrix for building
    call MatAssemblyBegin(F,MAT_FLUSH_ASSEMBLY,ierr)

    ! get maximum indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! begin loop around energy groups and spatial indices
    GROUP: do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            NFISS: do h = 1,ng

              ! get cell data
              nfissxs = cmfd%nfissxs(i,j,k,h,g)

              ! get matrix location
              gmat_idx = get_matrix_idx(i,j,k,g,nx,ny,nz)
              hmat_idx = get_matrix_idx(i,j,k,h,nx,ny,nz)

              ! reocrd value in matrix
              val = nfissxs
              call MatSetValue(F,gmat_idx-1,hmat_idx-1,val,INSERT_VALUES,ierr)

            end do NFISS

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

    ! finalize matrix assembly
    call MatAssemblyEnd(F,MAT_FINAL_ASSEMBLY,ierr)

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'prodmat.bin',FILE_MODE_WRITE, &
                               viewer,ierr)
    call MatView(F,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine prod_matrix

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine convergence(phi_n,phi_o,S_n,S_o,k_o,k_n,iconv)

#include <finclude/petsc.h90>

    ! arguments
    Vec         :: phi_n  ! new flux eigenvector
    Vec         :: phi_o  ! old flux eigenvector
    Vec         :: S_n    ! new source vector
    Vec         :: S_o    ! old source vector
    real(8)     :: k_n    ! new k-eigenvalue
    real(8)     :: k_o    ! old k-eigenvalue
    logical     :: iconv  ! is the problem converged

    ! local variables
    Vec         :: phi_v          ! flux temp vector 
    Vec         :: S_v            ! source temp vector
    real(8)     :: ktol = 1.e-6   ! tolerance on keff
    real(8)     :: ftol = 1.e-4   ! tolerance on flux
    real(8)     :: stol = 1.e-5   ! tolerance on source
    real(8)     :: kerr           ! error in keff
    real(8)     :: ferr           ! error in flux
    real(8)     :: serr           ! error in source
    real(8)     :: one = -1.0     ! one
    integer     :: floc           ! location of max error in flux
    integer     :: sloc           ! location of max error in source
    integer     :: ierr           ! petsc error code
    integer     :: n              ! vector size

    ! reset convergence flag
    iconv = .FALSE.

    ! initialize temp vectors
    call VecGetLocalSize(phi_n,n,ierr)
    call VecCreate(PETSC_COMM_WORLD,phi_v,ierr)
    call VecSetSizes(phi_v,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi_v,ierr)
    call VecCreate(PETSC_COMM_WORLD,S_v,ierr)
    call VecSetSizes(S_v,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(S_v,ierr)

    ! calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! calculate max error in flux
    call VecWAXPY(phi_v,one,phi_n,phi_o,ierr)
    call VecPointwiseDivide(phi_v,phi_v,phi_n,ierr)
    call VecAbs(phi_v,ierr)
    call VecMax(phi_v,floc,ferr,ierr)

    ! calculate max error in source
    call VecWAXPY(S_v,one,S_n,S_o,ierr)
    call VecPointwiseDivide(S_v,S_v,S_n,ierr)
    call VecAbs(S_v,ierr)
    call VecMax(S_v,sloc,serr,ierr)

    ! check for convergence
    if(kerr < ktol .and. ferr < ftol .and. serr < stol) iconv = .TRUE.

    ! destroy vectors
    call VecDestroy(phi_v,ierr)
    call VecDestroy(S_v,ierr)
 
  end subroutine convergence

!===============================================================================
! SOURCE_PDF calculates the probability distribution of the cmfd fission source
!===============================================================================

  subroutine source_pdf(source)

    use cmfd_utils,  only: get_matrix_idx

#include <finclude/petsc.h90>

    ! arguments
    Vec         :: source     ! new source vector
 
    ! local variables
    integer :: nx                          ! maximum number of cells in x direction
    integer :: ny                          ! maximum number of cells in y direction
    integer :: nz                          ! maximum number of cells in z direction
    integer :: ng                          ! maximum number of energy groups
    integer :: i                           ! iteration counter for x
    integer :: j                           ! iteration counter for y
    integer :: k                           ! iteration counter for z
    integer :: g                           ! iteration counter for groups
    integer :: ierr                        ! PETSC error code
    integer :: idx                         ! index in vector
    real(8) :: total                       ! sum of source vector
    real(8) :: hxyz(3)                     ! cell dimensions of current ijk cell
    PetscScalar, pointer :: source_ptr(:)  ! pointer to petsc vector for fortran

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! get source vector from petsc
    call VecGetArrayF90(source,source_ptr,ierr)

    ! loop around indices to map to cmfd object 
    GROUP:  do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! get dimensions of cell
            hxyz = cmfd%hxyz(i,j,k,:)

            ! get index
            idx = get_matrix_idx(i,j,k,g,nx,ny,nz)

            ! multiply source density by volume and record in object
            cmfd%sourcepdf(i,j,k,g) = source_ptr(idx)*hxyz(1)*hxyz(2)*hxyz(3)
            
          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

    ! normalize source such that it sums to 1.0
    cmfd%sourcepdf = cmfd%sourcepdf/sum(cmfd%sourcepdf)

    ! restore petsc vector
    call VecRestoreArrayF90(source,source_ptr,ierr)

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
