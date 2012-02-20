module cmfd_data

  implicit none
  private
  public :: set_up_cmfd

contains

!==============================================================================
! SET_UP_CMFD
!===============================================================================

  subroutine set_up_cmfd()

    use cmfd_input

    ! read the input file
    call read_input()

    ! compute dtilde terms
    call compute_dtilde()

  end subroutine set_up_cmfd

!===============================================================================
! COMPUTE_DTILDE computes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_dtilde()

    use global, only: cmfd

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
    real(8) :: dtilde             ! finite difference coupling parameter 
    real(8) :: ref_albedo         ! albedo to reflector

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
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUP: do g = 1,ng

            ! get cell data
            cell_dc = cmfd%diffcof(g,i,j,k)
            cell_hxyz = cmfd%hxyz(:,i,j,k)

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

                ! compute dtilde
                dtilde = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) +  &
               &         (1-albedo(l))*cell_hxyz(xyz_idx))

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor cell data
                neig_dc = cmfd%diffcof(g,neig_idx(1),neig_idx(2),neig_idx(3))
                neig_hxyz = cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3))

                ! compute dtilde
                dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc +     &
               &         cell_hxyz(xyz_idx)*neig_dc)

              end if

              ! record dtilde in cmfd object
              cmfd%dtilde(l,g,i,j,k) = dtilde

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_dtilde

end module cmfd_data 
