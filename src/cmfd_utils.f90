module cmfd_utils

  use global

 implicit none

contains

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

  end subroutine read_input

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(i,j,k,g,nx,ny,nz)

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: nx              ! maximum cells in x direction
    integer :: ny              ! maximum cells in y direction
    integer :: nz              ! maximum cells in z direction

    ! local variables
    integer :: nidx            ! index in matrix

    ! compute index
    nidx = i + nx*(j - 1) + nx*ny*(k - 1) + nx*ny*nz*(g - 1)

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

end module cmfd_utils
