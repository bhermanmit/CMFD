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
    integer :: nx      ! number of volumes in x-direction
    integer :: ny      ! number of volumes in y-direction
    integer :: nz      ! number of volumes in z-direction
    integer :: ng      ! number of energy groups
    integer :: i       ! x iteration counter
    integer :: j       ! y iteration counter
    integer :: k       ! z iteration counter
    integer :: g       ! group iteration counter
    integer :: h       ! group iteration counter
    integer :: map_idx ! vector location for core map
    integer :: matid   ! material id number
    integer :: hg_idx  ! energy h-->g vector location 
    integer :: dim_idx ! vector location for dimensions

    ! read xml input file
    call read_xml_file_cmfd_t('cmfd.xml')

    ! get mesh and group indices
    nx = geometry%nx
    ny = geometry%ny
    nz = geometry%nz
    ng = geometry%ng

    ! allocate cmfd object
    allocate(cmfd%totalxs(nx,ny,nz,ng))
    allocate(cmfd%scattxs(nx,ny,nz,ng,ng))
    allocate(cmfd%nfissxs(nx,ny,nz,ng,ng))
    allocate(cmfd%diffcof(nx,ny,nz,ng))
    allocate(cmfd%dtilda(nx,ny,nz,ng,6))
    allocate(cmfd%dhat(nx,ny,nz,ng,6))
    allocate(cmfd%hxyz(nx,ny,nz,3))
    allocate(cmfd%coremap(nx,ny,nz))
    allocate(cmfd%sourcepdf(nx,ny,nz,ng))  ! take this out when interface with openMC

    ! record indices in object
    cmfd%indices(1) = nx
    cmfd%indices(2) = ny
    cmfd%indices(3) = nz
    cmfd%indices(4) = ng

    ! set boundary conditions
    cmfd%albedo = geometry%bc

    ! set dhat to 0.0
    cmfd%dhat = 0.0

    ! check core map dimensions
    if (size(geometry%mesh,1) /= nx*ny*nz) then
    
      ! write out fatal error
      print *,'FATAL ===> core map dimensions not consistent'
      STOP

    end if

    ! read in core map and xs
    GROUP: do g = 1,ng

      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! get vector idx for core map
            map_idx = get_matrix_idx(i,j,k,1,nx,ny,nz)

            ! extract material identifier
            matid = geometry%mesh(map_idx)

            ! record to core map
            cmfd%coremap(i,j,k) = matid

            ! check to see if matid is there
            if (matid > size(mat,1)) then

              ! write out fatal error
              print *, 'Fatal Error ===> MATERIAL ID',matid,' NOT SET!'
              STOP

            end if

            ! set tot xs and diff coef
            cmfd%totalxs(i,j,k,g) = mat(matid)%totxs(g)
            cmfd%diffcof(i,j,k,g) = mat(matid)%diffcoef(g)

            ! loop around outgoing energy groups 
            ELOOP: do h = 1,ng

              ! get vector h-->g index
              hg_idx = g + ng*(h - 1)
              cmfd%scattxs(i,j,k,h,g) = mat(matid)%scattxs(hg_idx)
              cmfd%nfissxs(i,j,k,h,g) = mat(matid)%nfissxs(hg_idx)

            end do ELOOP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do GROUP

    ! get dimensions
    if (associated(geometry%uniform)) then

      ! record uniform dimensions
      cmfd%hxyz(:,:,:,1) = geometry%uniform(1)
      cmfd%hxyz(:,:,:,2) = geometry%uniform(2)
      cmfd%hxyz(:,:,:,3) = geometry%uniform(3)

    else if (associated(geometry%dx)) then

      ! loop through to get nonuniform dimensions
      ZLOOP2: do k = 1,nz

        YLOOP2: do j = 1,ny

          XLOOP2: do i = 1,nx

            ! get vector idx for dimension 
            dim_idx = get_matrix_idx(i,j,k,1,nx,ny,nz)

            ! record dimension
            cmfd%hxyz(i,j,k,1) = geometry%dx(dim_idx)
            cmfd%hxyz(i,j,k,2) = geometry%dy(dim_idx)
            cmfd%hxyz(i,j,k,3) = geometry%dz(dim_idx)

          end do XLOOP2

        end do YLOOP2

      end do ZLOOP2

    else

      ! user forgot dimensions
      print *,'Fatal Error ===> Dimensions not entered correctly!'
      STOP

    end if

    ! echo input
    print *, 'Dimensions:'
    print *,cmfd%indices
    print *, 'CORE MAP:'
    print *,cmfd%coremap
    print *, 'TOTAL XS:'
    print *,cmfd%totalxs
    print *, 'SCATTERING XS:'
    print *,cmfd%scattxs
    print *, 'Nu-FISSION XS:'
    print *,cmfd%nfissxs
    print *, 'DIFFUSION COEFFICIENT:'
    print *,cmfd%diffcof
    print *, 'BOUNDARY CONDITIONS:'
    print *,cmfd%albedo
    print *, 'CORE CELL DIMENSIONS X:'
    print *,cmfd%hxyz(:,:,:,1)
    print *, 'CORE CELL DIMENSIONS Y:'
    print *,cmfd%hxyz(:,:,:,2)
    print *, 'CORE CELL DIMENSIONS Z:'
    print *,cmfd%hxyz(:,:,:,3)

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
