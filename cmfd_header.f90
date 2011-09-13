module cmfd_header

  implicit none

!===============================================================================
! cmfd is used to store diffusion parameters and other information for CMFD
! analysis.  
!===============================================================================

  type cmfd_obj

    ! array indices([1-x,2-y,3-z,4-g],upper bound)
    real              :: indices(4,1)

    ! cross sections
    real, allocatable :: totalxs(:,:,:,:)
    real, allocatable :: scattxs(:,:,:,:,:)
    real, allocatable :: fissxs(:,:,:,:,:)

    ! currents
    real, allocatable :: currentX(:,:,:,:)
    real, allocatable :: currentY(:,:,:,:)
    real, allocatable :: currentZ(:,:,:,:)

    ! coupling coefficients
    real, allocatable :: dtilda(:,:,:,:,:)
    real, allocatable :: dhat(:,:,:,:,:)

    ! core albedo boundary conditions
    real              :: albedo(6)

    ! dimensions of mesh cells (xloc,yloc,zloc,[hu,hv,hw])
    real, allocatable :: h(:,:,:,:)

    ! source probability distribution
    real, allocatable :: sourcepdf(:,:,:,:)

    ! source sites in each mesh box
    real, allocatable :: sourcecounts(:,:,:,:)

    ! weight adjustment factors
    real, allocatable :: weightfactors(:,:,:,:)

    ! we may need to add the mesh object

  end type cmfd_obj

end module cmfd_header
