module cmfd_execute

  use cmfd_header, only: cmfd_obj

  implicit none

  type(cmfd_obj) :: cmfd  ! object for cmfd calculation

! may just need the object from global data then get rid of above

contains

!===============================================================================
! ALLOCATE_CMFD allocates all of the space for the cmfd object based on tallies.
!===============================================================================

  subroutine allocate_cmfd

  end subroutine allocate_cmfd

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections.
!===============================================================================

  subroutine compute_xs

  end subroutine compute_xs

!===============================================================================
! COMPUTE_DIFFCOEF computes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_diffcoef

  end subroutine compute_diffcoef

!===============================================================================
! COMPUTE_DHAT computes the nonlinear coupling coefficient
!===============================================================================

  subroutine compute_dhat

  end subroutine compute_dhat

!===============================================================================
! ACCUMULATE accumulates cycle to cycle diffusion parameters
!===============================================================================

  subroutine accumulate

  end subroutine accumulate

!===============================================================================
! CMFD_SOLVER in the main power iteration routine for the cmfd calculation
!===============================================================================

  subroutine cmfd_solver

  end subroutine cmfd_solver

!===============================================================================
! LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine loss_matrix

  end subroutine loss_matrix

!===============================================================================
! PROD_MATRIX creates the matrix representing production of neutrons
!===============================================================================

  subroutine prod_matrix

  end subroutine prod_matrix

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine convergence

  end subroutine convergence

!===============================================================================
! SOURCE_PDF calculates the probability distribution of the cmfd fission source
!===============================================================================

  subroutine source_pdf

  end subroutine source_pdf

!===============================================================================
! COUNT_SOURCE determines the number of source sites in each mesh box
!===============================================================================

  subroutine count_source

  end subroutine count_source

!===============================================================================
! WEIGHT_FACTORS calculates the weight adjustment factors for next MC cycle
!===============================================================================

  subroutine weight_factors

  end subroutine weight_factors

!===============================================================================
! ADJUST_WEIGHT adjusts the initial weight of the particle
!===============================================================================

  subroutine adjust_weight

! should this do all source particles at once or called each time before a
! neutron is born in the MC cycle?  MC21 is the latter

  end subroutine adjust_weight

end module cmfd_execute
