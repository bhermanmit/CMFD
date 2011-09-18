module xml_data_cmfd_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   integer, private :: lurep_
   logical, private :: strict_

type geometry_xml
   integer                                         :: nx
   integer                                         :: ny
   integer                                         :: nz
   integer                                         :: ng
   integer, dimension(:), pointer                  :: mesh => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: uniform => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: dx => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: dy => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: dz => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: bc => null()
end type geometry_xml

type mat_xml
   integer                                         :: uid
   real(kind=kind(1.0d0)), dimension(:), pointer   :: totxs => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: scattxs => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: nfissxs => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: diffcoef => null()
end type mat_xml
   type(geometry_xml)                              :: geometry
   type(mat_xml), dimension(:), pointer            :: mat => null()
contains
subroutine read_xml_type_geometry_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(geometry_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(geometry_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_geometry_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_geometry_xml_array

subroutine read_xml_type_geometry_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(geometry_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_nx
   logical                                         :: has_ny
   logical                                         :: has_nz
   logical                                         :: has_ng
   logical                                         :: has_mesh
   logical                                         :: has_uniform
   logical                                         :: has_dx
   logical                                         :: has_dy
   logical                                         :: has_dz
   logical                                         :: has_bc
   has_nx                               = .false.
   has_ny                               = .false.
   has_nz                               = .false.
   has_ng                               = .false.
   has_mesh                             = .false.
   has_uniform                          = .false.
   has_dx                               = .false.
   has_dy                               = .false.
   has_dz                               = .false.
   has_bc                               = .false.
   call init_xml_type_geometry_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('nx')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nx, has_nx )
      case('ny')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%ny, has_ny )
      case('nz')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nz, has_nz )
      case('ng')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%ng, has_ng )
      case('mesh')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%mesh, has_mesh )
      case('uniform')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%uniform, has_uniform )
      case('dx')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dx, has_dx )
      case('dy')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dy, has_dy )
      case('dz')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dz, has_dz )
      case('bc')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%bc, has_bc )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_nx ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on nx')
   endif
   if ( .not. has_ny ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on ny')
   endif
   if ( .not. has_nz ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on nz')
   endif
   if ( .not. has_ng ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on ng')
   endif
   if ( .not. has_mesh ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on mesh')
   endif
   if ( .not. has_uniform ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on uniform')
   endif
   if ( .not. has_dx ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on dx')
   endif
   if ( .not. has_dy ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on dy')
   endif
   if ( .not. has_dz ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on dz')
   endif
   if ( .not. has_bc ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on bc')
   endif
end subroutine read_xml_type_geometry_xml
subroutine init_xml_type_geometry_xml_array( dvar )
   type(geometry_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_geometry_xml_array
subroutine init_xml_type_geometry_xml(dvar)
   type(geometry_xml) :: dvar
end subroutine init_xml_type_geometry_xml
subroutine write_xml_type_geometry_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(geometry_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_geometry_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_geometry_xml_array

subroutine write_xml_type_geometry_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(geometry_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'nx', indent+3, dvar%nx)
   call write_to_xml_integer( info, 'ny', indent+3, dvar%ny)
   call write_to_xml_integer( info, 'nz', indent+3, dvar%nz)
   call write_to_xml_integer( info, 'ng', indent+3, dvar%ng)
   call write_to_xml_integer_array( info, 'mesh', indent+3, dvar%mesh)
   call write_to_xml_double_array( info, 'uniform', indent+3, dvar%uniform)
   call write_to_xml_double_array( info, 'dx', indent+3, dvar%dx)
   call write_to_xml_double_array( info, 'dy', indent+3, dvar%dy)
   call write_to_xml_double_array( info, 'dz', indent+3, dvar%dz)
   call write_to_xml_double_array( info, 'bc', indent+3, dvar%bc)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_geometry_xml

subroutine read_xml_type_mat_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mat_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(mat_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_mat_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_mat_xml_array

subroutine read_xml_type_mat_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mat_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_uid
   logical                                         :: has_totxs
   logical                                         :: has_scattxs
   logical                                         :: has_nfissxs
   logical                                         :: has_diffcoef
   has_uid                              = .false.
   has_totxs                            = .false.
   has_scattxs                          = .false.
   has_nfissxs                          = .false.
   has_diffcoef                         = .false.
   call init_xml_type_mat_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('uid')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%uid, has_uid )
      case('totxs')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%totxs, has_totxs )
      case('scattxs')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scattxs, has_scattxs )
      case('nfissxs')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nfissxs, has_nfissxs )
      case('diffcoef')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%diffcoef, has_diffcoef )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_uid ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on uid')
   endif
   if ( .not. has_totxs ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on totxs')
   endif
   if ( .not. has_scattxs ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scattxs')
   endif
   if ( .not. has_nfissxs ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on nfissxs')
   endif
   if ( .not. has_diffcoef ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on diffcoef')
   endif
end subroutine read_xml_type_mat_xml
subroutine init_xml_type_mat_xml_array( dvar )
   type(mat_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_mat_xml_array
subroutine init_xml_type_mat_xml(dvar)
   type(mat_xml) :: dvar
end subroutine init_xml_type_mat_xml
subroutine write_xml_type_mat_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mat_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_mat_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_mat_xml_array

subroutine write_xml_type_mat_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mat_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'uid', indent+3, dvar%uid)
   call write_to_xml_double_array( info, 'totxs', indent+3, dvar%totxs)
   call write_to_xml_double_array( info, 'scattxs', indent+3, dvar%scattxs)
   call write_to_xml_double_array( info, 'nfissxs', indent+3, dvar%nfissxs)
   call write_to_xml_double_array( info, 'diffcoef', indent+3, dvar%diffcoef)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_mat_xml

subroutine read_xml_file_cmfd_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=200), dimension(1:100)   :: data
   integer                                :: nodata
   logical                                         :: has_geometry
   logical                                         :: has_mat
   has_geometry                         = .false.
   has_mat                              = .false.
   allocate(mat(0))

   call init_xml_file_cmfd_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag /= '!--' ) exit
   enddo
   if ( starttag /= "cmfd" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "cmfd"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .false.
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('geometry')
         call read_xml_type_geometry_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            geometry, has_geometry )
      case('mat')
         call read_xml_type_mat_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            mat, has_mat )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_geometry ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on geometry')
   endif
   if ( .not. has_mat ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on mat')
   endif
   if ( present(errout) ) errout = error
end subroutine

subroutine write_xml_file_cmfd_t(fname, lurep)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep

   type(XML_PARSE)                        :: info
   integer                                :: indent = 0

   call xml_open( info, fname, .false. )
   call xml_options( info, report_errors=.true.)
   if ( present(lurep) ) then
       call xml_options( info, report_errors=.true.)
   endif
   write(info%lun,'(a)') &
      '<cmfd>'
   call write_xml_type_geometry_xml( info, 'geometry', indent+3, geometry)
   call write_xml_type_mat_xml_array( info, 'mat', indent+3, mat)
   write(info%lun,'(a)') '</cmfd>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_cmfd_t

end subroutine

end module
