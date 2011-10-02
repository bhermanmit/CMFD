   if ( present(errout) ) errout = error
end subroutine

   write(info%lun,'(a)') '</cmfd_t.xml>'
   call xml_close(info)
end subroutine


end subroutine

end module
