subroutine read_parameters()

    implicit none

    integer::namelist_err=0

    ! for now, there is only one parameter to read
    namelist /indi_star/ mist_sample

    rewind(1)
    read(1, nml=indi_star, end=111)
    rewind(1)

    ! why it becomes -1? Error happens?
    ! write(*,*) 'namelist_err: ', namelist_err
    ! if (namelist_err < 0) then
    !     namelist_ok = .false.
    ! else if (namelist_err > 0) then
    !     namelist_ok = .false.
    ! end if

    ! why re-reading makes error????
    ! close(namelist_unit)
111 return

end subroutine read_parameters
