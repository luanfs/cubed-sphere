module miscellaneous
  !===============================================================================================
  !   This module contains all the miscellaneous routines
  !===============================================================================================

  !Global constants
  use constants, only: &
        i4

 implicit none

contains 

  !===============================================================================================
  !    MISCELLANEOUS
  !===============================================================================================

  subroutine getunit ( iunit )
    !----------------------------------------------------------
    ! GETUNIT returns a free FORTRAN unit number.
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !    John Burkardt
    !    18 September 2005
    !----------------------------------------------------------------------------
    integer ( i4 ):: i
    integer ( i4 ):: ios
    integer ( i4 ):: iunit
    logical:: lopen

    iunit = 0
    do i = 11, 99
       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
          inquire ( unit = i, opened = lopen, iostat = ios )
          if ( ios == 0 ) then
             if ( .not. lopen ) then
                iunit = i
                return
             end if
          end if
       end if
    end do

    return
  end subroutine getunit

end module miscellaneous
