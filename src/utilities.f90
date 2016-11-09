! ------------------------------------------------------------------------------------------------------
! General program utilities 
! ------------------------------------------------------------------------------------------------------
! daan.camps [at] cs.kuleuven.be
! Oct. 2015
! v0.1
! ------------------------------------------------------------------------------------------------------
module utilities
  use units
  implicit none
  private
  save
  public ::error, warning, progress
  
  contains
  
  ! MESSAGING
  ! -------------------------------------------------------------------------------------------
  subroutine error(msg)
    ! -------------------------------------------------------------------------------------------
    ! error
    ! -------------------------------------------------------------------------------------------
    ! Displays the error passed in msg and stops the program execution
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in)	::	msg
    
    write(PRINT_UNIT, *) ''
    write(PRINT_UNIT, *) 'ERROR'
    write(PRINT_UNIT, *) ' ', msg
    write(PRINT_UNIT, *) ''
    write(PRINT_UNIT, *) ' Program execution terminated.'
    write(PRINT_UNIT, *) '----------------------------------------------------------------------------'
    write(PRINT_UNIT, *) ''
    stop
  end subroutine error
  
  subroutine warning(msg)
    ! -------------------------------------------------------------------------------------------
    ! warning
    ! -------------------------------------------------------------------------------------------
    ! Displays the warning passed in msg. Does not stop the program execution
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in)	::	msg
    
    write(PRINT_UNIT, *) ''
    write(PRINT_UNIT, *) 'WARNING'
    write(PRINT_UNIT, *) ' ', msg
    write(PRINT_UNIT, *) '----------------------------------------------------------------------------'
    write(PRINT_UNIT, *) ''
  end subroutine warning
  
  subroutine progress(msg)
    ! -------------------------------------------------------------------------------------------
    ! progress
    ! -------------------------------------------------------------------------------------------
    ! Displays a progress message
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in)	::	msg
    
    write(PRINT_UNIT, *) ''
    write(PRINT_UNIT, *) '... ', msg  
  end subroutine progress
end module utilities
