! ------------------------------------------------------------------------------------------------------
! General program utilities that do not fit in other modules
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
  ! DERIVED TYPES
  ! -------------------------------------------------------------------------------------------
  type sparseMat
    real(wp), pointer :: A(:)
    integer, pointer :: irow(:)
    integer, pointer :: jcol(:)
    integer :: m, n
    integer :: nnz
  endtype sparseMat
  
  type seqarray
    integer, allocatable :: el(:)
  endtype seqarray
  public :: seqarray, error, warning, progress, configread
  
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
  
  ! CONFIGURATION FILE
  ! -------------------------------------------------------------------------------------------
  subroutine configread(filename, imat, ivec, N, prec, dig, theta, alg, dt, nbt, timing, ofile, matfile)
    ! -------------------------------------------------------------------------------------------
    ! configread
    ! -------------------------------------------------------------------------------------------
    ! Reads in the configuration file and returns the parameters accordingly
    ! Parameters: INPUTMAT, INPUTVEC, PRECISION, DIGITS, THETA, TIMESTEP, NUMBER, OUTPUTFILE, 
    ! MATRIXFILE
    !
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in)	::	filename
    character(len=*), intent(out)	::	imat, ivec, prec, alg, ofile, matfile
    integer, intent(out)		::	dig, dt, nbt, N
    real(wp), intent(out)		::	theta
    logical, intent(out)		::	timing
    
    integer, parameter			::	fh=15
    integer				::	ios=0, idx
    character(len=100)			::	buffer, arg1, arg2
    character				::	tab
    
    tab = char(9)
    ! Open the file
    open(unit=fh, file=filename, action='read', status='old')
    
    ! Set default values for the outputs
    imat = 'uninitialized'
    ivec = 'uninitialized'
    prec = 'uninitialized'
    alg = 'uninitialized'
    ofile = 'uninitialized'
    matfile ='uninitialized'
    N = -1
    dig = -1
    dt = -1
    nbt = -1
    timing = .false.
    theta = -1.0
    
    ! Read the file
    do while (ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if ( buffer(1:1) .ne. '#') then	! it's not a comment line, we interpret it
	buffer = trim(buffer)
	idx = scan(buffer,tab)
	arg1 = trim(buffer(1:idx-1))
	arg2 = trim(buffer(idx+1:))
	select case( arg1 )
	  case('INPUTMAT')	! Input file Bateman matrix
	    imat = arg2
	  case('INPUTVEC')	! Input file initial vector
	    ivec = arg2
	  case('PROBSIZE')	! Problem size
	    read(arg2, *) N
	  case('PRECISION')	! Requested precision
	    prec = arg2
	  case('DIGITS')	! Requested digits for multiprecision
	    read(arg2, *) dig
	  case('THETA')		! Threshold on half-life time
	    if(arg2 .eq. 'none' .or. arg2 .eq. 'NONE' .or. arg2 .eq. 'None') then
	      theta = 0.0
	    else
	      read(arg2,*) theta
	    endif
	  case('ALGORITHM')	! Requested Algorithm
	    alg = arg2
	  case('TIMESTEP')	! Requested timestep (t in exp(A*t))
	    read(arg2, *) dt
	  case('NUMBER')	! Number of timesteps
	    read(arg2, *) nbt
	  case('TIMECALCS')	! Display timing or not
	    read(arg2, *) timing
	  case('OUTPUTFILE')	! Time series output file
	    ofile = arg2
	  case('MATRIXFILE')	! Matrix exponential output file
	    matfile = arg2
	  case default
	    call error('in matrixrw::configread - unknown argument found in config file')
	end select
      endif
    end do
  end subroutine configread
end module utilities
