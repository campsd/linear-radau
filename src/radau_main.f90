! ------------------------------------------------------------------------------------------------------
! Main file
! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------
! daan.camps [at] cs.kuleuven.be
! Nov. 2016
! v0.1
! ------------------------------------------------------------------------------------------------------
program radau_main
  use units
  use utilities
  use matrixrw
  use radau_common
  use radau_ad_mod
  use radau_sa_mod
  use radau_fx_mod

  implicit none
  
  include 'mpif.h'
  include 'dmumps_struc.h'
  include 'zmumps_struc.h'
  
  ! Radau variables
  ! ----------------------------------------------------------------------------------------------------
  integer		::	N, NNZ, NVEC, NS, HEVALS
  integer, allocatable	::	RI(:), CJ(:)
  real(wp)		::	HFX, HINIT, TEND, TFX, TSA, TAD
  real(wp), allocatable	::	A(:), YZERO(:), YFINAL_FX(:), YFINAL_SA(:), YFINAL_AD(:), YINT(:,:)
  integer, dimension(4)		::	STEPS
  

  ! OpenMPI variables
  integer ( kind = 4 ) 	::	ierr, color
  real ( kind = 8 ) 	::	wtime
  integer ( kind = 4 ) 	::	world_id, mumps_id
  integer ( kind = 4 ) 	::	world_size, mumps_rank
  integer		::	MUMPS_COMM
  integer i

  call MPI_INIT(ierr)
  call MPI_Comm_size ( MPI_COMM_WORLD, world_size, ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD, world_id, ierr )

  ! Load Data
  ! ----------------------------------------------------------------------------------------------------	
  call matrixread_tosparse('PROVIDE MATRIX FILENAME', N, NNZ, RI, CJ, A)
  call vectorread_todense('PROVIDE VECTOR FILENAME', NVEC, YZERO)
  if( N .ne. NVEC ) then
	call error('matvec dimensions incompatible')
  endif

  ! Integration parameters
  ! ----------------------------------------------------------------------------------------------------
  HINIT = 0.25
  HFX = 200
  TEND = 2592000
  HEVALS = 12

  ! Run methods
  ! ----------------------------------------------------------------------------------------------------
  allocate(YFINAL_FX(N))
  allocate(YFINAL_SA(N))
  allocate(YFINAL_AD(N))
  allocate(YINT(N,NS-1))
  write(*,*) ''
  write(*,*) '--------------------------------------------------------'
  write(*,*) ''
  call progress('CALLING FIXED STEP RADAU')
  wtime = MPI_Wtime ( )
  call radau_fx( A, RI, CJ, N, NNZ, YZERO, TEND, HFX, YFINAL_FX )
  TFX = MPI_Wtime ( ) - wtime

  write(*,*) ''
  write(*,*) '--------------------------------------------------------'
  write(*,*) ''
  call progress('CALLING SEMI ADAPTIVE STEP RADAU')
  wtime = MPI_Wtime()
  call radau_sa( A, RI, CJ, N, NNZ, YZERO, TEND, HINIT, HEVALS, YFINAL_SA)
  TSA = MPI_Wtime ( ) - wtime

  write(*,*) ''
  write(*,*) '--------------------------------------------------------'
  write(*,*) ''
  call progress('CALLING ADAPTIVE STEP RADAU')
  wtime = MPI_Wtime ( )
  call radau_ad( A, RI, CJ, N, NNZ, YZERO, TEND, HINIT, YFINAL_AD)
  TAD = MPI_Wtime ( ) - wtime


  write(*,*) ''
  write(*,*) '--------------------------------------------------------'
  write(*,*) ''
  ! FX WITH INTERMEDIATE RESULTS
  !NS = 4
  !STEPS(1:4) =  (/ 2, 5, 10, 50 /)
  !call radau_fx_int(A, RI, CJ, N, NNZ, YZERO, H, STEPS, NS, YINT, YFINAL_FX)

  ! Analyse results
  ! ----------------------------------------------------------------------------------------------------
  call progress('COMPARE RESULTS')
  write ( *, '(a,g14.6,a)' ) '  Elapsed wall clock time FX = ', TFX, ' seconds.'
  write ( *, '(a,g14.6,a)' ) '  Elapsed wall clock time SA = ', TSA, ' seconds.'
  write ( *, '(a,g14.6,a)' ) '  Elapsed wall clock time AD = ', TAD, ' seconds.'
  write(*,*) 'ABSOLUTE DIFFERENCE FX VS SA: ', NORM2(YFINAL_SA - YFINAL_FX)
  write(*,*) 'RELATIVE DIFFERENCE FX VS SA: ', NORM2(YFINAL_SA - YFINAL_FX)/NORM2(YFINAL_SA)
  write(*,*) 'ABSOLUTE DIFFERENCE FX VS AD: ', NORM2(YFINAL_AD - YFINAL_FX)
  write(*,*) 'RELATIVE DIFFERENCE FX VS AD: ', NORM2(YFINAL_AD - YFINAL_FX)/NORM2(YFINAL_AD)
  write(*,*) 'ABSOLUTE DIFFERENCE SA VS AD: ', NORM2(YFINAL_AD - YFINAL_SA)
  write(*,*) 'RELATIVE DIFFERENCE SA VS AD: ', NORM2(YFINAL_AD - YFINAL_SA)/NORM2(YFINAL_AD)

  write(*,*) ''
  write(*,*) '--------------------------------------------------------'
  write(*,*) ''
  call progress('WRITE RESULTS TO FILE')
  call matrixwrite('outputFX.dat', YFINAL_FX, N, 1)
  call matrixwrite('outputSA.dat', YFINAL_SA, N, 1)
  call matrixwrite('outputAD.dat', YFINAL_AD, N, 1)
      
  call MPI_FINALIZE(IERR)

end program
