! ------------------------------------------------------------------------------------------------------
! daan.camps [at] cs.kuleuven.be
! Nov. 2016
! v0.1
! ------------------------------------------------------------------------------------------------------
module units
  implicit none
  !integer, parameter 	::	wp = selected_real_kind(33, 4931) ! quadruple precision (https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format)
  integer, parameter	::	wp = selected_real_kind(15, 307) ! double precision
  ! IF YOU CHANGE THIS PARAMETER TO QP, WE NEED TO UPDATE THE MPI SEND !!
  integer, parameter	::	PRINT_UNIT = 0	! Used for writing and reading
! RADAU integration parameters	
	! ----------------------------------------------------------------------------------------------
	! Similarity matrices Butcher tableau

	! Last row of similarity matrix T
	real(wp), parameter ::	T31 = 0.96604818261509293619D0
	! Rowsums of inverse similarity matrix TI	
	real(wp), parameter :: TI1 = 5.2065496818148400137715725577436D0
	real(wp), parameter :: TI2 = -4.0297778578124168191720855247695D0
	real(wp), parameter :: TI3 = 1.4730151100815938036703300895169D0
	! Eigenvalues inverse Butcher tableau	
	real(wp), parameter :: GAMM = 3.6378342527444960731080210111D0
	real(wp), parameter :: ALPH = 2.68108287362775247213855488906D0
	real(wp), parameter :: BETA = 3.05043019924741058363858025365D0
	! Error estimation
	real(wp), parameter :: GAMM0 = 0.274888829595677341988135822248D0
	real(wp), parameter :: E1 = -2.76230545474859913949568563737D0
	real(wp), parameter :: E2 = 0.379935598252728842265175177889D0
	real(wp), parameter :: E3 = -0.0916296098652257806627119407494D0
	real(wp), parameter :: ATOL = 1e2
	real(wp), parameter :: RTOL = 1e-4
	real(wp), parameter :: FAC = 0.9
	real(wp), parameter :: C1 = 1
	real(wp), parameter :: C2 = 4
	real(wp), parameter ::	T11 = 0.091232394870892942792D0
	real(wp), parameter ::	T12 = -0.14125529502095420843D0
	real(wp), parameter ::	T13 = -0.030029194105147424492D0
	real(wp), parameter ::	T21 = 0.24171793270710701896D0
	real(wp), parameter ::	T22 = 0.20412935229379993199D0
	real(wp), parameter ::	T23 = 0.38294211275726193779D0
end module units
