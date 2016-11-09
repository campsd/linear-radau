! ------------------------------------------------------------------------------------------------------
! radaucommon
! ------------------------------------------------------------------------------------------------------
! Common functionality for all radau versions
! ------------------------------------------------------------------------------------------------------
! daan.camps@cs.kuleuven.be
! Aug. 2016
! v0.1
! ------------------------------------------------------------------------------------------------------
module radau_common
	use units
	use utilities
	implicit none
	private
	save
	public	:: constructRealN, constructComplexN, constructRHS, estimateError, estimateErrorTilde
	integer	:: i,j
	contains
	
	! RADAU common subroutines
	! ----------------------------------------------------------------------------------------------
	subroutine constructRealN(A, RI, CJ, N, NNZ, H, ANR, RINR, CJNR, NNZNR)
		! we will assume that ordering of sparse matrix is row per row
		! we will use a temp array of size NNZ + N (max growth in NNZ)
		integer, intent(in) 	::	N, NNZ		
		real(wp), intent(in)	::	A(NNZ), H
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)
		real(wp), allocatable, intent(out)	::	ANR(:)
		integer, allocatable, intent(out)	::	RINR(:), CJNR(:)
		integer, intent(out)	::	NNZNR

		real(wp)	::	TMPA(NNZ+N)
		integer		::	TMPRI(NNZ+N), TMPCJ(NNZ+N)
		integer		::	d 	! keeps track of the current diagonal element we're handling
		real(wp)	::	par

		j = 1
		d = 1
		par = GAMM/H	! diagonal parameter
		do i = 1, NNZ
			if (RI(i) > CJ(i)) then	! left of diagonal
				if (RI(i)  > d) then	! we need to include previous diagonal elements
					do while (RI(i) > d)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! add the current element
				TMPA(j) = A(i)
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
			elseif (RI(i) < CJ(i)) then ! right of diagonal
				if (RI(i) > d - 1) then	! we need to include previous diagonal elements
					do while (RI(i) > d - 1)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! add the current element
				TMPA(j) = A(i)
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
			else	! on diagonal
				if (RI(i) > d) then ! we need to include previous diagonal elements
					do while (RI(i) > d)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! modify the current diagonal element
				TMPA(j) = A(i) - par
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
				d = d + 1
			endif
		end do

		! potentially some last diagonal entries need to be included
		if (d < N + 1) then
			do while (d < N + 1)
				TMPA(j) = -par
				TMPRI(j) = d
				TMPCJ(j) = d
				d = d + 1
				j = j + 1
			enddo
		endif

		! allocate and create final arrays
		j = j - 1		
		allocate(ANR(j))
		allocate(RINR(j))
		allocate(CJNR(j))
		ANR = TMPA(:j)
		RINR = TMPRI(:j)
		CJNR = TMPCJ(:j)
		NNZNR = j

	end subroutine constructRealN

	subroutine constructComplexN(A, RI, CJ, N, NNZ, H, ANC, RINC, CJNC, NNZNC)
		integer, intent(in) 	::	N, NNZ		
		real(wp), intent(in)	::	A(NNZ), H
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)
		complex(wp), allocatable, intent(out)	::	ANC(:)
		integer, allocatable, intent(out)	::	RINC(:), CJNC(:)
		integer, intent(out)	::	NNZNC

		complex(wp)	::	TMPA(NNZ+N)
		integer		::	TMPRI(NNZ+N), TMPCJ(NNZ+N)
		integer		::	d
		complex(wp)	::	par
		
		d = 1
		j=1
		par = complex(ALPH/H, BETA/H)	! diagonal parameter

		do i = 1, NNZ
			if (RI(i) > CJ(i)) then ! left of diagonal
				if (RI(i)  > d) then	! we need to include previous diagonal elements
					do while (RI(i) > d)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! add the current element
				TMPA(j) = complex(A(i), 0.0)
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
			elseif (RI(i) < CJ(i)) then ! right of diagonal
				if (RI(i) > d - 1) then	! we need to include previous diagonal elements
					do while (RI(i) > d - 1)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! add the current element
				TMPA(j) = complex(A(i), 0.0)
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
			else	! on diagonal
				if (RI(i) > d) then ! we need to include previous diagonal elements
					do while (RI(i) > d)
						TMPA(j) = - par
						TMPRI(j) = d
						TMPCJ(j) = d
						d = d + 1
						j = j + 1
					enddo
				endif
				! modify the current diagonal element
				TMPA(j) = complex(A(i), 0.0) - par
				TMPRI(j) = RI(i)
				TMPCJ(j) = CJ(i)
				j = j + 1
				d = d + 1
			endif
		end do

		! potentially some last diagonal entries need to be included
		if (d < N + 1) then
			do while (d < N + 1)
				TMPA(j) = -par
				TMPRI(j) = d
				TMPCJ(j) = d
				d = d + 1
				j = j + 1
			enddo
		endif

		! allocate and create final arrays
		j = j - 1		
		allocate(ANC(j))
		allocate(RINC(j))
		allocate(CJNC(j))
		ANC = TMPA(:j)
		RINC = TMPRI(:j)
		CJNC = TMPCJ(:j)
		NNZNC = j
		
	end subroutine constructComplexN

	subroutine constructRHS(A, RI, CJ, N, NNZ, Y, R1, R2, R3)
		integer, intent(in) 	::	N, NNZ		
		real(wp), intent(in)	::	A(NNZ), Y(N)
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)
		real(wp), intent(out)	::	R1(N), R2(N), R3(N)
		
		real(wp)	:: SPMV(N)
		integer		:: i		
		! Compute SpMV -A*Y
		call SparseMV(A, RI, CJ, N, NNZ, Y, SPMV)
		SPMV = -SPMV

		! Deduce R1, R2, R3
		R1 = TI1 * SPMV
		R2 = TI2 * SPMV
		R3 = TI3 * SPMV

	end subroutine constructRHS

	subroutine convertSolution(W1, W2, Z3, N)
		integer, intent(in) ::		N
		real(wp), intent(in) ::		W1(N), W2(N)
		real(wp), intent(out) ::	Z3(N)

		Z3 = T31 * W1 + W2
	end subroutine

	subroutine estimateError(A, RI, CJ, N, NNZ, YCURR, YPREV, W1, W2, W3, HCURR, HPREV, AERRPREV, CHANGE, mumps_R)
		! Estimate the error based on Eq. (8.19)
		include 'mpif.h'		
		include 'dmumps_struc.h'

		!I/O
		integer, intent(in)	::	N, NNZ
		real(wp), intent(in)	::	A(NNZ), YCURR(N), YPREV(N), W1(N), W2(N), W3(N)
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)		
		real(wp), intent(inout)	::	AERRPREV, HCURR, HPREV
		logical, intent(out)	::	CHANGE
		type (DMUMPS_STRUC)	::	mumps_R

		!Subroutine variables
		real(wp)	::	Z1(N), Z2(N), Z3(N), ERRX(N), SPMV(N), XK(N), ERR1(N), AERR1, ERR2(N), AERR2, AERR, H1, H2, HIDEAL
		AERR1 = 0.0		
		AERR2 = 0.0
		!write(*,*) 'ERROR EST START', AERRPREV
		Z1 = T11 * W1 + T12 * W2 + T13 * W3
		Z2 = T21 * W1 + T22 * W2 + T23 * W3
		Z3 = T31 * W1 + W2
		ERRX = E1 * Z1 + E2 * Z2 + E3 * Z3
		call SparseMV(A, RI, CJ, N, NNZ, YPREV, SPMV)
		XK = GAMM0 * HCURR * SPMV + ERRX

		mumps_R%RHS = XK
		mumps_R%JOB = 3
		call DMUMPS(mumps_R)
		ERR1 = - GAMM0 * HCURR * mumps_R%RHS

		
		do i=1, N
			AERR1 = AERR1 + (ERR1(i) / (ATOL + RTOL * max(YCURR(i), YPREV(i))))**2
		enddo
		AERR1 = sqrt(AERR1/N)
		AERR = AERR1


		if (AERR1 < -1) then
			call SparseMV(A, RI, CJ, N, NNZ, YPREV + ERR1, SPMV)
			XK = GAMM0 * HCURR * SPMV + ERRX
			
			mumps_R%RHS = XK
			mumps_R%JOB = 3
			call DMUMPS(mumps_R)
			ERR2 = - GAMM0 * HCURR * mumps_R%RHS

			
			do i=1, N
				AERR2 = AERR2 + (ERR2(i) / (ATOL + RTOL * max(YCURR(i), YPREV(i))))**2
			enddo
			AERR2 = sqrt(AERR2/N)
			AERR = AERR2
		endif

		H1 = FAC * HCURR * AERR**(-0.25)
		H2 = H1 * HCURR / HPREV * ( (AERRPREV / AERR)**0.25 )
		HIDEAL = min(H1,H2)


		AERRPREV = AERR
		HPREV = HCURR
		CHANGE = .false.

		if (C1*HCURR > HIDEAL .or. HIDEAL > C2*HCURR ) then
			HCURR = HIDEAL
			CHANGE = .true.
		endif
		!write(*,*) 'ERROR EST END',	AERR, AERR1, AERR2, H1, H2
	end subroutine

	subroutine estimateErrorTilde(A, RI, CJ, N, NNZ, YCURR, YPREV, W1, W2, W3, HCURR, HPREV, AERRPREV, CHANGE, mumps_R)
		! Estimate the error based on Eq. (8.20)
		include 'mpif.h'		
		include 'dmumps_struc.h'

		!I/O
		integer, intent(in)	::	N, NNZ
		real(wp), intent(in)	::	A(NNZ), YCURR(N), YPREV(N), W1(N), W2(N), W3(N)
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)		
		real(wp), intent(inout)	::	AERRPREV, HCURR, HPREV
		logical, intent(out)	::	CHANGE
		type (DMUMPS_STRUC)	::	mumps_R

		!Subroutine variables
		real(wp)	::	Z1(N), Z2(N), Z3(N), ERRX(N), SPMV(N), XK(N), ERR1(N), AERR1, ERR2(N), AERR2, AERR, H1, H2, HIDEAL
		AERR1 = 0.0		
		AERR2 = 0.0
		!write(*,*) 'ERROR EST START', AERRPREV
		Z1 = T11 * W1 + T12 * W2 + T13 * W3
		Z2 = T21 * W1 + T22 * W2 + T23 * W3
		Z3 = T31 * W1 + W2
		ERRX = E1 * Z1 + E2 * Z2 + E3 * Z3
		call SparseMV(A, RI, CJ, N, NNZ, YPREV, SPMV)
		XK = GAMM0 * HCURR * SPMV + ERRX

		mumps_R%RHS = XK
		mumps_R%JOB = 3
		call DMUMPS(mumps_R)
		ERR1 = - GAMM0 * HCURR * mumps_R%RHS

		call SparseMV(A, RI, CJ, N, NNZ, YPREV + ERR1, SPMV)
		XK = GAMM0 * HCURR * SPMV + ERRX
			
		mumps_R%RHS = XK
		mumps_R%JOB = 3
		call DMUMPS(mumps_R)
		ERR2 = - GAMM0 * HCURR * mumps_R%RHS
			
		do i=1, N
			AERR2 = AERR2 + (ERR2(i) / (ATOL + RTOL * max(YCURR(i), YPREV(i))))**2
		enddo
		AERR2 = sqrt(AERR2/N)
		AERR = AERR2

		H1 = FAC * HCURR * AERR**(-0.25)
		H2 = H1 * HCURR / HPREV * ( (AERRPREV / AERR)**0.25 )
		HIDEAL = min(H1,H2)

		AERRPREV = AERR
		HPREV = HCURR
		CHANGE = .false.

		if (C1*HCURR > HIDEAL .or. HIDEAL > C2*HCURR ) then
			HCURR = HIDEAL
			CHANGE = .true.
		endif
		!write(*,*) 'ERROR EST END',	AERR, AERR1, AERR2, H1, H2
	end subroutine

	subroutine SparseMV(A, RI, CJ, N, NNZ, Y, SPMV)
		integer, intent(in)	::	N, NNZ
		real(wp), intent(in)	::	A(NNZ), Y(N) 
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)
		real(wp), intent(out)	::	SPMV(N)
		
		SPMV(:) = 0.0
		do i=1, NNZ
			SPMV(RI(i)) = SPMV(RI(i)) + A(i)*Y(CJ(i))
		enddo
	end subroutine

end module radau_common
