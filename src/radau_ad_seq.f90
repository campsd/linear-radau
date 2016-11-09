! ------------------------------------------------------------------------------------------------------
! radau_ad_seq
! ------------------------------------------------------------------------------------------------------
! Sequential implementation of IRK integrator radauIIA with fully adaptive stepsize
! ------------------------------------------------------------------------------------------------------
! Nov. 2016
! ------------------------------------------------------------------------------------------------------

module radau_ad_mod
	use units
	use utilities
	use radau_common
	use matrixrw

	implicit none
	private
	public	::	radau_ad
	contains

		subroutine radau_ad( A, RI, CJ, N, NNZ, YZERO, TEND, HINIT, YFINAL)
		
		include 'mpif.h'
		include 'dmumps_struc.h'
		include 'zmumps_struc.h'

		! declare all variables
		! I/O
		integer, intent(in)	::	N, NNZ	
		real(wp), intent(in)	::	A(NNZ), YZERO(N)
		integer, intent(in)	::	RI(NNZ), CJ(NNZ)
		real(wp), intent(in)	::	TEND, HINIT
		real(wp), intent(out)	::	YFINAL(N)

		! subroutine specific
		integer 			::	IERR, NNZNR, NNZNC, NSTEPS, i, j, k, world_id, world_size
		real(wp), allocatable		::	ANR(:)
		complex(wp), allocatable	::	ANC(:)
		integer, allocatable		::	RINR(:), CJNR(:), RINC(:), CJNC(:)	
		real(wp)			::	YCURR(N), YPREV(N), R1(N), R2(N), R3(N), W1(N), W2(N), TCURR, HCURR, HPREV, AERRPREV
                integer p, id, cnt
		logical	CHANGE

		! MUMPS variables
		type (DMUMPS_STRUC)	::	mumps_R
		type (ZMUMPS_STRUC)	:: 	mumps_C
		
		! INITIALIZE REAL SOLVER
		mumps_R%COMM = MPI_COMM_WORLD
		mumps_R%JOB = -1
		mumps_R%SYM = 0
		mumps_R%PAR = 1
		call DMUMPS(mumps_R)

		! INITIALIZE COMPLEX SOLVER
		mumps_C%COMM = MPI_COMM_WORLD
		mumps_C%JOB = -1
		mumps_C%SYM = 0
		mumps_C%PAR = 1
		call ZMUMPS(mumps_C)

		! (sequential - no processor check)
		! REAL
		call constructRealN(A, RI, CJ, N, NNZ, HINIT, ANR, RINR, CJNR, NNZNR)
		mumps_R%N = N
		mumps_R%NZ = NNZNR
		allocate( mumps_R%IRN( NNZNR ) )
		mumps_R%IRN = RINR
		allocate( mumps_R%JCN( NNZNR ) )
		mumps_R%JCN = CJNR
		allocate( mumps_R%A( NNZNR ) )
		mumps_R%A = ANR
		allocate( mumps_R%RHS ( N  ) )

		! COMPLEX
		call constructComplexN(A, RI, CJ, N, NNZ, HINIT, ANC, RINC, CJNC, NNZNC)
		mumps_C%N = N
		mumps_C%NZ = NNZNC
		allocate( mumps_C%IRN( NNZNC ) )
		mumps_C%IRN = RINC
		allocate( mumps_C%JCN( NNZNC ) )
		mumps_C%JCN = CJNC
		allocate( mumps_C%A( NNZNC ) )
		mumps_C%A = ANC
		allocate( mumps_C%RHS ( N  ) )

		! Factorise the problems
		mumps_R%JOB = 4
		mumps_R%ICNTL(4) = 0 !no printing
		mumps_R%KEEP(1) = 45 !speeds up fwd/bwd step
		call DMUMPS(mumps_R)

		mumps_C%JOB = 4
		mumps_C%ICNTL(4) = 0 !no printing	
		mumps_C%KEEP(1) = 40 !speeds up fwd/bwd step		
		call ZMUMPS(mumps_C)

		! Main loop
		YPREV = YZERO
		YCURR = YZERO
		TCURR = 0.0
		HCURR = HINIT
		HPREV = HINIT
		AERRPREV = 1.0
		k = 1
		CHANGE = .false.

		do !while (TCURR < TEND)
			! Set RHS	
			write (*,*) 'Current time : ', TCURR, ' current step : ', HCURR				
			call constructRHS(A, RI, CJ, N, NNZ, YPREV, R1, R2, R3)
			mumps_R%RHS = R1
			do j = 1, N
				mumps_C%RHS(j) = complex(R2(j), R3(j))
			enddo

			! Solve
			mumps_R%JOB = 3
			mumps_C%JOB = 3
			call DMUMPS(mumps_R)			
			call ZMUMPS(mumps_C)

			! Update timestep
			YCURR = YCURR + T31 * mumps_R%RHS + real(mumps_C%RHS)
			TCURR = TCURR + HCURR

			! Estimate error
			if ((k .eq. 1) .or. ((CHANGE .eqv. .true.) .and. (AERRPREV > 10))) then
				write(*,*) 'USING TILDE ESTIMATE  AERRPREV: ', AERRPREV
				call estimateErrorTilde(A, RI, CJ, N, NNZ, YCURR, YPREV, mumps_R%RHS, real(mumps_C%RHS), dimag(mumps_C%RHS), HCURR, &
					   HPREV, AERRPREV, CHANGE, mumps_R)
			else
				call estimateError(A, RI, CJ, N, NNZ, YCURR, YPREV, mumps_R%RHS, real(mumps_C%RHS), dimag(mumps_C%RHS), HCURR, &
					   HPREV, AERRPREV, CHANGE, mumps_R)
			endif

			if (CHANGE) then
				! Refactor
				write(*,*) 'CHANGE OF STEPSIZE'
				call constructRealN(A, RI, CJ, N, NNZ, HCURR, ANR, RINR, CJNR, NNZNR)				
				mumps_R%NZ = NNZNR
				mumps_R%IRN = RINR
				mumps_R%JCN = CJNR
				mumps_R%A = ANR
				mumps_R%JOB = 2
				call DMUMPS(mumps_R)
				call constructComplexN(A, RI, CJ, N, NNZ, HCURR, ANC, RINC, CJNC, NNZNC)
				mumps_C%NZ = NNZNC
				mumps_C%IRN = RINC
				mumps_C%JCN = CJNC
				mumps_C%A = ANC
				mumps_C%JOB = 2
				call ZMUMPS(mumps_C)

				! Step size is decreased (discard current step)
				if (HCURR < HPREV) then 
					write(*,*) 'DECREASING STEPSIZE, DISCARD LAST STEP'
					YCURR = YPREV
					TCURR = TCURR - HPREV
				endif
			endif

			YPREV = YCURR
			! Check if next step is too far
			if (TCURR + HCURR > TEND) then
				EXIT
			endif

			! Check if current step is 0
			if (HCURR .lt. 1e-20) then
				call error('ZERO STEP')
			endif

			k = k + 1
		enddo

		! Last step to get at TEND
		if (TEND - TCURR > 0.0) then
			HCURR = TEND - TCURR
			write (*,*) 'Current time : ', TCURR, ' final step : ', HCURR	
			! Refactor
			call constructRealN(A, RI, CJ, N, NNZ, HCURR, ANR, RINR, CJNR, NNZNR)
			mumps_R%NZ = NNZNR
			mumps_R%IRN = RINR
			mumps_R%JCN = CJNR
			mumps_R%A = ANR
			mumps_R%JOB = 2
			call DMUMPS(mumps_R)
			call constructComplexN(A, RI, CJ, N, NNZ, HCURR, ANC, RINC, CJNC, NNZNC)
			mumps_C%NZ = NNZNC
			mumps_C%IRN = RINC
			mumps_C%JCN = CJNC
			mumps_C%A = ANC
			mumps_C%JOB = 2
			call ZMUMPS(mumps_C)

			! Set RHS					
			call constructRHS(A, RI, CJ, N, NNZ, YPREV, R1, R2, R3)
			mumps_R%RHS = R1
			do j = 1, N
				mumps_C%RHS(j) = complex(R2(j), R3(j))
			enddo

			! Solve
			mumps_R%JOB = 3
			mumps_C%JOB = 3
			call DMUMPS(mumps_R)			
			call ZMUMPS(mumps_C)

			! Update timestep
			YCURR = YCURR + T31 * mumps_R%RHS + real(mumps_C%RHS)
			TCURR = TCURR + HCURR
		endif
		
		! Set final solution                        
		YFINAL = YCURR	 

		! WRAP UP SOLVERS
		deallocate( mumps_R%IRN )
		deallocate( mumps_R%JCN )
		deallocate( mumps_R%A   )
		deallocate( mumps_R%RHS )
		mumps_R%JOB = -2
		call DMUMPS(mumps_R)
		
		deallocate( mumps_C%IRN )
		deallocate( mumps_C%JCN )
		deallocate( mumps_C%A   )
		deallocate( mumps_C%RHS )
		mumps_C%JOB = -2
		call ZMUMPS(mumps_C)
	 
		end subroutine radau_ad

end module
