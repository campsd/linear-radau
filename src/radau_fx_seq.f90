! ------------------------------------------------------------------------------------------------------
! radau_fx
! ------------------------------------------------------------------------------------------------------
! Implementation of IRK integrator radauIIA with a fixed stepsize and sequential approach
! ------------------------------------------------------------------------------------------------------
! daan.camps@cs.kuleuven.be
! Nov. 2016
! v0.1
! ------------------------------------------------------------------------------------------------------
module radau_fx_mod
	use units
	use utilities
	use radau_common

	implicit none
	private
	save
	public	:: radau_fx, radau_fx_int
	contains
		subroutine radau_fx( A, RI, CJ, N, NNZ, YZERO, TEND, H, YFINAL)
		! ---------------------------------------------------------------------------------------
		! radau_fx( A, RI, CJ, N, NNZ, YZERO, TEND, H, YFINAL )
		! ---------------------------------------------------------------------------------------
		! IN:
		! A(NNZ)	nonzero entries of the Bateman matrix
		! RI(NNZ)	row indices of nonzero entries
		! CJ(NNZ)	col indices of nonzero entries
		! N		size of the system
		! NNZ		number of nonzero entries
		! YZERO(N)	initial conditions of system
		! TEND		end time of integration	(implicitly assumed TSTART=0 & TEND%H=0)
		! H		timestep
		! OUT:
		! YFINAL(N)	result of the integration at TEND
		! 
		! NOTE:
		! IT IS ASSUMED TEND IS DIVISBLE BY H, IF NOT, YFINAL WILL NOT EXACTLY BE AT TEND BUT AT
		! THE CLOSEST SMALLER MULTIPLE OF H
		! ---------------------------------------------------------------------------------------
			include 'mpif.h'
      			include 'dmumps_struc.h'
			include 'zmumps_struc.h'
			
			! declare all variables
			! I/O
			integer, intent(in)	::	N, NNZ			
			real(wp), intent(in)	::	A(NNZ), YZERO(N)
			integer, intent(in)	::	RI(NNZ), CJ(NNZ)
			real(wp), intent(in)	::	TEND, H
			real(wp), intent(out)	::	YFINAL(N)
			! subroutine specific
			integer 			::	IERR, NNZNR, NNZNC, NSTEPS, i, j, world_id, world_size
			real(wp), allocatable		::	ANR(:)
			complex(wp), allocatable	::	ANC(:)
			integer, allocatable		::	RINR(:), CJNR(:), RINC(:), CJNC(:)	
			real(wp)			::	YPREV(N), R1(N), R2(N), R3(N), W1(N), W2(N)
                        integer p, id, cnt

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
			call constructRealN(A, RI, CJ, N, NNZ, H, ANR, RINR, CJNR, NNZNR)
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
			call constructComplexN(A, RI, CJ, N, NNZ, H, ANC, RINC, CJNC, NNZNC)
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
			NSTEPS = int(TEND/H)
			
			do i = 1, NSTEPS
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
				YPREV = YPREV + T31 * mumps_R%RHS + real(mumps_C%RHS)
			enddo
			! Set final solution                        
			YFINAL = YPREV	 

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

		end subroutine radau_fx

		subroutine radau_fx_int(A, RI, CJ, N, NNZ, YZERO, H, STEPS, NS, YINT, YFINAL)
		! ---------------------------------------------------------------------------------------
		! radau_fx_int(A, RI, CJ, N, NNZ, YZERO, STEPS, H, Y, YFINAL)
		! ---------------------------------------------------------------------------------------
		! IN:
		! A(NNZ)	nonzero entries of the Bateman matrix
		! RI(NNZ)	row indices of nonzero entries
		! CJ(NNZ)	col indices of nonzero entries
		! N		size of the system
		! NNZ		number of nonzero entries
		! YZERO(N)	initial conditions of system at t=0
		! H		stepsize
		! STEPS(NS)	(ordered) list of timesteps for which intermediate results are requested
		!		the last entry is the final integration stepsize
		! NS 		number of results requested
		!
		! OUT:
		! YINT(N, NS-1)	intermediate results
		! YFINAL(N)	final result
		! 
		! NOTE:
		! IN CONTRAST TO THE FIRST VERSION WE COUNT BY STEPS AND NOT BY TIME
		! ---------------------------------------------------------------------------------------

			include 'mpif.h'
      			include 'dmumps_struc.h'
			include 'zmumps_struc.h'

			! declare all variables
			! I/O
			integer, intent(in)	::	N, NNZ, NS
			real(wp), intent(in)	::	A(NNZ), YZERO(N)
			integer, intent(in)	::	RI(NNZ), CJ(NNZ), STEPS(NS)
			real(wp), intent(in)	::	H
			real(wp), intent(out)	::	YINT(N, NS-1)
			real(wp), intent(out)	::	YFINAL(N)
			
			! subroutine specific
			integer 			::	IERR, NNZNR, NNZNC, i, j, k, world_id, world_size
			real(wp), allocatable		::	ANR(:)
			complex(wp), allocatable	::	ANC(:)
			integer, allocatable		::	RINR(:), CJNR(:), RINC(:), CJNC(:)	
			real(wp)			::	YPREV(N), R1(N), R2(N), R3(N), W1(N), W2(N)
                        integer p, id, cnt

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
			call constructRealN(A, RI, CJ, N, NNZ, H, ANR, RINR, CJNR, NNZNR)
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
			call constructComplexN(A, RI, CJ, N, NNZ, H, ANC, RINC, CJNC, NNZNC)
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
			k = 1
			do i = 1, STEPS(NS)
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
				YPREV = YPREV + T31 * mumps_R%RHS + real(mumps_C%RHS)

				! Check if result needs to be stored
				if (i .eq. STEPS(k) .and. i .lt. STEPS(NS)) then
					YINT(1:N,k) = YPREV
					k = k + 1
				endif
			enddo
			! Set final solution                        
			YFINAL = YPREV

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

		end subroutine radau_fx_int
end module
		
