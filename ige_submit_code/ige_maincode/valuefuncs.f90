MODULE valuefuncs
USE GLOBAL
USE old_parent
USE mid_parent
USE sec_parent
USE pri_parent
USE new_parent
USE not_parent
!USE printresults
IMPLICIT NONE

	CONTAINS

!------------------------------------------------------------------
SUBROUTINE vfiter
IMPLICIT NONE

	INTEGER,PARAMETER :: vnmax=20,vn=10!vnmax=20,vn=10!
	INTEGER :: viter,vmiter

	REAL(dp) :: vdist,vdist0
	REAL(dp),DIMENSION(sn,hn,an,cn) :: v4fun1,v4funn1

	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn,cn) :: v9funn,s10funn,n9funn,s4funn
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn,cn) :: w8funn,s9funn, n8funn,n3funn
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn) :: v6funn,v7funn,v8funn
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn) :: s7funn,s8funn,n6funn,n7funn,l1funn,l2funn
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn) :: c3funn
	REAL(dp),DIMENSION(sn,hn,an,an,cn) :: v5funn,s6funn,n5funn,l0funn
	REAL(dp),DIMENSION(sn,hn,an,cn) :: s5funn,n4funn

! retirement is deterministic
	CALL age10retirement

! fix first guess
	CALL age4init

!------------------------------------------------------------------
! now iterate
!------------------------------------------------------------------
	viter = 0
	vmiter= 0
	vdist0 = infty**2	!makes 1st iter always successful
	DO

		CALL age9bequests
		CALL age8after
		CALL age8college
		CALL age7secondary
		CALL age6primary
		CALL age5infant
		CALL age4young
!------------------------------------------------------------------
		vdist = MAXVAL(ABS(v4fun-v4funn))
		IF (myid==0) THEN
			WRITE(*,'(4I4,F16.8)'),caliter,piter,vmiter,viter,vdist
		ENDIF

! IF SUCCESSFUL, COMPLETELY REPLACE V4: first step is always successful
		IF (vdist<vdist0) THEN

			v4fun1 = v4fun
			v4fun  = v4funn
			vdist0 = vdist
			IF (vdist0<vtol .OR. vmiter>=vnmax .OR. viter>=vn) RETURN

			v4funn1= v4funn
							n4funn = n4fun; s5funn = s5fun;
			v5funn = v5fun; n5funn = n5fun; s6funn = s6fun; l0funn = l0fun
			v6funn = v6fun; n6funn = n6fun; s7funn = s7fun; l1funn = l1fun
			v7funn = v7fun; n7funn = n7fun; s8funn = s8fun; l2funn = l2fun
			v8funn = v8fun; n8funn = n8fun; s9funn = s9fun; n3funn = n3fun
			w8funn = w8fun; c3funn = c3fun
			v9funn = v9fun; n9funn = n9fun; s10funn=s10fun; s4funn = s4fun

			viter = viter+1

! OTHERWISE, BACKTRACK
		ELSEIF (vmiter<vnmax) THEN
			v4fun = adjv*v4fun1 + (1d0-adjv)*v4fun
		ELSE

! AT FINAL ITER, SWAP BACK TO PREVIUS BEST ITER IF NOT CONVERGED
			v4fun = v4funn1
							n4fun = n4funn; s5fun = s5funn;
			v5fun = v5funn; n5fun = n5funn; s6fun = s6funn; l0fun = l0funn
			v6fun = v6funn; n6fun = n6funn; s7fun = s7funn; l1fun = l1funn
			v7fun = v7funn; n7fun = n7funn; s8fun = s8funn; l2fun = l2funn
			v8fun = v8funn; n8fun = n8funn; s9fun = s9funn; n3fun = n3funn
			w8fun = w8funn; c3fun = c3funn
			v9fun = v9funn; n9fun = n9funn; s10fun=s10funn; s4fun = s4funn

			RETURN
		ENDIF

		vmiter = vmiter+1

	ENDDO
!------------------------------------------------------------------
RETURN
END SUBROUTINE
END MODULE
