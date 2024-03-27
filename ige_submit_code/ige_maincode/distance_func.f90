MODULE distance_func
USE grids 			! change to equilibrium later
USE simulation, ONLY:fix_proc,mcmc 		! change to equilibrium later
USE valuefuncs	 	! change to equilibrium later
USE compute_iter
USE printresults	!,ONLY: printline,printmoments
IMPLICIT NONE

!	REAL(dp),DIMENSION(7),PARAMETER :: avgweights=(/1d0,1d0,1d0,1d0,1d0,0d0,1d0/)
!	REAL(dp),DIMENSION(7),PARAMETER :: colweights=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
!	REAL(dp),DIMENSION(7),PARAMETER :: stdweights=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0/)

	REAL(dp),DIMENSION(7),PARAMETER :: lowweights=(/1d0,0d0,0d0,0d0,0d0,0d0,1d0/)
	REAL(dp),DIMENSION(7),PARAMETER :: hghweights=(/1d0,0d0,0d0,0d0,0d0,0d0,1d0/)
	REAL(dp),DIMENSION(7),PARAMETER :: stdweights=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0/)

	REAL(dp),DIMENSION(5),PARAMETER :: igcweights=(/1d0,0d0,0d0,0d0,0d0/)
	REAL(dp),DIMENSION(3),PARAMETER :: tmeweights=(/1d0,1d0,1d0/)
!	REAL(dp),DIMENSION(6),PARAMETER :: tstweights=(/0d0,0d0,0d0,0d0,0d0,0d0/)

	CONTAINS
!----------------------------------------------------
FUNCTION func(xx)
IMPLICIT NONE
	REAL(dp),DIMENSION(:),INTENT(IN) :: xx
	REAL(dp),DIMENSION(nop) :: intemp
	REAL(dp),DIMENSION(nout) :: outtemp
	REAL(dp) :: func

	INTEGER :: i
!	REAL(dp) :: dumm,dum3(3)
	CHARACTER(LEN=20),PARAMETER :: afmt = '(20A8)'
	CHARACTER(LEN=20),PARAMETER :: ffmt = '(20F8.3)'

!----------------------------------------------------
	caliter = caliter+1;
	intemp = (pub-plb)*( EXP(xx)/(1d0+EXP(xx)) ) + plb
!----------------------------------------------------
	IF (myid==0) THEN
		CALL printline(6)
		WRITE(*  ,'(A10,I10)'), 'CALITER:',caliter
		WRITE(*  ,afmt),'thet','rho_a','zeta','omeg1','omeg2','bp','zcols','zcolm','mu_a','sig_a','bta','ww','tfp'
		WRITE(*  ,ffmt), intemp
		CALL printline(6)

		CALL printline(100)
		WRITE(100,'(A10,I10)'), 'CALITER:',caliter
		WRITE(100,afmt),'thet','rho_a','zeta','omeg1','omeg2','bp','zcols','zcolm','mu_a','sig_a','bta','ww','tfp'
		WRITE(100,ffmt), intemp
		CALL printline(100)
	ENDIF
!----------------------------------------------------
! set params and other implied params
!---------------------

!	intemp = FLOOR(intemp*1D3+.5d0)/1D3

	thet =	intemp(1)
	rh_a =	intemp(2)
	zeta = intemp(3)
	omega = intemp(4:5)

	gamm =	intemp(6)
!	gamm(1) = gamm(1)-intemp(7)
!	gamm(2) = gamm(2)+intemp(7)

	zcols = intemp(7)
	zcolm =	intemp(8)

	mu_a =  intemp(9)
	sig_a = intemp(10)

!	nu = 1d0
	bta = 	intemp(11)**per
	ww  =	intemp(12)
	tfp =	intemp(13)

! ADDED
	Windx = Windx0 * tfp
	hmax = nmdlavgearn/Windx*hfrac
	CALL gen_hgrid(hmax)
	hkmax = (hmax*hkfrac)	!**nu
	CALL gen_hkgrid(hkmax)
! ADDED


! IMPLIED PARAMS
	qq = (1d0+thet**(1d0/chi))**chi
	sig_eta = sig_a*SQRT(1d0-rh_a**2)

	ups = colshare/(1d0-colshare) *eprem
	ups = ww**ces *ups**(1d0-ces)
	ups = 1d0/(1d0+ups)

! IMPLIED PRICES
	perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0)
	perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0) *perm_cons_denom

!	ww = eprem*colshare/(1d0-colshare)
!	ww = ww**((ces-1d0)/ces) * ((1d0-ups)/ups)**(1d0/ces)

	wskill(1) = (1d0-ups)**(1d0/(1d0-ces)) *ww**(ces/(ces-1d0))
	wskill(1) = ups**(1d0/(1d0-ces)) + wskill(1)
	wskill(1) = Windx *wskill(1)**((1d0-ces)/ces)
	wskill(2) = wskill(1) *ww

	DO i = 0,2
		lamb(i,:) = lamb0(i)/wskill**kamm(i)
	ENDDO
	
!---------------------
! fix grids
!---------------------
	CALL gen_xagrids		! NEEDS MU_A AND SIG_A (I.E. SIG_ETA)
	CALL fix_proc			! NEEDS AGRID,EGRID

!---------------------
! now iterate val func, and solve eqm
!---------------------
	CALL vfiter
	CALL mcmc
	CALL compute_moments

! implied moments
	moments = 1d2*(/ beq,earn12,cinv,eprem1,prof,logvar,r1,enroll,avgearn/nmdlavgearn /)
!---------------------------------------------------------------------------------------
	moments = MAX(zero,moments)

	mavgprof = MAX(zero,mavgprof)
	mpremprof = MAX(zero,mpremprof)
	
	mhscprof = MAX(zero,mhscprof)
	mcolprof = MAX(zero,mcolprof)

	mstdprof = MAX(zero,mstdprof)
	mtimeprof = MAX(zero,mtimeprof)
	migcprof = MAX(zero,migcprof)
!----------------------------
!	outtemp(1) = DOT_PRODUCT(lowweights, LOG(mhscprof/hscprof)**2)
!	outtemp(2) = DOT_PRODUCT(hghweights, LOG(mcolprof/colprof)**2)

	outtemp(1) = DOT_PRODUCT(lowweights, LOG(mavgprof/avgprof)**2)
	outtemp(2) = DOT_PRODUCT(hghweights, LOG(mpremprof/premprof)**2)
	outtemp(3) = DOT_PRODUCT(stdweights, LOG(mstdprof/stdprof)**2)
!----------------------------
!	dumm = SUM(testprof(5:6))/SUM(mtestprof(5:6))
!	mtestprof = mtestprof * dumm
!	dum3 = testprof(1:3)-mtestprof(1:3)
!	mtestprof(1:3) = mtestprof(1:3) + dum3
!	mtestprof = MAX(zero,mtestprof)
!----------------------------
	outtemp(4) = DOT_PRODUCT(tmeweights, LOG(mtimeprof/timeprof)**2)
!	outtemp(4) = DOT_PRODUCT(tstweights, LOG(mtestprof/testprof)**2) + outtemp(4)
	outtemp(5) = DOT_PRODUCT(igcweights, LOG(migcprof/igcprof)**2)
!----------------------------
	outtemp(6) = 	LOG(moments(4)/targets(4))**2
	outtemp(7:9) = 	LOG(moments(9:11)/targets(9:11))**2
!----------------------------
	func = DOT_PRODUCT(outweights,outtemp)/SUM(outweights)	!SQRT(func)*1d2
	func = SQRT(func)*1d2

	IF (func/=func) func = infty

! write out results  !------------------------------------------------------------------
	IF (myid==0) THEN
		CALL printmoments(6)
		CALL printline(6)
		CALL printmoments(100)
		CALL printline(100)

		finish = MPI_wtime()
		WRITE(*  ,'(I10,A10,F16.2,A10,F8.2,A5)'),caliter,'DISTANCE:',func,'TIME:',(finish-start)/60,'MIN'
		WRITE(100,'(I10,A10,F16.2,A10,F8.2,A5)'),caliter,'DISTANCE:',func,'TIME:',(finish-start)/60,'MIN'
	ENDIF
!	func = 0d0

RETURN
END FUNCTION

END MODULE
