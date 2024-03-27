MODULE compute_iter
USE GLOBAL
IMPLICIT NONE

	CONTAINS

!----------------------------------------------------
SUBROUTINE compute_prices
IMPLICIT NONE
!----------------------------------------------------

	INTEGER :: ci,ind
	REAL(dp) :: tlab

! WE HAVE ALL REQUIRED INFO IN HHIDPROC'S
! THE MATRIX IS HUGE; ONLY ALLGATHER THEN WHEN WE ARE DONE
! HERE JUST COMPUTE AGGREGATE MOMENTS BY ALLREDUCE

! TOTAL P. AND H. CAPITAL
	cap = 0d0
	lab = 0d0
	DO ind=1,NI
		cap = SUM(hhid(ind)%ss_k(5:12)) + cap

		IF (hhid(ind)%col_k==1) THEN
			lab(1) = lab(1) + SUM(hhid(ind)%ee_k(3:10))
		ELSE
			lab(2) = lab(2) + SUM(hhid(ind)%ee_k(3:10))
		ENDIF
	ENDDO
	cap = cap/DBLE(NI)
	lab = lab/DBLE(NI)/wskill

! PRICES
	cap = MAX(zero,cap)
	lab = MAX(zero,lab)

	tlab = ups*lab(1)**ces + (1d0-ups)*lab(2)**ces
	tlab = tlab**(1d0/ces)
	gdp = cap**alph * (tfp*tlab)**(1d0-alph)

	r1 = alph*gdp/cap
	r1 = (1d0+r1)**(1d0/per) -1d0 -delk

	ww1 = (1d0-ups)/ups *(lab(2)/lab(1))**(ces-1d0)
	ups1 = colshare/(1d0-colshare) *eprem
	ups1 = ww1**ces *ups1**(1d0-ces)
	ups1 = 1d0/(1d0+ups1)

RETURN
END SUBROUTINE
!----------------------------------------------------
SUBROUTINE compute_moments
IMPLICIT NONE


	REAL(dp) :: bzero,azero,lbzero,lazero

	REAL(dp),DIMENSION(NI,4:10) :: earnmat

	INTEGER :: ci,cn,ind,agei,avgpos

	REAL(dp),DIMENSION(NI,2) :: earnall

	REAL(dp),DIMENSION(2) :: mvar
	INTEGER,DIMENSION(2) :: mcount
	REAL(dp),DIMENSION(3) :: tvar

!----------------------------------------------------

! 1. interest rate
! 2. wage ratio
	CALL compute_prices

! 3. BEQUESTS
! 4. COLLEGE ENROLLMENT
	beq =    SUM(hhid%ss_g(4))/DBLE(NI)/cap

	enroll = SUM(colfk)/DBLE(NI)
if (myid==0) print*,enroll, &
	SUM(colfk,hhid%col_p>1.5d0) / DBLE(COUNT(hhid%col_p>1.5d0)) -enroll

	enroll = DBLE(COUNT(hhid%col_k>1.5d0))/DBLE(NI)

	migcprof(3) = DBLE(COUNT(hhid%col_k>1.5d0 .and. hhid%col_p>1.5d0))
!	migcprof(3) = SUM(colfk,hhid%col_p>1.5d0)
	migcprof(3) = migcprof(3) / DBLE(COUNT(hhid%col_p>1.5d0)) - enroll

!	migcprof(4) = DBLE(COUNT(hhid%col_k>1.5d0 .and. hhid%col_p<1.5d0))
!	migcprof(4) = migcprof(4) / DBLE(COUNT(hhid%col_p<1.5d0))
!	migcprof(3) = migcprof(3)-migcprof(4)

	migcprof(4) = beq
	migcprof(5) = DBLE(COUNT(hhid%ss_g(4)<zbeq))/DBLE(NI)

!------------------
! DROP ZERO OBS
! 5. TIME INVESTMENT
! 6. AVG EARNINGS

!	cinv = SUM(hhid%nn_p(0),hhid%nn_p(0)>zero)/DBLE(COUNT(hhid%nn_p(0)>zero))
!	mtimeprof(1) = cinv
!	DO ci=1,2
!		mtimeprof(ci+1) = SUM(hhid%nn_p(ci),hhid%nn_p(ci)>zero) / DBLE(COUNT(hhid%nn_p(ci)>zero))
!	ENDDO

	avgearn= SUM(hhid%ee_k(10),hhid%ee_k(10)>=cuta30)
	avgpos = COUNT(hhid%ee_k(10)>=cuta30)
	DO ind = 1,NI
		avgearn = avgearn + MERGE(hhid(ind)%ee_k(4 ),0D0,hhid(ind)%ee_k(4 )>=cutb30 .and. hhid(ind)%nn_k(4)<1d0-timeb30)
!		avgearn = avgearn + MERGE(hhid(ind)%ee_k(10),0D0,hhid(ind)%ee_k(10)>=cuta30)

		avgearn = avgearn + SUM(hhid(ind)%ee_k(5:7),hhid(ind)%ee_k(5:7)>=cuta30 .and. hhid(ind)%nn_k(5:7) + hhid(ind)%nn_k(0:2) < 1d0-timea30)
		avgearn = avgearn + SUM(hhid(ind)%ee_k(8:9),hhid(ind)%ee_k(8:9)>=cuta30 .and. hhid(ind)%nn_k(8:9) < 1d0-timea30)

		avgpos = avgpos + MERGE(1,0,hhid(ind)%ee_k(4 )>=cutb30 .and. hhid(ind)%nn_k(4)<1d0-timeb30)
!		avgpos = avgpos + MERGE(1,0,hhid(ind)%ee_k(10)>=cuta30)

		avgpos = avgpos + COUNT(hhid(ind)%ee_k(5:7)>=cuta30 .and. hhid(ind)%nn_k(5:7) + hhid(ind)%nn_k(0:2) < 1d0-timea30)
		avgpos = avgpos + COUNT(hhid(ind)%ee_k(8:9)>=cuta30 .and. hhid(ind)%nn_k(8:9) < 1d0-timea30)
	ENDDO
	avgearn = avgearn/MAX(1d0,DBLE(avgpos))
!------------------
	bzero = cutb30 *avgearn
	azero = cuta30 *avgearn

	lbzero = LOG(bzero)
	lazero = LOG(azero)
!------------------

!----
!----
! DON'T DROP
	cinv = SUM(hhid%nn_p(0))/DBLE(NI)
	mtimeprof(1) = cinv
	DO ci=1,2
		mtimeprof(ci+1) = SUM(hhid%nn_p(ci))/DBLE(NI)
	ENDDO
!	mtimeprof(1) = mtimeprof(1)/mtimeprof(2)
!	mtimeprof(3) = mtimeprof(3)/mtimeprof(2)

!	avgearn = 0d0
!	DO ind = 1,NI
!		avgearn = avgearn + SUM(hhid(ind)%ee_p(3:10))/8d0
!	ENDDO
!	avgearn = avgearn/DBLE(NI)
!
!	bzero = zero
!	azero = zero
!
!	lbzero = LOG(bzero)
!	lazero = LOG(azero)
!----
!----

!------------------
! 7. COLLEGE PREMIUM: MUST USE BENCHMARK, NOT CALIBRATED WAGE RAT!
! ALSO KEEP ZERO OBS, SINCE THIS IS EQM VALUE

!	eprem1 = lab(2)/lab(1) *ww *(1d0-enroll)/enroll
!-----------------
! NOT! DON'T KEEP ZERO OBS?

	DO ind = 1,NI
		earnmat(ind,4:10) = LOG(MAX(zero,hhid(ind)%ee_k(4:10)))

		earnmat(ind,4) = MERGE(earnmat(ind,4),LOG(zero),hhid(ind)%nn_k(4)<1d0-timeb30)
		DO agei = 5,7
			earnmat(ind,agei) = MERGE(earnmat(ind,agei),LOG(zero),hhid(ind)%nn_k(agei)+hhid(ind)%nn_k(agei-5)<1d0-timea30)
		ENDDO
		DO agei = 8,9
			earnmat(ind,agei) = MERGE(earnmat(ind,agei),LOG(zero),hhid(ind)%nn_k(agei)<1d0-timea30)
		ENDDO
	ENDDO
!-----------------

! 8. EDUC PROFILES
! 9.,10. LOG RESIDUAL EARNINGS VARIANCE AT STAGES 4/10


	mcount(1) = COUNT(hhid%col_k==1 .and. earnmat(:,4)>=lbzero)
	mcount(2) = COUNT(hhid%col_k==2 .and. earnmat(:,4)>=lbzero)

	mhscprof(1) = SUM( earnmat(:,4),hhid%col_k==1 .and. earnmat(:,4)>=lbzero) / MAX(1d0,DBLE(mcount(1)))
	mcolprof(1) = SUM( earnmat(:,4),hhid%col_k==2 .and. earnmat(:,4)>=lbzero) / MAX(1d0,DBLE(mcount(2)))

	mvar(1) = SUM( (earnmat(:,4)-mhscprof(1))**2,hhid%col_k==1 .and. earnmat(:,4)>=lbzero)
	mvar(2) = SUM( (earnmat(:,4)-mcolprof(1))**2,hhid%col_k==2 .and. earnmat(:,4)>=lbzero)

	mstdprof(1) = SUM(mvar)/MAX(1d0,DBLE(SUM(mcount)))
	DO ci = 5,10
		mcount(1) = COUNT(hhid%col_k==1 .and. earnmat(:,ci)>=lazero)
		mcount(2) = COUNT(hhid%col_k==2 .and. earnmat(:,ci)>=lazero)

		mhscprof(ci-3) = SUM( earnmat(:,ci),hhid%col_k==1 .and. earnmat(:,ci)>=lazero) / MAX(1d0,DBLE(mcount(1)))
		mcolprof(ci-3) = SUM( earnmat(:,ci),hhid%col_k==2 .and. earnmat(:,ci)>=lazero) / MAX(1d0,DBLE(mcount(2)))

		mvar(1) = SUM( (earnmat(:,ci)-mhscprof(ci-3))**2,hhid%col_k==1 .and. earnmat(:,ci)>=lazero)
		mvar(2) = SUM( (earnmat(:,ci)-mcolprof(ci-3))**2,hhid%col_k==2 .and. earnmat(:,ci)>=lazero)

		mstdprof(ci-3) = SUM(mvar)/MAX(1d0,DBLE(SUM(mcount)))
	ENDDO
	mstdprof = SQRT(mstdprof)

	eprem1 = EXP(SUM(mcolprof-mhscprof)/7)
	prof(1) = EXP(mhscprof(6)-mhscprof(1))
	prof(2) = EXP(mcolprof(6)-mcolprof(1))

	logvar(1) = mstdprof(1)
	logvar(2) = mstdprof(7)

	mavgprof = enroll*mcolprof + (1d0-enroll)*mhscprof
	mavgprof = EXP(mavgprof-mavgprof(6))

	mpremprof = EXP(mcolprof-mhscprof)
	mcolprof = EXP(mcolprof-mhscprof(6))
	mhscprof = EXP(mhscprof-mhscprof(6))

! 11. IGE: rank-rank

!!!! CHECK THIS CODE !!!!

!i. lifetime average
	DO ind = 1,NI
		earnall(ind,1) = SUM(hhid(ind)%ee_p(4:10))/7d0 !+ rr*SUM(hhid(ind)%ss_p(7:8))/2d0
		earnall(ind,2) = SUM(hhid(ind)%ee_k(4:10))/7d0 !+ rr*SUM(hhid(ind)%ss_p(7:8))/2d0
	ENDDO
	CALL compute_ige(earnall,migcprof(1),migcprof(2),1d2/psidavgearn)

!	migcprof(2:5) = 0d0
!!ii. parent age 42-53, child age 30-41
!	DO ind = 1,NI
!		earnall(ind,1) = SUM(hhid(ind)%ee_p(7:8))/2d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!		earnall(ind,2) = SUM(hhid(ind)%ee_k(5:5))/2d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!	ENDDO
!	CALL compute_ige(earnall,migcprof(2))
!
!!iii. parent age 42-65, child age 36
!	DO ind = 1,NI
!		earnall(ind,1) = SUM(hhid(ind)%ee_p(7:9))/3d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!	ENDDO
!	CALL compute_ige(earnall,migcprof(3))
!
!!iv. parent age 42-53, child age 42
!	DO ind = 1,NI
!		earnall(ind,1) = SUM(hhid(ind)%ee_p(7:10))/4d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!	ENDDO
!	CALL compute_ige(earnall,migcprof(4))
!
!!iv. parent age 42-53, child age 36-47
!	DO ind = 1,NI
!		earnall(ind,1) = SUM(hhid(ind)%ee_p(7:8))/2d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!		earnall(ind,2) = SUM(hhid(ind)%ee_k(6:6))/2d0 !+ rr*SUM(hhid(ind)%ss_p(7:9))/3d0
!	ENDDO
!	CALL compute_ige(earnall,migcprof(5))

! 12 . testscore stats
!
!	DO ci = 1,3
!		tvar(ci) = SUM( LOG(hhid%hh_k(ci)) )/NI
!	ENDDO
!	mtestprof(1:3) = tvar
!
!	DO ci = 1,3
!		mtestprof(ci+3) = SUM( (LOG(hhid%hh_k(ci))-tvar(ci))**2 )/NI
!	ENDDO
!	mtestprof(4:6) = SQRT(mtestprof(4:6))
!
!	tvar = 0d0
!	DO ci=1,3
!		DO ind = 1,NI
!			tvar(ci) = hhid(ind)%nn_p(ci-1)*hhid(ind)%hh_p(ci+4) *wskill(hhid(ind)%col_p) + tvar(ci)
!		ENDDO
!	ENDDO
!	tvar = tvar/DBLE(NI)
!
!------------------
!if (myid==0) then
!print*,"mean time cost"
!print*,tvar
!endif
!------------------

RETURN
END SUBROUTINE
!----------------------------------------------------
SUBROUTINE compute_check
IMPLICIT NONE
	INTEGER :: agei,hi,si
	REAL(dp) :: hbinsproc(hn+1),hkbinsproc(hn+1)
	REAL(dp) :: sbinsproc(sn+1),skbinsproc(sn+1),s11,s12

!-----------------
! count hitting max

	hkbinsproc= 0d0
	hbinsproc = 0d0
	DO hi=1,hn
		DO agei=1,2
			hkbinsproc(hi)= COUNT(hhidproc%hh_k(agei)<hkgrid(hi))+ hkbinsproc(hi)
		ENDDO
		DO agei=3,10
			hbinsproc(hi) = COUNT(hhidproc%hh_k(agei)<hgrid(hi)) + hbinsproc(hi)
		ENDDO
	ENDDO
	hkbinsproc= hkbinsproc/DBLE(NI)/2d0
	hbinsproc = hbinsproc /DBLE(NI)/8d0

	hkbinsproc(hn+1)= 1d0/nprocs
	hbinsproc(hn+1) = 1d0/nprocs
	hkbinsproc(2:hn+1) =hkbinsproc(2:hn+1)-hkbinsproc(1:hn)
	hbinsproc(2:hn+1) = hbinsproc(2:hn+1) -hbinsproc(1:hn)

	skbinsproc= 0d0
	sbinsproc = 0d0
	DO si=1,sn
			skbinsproc(si)= COUNT(hhidproc%ss_k(4)   <skgrid(si))+ skbinsproc(si)
		DO agei=5,10
			sbinsproc(si) = COUNT(hhidproc%ss_k(agei)<sgrid(si)) + sbinsproc(si)
		ENDDO
	ENDDO
	skbinsproc= skbinsproc/DBLE(NI)
	sbinsproc = sbinsproc /DBLE(NI)/6d0

	skbinsproc(sn+1)= 1d0/nprocs
	sbinsproc(sn+1) = 1d0/nprocs
	skbinsproc(2:sn+1) =skbinsproc(2:sn+1)-skbinsproc(1:sn)
	sbinsproc(2:sn+1) = sbinsproc(2:sn+1) -sbinsproc(1:sn)

!-----------------
	CALL MPI_ALLREDUCE(hbinsproc,hbins,hn+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLREDUCE(sbinsproc,sbins,sn+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLREDUCE(hkbinsproc,hkbins,hn+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLREDUCE(skbinsproc,skbins,sn+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!-----------------
END SUBROUTINE
!----------------------------------------------------
SUBROUTINE compute_ige(earnall,ppigc,earn12,mzero)
USE qsort_module
USE ols
IMPLICIT NONE
	REAL(dp),DIMENSION(NI,2),INTENT(IN) :: earnall
	REAL(dp),INTENT(OUT) :: ppigc
	REAL(dp),OPTIONAL,INTENT(INOUT) :: earn12
	REAL(dp),OPTIONAL,INTENT(IN)  :: mzero
	REAL(dp),DIMENSION(NI,2) :: earngen

	REAL(dp),DIMENSION(100,2) :: Ek_slope

	INTEGER :: i,ind,ind0	,ipctl1
	REAL(dp) :: pctl1

! i indexes generations
	earngen = earnall
	DO i = 2,1,-1
		CALL qsort(earngen,i)
		DO ind0 = 1,NI
			IF (earngen(ind0,i)>=zero) EXIT
		ENDDO
		earngen(:,i) = (/ (DBLE(ind),ind=1,NI) /) / DBLE(NI)*1d2
		earngen(:,i) = DBLE(CEILING(earngen(:,i)))

		pctl1 = DBLE(ind0)/DBLE(NI)*1d2
		IF (pctl1>1d0) earngen(1:ind0,i) = pctl1/2d0
	ENDDO

!	ppigc = regress(earngen(:,1),earngen(:,2))

!------------------------
	Ek_slope(:,1) = (/ (DBLE(i),i=1,100) /)
	DO i = 1,100
		Ek_slope(i,2) = SUM(earngen(:,2),earngen(:,1)==DBLE(i))/DBLE(COUNT(earngen(:,1)==DBLE(i)))
	ENDDO

	ipctl1 = MAX0(1,FLOOR(pctl1))
	IF (pctl1>1d0) THEN
		Ek_slope(1:ipctl1,1) = pctl1/2d0
		Ek_slope(1:ipctl1,2) = SUM(earngen(:,2),earngen(:,1)==pctl1/2d0)/DBLE(COUNT(earngen(:,1)==pctl1/2d0))
	ENDIF

	ppigc = regress(Ek_slope(MAX0(1,ipctl1):,1),Ek_slope(MAX0(1,ipctl1):,2))
!------------------------

	IF (.NOT. PRESENT(earn12)) RETURN

! schooling slope
	earngen(:,1) = earnall(:,1)
	earngen(:,2) = DBLE(hhid%col_k-1)*1d2
	CALL qsort(earngen,1)
	earngen(:,1) = (/ (DBLE(ind),ind=1,NI) /) / DBLE(NI)*1d2
	earngen(:,1) = DBLE(CEILING(earngen(:,1)))

	IF (pctl1>1d0) earngen(1:ind0,1) = pctl1/2d0

!	earn12 = regress(earngen(:,1),earngen(:,2))

!------------------------
	DO i = 1,100
		Ek_slope(i,2) = SUM(earngen(:,2),earngen(:,1)==DBLE(i))/DBLE(COUNT(earngen(:,1)==DBLE(i)))
	ENDDO
	IF (pctl1>1d0) THEN
		Ek_slope(1:ipctl1,2) = SUM(earngen(:,2),earngen(:,1)==pctl1/2d0)/DBLE(COUNT(earngen(:,1)==pctl1/2d0))
	ENDIF

	earn12 = regress(Ek_slope(MAX0(1,ipctl1):,1),Ek_slope(MAX0(1,ipctl1):,2))
!------------------------


!	earngen = earnall
!	DO ind=1,NI
!		IF (earngen(ind,1)<=mzero .OR. earngen(ind,2)<=mzero) earngen(ind,:) = 0d0
!	ENDDO
!
!	CALL qsort(earngen,1)
!	DO ind0 = 1,NI
!		IF (earngen(ind0,1)>mzero) EXIT
!	ENDDO
!	ind0 = MAX(ind0,NI/10)
!	earn12 = regress(LOG(earngen(ind0+1:NI,2)),LOG(earngen(ind0+1:NI,1)))

	RETURN


!-----------------
END SUBROUTINE
!----------------------------------------------------
END MODULE
