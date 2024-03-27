MODULE bt_module
USE GLOBAL
IMPLICIT NONE

	PRIVATE

	CHARACTER(LEN=20),PARAMETER :: pfmt = '(20F9.4)'
	INTEGER,PARAMETER :: bt_nop = 3, calmax = 1000
	REAL(dp) :: bt_plb(bt_nop)=(/1d-2,1d-1,1d-1/)
	REAL(dp) :: bt_pub(bt_nop)=(/0.9d0,5d0,.5d0/)

	REAL(dp),PARAMETER :: ltax=0.118d0
!	REAL(dp),PARAMETER :: gbar = 0.81375541962086451d0
!	REAL(dp),PARAMETER :: thbar = 0.32939718057342166d0

	REAL(dp) :: gbar,zbar,thbar
	REAL(dp) :: xmean,hmean,smean
	REAL(dp) :: esub,budget,ak,hstar

	REAL(dp) :: resource(NI)
	REAL(dp) :: hparent(NI),sparent(NI)
	REAL(dp) :: xchild(NI),hchild(NI),schild(NI)

	INTEGER :: bt_iter

	PUBLIC :: bt_calibrate

CONTAINS
!-----------------------------------------------------------
SUBROUTINE bt_calibrate
USE random
USE simplex
IMPLICIT NONE
	REAL(dp),DIMENSION(bt_nop) :: bt_params,bt_guess,bt_web,step
	REAL(dp) :: bt_distance,bt0,add
	INTEGER :: agei,ind,iflg

	INTEGER :: ii,jj,kk
	REAL(dp),PARAMETER :: pwgts(3) = (/0.25d0,0.5d0,0.75d0/)


! targets & inputs for simple model
! discounted educ subs

	DO ind=1,NI
		hparent(ind) = hhid(ind)%hh_p(4)*(1d0-hhid(ind)%nn_p(4))

		DO agei = 5,9
			hparent(ind) = hhid(ind)%hh_p(agei)*(1d0-hhid(ind)%nn_p(agei))/(1d0+rr)**(agei-4) + hparent(ind)
		ENDDO
		hparent(ind) = hhid(ind)%hh_p(10)/(1d0+rr)**(10-4) + hparent(ind)
		hparent(ind) = hparent(ind)*wskill(hhid(ind)%col_p)
	ENDDO
	hmean = SUM(hparent)/NI

	esub = d012(1)/(1d0+rr) + d012(2)/(1d0+rr)**2 + d012(3)/(1d0+rr)**3
	esub = esub*hmean

! absolute bequest amounts
	xmean = 0d0
	DO ind = 1,NI
!		xmean = wskill(hhid(ind)%col_k)*hhid(ind)%nn_k(3) + xmean
!		xmean = MERGE(0D0,kappa,hhid(ind)%col_k==1) + xmean
		DO jj = 0,2
			add = hhid(ind)%nn_p(jj) * hhid(ind)%hh_p(jj+5)/kamm(jj)/(1d0+rr)**jj
			add = wskill(hhid(ind)%col_p)* add

			xmean = add+xmean
		ENDDO
	ENDDO
	xmean = xmean/(1d0+rr)/NI
	smean = SUM(hhid%ss_k(4))/NI

	resource = (1d0-ltax)*hparent + hhid%ss_p(4)

! now calibrate

	OPEN(210,FILE="results/bt_cal.out")
		bt_distance = infty

		bt_iter = 0
		DO ii=1,3
			bt_web(1) = bt_plb(1)*pwgts(4-ii) + bt_pub(1)*pwgts(ii)
		DO jj=1,3
			bt_web(2) = bt_plb(2)*pwgts(4-jj) + bt_pub(2)*pwgts(jj)
		DO kk=1,3
			bt_web(3) = bt_plb(3)*pwgts(4-kk) + bt_pub(3)*pwgts(kk)

			bt_guess = LOG(bt_web-bt_plb) - LOG(bt_pub-bt_web)
			bt0 = bt_model(bt_guess)
			IF (bt0<bt_distance) THEN
				bt_distance = bt0
				bt_params = bt_guess
			ENDIF
		ENDDO;ENDDO;ENDDO

		bt_iter = 0
		step = 1d0
		CALL NMsimplex(bt_params,bt_distance,step,tol,bt_model,itmax,iflg)

		bt_iter = 0
		bt_distance = bt_model(bt_params)
		bt_params = (bt_pub-bt_plb)*( EXP(bt_params)/(1d0+EXP(bt_params)) ) + bt_plb

	WRITE(210,*),bt_distance
	CLOSE(210)

	OPEN(1000,FILE="results/bt_dist.out",STATUS='replace')
		WRITE(1000,'(12A9)'),'xchild','hchild','schild','hparent','resource'
		DO ind=1,NI
			WRITE(1000,pfmt),xchild(ind),hchild(ind),schild(ind),hparent(ind),resource(ind)
		ENDDO
	CLOSE(1000)

	CALL bt_ss

RETURN
END SUBROUTINE
!-----------------------------------------------------------
SUBROUTINE bt_ss
USE random
USE brent, ONLY: mybrent
IMPLICIT NONE
	INTEGER :: akivec(NI),apivec(NI)

	REAL(dp) :: cons1,hmean,hmean1,smean,smean1
	REAL(dp) :: xstar,sstar	!,esub0

	INTEGER,PARAMETER :: maxit=100
	INTEGER :: ind,iflt,sim
	INTEGER :: aii,aki

	cons1 = thbar*(1d0+rr)**5
	cons1 = cons1**(1d0/chi)

	CALL Set_Seed(123456789+myid,362436069+myid,521288629+myid,916191069+myid)
	! SOLVE ALL FAMILIES' EULER EQUATIONS

!	esub0 = esub
	hmean1 = SUM(hchild)/NI
	smean1 = SUM(schild)/NI

	akivec = hhid%ability_k

	DO sim=1,TT
		hmean = hmean1
		smean = smean1
!		esub = esub0*hmean

		apivec = akivec
		hparent = hchild
		sparent = schild
		resource = (1d0-ltax)*hparent + sparent
		DO ind=1,NI

			ak = sample_uniform(0d0,1d0)
			DO aii=1,an
				IF ( ak<=SUM(atrans(apivec(ind),:aii)) ) THEN
					akivec(ind) = aii
					EXIT
				ENDIF
			ENDDO
			ak = agrid(akivec(ind))

			xstar = gbar*zbar*ak / (1d0+rr)**5
			xstar = xstar**(1d0/(1d0-gbar))
			xstar = MAX(0d0,xstar-esub)

			budget = resource(ind)-(1d0-ltax)*xstar
			hstar = zbar*ak*(xstar+esub)**gbar

			sstar = cons1*budget - (1d0-ltax)*hstar
			sstar = sstar / (1d0 + cons1/(1d0+rr)**5)

			IF (sstar<=0d0) THEN
				budget = resource(ind)
				xstar = mybrent(euler,zero,budget,tol,maxit,iflt)
				hstar = zbar*ak*(xstar+esub)**gbar
				sstar = 0d0
			ENDIF

			xchild(ind) = xstar
			hchild(ind) = hstar
			schild(ind) = sstar

		ENDDO

		hmean1 = SUM(hchild)/NI
		smean1 = SUM(schild)/NI

		IF (MOD(sim,10)==0) WRITE(*,'(I4,5F10.4)'),sim,hmean,hmean1,smean,smean1
		IF (MAX(ABS(smean-smean1),ABS(hmean-hmean1))<tol) EXIT

	ENDDO
	PRINT*,'BT_SS DONE!'

	OPEN(1000,FILE="results/bt_ss.out",STATUS='replace')
		WRITE(1000,'(12A9)'),'abp_ss','hparent','sparent','abk_ss','xchild','hchild','schild'
		DO ind=1,NI
			WRITE(1000,pfmt),agrid(apivec(ind)),hparent(ind),sparent(ind),agrid(akivec(ind)),xchild(ind),hchild(ind),schild(ind)
		ENDDO
	CLOSE(1000)


RETURN
END SUBROUTINE
!-----------------------------------------------------------
FUNCTION bt_model(bt_inparm)
USE brent, ONLY: mybrent
IMPLICIT NONE
	REAL(dp),DIMENSION(:),INTENT(IN) :: bt_inparm
	REAL(dp),DIMENSION(bt_nop) :: bt_rparm
	REAL(dp) :: bt_model

	REAL(dp) :: cons1
	REAL(dp) :: xstar,sstar
	REAL(dp) :: avgx,avgh,avgs

	INTEGER,PARAMETER :: maxit=100
	INTEGER :: ind,iflt

	bt_rparm = (bt_pub-bt_plb)*( EXP(bt_inparm)/(1d0+EXP(bt_inparm)) ) + bt_plb

	gbar = bt_rparm(1)
	zbar = bt_rparm(2)
	thbar = bt_rparm(3)

	cons1 = thbar*(1d0+rr)**5
	cons1 = cons1**(1d0/chi)

	! SOLVE ALL FAMILIES' EULER EQUATIONS
	DO ind=1,NI
		ak = agrid(hhid(ind)%ability_k)

		xstar = gbar*zbar*ak / (1d0+rr)**5
		xstar = xstar**(1d0/(1d0-gbar))
		xstar = MAX(0d0,xstar-esub)

		budget = resource(ind)-(1d0-ltax)*xstar
		hstar = zbar*ak*(xstar+esub)**gbar

		sstar = cons1*budget - (1d0-ltax)*hstar
		sstar = sstar / (1d0 + cons1/(1d0+rr)**5)

		IF (sstar<=0d0) THEN
			budget = resource(ind)
			xstar = mybrent(euler,zero,budget,tol,maxit,iflt)
			hstar = zbar*ak*(xstar+esub)**gbar
			sstar = 0d0
		ENDIF

		xchild(ind) = xstar
		hchild(ind) = hstar
		schild(ind) = sstar
	ENDDO

	! MOMENTS

	avgx = SUM(xchild)/NI
	avgh = SUM(hchild)/NI
	avgs = SUM(schild)/NI

	bt_model = 0d0
	bt_model = (avgx/xmean-1d0)**2 + bt_model
	bt_model = (avgh/hmean-1d0)**2 + bt_model
	bt_model = (avgs/smean-1d0)**2 + bt_model

	bt_iter = bt_iter+1

	WRITE(210,'(I4,3F8.3)'),bt_iter,bt_rparm
	WRITE(210,'(I4,6F8.3)'),bt_iter,avgx,avgh,avgs
	WRITE(210,'(I4,6F8.3)'),bt_iter,xmean,hmean,smean,migcprof(1)
	WRITE(210,*),'-------------------'

!	WRITE(*,'(I4,3F8.3)'),bt_iter,bt_rparm
!	WRITE(*,'(I4,6F8.3)'),bt_iter,avgk,ppik,beqk
!	WRITE(*,'(I4,6F8.3)'),bt_iter,kmean,migcprof(1),smean
!	WRITE(*,*),'-------------------'

RETURN
END FUNCTION
!-----------------------------------------------------------
FUNCTION euler(xx)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: xx
	REAL(dp) :: euler,lfs,cp,ck,mph

	cp = budget - (1d0-ltax)*xx
	ck = (1d0-ltax)*zbar*ak*(xx+esub)**gbar
	mph = gbar*zbar*ak/(xx+esub)**(1d0-gbar)

	euler = cp**(-chi) - thbar*mph*ck**(-chi)

RETURN
END FUNCTION
!-----------------------------------------------------------
END MODULE






!!-----------------------------------------------------------
!FUNCTION bt_ige(earnall)
!USE qsort_module
!USE ols
!IMPLICIT NONE
!	REAL(dp),DIMENSION(NI,2),INTENT(IN) :: earnall
!	REAL(dp) :: bt_ige
!
!	REAL(dp),DIMENSION(NI,2) :: earngen
!	REAL(dp),DIMENSION(100,2) :: Ek_slope
!
!	INTEGER :: i,ind,ind0	,ipctl1
!	REAL(dp) :: pctl1
!
!! i indexes generations
!	earngen = earnall
!	DO i = 2,1,-1
!		CALL qsort(earngen,i)
!		DO ind0 = 1,NI
!			IF (earngen(ind0,i)>=zero) EXIT
!		ENDDO
!		earngen(:,i) = (/ (DBLE(ind),ind=1,NI) /) / DBLE(NI)*1d2
!		earngen(:,i) = DBLE(CEILING(earngen(:,i)))
!
!		pctl1 = DBLE(ind0)/DBLE(NI)*1d2
!		IF (pctl1>1d0) earngen(1:ind0,i) = pctl1/2d0
!	ENDDO
!
!!	ppigc = regress(earngen(:,1),earngen(:,2))
!
!!------------------------
!	Ek_slope(:,1) = (/ (DBLE(i),i=1,100) /)
!	DO i = 1,100
!		Ek_slope(i,2) = SUM(earngen(:,2),earngen(:,1)==DBLE(i))/DBLE(COUNT(earngen(:,1)==DBLE!(i)))
!	ENDDO
!
!	ipctl1 = MAX0(1,FLOOR(pctl1))
!	IF (pctl1>1d0) THEN
!		Ek_slope(1:ipctl1,1) = pctl1/2d0
!		Ek_slope(1:ipctl1,2) = SUM(earngen(:,2),earngen(:,1)==pctl1/2d0)/DBLE(COUNT(earngen!(:,1)==pctl1/2d0))
!	ENDIF
!
!	bt_ige = regress(Ek_slope(MAX0(1,ipctl1):,1),Ek_slope(MAX0(1,ipctl1):,2))
!
!RETURN
!END FUNCTION
