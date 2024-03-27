MODULE grids
USE GLOBAL
IMPLICIT NONE

	INTEGER,PRIVATE :: ci,ai,ei,si,hi

CONTAINS
!------------------------------------------------------------------
SUBROUTINE gen_xegrids ! EXOGENOUS GRIDS
USE random, ONLY: P_CDF_Normal_Ran, P_CDF_Normal_Inverse_Ran
IMPLICIT NONE

!--------------------------
! exogenous shocks: earnings
	DO ci = 1,cn
	DO ei = 1,en
		egrid(ei,ci) =  P_CDF_Normal_Inverse_Ran((ei-0.5d0)/en,mu_e(ci),sig_e(ci)**2)
	ENDDO
	egrid(:,ci) = EXP(egrid(:,ci))	!-sig_e(ci)**2/2d0)
	ENDDO

RETURN
END SUBROUTINE
!------------------------------------------------------------------
SUBROUTINE gen_xagrids ! EXOGENOUS GRIDS
USE rouwenhorst
IMPLICIT NONE

!--------------------------
! exogenous shocks: abilities
	CALL ar1(rh_a,sig_eta,agrid,atrans,astat)
	agrid = EXP(agrid)

! normalize ability shocks so that mean is mu_a
	agrid = agrid/DOT_PRODUCT(astat,agrid)
	agrid = agrid *mu_a

! mean of log normalized agrid
	mloga = DOT_PRODUCT(astat,LOG(agrid))

RETURN
END SUBROUTINE
!------------------------------------------------------------------
! ENDONGENOUS GRIDS
!------------------------------------------------------------------
SUBROUTINE gen_hgrid(hmax)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: hmax

	hmin = SQRT(zero);
	hstep = (hmax-hmin) / (hn-1)**powh
	FORALL (hi=1:hn) hgrid(hi) = hstep*(hi-1)**powh + hmin

END SUBROUTINE
!------------------------------------------------------------------
SUBROUTINE gen_hkgrid(hkmax)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: hkmax

	hkmin = SQRT(zero);
	hkstep = (hkmax-hkmin) / (hn-1)**powh
	FORALL (hi=1:hn) hkgrid(hi) = hkstep*(hi-1)**powh + hkmin

END SUBROUTINE
!------------------------------------------------------------------
SUBROUTINE gen_sgrid(smin,smax)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: smin,smax

	sstep = (smax-smin) / (sn-1)**pows
	FORALL (si=1:sn) sgrid(si) = sstep*(si-1)**pows + smin

END SUBROUTINE
!------------------------------------------------------------------
SUBROUTINE gen_skgrid(skmin,skmax)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: skmin,skmax

	skstep = (skmax-skmin) / (sn-1)**pows
	FORALL (si=1:sn) skgrid(si) = skstep*(si-1)**pows + skmin

END SUBROUTINE
!------------------------------------------------------------------
!------------------------------------------------------------------
SUBROUTINE printgrids
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: ffmt0 = '(10F8.4)'
	CHARACTER(LEN=20),PARAMETER :: ffmt1 = '(10F16.3)'

	IF (myid==0) THEN
		OPEN(22,FILE='output/agrid.out')
		OPEN(33,FILE='output/egrid.out')

		DO ci=1,2
			WRITE(33,ffmt0),egrid(:,ci)
		ENDDO
		DO ai = 1,an
			WRITE(22,ffmt0),agrid(ai),astat(ai),atrans(ai,:),SUM(atrans(ai,:))
		ENDDO
		WRITE(22,ffmt0),DOT_PRODUCT(agrid,astat),SUM(astat)
		CLOSE(22);CLOSE(33)

		OPEN(22,FILE="output/hgrid.out")
		OPEN(33,FILE="output/sgrid.out")
		DO hi = 1,hn
			WRITE(22,ffmt1),hkgrid(hi),hgrid(hi)
		ENDDO
		DO si = 1,sn
			WRITE(33,ffmt1),skgrid(si),sgrid(si)
		ENDDO
		CLOSE(22);CLOSE(33)
	ENDIF
END SUBROUTINE
!------------------------------------------------------------------
END MODULE
