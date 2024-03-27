MODULE functions
USE GLOBAL
IMPLICIT NONE

CONTAINS
!-----------------------------------
REAL(dp) FUNCTION du(cx)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: cx

	IF (cx < zero) THEN
		du = zero**(-chi) - (cx-zero)*chi*zero**(-1d0-chi)
		RETURN
	ENDIF
	du = cx**(-chi)
	du = du *(1d0-bta) *(1d0-bta**5*thet) /(1d0-bta**9)

RETURN
END FUNCTION
!-----------------------------------
REAL(dp) FUNCTION u(cx)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: cx

	IF (cx < zero) THEN
		u = (zero**(1d0-chi))/(1d0-chi) + (cx-zero)*du(zero)
		RETURN
	ENDIF
	u = cx**(1d0-chi)/(1d0-chi)
	u = u *(1d0-bta) *(1d0-bta**5*thet) /(1d0-bta**9)

RETURN
END FUNCTION
!-----------------------------------
REAL(dp) FUNCTION taxrate(yx)
! relative earnings: avgearn assumed to be 100
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: yx

	IF (yx<=0d0) THEN
		taxrate=0d0
		RETURN
	ENDIF

	taxrate = tau0 + tau1*LOG(yx/nmdlavginc)	!/avgearn)
	taxrate = MAX(0d0,MIN(1d0,taxrate))

RETURN
END FUNCTION
!-----------------------------------
REAL(dp) FUNCTION atinc(ex,sx,agex)
! relative wages
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: ex,sx
	INTEGER, INTENT(IN) :: agex
	REAL(dp) :: yx

	yx = ex + rr*sx
	atinc = (1d0-taxrate(yx))*yx - tauS*ex

	yx = gg*MERGE(aeq,1d0,agex==5 .or. agex==6 .or. agex==7)
	atinc = MAX(0d0,atinc) + yx
!	atinc = MAX(yx,atinc)

RETURN
END FUNCTION
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE locates(sx,swgts,si1)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: sx
	REAL(dp),INTENT(OUT):: swgts(2)
	INTEGER,INTENT(OUT) :: si1

	si1 = MAX0(1,MIN0(sn-1,INT(((sx-smin)/sstep)**(1d0/pows))+1))
	swgts(1) = ( sgrid(si1+1)-sx ) / (sgrid(si1+1)-sgrid(si1))
	swgts(2) = 1d0-swgts(1)

	swgts = MAX(0d0,MIN(swgts,1d0))
RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE locateh(hx,hwgts,hi1)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: hx
	REAL(dp),INTENT(OUT):: hwgts(2)
	INTEGER,INTENT(OUT) :: hi1

	hi1 = MAX0(1,MIN0(hn-1,INT(((hx-hmin)/hstep)**(1d0/powh))+1))
	hwgts(1) = ( hgrid(hi1+1)-hx ) / (hgrid(hi1+1)-hgrid(hi1))
	hwgts(2) = 1d0-hwgts(1)

	hwgts = MAX(0d0,MIN(hwgts,1d0))
RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE locatesk(skx,skwgts,ski1)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: skx
	REAL(dp),INTENT(OUT):: skwgts(2)
	INTEGER,INTENT(OUT) :: ski1

	ski1 = MAX0(1,MIN0(sn-1,INT(((skx-skmin)/skstep)**(1d0/pows))+1))
	skwgts(1) = ( skgrid(ski1+1)-skx ) / (skgrid(ski1+1)-skgrid(ski1))
	skwgts(2) = 1d0-skwgts(1)

	skwgts = MAX(0d0,MIN(skwgts,1d0))
RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE locatehk(hkx,hkwgts,hki1)
IMPLICIT NONE
	REAL(dp),INTENT(IN) :: hkx
	REAL(dp),INTENT(OUT):: hkwgts(2)
	INTEGER,INTENT(OUT) :: hki1

	hki1 = MAX0(1,MIN0(hn-1,INT(((hkx-hkmin)/hkstep)**(1d0/powh))+1))
	hkwgts(1) = ( hkgrid(hki1+1)-hkx ) / (hkgrid(hki1+1)-hkgrid(hki1))
	hkwgts(2) = 1d0-hkwgts(1)

	hkwgts = MAX(0d0,MIN(hkwgts,1d0))
RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
END MODULE



!-----------------------------------
!SUBROUTINE apf(capx,labx,tlabx,gdpx)
!IMPLICIT NONE
!	REAL(dp),INTENT(IN) :: capx,labx(2)
!	REAL(dp),INTENT(OUT) :: tlabx,gdpx
!	REAL(dp) :: col,hs
!
!	hs = ups*labx(1)**ces
!	col = (1d0-ups)*labx(2)**ces
!
!	tlabx = (col + hs)**(1d0/ces)
!
!	gdpx = capx**alph * tlabx**(1d0-alph)
!RETURN
!END SUBROUTINE
