MODULE new_parent
USE functions
IMPLICIT NONE

! -----------------------------------------------------------
! states dummies
	PRIVATE
	REAL(dp) :: kamk,ddk,lamk	!,phik
	REAL(dp) :: ap,ak
	REAL(dp) :: wp,gamp
	REAL(dp) :: ss,hh

	INTEGER :: ci,ai,aii,ei,hi,si
	INTEGER :: si1,hi1,hki1
	REAL(dp),DIMENSION(2) :: swgts,hwgts,hkwgts
! -----------------------------------------------------------
	PUBLIC :: age5infant

CONTAINS

! NEW PARENT'S PROBLEM
! -----------------------------------------------------------
! 5: new born
! -----------------------------------------------------------
SUBROUTINE age5infant
USE BFGS_PQN
IMPLICIT NONE
	INTEGER :: ii,jj,kk,bounds(3)
	REAL(dp),PARAMETER :: pwgts(3) = (/0.1d0,0.5d0,0.9d0/)
	REAL(dp),DIMENSION(3) :: pl0,pu0,gue0,pl,pu,gue
	REAL(dp) :: vt,vt0
	INTEGER :: iflg

	kamk = kamm(0)
!	phik = phi(0)
	ddk = dd(0)

	DO ti = ts5,tf5,sn
		hi = MOD(ti/sn,hn);			hi =hi +1
		ai = MOD(ti/sn/hn,an);		ai =ai +1
		aii= MOD(ti/sn/hn/an,an);	aii=aii+1
		ci = ti/sn/hn/an/an;		ci =ci +1

		hh = hgrid(hi);  ap = agrid(ai); wp = wskill(ci)*hh; gamp = gamm(ci);
		ak = agrid(aii)

		lamk = lamb(0,ci)

! starts here-----------------------------------------------------------
! we need to choose n5,x0,s6

	pl0 = (/0d0,0d0,0d0/)
	pu0 = (/1d0,1d0,1d0/)
	DO si = 1,sn!,1,-1 !

		ss = sgrid(si)

!		IF (si==sn) THEN
!		IF (si==1) THEN
			vt0 = infty
			DO ii = 1,3
				gue(1) = pl0(1)*pwgts(4-ii) + pu0(1)*pwgts(ii)
				DO jj = 1,3
					gue(2) = pl0(2)*pwgts(4-jj) + pu0(2)*pwgts(jj)
					DO kk = 1,3
						gue(3) = pl0(3)*pwgts(4-kk) + pu0(3)*pwgts(kk)

						vt = obj5(gue)
						IF (vt<vt0) THEN
							bounds = (/ii,jj,kk/)
							vt0 = vt
							gue0 = gue
						ENDIF
					ENDDO
				ENDDO
			ENDDO

			pl = pl0; pu = pu0
			DO ii=1,3
				IF (bounds(ii)==1) THEN
					pu(ii) = (pl0(ii)+pu0(ii))/2d0
					CYCLE
				ELSEIF (bounds(ii)==3) THEN
					pl(ii) = (pl0(ii)+pu0(ii))/2d0
					CYCLE
				ENDIF
				pu(ii) = pl0(ii)*pwgts(1) + pu0(ii)*pwgts(3)
				pl(ii) = pl0(ii)*pwgts(3) + pu0(ii)*pwgts(1)
			ENDDO

!		ELSE
!			gue0 = gue0*9d-1 + pu0*1d-1
!			vt0 = obj5(gue0)
!		ENDIF

		CALL dfpmin_pqn(gue0,vt0,pl,pu,obj5,tol,maxfn,iflg)
!		pu0 = gue0
		pl0 = gue0

		v5split(ti-ts5+si) = -vt0
		n5split(ti-ts5+si) = gue0(1)
		l0split(ti-ts5+si) = gue0(2)
		s6split(ti-ts5+si) = gue0(3)

	ENDDO;ENDDO

	CALL MPI_ALLGATHER(v5split,nwp5,MPI_REAL8,v5fun,nwp5,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n5split,nwp5,MPI_REAL8,n5fun,nwp5,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(l0split,nwp5,MPI_REAL8,l0fun,nwp5,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(s6split,nwp5,MPI_REAL8,s6fun,nwp5,MPI_REAL8,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE
! -----------------------------------------------------------
! -----------------------------------------------------------
REAL(dp) FUNCTION obj5(choices)
IMPLICIT NONE
	REAL(dp), INTENT(IN) :: choices(:)
	REAL(dp) :: hpx,nx,sx
	REAL(dp) :: cx,v6x,budget
	REAL(dp) :: lx,mx,xx,hkx

	nx = choices(1)
	lx = choices(2)
!	nx = nx*(1d0-lx)
	lx = lx*(1d0-nx)

	hpx = ap*(nx*hh)**gamp + hh

	mx = wp*lx *(1d0-kamk)/kamk
	xx = lamk*(wp*lx + mx + ddk)

	hkx = zeta* xx	!**phi0

	budget = atinc(wp*(1d0-nx-lx)-mx,ss,5) +ss -smin
	budget = MAX(0d0,budget)

	sx = choices(3)*budget
	cx = budget - sx

	sx = sx +smin

	! interpolate v6
	! ipolate needs locatex always: need xi1,xwgts
	CALL locates(sx,swgts,si1)
	CALL locatehk(hkx,hkwgts,hki1)
	v6x = ipolatev6(hpx,hkwgts,hki1,swgts,si1)

	obj5 = -qq*u(cx) -bta*v6x

RETURN
END FUNCTION
! -----------------------------------------------------------
REAL(dp) FUNCTION ipolatev6(hpx,hkwgts,hki1,swgts,si1)
IMPLICIT NONE
	REAL(dp),INTENT(IN):: hpx,hkwgts(2),swgts(2)
	INTEGER,INTENT(IN) :: hki1,si1
	REAL(dp) :: hpxx,temp(en,3),temk(2)

	DO ei=1,en
		hpxx = egrid(ei,ci)*hpx
		CALL locateh(hpxx,hwgts,hi1)

		temp(ei,1) = DOT_PRODUCT(swgts,v6fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(ei,2) = DOT_PRODUCT(swgts,v6fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		temp(ei,1) = DOT_PRODUCT(hwgts,temp(ei,1:2))

		temp(ei,2) = DOT_PRODUCT(swgts,v6fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(ei,3) = DOT_PRODUCT(swgts,v6fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		temp(ei,2) = DOT_PRODUCT(hwgts,temp(ei,2:3))
	ENDDO
	temk(1) = SUM(temp(:,1))/en
	temk(2) = SUM(temp(:,2))/en

	ipolatev6 = DOT_PRODUCT(hkwgts,temk)

RETURN
END FUNCTION
! -----------------------------------------------------------
END MODULE

!		gue0 = 0.9d0*pl0 + 0.1d0*pu0
!		gue1 = 0.5d0*pl0 + 0.5d0*pu0
!		gue2 = 0.1d0*pl0 + 0.9d0*pu0
!
!		vt0 = obj5(gue0)
!		vt1 = obj5(gue1)
!		vt2 = obj5(gue2)
!
!		IF (vt1<vt0) THEN
!			vt0 = vt1
!			gue0= gue1
!		ENDIF
!		IF (vt2<vt0) THEN
!			vt0 = vt2
!			gue0= gue2
!		ENDIF
