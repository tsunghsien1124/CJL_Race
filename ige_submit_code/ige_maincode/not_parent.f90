MODULE not_parent
USE functions
IMPLICIT NONE

! -----------------------------------------------------------
! states dummies
	PRIVATE
	REAL(dp) :: ap
	REAL(dp) :: wp,gamp
	REAL(dp) :: ss,hh

	INTEGER :: ci,ai,aii,ei,hi,si
	INTEGER :: si1,hi1
	REAL(dp),DIMENSION(2) :: swgts,hwgts
! -----------------------------------------------------------
	REAL(dp),DIMENSION(sn,hn,an,cn) :: v4funn
	PUBLIC :: v4funn,age4init,age4young


CONTAINS
! -----------------------------------------------------------
SUBROUTINE age4init
IMPLICIT NONE

	REAL(dp) :: earn,cc

	DO ti = ts4,tf4,sn
		hi = MOD(ti/sn,hn);		hi = hi+1
		ai = MOD(ti/sn/hn,an);	ai = ai+1
		ci = ti/sn/hn/an;		ci = ci+1

		earn = wskill(ci)*hgrid(hi)
		DO si = 1,sn
			cc = atinc(earn,0d0,4) +skgrid(si) -smin
			v4split(ti-ts4+si) = u(cc)
		ENDDO
	ENDDO
	CALL MPI_ALLGATHER(v4split,nwp4,MPI_REAL8,v4fun,nwp4,MPI_REAL8,MPI_COMM_WORLD,ierr)
!	v4fun = -zero
END SUBROUTINE
! -----------------------------------------------------------
SUBROUTINE age4young
USE BFGS_PQN
IMPLICIT NONE
	INTEGER :: ii,jj,bounds(2)
	REAL(dp),PARAMETER :: pwgts(3) = (/0.1d0,0.5d0,0.9d0/)
	REAL(dp),DIMENSION(2) :: pl0,pu0,gue0,pl,pu,gue
	REAL(dp) :: vt,vt0
	INTEGER :: iflg

	DO ti = ts4,tf4,sn
		hi = MOD(ti/sn,hn);		hi =hi +1
		ai = MOD(ti/sn/hn,an);	ai =ai +1
		ci = ti/sn/hn/an;		ci =ci +1

		hh = hgrid(hi);  ap = agrid(ai); wp = wskill(ci)*hh; gamp = gamm(ci);

! starts here-----------------------------------------------------------
! we need to choose n4,s5

	pl0 = (/0d0,0d0/)
	pu0 = (/1d0,1d0/)
	DO si = 1,sn!,1,-1 !

		ss = skgrid(si)

!		IF (si==sn) THEN
!		IF (si==1) THEN
			vt0 = infty
			DO ii = 1,3
				gue(1) = pl0(1)*pwgts(4-ii) + pu0(1)*pwgts(ii)
				DO jj = 1,3
					gue(2) = pl0(2)*pwgts(4-jj) + pu0(2)*pwgts(jj)

					vt = obj4(gue)
					IF (vt<vt0) THEN
						bounds = (/ii,jj/)
						vt0 = vt
						gue0 = gue
					ENDIF
				ENDDO
			ENDDO

			pl = pl0; pu = pu0
			DO ii=1,2
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
!			vt0 = obj4(gue0)
!		ENDIF

		CALL dfpmin_pqn(gue0,vt0,pl,pu,obj4,tol,maxfn,iflg)
!		pu0 = gue0
		pl0 = gue0

		v4split(ti-ts4+si) = -vt0
		n4split(ti-ts4+si) = gue0(1)
		s5split(ti-ts4+si) = gue0(2)

	ENDDO;ENDDO

	CALL MPI_ALLGATHER(v4split,nwp4,MPI_REAL8,v4funn,nwp4,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n4split,nwp4,MPI_REAL8,n4fun, nwp4,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(s5split,nwp4,MPI_REAL8,s5fun, nwp4,MPI_REAL8,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE
! -----------------------------------------------------------
! -----------------------------------------------------------
REAL(dp) FUNCTION obj4(choices)
IMPLICIT NONE
	REAL(dp), INTENT(IN) :: choices(:)
	REAL(dp) :: hpx,nx,sx
	REAL(dp) :: cx,v5x,budget

	nx = choices(1)
	hpx = ap*(nx*hh)**gamp + hh

	budget = atinc(wp*(1d0-nx),0d0,4) +ss -smin
	budget = MAX(0d0,budget)

	sx = choices(2)*budget
	cx = budget - sx

	sx = sx +smin

	! interpolate v5
	! ipolate needs locatex always: need xi1,xwgts
	CALL locates(sx,swgts,si1)
	v5x = ipolatev5(hpx,swgts,si1)

	obj4 = -u(cx) -bta*v5x

RETURN
END FUNCTION
! -----------------------------------------------------------
REAL(dp) FUNCTION ipolatev5(hpx,swgts,si1)
IMPLICIT NONE
	REAL(dp),INTENT(IN):: hpx,swgts(2)
	INTEGER,INTENT(IN) :: si1
	REAL(dp) :: hpxx,temp(en,2),temk(an)

	DO aii=1,an
		DO ei=1,en
			hpxx = egrid(ei,ci)*hpx
			CALL locateh(hpxx,hwgts,hi1)

			temp(ei,1) = DOT_PRODUCT(swgts,v5fun(si1:si1+1,hi1,ai,aii,ci))
			temp(ei,2) = DOT_PRODUCT(swgts,v5fun(si1:si1+1,hi1+1,ai,aii,ci))
			temp(ei,1) = DOT_PRODUCT(hwgts,temp(ei,1:2))
		ENDDO
		temk(aii) = SUM(temp(:,1))/en
	ENDDO

	ipolatev5 = DOT_PRODUCT(atrans(ai,:),temk)

RETURN
END FUNCTION
! -----------------------------------------------------------
END MODULE

!! V4 INIT WITHOUT MPI
!	DO ci=1,cn;
!		wp = wskill(ci)
!		DO hi=1,hn;
!			earn = wp*hgrid(hi)
!			DO si=1,sn
!
!				cc = atinc(earn,skgrid(si),0d0,4)
!				v4fun(si,hi,:,ci) = u(cc)
!
!			ENDDO
!		ENDDO
!	ENDDO

!! OLD INITIALIZATION
!		gue0 = 0.9d0*pl0 + 0.1d0*pu0
!		gue1 = 0.5d0*pl0 + 0.5d0*pu0
!		gue2 = 0.1d0*pl0 + 0.9d0*pu0
!
!		vt0 = obj4(gue0)
!		vt1 = obj4(gue1)
!		vt2 = obj4(gue2)
!
!		IF (vt1<vt0) THEN
!			vt0 = vt1
!			gue0= gue1
!		ENDIF
!		IF (vt2<vt0) THEN
!			vt0 = vt2
!			gue0= gue2
!		ENDIF
