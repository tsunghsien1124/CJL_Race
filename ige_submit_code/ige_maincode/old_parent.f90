MODULE old_parent
USE functions
IMPLICIT NONE

! -----------------------------------------------------------
! states dummies
	PRIVATE
	REAL(dp) :: ap
	REAL(dp) :: wp,gamp
	REAL(dp) :: ss,hh

	INTEGER :: ci,cii,ai,aii,ei,hi,hii,si
	INTEGER :: si1,hi1,ski1
	REAL(dp),DIMENSION(2) :: swgts,hwgts,skwgts
! -----------------------------------------------------------
	PUBLIC :: age10retirement,age9bequests

CONTAINS

! OLD PARENT'S PROBLEM
! old parent has nothing to backward induct

! -----------------------------------------------------------
! 10-12: retirement is deterministic
! -----------------------------------------------------------
SUBROUTINE age10retirement
IMPLICIT NONE

	REAL(dp) :: earn,cc

	DO ti = ts10,tf10,sn
		hi = MOD(ti/sn,hn);	hi = hi+1
		ci = ti/sn/hn;		ci = ci+1

! starts here-----------------------------------------------------------

	earn =  wskill(ci)*hgrid(hi)
	DO si = 1,sn

		cc = (2d0+(1d0-rtax)*rr)/(1d0+(1d0-rtax)*rr)**2 *(p0+p1*earn+gg)
		cc = atinc(earn,sgrid(si),10) +sgrid(si) +cc
		cc = cc / perm_cons_denom

		v10split(ti-ts10+si) = u(cc) * perm_cons_denom

	ENDDO;ENDDO

	CALL MPI_ALLGATHER(v10split,nwp10,MPI_REAL8,v10fun,nwp10,MPI_REAL8,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE

! -----------------------------------------------------------
! 9: bequests (inter-vivos)
! -----------------------------------------------------------
SUBROUTINE age9bequests
USE BFGS_PQN
IMPLICIT NONE
	INTEGER :: ii,jj,kk,bounds(3)
	REAL(dp),PARAMETER :: pwgts(3) = (/0.1d0,0.5d0,0.9d0/)
	REAL(dp),DIMENSION(3) :: pl0,pu0,gue0,pl,pu,gue
	REAL(dp) :: vt,vt0
	INTEGER :: iflg

	DO ti = ts89,tf89,sn
		hi = MOD(ti/sn,hn);				hi =hi +1
		hii= MOD(ti/sn/hn,hn);			hii=hii+1
		ai = MOD(ti/sn/hn/hn,an);		ai =ai +1
		aii= MOD(ti/sn/hn/hn/an,an);	aii=aii+1
		ci = MOD(ti/sn/hn/hn/an/an,cn);	ci =ci +1
		cii= ti/sn/hn/hn/an/an/cn;		cii=cii+1

		hh = hgrid(hi); ap = agrid(ai); wp = wskill(ci) *hh; gamp = gamm(ci);

! starts here-----------------------------------------------------------
! we need to choose n9,s10,s4

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

						vt = obj9(gue)
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
!			vt0 = obj9(gue0)
!		ENDIF

		CALL dfpmin_pqn(gue0,vt0,pl,pu,obj9,tol,maxfn,iflg)
!		pu0 = gue0
		pl0 = gue0

		w8v9split(ti-ts89+si) = -vt0
		n89split(ti-ts89+si) = gue0(1)
		s910split(ti-ts89+si) = gue0(2)
		n3s4split(ti-ts89+si) = gue0(3)

	ENDDO;ENDDO

	CALL MPI_ALLGATHER(w8v9split,nwp89,MPI_REAL8,v9fun ,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n89split ,nwp89,MPI_REAL8,n9fun ,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(s910split,nwp89,MPI_REAL8,s10fun,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n3s4split,nwp89,MPI_REAL8,s4fun ,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE
! -----------------------------------------------------------
! -----------------------------------------------------------
REAL(dp) FUNCTION obj9(choices)
IMPLICIT NONE
	REAL(dp), INTENT(IN) :: choices(:)
	REAL(dp) :: hpx,nx,spx,skx
	REAL(dp) :: cx,v4x,v10x,budget

	nx = choices(1)
	hpx = ap*(nx*hh)**gamp + hh
	budget = atinc(wp*(1d0-nx),ss,9) +ss -smin -skmin
	budget = MAX(0d0,budget)

	spx = choices(2)*budget
	skx = choices(3)*(budget-spx)
!	skx = MAX(0d0,skx)		! roundoff problem for bequest constraint

	cx = budget - spx - skx

	spx = spx +smin
	skx = skx +skmin

	! interpolate v4 and v10
	! ipolate needs locatex always: need xi1,xwgts
	CALL locatesk(skx,skwgts,ski1)
	v4x = DOT_PRODUCT(skwgts,v4fun(ski1:ski1+1,hii,aii,cii))

	CALL locates(spx,swgts,si1)
	v10x = ipolatev10(hpx,swgts,si1)

	obj9 = -u(cx) -thet*v4x -bta*v10x

RETURN
END FUNCTION
! -----------------------------------------------------------
REAL(dp) FUNCTION ipolatev10(hpx,swgts,si1)
IMPLICIT NONE
	REAL(dp),INTENT(IN):: hpx,swgts(2)
	INTEGER,INTENT(IN) :: si1
	REAL(dp) :: hpxx,temp(en,2)

	DO ei=1,en
		hpxx = egrid(ei,ci)*hpx
		CALL locateh(hpxx,hwgts,hi1)
		temp(ei,1) = DOT_PRODUCT(swgts,v10fun(si1:si1+1,hi1,ci))
		temp(ei,2) = DOT_PRODUCT(swgts,v10fun(si1:si1+1,hi1+1,ci))
		temp(ei,1) = DOT_PRODUCT(hwgts,temp(ei,:))
	ENDDO

	ipolatev10 = SUM(temp(:,1))/en

RETURN
END FUNCTION
! -----------------------------------------------------------
END MODULE

!		gue0 = 0.9d0*pl0 + 0.1d0*pu0
!		gue1 = 0.5d0*pl0 + 0.5d0*pu0
!		gue2 = 0.1d0*pl0 + 0.9d0*pu0
!
!		vt0 = obj9(gue0)
!		vt1 = obj9(gue1)
!		vt2 = obj9(gue2)
!
!		IF (vt1<vt0) THEN
!			vt0 = vt1
!			gue0= gue1
!		ENDIF
!		IF (vt2<vt0) THEN
!			vt0 = vt2
!			gue0= gue2
!		ENDIF
