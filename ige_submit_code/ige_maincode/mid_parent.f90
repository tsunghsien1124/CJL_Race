MODULE mid_parent
USE functions
IMPLICIT NONE

! -----------------------------------------------------------
! states dummies
	PRIVATE
	REAL(dp) :: ap,ak
	REAL(dp) :: wp,wk,gamp,gamk
	REAL(dp) :: ss,hh,hk

	INTEGER :: ci,cii,ai,aii,ei,hi,hii,si
	INTEGER :: si1,hi1,hki1
	REAL(dp),DIMENSION(2) :: swgts,hwgts,hkwgts
! -----------------------------------------------------------
	PUBLIC :: age8after,age8college

CONTAINS

! COLLEGE PARENT'S PROBLEM
! -----------------------------------------------------------
! 8.1: after college choice
! -----------------------------------------------------------
SUBROUTINE age8after
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
		hk = hgrid(hii);ak = agrid(aii);wk = wskill(cii)*hk; gamk = gamm(cii);

! starts here-----------------------------------------------------------
! we need to choose n8,n3,s9

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

						vt = obj8(gue)
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
!			vt0 = obj8(gue0)
!		ENDIF

		CALL dfpmin_pqn(gue0,vt0,pl,pu,obj8,tol,maxfn,iflg)
!		pu0 = gue0
		pl0 = gue0

		w8v9split(ti-ts89+si) = -vt0
		n89split(ti-ts89+si) = gue0(1)
		n3s4split(ti-ts89+si) = gue0(2)
		s910split(ti-ts89+si) = gue0(3)

	ENDDO;ENDDO

	CALL MPI_ALLGATHER(w8v9split,nwp89,MPI_REAL8,w8fun,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n89split ,nwp89,MPI_REAL8,n8fun,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(n3s4split,nwp89,MPI_REAL8,n3fun,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(s910split,nwp89,MPI_REAL8,s9fun,nwp89,MPI_REAL8,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE

! -----------------------------------------------------------
! 8: college choice
! -----------------------------------------------------------
SUBROUTINE age8college
IMPLICIT NONE
	REAL(DP) :: vht,vct,pref
	INTEGER :: ii

	DO ti = ts678,tf678,sn
		hi = MOD(ti/sn,hn);				hi =hi +1
		hii= MOD(ti/sn/hn,hn);			hii=hii+1
		ai = MOD(ti/sn/hn/hn,an);		ai =ai +1
		aii= MOD(ti/sn/hn/hn/an,an);	aii=aii+1
		ci = ti/sn/hn/hn/an/an;			ci =ci +1

		hh = hgrid(hi); ap = agrid(ai); wp = wskill(ci)*hh;
		hk = hgrid(hii);ak = agrid(aii)

! starts here-----------------------------------------------------------
! just choose the max

!	pref = zcolm + zcols*(LOG(ak)-mloga)
	pref = zcolm + MERGE(0d0,zcols,ci==1)

!	ii = sn
	ii = 1
	DO si = 1,sn!,1,-1

		vht = w8fun(si,hi,hii,ai,aii,ci,1)
		vct = w8fun(si,hi,hii,ai,aii,ci,2) + pref
		IF (vht < vct) EXIT

		c3split(ti-ts678+si) = 0
		v678split(ti-ts678+si) = vht

		ii = si+1

	ENDDO

	IF (ii<sn+1) THEN
		c3split(ti-ts678+si:ti-ts678+sn) = 1
		v678split(ti-ts678+si:ti-ts678+sn) = w8fun(si:sn,hi,hii,ai,aii,ci,2) + pref
	ENDIF

	ENDDO

	CALL MPI_ALLGATHER(v678split,nwp678,MPI_REAL8,  v8fun,nwp678,MPI_REAL8,  MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(c3split,  nwp678,MPI_INTEGER,c3fun,nwp678,MPI_INTEGER,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE
! -----------------------------------------------------------
! -----------------------------------------------------------
REAL(dp) FUNCTION obj8(choices)
IMPLICIT NONE
	REAL(dp), INTENT(IN) :: choices(:)
	REAL(dp) :: hpx,npx,sx
	REAL(dp) :: hkx,nkx
	REAL(dp) :: cx,v9x,budget

	npx = choices(1)
	nkx = ctime(1,cii) + ctime(2,cii)*choices(2)

	hpx = ap*(npx*hh)**gamp + hh
	hkx = ak*(nkx*hk)**gamk + hk

	budget = atinc(wp*(1d0-npx),ss, 8) +ss -smin
	budget = atinc(wk*(1d0-nkx),0d0,3) +budget
	budget = MAX(0d0,budget - MERGE(0d0,kappa,cii==1))

	sx = choices(3)*budget
	cx = budget - sx

	sx = sx +smin

	! interpolate v9
	! ipolate needs locatex always: need xi1,xwgts
	CALL locates(sx,swgts,si1)
	v9x = ipolatev9(hpx,hkx,swgts,si1)

	obj8 = -qq*u(cx) -bta*v9x

RETURN
END FUNCTION
! -----------------------------------------------------------
REAL(dp) FUNCTION ipolatev9(hpx,hkx,swgts,si1)
IMPLICIT NONE
	REAL(dp),INTENT(IN):: hpx,hkx,swgts(2)
	INTEGER,INTENT(IN) :: si1
	REAL(dp) :: hpxx,hkxx,temp(en,3),temk(en,2)
	INTEGER :: eii

	DO eii=1,en

		hkxx = egrid(eii,cii)*hkx
		CALL locateh(hkxx,hkwgts,hki1)

		DO ei=1,en
			hpxx = egrid(ei,ci)*hpx
			CALL locateh(hpxx,hwgts,hi1)

			temp(ei,1) = DOT_PRODUCT(swgts,v9fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
			temp(ei,2) = DOT_PRODUCT(swgts,v9fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
			temp(ei,1) = DOT_PRODUCT(hwgts,temp(ei,1:2))

			temp(ei,2) = DOT_PRODUCT(swgts,v9fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
			temp(ei,3) = DOT_PRODUCT(swgts,v9fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
			temp(ei,2) = DOT_PRODUCT(hwgts,temp(ei,2:3))
		ENDDO
		temk(eii,1) = SUM(temp(:,1))/en
		temk(eii,2) = SUM(temp(:,2))/en

		temk(eii,1) = DOT_PRODUCT(hkwgts,temk(eii,1:2))
	ENDDO
!	ipolatev9 = DOT_PRODUCT(etrans(aii,:),temk(:,1))
	ipolatev9 = SUM(temk(:,1))/en

RETURN
END FUNCTION
! -----------------------------------------------------------
END MODULE

!		gue0 = 0.9d0*pl0 + 0.1d0*pu0
!		gue1 = 0.5d0*pl0 + 0.5d0*pu0
!		gue2 = 0.1d0*pl0 + 0.9d0*pu0
!
!		vt0 = obj8(gue0)
!		vt1 = obj8(gue1)
!		vt2 = obj8(gue2)
!
!		IF (vt1<vt0) THEN
!			vt0 = vt1
!			gue0= gue1
!		ENDIF
!		IF (vt2<vt0) THEN
!			vt0 = vt2
!			gue0= gue2
!		ENDIF
