MODULE simulation
USE functions
IMPLICIT NONE

	REAL(dp) :: ap,ak,ep
	INTEGER :: ai,aii,ei,agei,ind

	REAL(dp) :: l0,l1,l2,n3,n4,n5,n6,n7,n8,n9
	REAL(dp) :: h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
	REAL(dp) :: s4,s5,s6,s7,s8,s9,s10,c3
	REAL(dp) :: e3,e4,e5,e6,e7,e8,e9,e10
	REAL(dp) :: i4,i5,i6,i7,i8,i9,i10,i11,i12

	REAL(dp) :: wp,gamp,budget,temp(2),pmet(2)
	REAL(dp) :: wk,gamk,mt,xt

	INTEGER :: si,hi,ci,cii
	INTEGER :: si1,hi1,hki1
	REAL(dp),DIMENSION(2) :: swgts,hwgts,hkwgts

CONTAINS
!---------------------------------------------------------------------
SUBROUTINE fix_proc
USE random
IMPLICIT NONE

	INTEGER :: sim

	CALL Set_Seed(123456789+myid,362436069+myid,521288629+myid,916191069+myid)

	DO ind=1,NIproc
!---------------------------------------------------------------------
! draw first ability
		ap = sample_uniform(0d0,1d0)
		DO ai=1,an
			IF ( ap<=SUM(astat(:ai)) ) THEN
				abi(ind,1) = ai
				EXIT
			ENDIF
		ENDDO

! draw luck shocks, then son's ability
		DO sim=1,TT+2
!			ep = sample_uniform(0d0,1d0)
!			DO ei=1,en
!				IF ( ep<=SUM(etrans(abi(ind,sim),:ei)) ) THEN
!					epi(ind,sim,4)=ei
!					EXIT
!				ENDIF
!			ENDDO
			DO agei = 4,10
				ep = sample_uniform(0d0,1d0)
				DO ei=1,en
					IF ( ep <= DBLE(ei)/DBLE(en) ) THEN
						epi(ind,sim,agei)=ei
						EXIT
					ENDIF
				ENDDO
			ENDDO
		!-------------------------------
			IF (sim==TT+2) EXIT
		!-------------------------------
			ak = sample_uniform(0d0,1d0)
			DO aii=1,an
				IF ( ak<=SUM(atrans(abi(ind,sim),:aii)) ) THEN
					abi(ind,sim+1)=aii
					EXIT
				ENDIF
			ENDDO
!---------------------------------------------------------------------
		ENDDO
		hhidproc(ind)%ability_p = abi(ind,TT-1)
		hhidproc(ind)%ability_k = abi(ind,TT)
		hhidproc(ind)%ability_g = abi(ind,TT+1)

	ENDDO
RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE mcmc
IMPLICIT NONE

	INTEGER :: sim

	DO ind = 1,NIproc

! pick different starting points for each YOUNG individual
! age 4
		si = MOD(ind,sn)+1;
		hi = MOD(ind,hn)+1;
		ai = abi(ind,1)
		ci = MOD(ind,cn)+1;

		h4 = hgrid(hi); ap = agrid(ai); wp = wskill(ci)*h4; gamp = gamm(ci)
		s4 = skgrid(si)

		! earnings choice
		n4 = n4fun(si,hi,ai,ci)
		n4 = MAX(0d0,MIN(n4,1d0))

		e4 = MAX(0d0,wp*(1d0-n4))

		h5 = ap*(n4*h4)**gamp + h4
		h5 = egrid(epi(ind,1,5),ci) *h5
		h5 = MAX(hmin,MIN(h5,hmax))

		! savings choice
		i4 = atinc(e4,0d0,4)
		budget = i4 +s4 -smin
		budget = MAX(0d0,budget)
		s5 = s5fun(si,hi,ai,ci) *budget
		s5 = MAX(0d0,MIN(s5,smax)) +smin

! --------------
		DO sim = 1,TT-2
			CALL timeiter(sim)
		ENDDO
		colfpproc(ind) = c3
		netincpproc(4,ind) = i4

		hhidproc(ind)%col_p = ci
		hhidproc(ind)%nn_p(3:4) = (/n3,n4/)
		hhidproc(ind)%ee_p(3:4) = (/e3,e4/)
		hhidproc(ind)%hh_p(1:5) = (/h1,h2,h3,h4,h5/)
		hhidproc(ind)%ss_p(4:5) = (/s4,s5/)

		CALL simul_end

	ENDDO

!------------------------
	CALL MPI_ALLGATHER(hhidproc,NIproc,MPI_AGENT,hhid,NIproc,MPI_AGENT,MPI_COMM_WORLD,ierr)
!------------------------
	CALL MPI_ALLGATHER(colfpproc,NIproc,MPI_REAL8,colfp,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(colfkproc,NIproc,MPI_REAL8,colfk,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(colfgproc,NIproc,MPI_REAL8,colfg,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!------------------------
	CALL MPI_ALLGATHER(netincpproc(:,4:12),NIproc*9,MPI_REAL8,netincp,NIproc*9,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(netinckproc(:,4:12),NIproc*9,MPI_REAL8,netinck,NIproc*9,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(netincgproc(:,4:12),NIproc*9,MPI_REAL8,netincg,NIproc*9,MPI_REAL8,MPI_COMM_WORLD,ierr)
!------------------------

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE timeiter(sim)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: sim
	INTEGER :: iii
	REAL(dp) :: vht,vct

! age 5: kid is born
		wp = wskill(ci)*h5
		aii = abi(ind,sim+1)
		ak = agrid(aii)

		CALL locates(s5,swgts,si1)
		CALL locateh(h5,hwgts,hi1)

		temp(1) = DOT_PRODUCT(swgts,n5fun(si1:si1+1,hi1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,n5fun(si1:si1+1,hi1+1,ai,aii,ci))

		n5 = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,l0fun(si1:si1+1,hi1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,l0fun(si1:si1+1,hi1+1,ai,aii,ci))

		l0 = DOT_PRODUCT(hwgts,temp)

!		l0 = MAX(0d0,MIN(l0,1d0))
!		n5 = n5*(1d0-l0)
!		n5 = MAX(0d0,MIN(n5,1d0-l0))

		n5 = MAX(0d0,MIN(n5,1d0))
		l0 = l0*(1d0-n5)
		l0 = MAX(0d0,MIN(l0,1d0-n5))

		! earnings choice
		e5 = MAX(0d0,wp*(1d0-n5-l0))

		h6 = ap*(n5*h5)**gamp + h5
		h6 = egrid(epi(ind,sim,6),ci) *h6
		h6 = MAX(hmin,MIN(h6,hmax))

		! educ choice
		mt = wp*l0 *(1d0-kamm(0))/kamm(0)
		xt = lamb(0,ci)*(wp*l0 + mt + dd(0))
		h1 = zeta* xt		!**phi0
		h1 = MAX(hkmin,MIN(h1,hkmax))

		! savings choice
		i5 = atinc(e5-mt,s5,5)
		budget = i5 +s5 -smin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s6fun(si1:si1+1,hi1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,s6fun(si1:si1+1,hi1+1,ai,aii,ci))
		s6 = DOT_PRODUCT(hwgts,temp) *budget
		s6 = MAX(0d0,MIN(s6,smax)) +smin

! age 6: primary school
		wp = wskill(ci)*h6

		CALL locates(s6,swgts,si1)
		CALL locateh(h6,hwgts,hi1)
		CALL locatehk(h1,hkwgts,hki1)

		temp(1) = DOT_PRODUCT(swgts,n6fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,n6fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,n6fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,n6fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		n6 = DOT_PRODUCT(hkwgts,pmet)

		temp(1) = DOT_PRODUCT(swgts,l1fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,l1fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,l1fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,l1fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		l1 = DOT_PRODUCT(hkwgts,pmet)

!		l1 = MAX(0d0,MIN(l1,1d0))
!		n6 = n6*(1d0-l1)
!		n6 = MAX(0d0,MIN(n6,1d0-l1))

		n6 = MAX(0d0,MIN(n6,1d0))
		l1 = l1*(1d0-n6)
		l1 = MAX(0d0,MIN(l1,1d0-n6))

		! earnings choice
		e6 = MAX(0d0,wp*(1d0-n6-l1))

		h7 = ap*(n6*h6)**gamp + h6
		h7 = egrid(epi(ind,sim,7),ci) *h7
		h7 = MAX(hmin,MIN(h7,hmax))

		! educ choice
		mt = wp*l1 *(1d0-kamm(1))/kamm(1)
		xt = lamb(1,ci)*(wp*l1 + mt + dd(1))
		xt = zeta* xt

!		h2 = omega(1)* xt**phi(1) + (1d0-omega(1))*h1**phi(1)
!		h2 = h2**(1d0/phi(1))

		h2 = xt**omega(1) * h1**(1d0-omega(1))
		h2 = MAX(hkmin,MIN(h2,hkmax))

		! savings choice
		i6 = atinc(e6-mt,s6,6)
		budget = i6 +s6 -smin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s7fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,s7fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,s7fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,s7fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		s7 = DOT_PRODUCT(hkwgts,pmet) *budget
		s7 = MAX(0d0,MIN(s7,smax)) +smin

! age 7: secondary school
		wp = wskill(ci)*h7

		CALL locates(s7,swgts,si1)
		CALL locateh(h7,hwgts,hi1)
		CALL locatehk(h2,hkwgts,hki1)

		temp(1) = DOT_PRODUCT(swgts,n7fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,n7fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,n7fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,n7fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		n7 = DOT_PRODUCT(hkwgts,pmet)

		temp(1) = DOT_PRODUCT(swgts,l2fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,l2fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,l2fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,l2fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		l2 = DOT_PRODUCT(hkwgts,pmet)

!		l2 = MAX(0d0,MIN(l2,1d0))
!		n7 = n7*(1d0-l2)
!		n7 = MAX(0d0,MIN(n7,1d0-l2))

		n7 = MAX(0d0,MIN(n7,1d0))
		l2 = l2*(1d0-n7)
		l2 = MAX(0d0,MIN(l2,1d0-n7))

		! earnings choice
		e7 = MAX(0d0,wp*(1d0-n7-l2))

		h8 = ap*(n7*h7)**gamp + h7
		h8 = egrid(epi(ind,sim,8),ci) *h8
		h8 = MAX(hmin,MIN(h8,hmax))

		! educ choice
		mt = wp*l2 *(1d0-kamm(2))/kamm(2)
		xt = lamb(2,ci)*(wp*l2 + mt + dd(2))
		xt = zeta* xt

!		h3 = omega(2)* xt**phi(2) + (1d0-omega(2))*h2**phi(2)
!		h3 = h3**(1d0/phi(2))
!		h3 = h3**(nu/phi(2))

		h3 = xt**omega(2) * h2**(1d0-omega(2))
!		h3 = h3**nu
		h3 = MAX(hmin,MIN(h3,hmax))

		! savings choice
		i7 = atinc(e7-mt,s7,7)
		budget = i7 +s7 -smin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s8fun(si1:si1+1,hi1,hki1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,s8fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
		pmet(1) = DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,s8fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
		temp(2) = DOT_PRODUCT(swgts,s8fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
		pmet(2) = DOT_PRODUCT(hwgts,temp)

		s8 = DOT_PRODUCT(hkwgts,pmet) *budget
		s8 = MAX(0d0,MIN(s8,smax)) +smin

! age 8: kid's college choice
		wp = wskill(ci)*h8

		CALL locates(s8,swgts,si1)
		CALL locateh(h8,hwgts,hi1)
		CALL locateh(h3,hkwgts,hki1)

		!! college choice
		temp(1) = DOT_PRODUCT(swgts,DBLE(c3fun(si1:si1+1,hi1,hki1,ai,aii,ci)))
		temp(2) = DOT_PRODUCT(swgts,DBLE(c3fun(si1:si1+1,hi1+1,hki1,ai,aii,ci)))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,DBLE(c3fun(si1:si1+1,hi1,hki1+1,ai,aii,ci)))
		temp(2) = DOT_PRODUCT(swgts,DBLE(c3fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci)))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		c3 = DOT_PRODUCT(hkwgts,pmet)
		c3 = MAX(0d0,MIN(c3,1d0))

		! high school
		temp(1) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1,hki1,ai,aii,ci,1))
		temp(2) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,1))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,1))
		temp(2) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,1))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		vht = DOT_PRODUCT(hkwgts,pmet)

		! college
		temp(1) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1,hki1,ai,aii,ci,2))
		temp(2) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,2))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,2))
		temp(2) = DOT_PRODUCT(swgts,w8fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,2))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

!		vct = DOT_PRODUCT(hkwgts,pmet) + zcolm + zcols*(LOG(ak)-mloga)
		vct = DOT_PRODUCT(hkwgts,pmet) + zcolm + MERGE(0d0,zcols,ci==1)

		IF (vht<vct .OR. c3>0.5d0) THEN
			cii = 2
		ELSE
			cii = 1
		ENDIF

!		DO iii=1,2
			wk = wskill(cii)*h3
			gamk = gamm(cii)

	! age 8.1: kid's work choice

			! earnings choice
			temp(1) = DOT_PRODUCT(swgts,n8fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
			temp(2) = DOT_PRODUCT(swgts,n8fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
			pmet(1)= DOT_PRODUCT(hwgts,temp)

			temp(1) = DOT_PRODUCT(swgts,n8fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
			temp(2) = DOT_PRODUCT(swgts,n8fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
			pmet(2)= DOT_PRODUCT(hwgts,temp)

			n8 = DOT_PRODUCT(hkwgts,pmet)
			n8 = MAX(0d0,MIN(n8,1d0))

			e8 = MAX(0d0,wp*(1d0-n8))

			h9 = ap*(n8*h8)**gamp + h8
			h9 = egrid(epi(ind,sim,9),ci) *h9
			h9 = MAX(hmin,MIN(h9,hmax))

			! educ choice
			temp(1) = DOT_PRODUCT(swgts,n3fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
			temp(2) = DOT_PRODUCT(swgts,n3fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
			pmet(1)= DOT_PRODUCT(hwgts,temp)

			temp(1) = DOT_PRODUCT(swgts,n3fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
			temp(2) = DOT_PRODUCT(swgts,n3fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
			pmet(2)= DOT_PRODUCT(hwgts,temp)

			n3 = ctime(1,cii) + ctime(2,cii)*DOT_PRODUCT(hkwgts,pmet)
			n3 = MAX(ctime(1,cii),MIN(n3,1d0))

			e3 = MAX(0d0,wk*(1d0-n3))

			h4 = ak*(n3*h3)**gamk + h3
			h4 = egrid(epi(ind,sim+1,4),cii) *h4
			h4 = MAX(hmin,MIN(h4,hmax))

			! savings choice
			i8 = atinc(e8,s8,8) + atinc(e3,0d0,3)
			budget = i8 + s8 -smin
			budget = MAX(0d0,budget - MERGE(0d0,kappa,cii==1))

!			IF (cii==2 .and. budget-kappa<=gg*2d0) THEN
!				cii=1
!				CONTINUE
!			ELSE
!				budget = budget - MERGE(gg*2d0-smin,kappa,cii==1)
!				EXIT
!			ENDIF
!
!		ENDDO

		temp(1) = DOT_PRODUCT(swgts,s9fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s9fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,s9fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s9fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		s9 = DOT_PRODUCT(hkwgts,pmet) *budget
		s9 = MAX(0d0,MIN(s9,smax)) +smin

! age 9: bequests
		wp = wskill(ci)*h9

		CALL locates(s9,swgts,si1)
		CALL locateh(h9,hwgts,hi1)
		CALL locateh(h4,hkwgts,hki1)

		! earnings choice
		temp(1) = DOT_PRODUCT(swgts,n9fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,n9fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,n9fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,n9fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		n9 = DOT_PRODUCT(hkwgts,pmet)
		n9 = MAX(0d0,MIN(n9,1d0))

		e9 = MAX(0d0,wp*(1d0-n9))

		h10 = ap*(n9*h9)**gamp + h9
		h10 = egrid(epi(ind,sim,10),ci) *h10
		h10 = MAX(hmin,MIN(h10,hmax))

		! savings choice
		i9 = atinc(e9,s9,9)
		budget = i9 +s9 -smin -skmin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s10fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s10fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,s10fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s10fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		s10 = DOT_PRODUCT(hkwgts,pmet) *budget
		s10 = MAX(0d0,MIN(s10,smax))

		! bequest choice
		temp(1) = DOT_PRODUCT(swgts,s4fun(si1:si1+1,hi1,hki1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s4fun(si1:si1+1,hi1+1,hki1,ai,aii,ci,cii))
		pmet(1)= DOT_PRODUCT(hwgts,temp)

		temp(1) = DOT_PRODUCT(swgts,s4fun(si1:si1+1,hi1,hki1+1,ai,aii,ci,cii))
		temp(2) = DOT_PRODUCT(swgts,s4fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci,cii))
		pmet(2)= DOT_PRODUCT(hwgts,temp)

		s4 = DOT_PRODUCT(hkwgts,pmet) *(budget-s10)
		s4 = MAX(0d0,MIN(s4,skmax)) + skmin

		s10 = s10+smin

!--------------------------
! NEXT GENERATION age 4: bequests
		ai = aii; ap = ak; ci = cii;
		wp = wskill(ci)*h4; gamp = gamk

		CALL locatesk(s4,swgts,si1)

		! earnings choice
		temp(1) = DOT_PRODUCT(swgts,n4fun(si1:si1+1,hki1,ai,ci))
		temp(2) = DOT_PRODUCT(swgts,n4fun(si1:si1+1,hki1+1,ai,ci))
		n4 = DOT_PRODUCT(hwgts,temp)
		n4 = MAX(0d0,MIN(n4,1d0))

		e4 = MAX(0d0,wp*(1d0-n4))

		h5 = ap*(n4*h4)**gamp + h4
		h5 = egrid(epi(ind,sim+1,5),ci) *h5
		h5 = MAX(hmin,MIN(h5,hmax))

		! savings choice
		i4 = atinc(e4,0d0,4)
		budget = i4 +s4 -smin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s5fun(si1:si1+1,hki1,ai,ci))
		temp(2) = DOT_PRODUCT(swgts,s5fun(si1:si1+1,hki1+1,ai,ci))
		s5 = DOT_PRODUCT(hwgts,temp) *budget
		s5 = MAX(0d0,MIN(s5,smax)) +smin

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE retirement(cix,h10,s10,e10,i11,i12,s11,s12)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: cix
	REAL(dp),INTENT(IN) :: h10,s10
	REAL(dp),INTENT(OUT) :: e10,i11,i12,s11,s12
	REAL(dp) :: b10,c10,b11,c11

		e10 = wskill(cix)*h10

		i10 = atinc(e10,s10,10)
		b10 = i10 +s10
		c10 = b10 + (2d0+(1d0-rtax)*rr)/(1d0+(1d0-rtax)*rr)**2 *(p0+p1*e10+gg)
		c10 = c10 / perm_cons_denom

		s11 = b10 - c10

		i11 = (1d0-rtax)*rr*s11 + p0+p1*e10 + gg
		b11 = i11 + s11
		c11 = c10*(bta*(1d0+(1d0-rtax)*rr))**(1d0/chi)

		s12 = b11 - c11

		i12 = (1d0-rtax)*rr*s12 + p0+p1*e10 + gg

END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE simul_end
IMPLICIT NONE
	REAL(dp) :: s11,s12

		CALL timeiter(TT-1)
! COMPLETE OLD AGE SAVINGS: ASSUME NOT-TAXED
		netincpproc(5:10,ind) = (/i5,i6,i7,i8,i9,i10/)

		hhidproc(ind)%nn_p(0:2) = (/l0,l1,l2/)
		hhidproc(ind)%nn_p(5:9) = (/n5,n6,n7,n8,n9/)
		hhidproc(ind)%ee_p(5:9) = (/e5,e6,e7,e8,e9/)
		hhidproc(ind)%hh_p(6:10) = (/h6,h7,h8,h9,h10/)
		hhidproc(ind)%ss_p(6:10) = (/s6,s7,s8,s9,s10/)

		CALL retirement(hhidproc(ind)%col_p,h10,s10,e10,i11,i12,s11,s12)
		netincpproc(11:12,ind) = (/i11,i12/)

		hhidproc(ind)%ee_p(10) = 	e10
		hhidproc(ind)%ss_p(11:12) = (/s11,s12/)
! --------------
! --------------
		colfkproc(ind) = c3
		netinckproc(4,ind) = i4

		hhidproc(ind)%col_k = ci
		hhidproc(ind)%nn_k(3:4) = (/n3,n4/)
		hhidproc(ind)%ee_k(3:4) = (/e3,e4/)
		hhidproc(ind)%ss_k(4:5) = (/s4,s5/)
		hhidproc(ind)%hh_k(1:5) = (/h1,h2,h3,h4,h5/)

		CALL timeiter(TT)
! COMPLETE OLD AGE SAVINGS: ASSUME NOT-TAXED
		netinckproc(5:10,ind) = (/i5,i6,i7,i8,i9,i10/)

		hhidproc(ind)%nn_k(0:2) = (/l0,l1,l2/)
		hhidproc(ind)%nn_k(5:9) = (/n5,n6,n7,n8,n9/)
		hhidproc(ind)%ee_k(5:9) = (/e5,e6,e7,e8,e9/)
		hhidproc(ind)%ss_k(6:10) = (/s6,s7,s8,s9,s10/)
		hhidproc(ind)%hh_k(6:10) = (/h6,h7,h8,h9,h10/)

		CALL retirement(hhidproc(ind)%col_k,h10,s10,e10,i11,i12,s11,s12)
		netinckproc(11:12,ind) = (/i11,i12/)

		hhidproc(ind)%ee_k(10) = 	e10
		hhidproc(ind)%ss_k(11:12) = (/s11,s12/)
! --------------
! --------------
		colfgproc(ind) = c3
		netincgproc(4,ind) = i4

		hhidproc(ind)%col_g = ci
		hhidproc(ind)%nn_g(3:4) = (/n3,n4/)
		hhidproc(ind)%ee_g(3:4) = (/e3,e4/)
		hhidproc(ind)%ss_g(4:5) = (/s4,s5/)
		hhidproc(ind)%hh_g(1:5) = (/h1,h2,h3,h4,h5/)

		CALL timeiter(TT+1)
! COMPLETE OLD AGE SAVINGS: ASSUME NOT-TAXED
		netincgproc(5:10,ind) = (/i5,i6,i7,i8,i9,i10/)

		hhidproc(ind)%nn_g(0:2) = (/l0,l1,l2/)
		hhidproc(ind)%nn_g(5:9) = (/n5,n6,n7,n8,n9/)
		hhidproc(ind)%ee_g(5:9) = (/e5,e6,e7,e8,e9/)
		hhidproc(ind)%ss_g(6:10) = (/s6,s7,s8,s9,s10/)
		hhidproc(ind)%hh_g(6:10) = (/h6,h7,h8,h9,h10/)

		CALL retirement(hhidproc(ind)%col_g,h10,s10,e10,i11,i12,s11,s12)
		netincgproc(11:12,ind) = (/i11,i12/)

		hhidproc(ind)%ee_g(10) = 	e10
		hhidproc(ind)%ss_g(11:12) = (/s11,s12/)

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
END MODULE


!		! college choice
!		temp(1) = DOT_PRODUCT(swgts,c3fun(si1:si1+1,hi1,hki1,ai,aii,ci))
!		temp(2) = DOT_PRODUCT(swgts,c3fun(si1:si1+1,hi1+1,hki1,ai,aii,ci))
!		pmet(1) = DOT_PRODUCT(hwgts,temp)
!
!		temp(1) = DOT_PRODUCT(swgts,c3fun(si1:si1+1,hi1,hki1+1,ai,aii,ci))
!		temp(2) = DOT_PRODUCT(swgts,c3fun(si1:si1+1,hi1+1,hki1+1,ai,aii,ci))
!		pmet(2) = DOT_PRODUCT(hwgts,temp)
!
!		xt = DOT_PRODUCT(hkwgts,pmet)
!		IF (xt<=0.5d0) THEN
!			cii = 1
!			sc = s8
!		ELSE
!			cii = 2
!			sc = MAX(s8-kappa,smin)
!			CALL locates(sc,swgts,si1)
!		ENDIF
