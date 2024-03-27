MODULE lifetime_stats
USE global
IMPLICIT NONE

	PRIVATE
	REAL(dp),ALLOCATABLE :: consproc(:,:)
	PUBLIC :: gen_lifetime_stats

CONTAINS
!---------------------------------------------------------------------
SUBROUTINE gen_lifetime_stats
USE functions, ONLY: u,locateh,locatesk
IMPLICIT NONE

	REAL(dp),ALLOCATABLE :: lfavgeproc(:),lfearnproc(:),lfwlthproc(:),lfutilproc(:)

	INTEGER :: agei,ind
	INTEGER :: ci,ai,hi,si
	REAL(dp) :: h4,s4
	REAL(dp),DIMENSION(2) :: swgts,hwgts,temp

	REAL(dp) :: cmean(4:10), cmeanproc(4:10)
	REAL(dp) :: cstd(4:10),  cstdproc(4:10)
	INTEGER ::  ccount(4:10),ccountproc(4:10)

	ALLOCATE( consproc(NIproc,4:12) )
	ALLOCATE( lfavgeproc(NIproc),lfearnproc(NIproc),lfwlthproc(NIproc),lfutilproc(NIproc) )

!----------------------------
! GET SON'S LIFETIME OUTCOMES
!----------------------------
! LIFETIME EARNINGS AND WEALTH

	lfavgeproc = hhidproc%ee_k(4)
	lfearnproc = hhidproc%ee_k(4)
	DO agei = 5,10
		lfavgeproc = lfavgeproc + hhidproc%ee_k(agei)
		lfearnproc = lfearnproc + hhidproc%ee_k(agei)/(1d0+rr)**(agei-4)
	ENDDO
	lfavgeproc = lfavgeproc/7d0
	lfwlthproc = lfearnproc + hhidproc%ss_k(4)
!----------------------------
	CALL MPI_ALLGATHER(lfavgeproc,NIproc,MPI_REAL8,lfkavge,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfearnproc,NIproc,MPI_REAL8,lfkearn,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfwlthproc,NIproc,MPI_REAL8,lfkwlth,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!----------------------------

! LIFETIME UTILITY: SIMULATE SON'S OUTCOMES STARTING FROM AGE 4
! GET CONSUMPTION PROC
	DO ind = 1,NIproc
		CALL onegen(ind)
	ENDDO

! LF UTIL
	DO ind = 1,NIproc

		lfutilproc(ind) = u(consproc(ind,4))
		DO agei = 5,8
			lfutilproc(ind) = lfutilproc(ind) + bta**(agei-4)*qq*u(consproc(ind,agei))
		ENDDO
		DO agei = 9,12
			lfutilproc(ind) = lfutilproc(ind) + bta**(agei-4)*u(consproc(ind,agei))
		ENDDO

! INTERPOLATE FOR GRANDKID'S VALUE

		ci = hhidproc(ind)%col_g
		ai = hhidproc(ind)%ability_g

		h4 = hhidproc(ind)%hh_g(4)
		s4 = hhidproc(ind)%ss_g(4)
		CALL locateh(h4,hwgts,hi)
		CALL locatesk(s4,swgts,si)

		temp(1) = DOT_PRODUCT(swgts,v4fun(si:si+1,hi,ai,ci))
		temp(2) = DOT_PRODUCT(swgts,v4fun(si:si+1,hi+1,ai,ci))
		lfutilproc(ind) = lfutilproc(ind) + bta**5*thet*DOT_PRODUCT(hwgts,temp)

	ENDDO
!----------------------------
	CALL MPI_ALLGATHER(lfutilproc,NIproc,MPI_REAL8,lfkutil,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!----------------------------
	DEALLOCATE( lfutilproc )
!----------------------------



!----------------------------
! GET PARENT'S LIFETIME OUTCOMES
!----------------------------
	lfavgeproc = hhidproc%ee_p(4)
	lfearnproc = hhidproc%ee_p(4)
	DO agei = 5,10
		lfavgeproc = lfavgeproc + hhidproc%ee_p(agei)
		lfearnproc = lfearnproc + hhidproc%ee_p(agei)/(1d0+rr)**(agei-4)
	ENDDO
	lfavgeproc = lfavgeproc/7d0
	lfwlthproc = lfearnproc + hhidproc%ss_p(4)
!----------------------------
	CALL MPI_ALLGATHER(lfavgeproc,NIproc,MPI_REAL8,lfpavge,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfearnproc,NIproc,MPI_REAL8,lfpearn,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfwlthproc,NIproc,MPI_REAL8,lfpwlth,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!----------------------------

!----------------------------
! GET GRANDKID'S LIFETIME OUTCOMES
!----------------------------
	lfavgeproc = hhidproc%ee_g(4)
	lfearnproc = hhidproc%ee_g(4)
	DO agei = 5,10
		lfavgeproc = lfavgeproc + hhidproc%ee_g(agei)
		lfearnproc = lfearnproc + hhidproc%ee_g(agei)/(1d0+rr)**(agei-4)
	ENDDO
	lfavgeproc = lfavgeproc/7d0
	lfwlthproc = lfearnproc + hhidproc%ss_g(4)
!----------------------------
	CALL MPI_ALLGATHER(lfavgeproc,NIproc,MPI_REAL8,lfgavge,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfearnproc,NIproc,MPI_REAL8,lfgearn,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(lfwlthproc,NIproc,MPI_REAL8,lfgwlth,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!----------------------------

	DEALLOCATE( lfavgeproc,lfearnproc,lfwlthproc )


!----------------------------
! INTERIM: SON'S CONSUMPTION INEQUALITY
!----------------------------
	consproc = LOG(MAX(zero,consproc))
	DO agei = 4,10
		ccountproc(agei) = COUNT(consproc(:,agei)>lzero)
		cmeanproc(agei) = SUM(consproc(:,agei),consproc(:,agei)>lzero)
	ENDDO

	CALL MPI_ALLREDUCE(ccountproc,ccount,7,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLREDUCE(cmeanproc,cmean,7,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
	cmean = cmean/DBLE(ccount)

	DO agei = 4,10
		cstdproc(agei) = SUM( (consproc(:,agei)-cmean(agei))**2,consproc(:,agei)>lzero )
	ENDDO

	CALL MPI_ALLREDUCE(cstdproc,cstd,7,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
	cstd = cstd/DBLE(ccount)
	cstd = SQRT(cstd)

	IF (myid==0) THEN
		OPEN(999,FILE="results/cons_prof.out")
		WRITE(999,'(20F8.3)'),cmean(4),cmean-cmean(4)
		WRITE(999,'(20F8.3)'),cstd(4),cstd-cstd(4)
		CLOSE(999)
	ENDIF
!----------------------------
	DEALLOCATE( consproc )

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE onegen(ind)
USE functions
IMPLICIT NONE
	INTEGER,INTENT(IN) :: ind

	INTEGER :: ai,aii,ei,agei
	INTEGER :: si,hi,ci,cii
	INTEGER :: si1,hi1,hki1
	REAL(dp),DIMENSION(2) :: swgts,hwgts,hkwgts

	REAL(dp) :: ap,ak,ep
	REAL(dp) :: l0,l1,l2,n3,n4,n5,n6,n7,n8,n9
	REAL(dp) :: h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
	REAL(dp) :: s4,s5,s6,s7,s8,s9,s10,sc
	REAL(dp) :: e3,e4,e5,e6,e7,e8,e9,e10

	REAL(dp) :: wp,gamp,budget,temp(2),pmet(2)
	REAL(dp) :: wk,gamk,mt,xt

	REAL(dp) :: b10,b11,c10,c11,c12,s11,s12

! SON'S STATES AT AGE 4
		ci = hhidproc(ind)%col_k
		ai = hhidproc(ind)%ability_k
		ap = agrid(ai)

		h4 = hhidproc(ind)%hh_k(4)
		s4 = hhidproc(ind)%ss_k(4)

		n4 = hhidproc(ind)%nn_k(4)
		e4 = hhidproc(ind)%ee_k(4)

		h5 = hhidproc(ind)%hh_k(5)
		s5 = hhidproc(ind)%ss_k(5)

		budget = atinc(e4,0d0,4) + s4-smin
		budget = MAX(0d0,budget)
		consproc(ind,4) = budget - s5

! age 5: kid is born

		wp = wskill(ci)*h5

		cii = hhidproc(ind)%col_g
		aii = hhidproc(ind)%ability_g
		ak = agrid(aii)

		n5 = hhidproc(ind)%nn_k(5)
		l0 = hhidproc(ind)%nn_k(0)
		e5 = hhidproc(ind)%ee_k(5)

		mt = wp*l0 *(1d0-kamm(0))/kamm(0)
		h1 = hhidproc(ind)%hh_g(1)

		h6 = hhidproc(ind)%hh_k(6)
		s6 = hhidproc(ind)%ss_k(6)

		budget = atinc(e5-mt,s5,5) +s5-smin
		budget = MAX(0d0,budget)
		consproc(ind,5) = budget - s6

! age 6: primary school
		wp = wskill(ci)*h6

		n6 = hhidproc(ind)%nn_k(6)
		l1 = hhidproc(ind)%nn_k(1)
		e6 = hhidproc(ind)%ee_k(6)

		mt = wp*l1 *(1d0-kamm(1))/kamm(1)
		h2 = hhidproc(ind)%hh_g(2)

		h7 = hhidproc(ind)%hh_k(7)
		s7 = hhidproc(ind)%ss_k(7)

		budget = atinc(e6-mt,s6,6) +s6 -smin
		budget = MAX(0d0,budget)
		consproc(ind,6) = budget - s7

! age 7: secondary school
		wp = wskill(ci)*h7

		n7 = hhidproc(ind)%nn_k(7)
		l2 = hhidproc(ind)%nn_k(2)
		e7 = hhidproc(ind)%ee_k(7)

		mt = wp*l1 *(1d0-kamm(1))/kamm(1)
		h3 = hhidproc(ind)%hh_g(3)

		h8 = hhidproc(ind)%hh_k(8)
		s8 = hhidproc(ind)%ss_k(8)

		budget = atinc(e7-mt,s7,7) +s7 -smin
		budget = MAX(0d0,budget)
		consproc(ind,7) = budget - s8

! age 8: kid's college choice

		sc = s8 - MERGE(0d0,kappa,cii==1)
		wk = wskill(cii)*h3
		gamk = gamm(cii)

! age 8.1: kid's work choice

		n8 = hhidproc(ind)%nn_k(8)
		e8 = hhidproc(ind)%ee_k(8)

		h9 = hhidproc(ind)%hh_k(9)
		s9 = hhidproc(ind)%ss_k(9)

		n3 = hhidproc(ind)%nn_g(3)
		e3 = hhidproc(ind)%ee_g(3)

		h4 = hhidproc(ind)%hh_g(4)

		budget = atinc(e8,sc,8) +sc -smin
		budget = atinc(e3,0d0,3) +budget
		budget = MAX(0d0,budget)
		consproc(ind,8) = budget - s9

! age 9: bequests

		n9 = hhidproc(ind)%nn_k(9)
		e9 = hhidproc(ind)%ee_k(9)

		h10 = hhidproc(ind)%hh_k(10)
		s10 = hhidproc(ind)%ss_k(10)

		s4 = hhidproc(ind)%ss_g(4)

		! savings choice
		budget = atinc(e9,s9,9) +s9 -smin -skmin
		budget = MAX(0d0,budget)
		consproc(ind,9) = budget -s10 -s4

!--------------------------
! RETIREMENT
		e10 = wskill(ci)*h10
		b10 = atinc(e10,s10,10) +s10
		c10 = b10 + (2d0+(1d0-rtax)*rr)/(1d0+(1d0-rtax)*rr)**2 *(p0+p1*e10+gg)
		c10 = c10 / perm_cons_denom

		s11 = b10 - c10
		b11 = (1d0+(1d0-rtax)*rr)*s11 + p0+p1*e10 + gg
		c11 = c10*(bta*(1d0+(1d0-rtax)*rr))**(1d0/chi)
		s12 = b11 - c11

		consproc(ind,10:11) = (/c10,c11/)
		consproc(ind,12) = c11*(bta*(1d0+(1d0-rtax)*rr))**(1d0/chi)

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
END MODULE
