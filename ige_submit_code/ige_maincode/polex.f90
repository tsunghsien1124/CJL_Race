MODULE polex
USE simulation
USE valuefuncs
USE compute_iter
USE printresults
IMPLICIT NONE

	PRIVATE
	REAL(dp) :: rstar
	PUBLIC :: shortrunpe, longrunpe, longrunge


!---------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------
SUBROUTINE printpolex(polindx,exptype)
IMPLICIT NONE

	INTEGER,INTENT(IN) :: polindx,exptype
	CHARACTER(LEN=30) :: filename

	IF (myid==0) THEN
		WRITE(*,*),'POLEX',polindx,exptype

		WRITE(filename,'(A14,I2,I2,A4)'),"results/polex_",polindx,exptype,".out"
		OPEN(33,FILE=filename,STATUS="REPLACE")
		CALL printmoments(33)
		WRITE(33,'(A8,2F8.3)'),'tax,dd',(/tau0,SUM(dd)/)
		IF (exptype==30) WRITE(33,'(A8,F8.3)'),'IR',1D2*rstar
		CLOSE(33)
	ENDIF

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE shortrunpe(polindx)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: polindx

	DO ind=1,NIproc

! fix distribution for each generation
		ai = hhidproc(ind)%ability_p
		ap = agrid(ai)

		h4 = hhidproc(ind)%hh_p(4)
		s4 = hhidproc(ind)%ss_p(4)
		CALL locatesk(s4,swgts,si1)
		CALL locateh(h4,hwgts,hi1)

		ci = hhidproc(ind)%col_p
		gamp = gamm(ci)
		wp = wskill(ci)*h4

		! earnings choice
		temp(1) = DOT_PRODUCT(swgts,n4fun(si1:si1+1,hki1,ai,ci))
		temp(2) = DOT_PRODUCT(swgts,n4fun(si1:si1+1,hki1+1,ai,ci))
		n4 = DOT_PRODUCT(hwgts,temp)
		n4 = MAX(0d0,MIN(n4,1d0))

		e4 = MAX(0d0,wp*(1d0-n4))

		h5 = ap*(n4*h4)**gamp + h4
		h5 = egrid(epi(ind,TT-1,5),ci) *h5
		h5 = MAX(hmin,MIN(h5,hmax))

		! savings choice
		i4 = atinc(e4,0d0,4)
		budget = i4 +s4 -smin
		budget = MAX(0d0,budget)

		temp(1) = DOT_PRODUCT(swgts,s5fun(si1:si1+1,hki1,ai,ci))
		temp(2) = DOT_PRODUCT(swgts,s5fun(si1:si1+1,hki1+1,ai,ci))
		s5 = DOT_PRODUCT(hwgts,temp) *budget
		s5 = MAX(0d0,MIN(s5,smax)) +smin

		CALL simul_end

	ENDDO
!------------------------
	CALL MPI_ALLGATHER(hhidproc,NIproc,MPI_AGENT,hhid,NIproc,MPI_AGENT,MPI_COMM_WORLD,ierr)
!------------------------
	CALL MPI_ALLGATHER(colfpproc,NIproc,MPI_REAL8,colfp,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(colfkproc,NIproc,MPI_REAL8,colfk,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
	CALL MPI_ALLGATHER(colfgproc,NIproc,MPI_REAL8,colfg,NIproc,MPI_REAL8,MPI_COMM_WORLD,ierr)
!------------------------
	CALL compute_moments

	CALL printpolex(polindx,10)

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE longrunpe(polindx)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: polindx

	CALL vfiter
	CALL mcmc
	CALL compute_moments

	CALL printpolex(polindx,20)

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
SUBROUTINE longrunge(polindx)
USE grids
IMPLICIT NONE
	INTEGER,INTENT(IN) :: polindx

	REAL(dp),PARAMETER :: ptol=1d-1
	INTEGER,PARAMETER :: rn=10,wn=10
	INTEGER :: riter,witer,i

	REAL(dp) :: rmin,rmax
	REAL(dp) :: wmin,wmax

	wmin = zero
	wmax = infty
	DO witer=1,wn;

		rstar = r0
		rmin = 0d0
		rmax = 1d-1
		DO riter=1,rn

			rr = (1d0+rstar)**per-1d0;

			Windx = (1d0+rstar+delk)**per-1d0
			Windx = (1d0-alph)*(alph/Windx)**(alph/(1d0-alph))
			Windx = Windx * tfp

			hmax = nmdlavgearn/Windx*hfrac
			CALL gen_hgrid(hmax)
			hkmax = (hmax*hkfrac)
			CALL gen_hkgrid(hkmax)

			smin = bc_lc/(1d0+rr);
			CALL gen_sgrid(smin,smax)

			perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0)
			perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0) *perm_cons_denom

			wskill(1) = (1d0-ups)**(1d0/(1d0-ces)) *ww**(ces/(ces-1d0))
			wskill(1) = ups**(1d0/(1d0-ces)) + wskill(1)
			wskill(1) = Windx *wskill(1)**((1d0-ces)/ces)
			wskill(2) = wskill(1) *ww

			DO i = 0,2
				lamb(i,:) = lamb0(i)/wskill**kamm(i)
			ENDDO

			CALL vfiter
			CALL mcmc
			CALL compute_prices

			IF (myid==0) THEN
				CALL printline(6)
				WRITE(*,'(2I4,A8,2F8.3,A8,2F8.3)'),witer,riter,'r01:', 1d2*(/rstar,r1/)
				CALL printline(6)
			ENDIF

			IF (ABS(r1/rstar-1d0)<ptol) EXIT

			rmin = MAX(rmin,MIN(r1,rstar))
			rmax = MIN(rmax,MAX(r1,rstar))

			rstar = (rmin+rmax)/2d0

		ENDDO

		IF (myid==0) THEN
			CALL printline(6)
			WRITE(*,'(I8,A4,2F8.3,A8,2F8.3)'),witer,'ww01:', 1d2*(/ww,ww1/)
			CALL printline(6)
		ENDIF

		IF (ABS(ww1/ww-1d0)<ptol) EXIT

		wmin = MAX(wmin,MIN(ww1,ww))
		wmax = MIN(wmax,MAX(ww1,ww))

		ww = (wmin+wmax)/2d0

! need to delete this block in revision?
!		ups = colshare/(1d0-colshare) *eprem
!		ups = ww**ces *ups**(1d0-ces)
!		ups = 1d0/(1d0+ups)

	ENDDO

	CALL compute_moments

	CALL printpolex(polindx,30)

RETURN
END SUBROUTINE
!---------------------------------------------------------------------
END MODULE
