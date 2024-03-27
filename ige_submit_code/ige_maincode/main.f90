! PHI COMES IN VALUE FUNCS AND SIMULATIONS ONLY

PROGRAM MAIN
!USE grids 		! change to equilibrium later
USE distance_func
!USE compute_iter, ONLY: compute_check
USE simplex
USE lifetime_stats
USE bt_module
!USE hm_module
USE polex
IMPLICIT NONE
! simplex variables
	REAL(dp) :: step(nop),rhobeg
	INTEGER :: iflg,ind,bigi,i,polindx

! experiments
	REAL(dp) :: tearn,tsave,temp

!-----------------------------------------------------------
	CALL parallelize
	CALL initialize
	
!-----------------------------------------------------------
! CALIBRATE: POLICY PARAMS CAN BE FIXED
! avgearn assumed to be 100
!-----------------------------------------------------------
! tax params
	tau0 = taues(1); 
	tau1 = taues(2); 
	tauS = taues(3);
! transfer params
	gg = gp01(1) *nmdlavgearn;
	p0 = gp01(2) *nmdlavgearn;
	p1 = gp01(3);
! educ subsidies
	dd(0:2) = d012(1:3) *nmdlavgearn;
	kappa = ccost 		*nmdlavgearn;

!-----------------------------------------------------------
! target fixed interest rate and premium; set average wage index
	rr = (1d0+r0)**per-1d0;
	Windx0 = (1d0+r0+delk)**per-1d0
	Windx0 = (1d0-alph)*(alph/Windx0)**(alph/(1d0-alph))

! fix value function sgrid (depends on rr),skgrids (0d0 for bequest bound)
	bc_lc =  -gg;
	smin = bc_lc/(1d0+rr)
	smax = nmdlavgearn*sfrac
	CALL gen_sgrid(smin,smax)

	skmin = 0d0
	skmax = smax*skfrac
	CALL gen_skgrid(skmin,skmax)
	
	CALL gen_xegrids

!-----------------------------------------------------------
DO bigi = 1,bign
!-----------------------------------------------------------
IF (bigi>1) THEN
	IF (myid==0) THEN
		WRITE(*,*),'AUTO',bigi
		WRITE(100,*),'AUTO',bigi
	ENDIF

	plb(1)  = params(1) -0.03d0; 	pub(1)  = params(1) +0.03d0
	plb(2)  = params(2) -0.03d0; 	pub(2)  = params(2) +0.03d0
	plb(3)  = params(3) -.125d0; 		pub(3)  = params(3) +.125d0
	plb(4)  = params(4) -0.05d0; 	pub(4)  = params(4) +0.05d0
	plb(5)  = params(5) -0.05d0; 	pub(5)  = params(5) +0.05d0
	plb(6)  = params(6) -0.03d0; 	pub(6)  = params(6) +0.03d0
	plb(7)  = params(7) -0.03d0; 	pub(7)  = params(7) +0.03d0
	plb(8)  = params(8) -0.05d0; 	pub(8)  = params(8) +0.05d0
	plb(9)  = params(9) -0.02d0; 	pub(9)  = params(9) +0.02d0
	plb(10) = params(10)-0.02d0; 	pub(10) = params(10)+0.02d0
	plb(11) = params(11)-0.001d0; 	pub(11) = params(11)+0.001d0
	plb(12) = params(12)-0.03d0; 	pub(12) = params(12)+0.03d0
	plb(13) = params(13)-0.5d0; 	pub(13) = params(13)+0.5d0

	plb(1) =  MAX(zero,plb(1));	pub(1) =  MIN(.4d0,pub(1))
	plb(2) =  MAX(zero,plb(2));	pub(2) =  MIN(.22d0,pub(2))
	plb(3) =  MAX(zero,plb(3));
	plb(4) =  MAX(zero,plb(4));	pub(4) =  MIN(1d0,pub(4))
	plb(5) =  MAX(zero,plb(5));	pub(5) =  MIN(1d0,pub(5))

	plb(6) =  MIN(.899d0,MAX(.7d0,plb(6)))	
	pub(6) =  MAX(.701d0,MIN(.9d0,pub(6)))

	plb(9) =  MAX(zero,plb(9));		pub(9) =  MIN(.9d0,pub(9))
	plb(10) = MAX(zero,plb(10));	pub(10) = MIN(.4d0,pub(10))
	plb(11) = MAX(.93d0,plb(11));	pub(11) = MIN(.995d0,pub(11))
	plb(12) = MAX(.5d0,plb(12));	pub(12) = MIN(1.5d0,pub(12));
	plb(13) = MAX(zero,plb(13));

ENDIF
!-----------------------------------------------------------
! search for parameters ::: CHANGE TO NEWUOA AFTER REGION FOUND
!-----------------------------------------------------------
	iflg = 0
	params = LOG(params-plb) - LOG(pub-params)
	IF (switch(1)==2) THEN
		caliter = -1; step = MERGE(1D0,EXP(1d0),bigi>1)

		CALL NMsimplex(params,distance,step,ctol,func,itmax,iflg)
!		rhobeg = MAX(EXP(1d0),ctol)
!		CALL newuoa(func,params,rhobeg,ctol,itmax,0,100)
	ENDIF
!----------------------------------------
	caliter = -1;
	distance = func(params)
	params = (pub-plb)*( EXP(params)/(1d0+EXP(params)) ) + plb

!-----------------------------------------------------------
ENDDO
!-----------------------------------------------------------
! BOOKKEEPING RESULTS
!----------------------------------------
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL compute_check
!----------------------------------------
!	CALL MPI_ALLGATHER(hhidproc,NIproc,MPI_AGENT,hhid,NIproc,MPI_AGENT,MPI_COMM_WORLD,ierr)
!----------------------------------------
	IF (myid==0) THEN
		IF (switch(1)/=1) CALL printguess
		IF (switch(1)>1) CALL printguessreal
!
!		CALL printearnings
		CALL printgrids
!		CALL printvalues
!		CALL printpolicy
	ENDIF
!-----------------------------------------------------------


!-----------------------------------------------------------

!-----------------------------------------------------------
! EXPERIMENTS
! EQM IS STORED IN HHID
!-----------------------------------------------------------
! 1. lifetime k and p, and consumption
!-----------------------

	IF (switch(2)>=1) THEN

! parallel
		CALL gen_lifetime_stats
		IF (myid==0) THEN
			CALL printlifetime
			CALL printsimulation
		ENDIF
	ENDIF

!-----------------------
! 2. lifetime profile, and bt-hm experiments
!-----------------------

	IF (switch(2)>=2 .AND. myid==0) THEN
!		CALL compute_netinc
		CALL printnetinc
		CALL bt_calibrate
	ENDIF

!	IF (switch(2)>=3 .AND. myid==0) THEN
!		CALL hm_calibrate
!	ENDIF

!--------------------------------------------------------------------------------------------
! 3. policy experiments: all are parallelized
!--------------------------------------------------------------------------------------------

! flat tax
	IF (switch(3)==1) THEN				! .OR. switch(3)==5) THEN

		IF (myid==1) THEN

			tearn = 0d0; 		tsave = 0d0
			DO ind = 1,NI
				tearn = SUM(hhid(ind)%ee_k(4:10)) + hhid(ind)%ee_g(3) + tearn
				tsave = SUM(hhid(ind)%ss_k(5:10)) + tsave
			ENDDO

			tau0 = SUM(netinck(4:10,:)) + tauS*tearn - gg*(5d0 + 3d0*aeq)*NI
			tau1 = tearn + rr*tsave
			tau0 = tau0/tau1
			tau0 = 1d0-tau0
		ENDIF
		CALL MPI_BCAST(tau0,1,MPI_REAL8,1,MPI_COMM_WORLD,ierr)		
		tau1 = 0d0
		
		CALL polex_compute_print(10)
		tau0 = taues(1); 
		tau1 = taues(2); 
	ENDIF

! all 0
	IF (switch(3)==2 .OR. switch(3)==5) THEN

		dd = 0d0
		dd(0) = SUM(d012) 	!*nmdlavgearn;

!		temp = 1d0 - 0.1d0*d012(1)/SUM(d012(2:3))
!		dd = d012*temp*nmdlavgearn
!		dd(0) = d012(1)*1.1d0 *nmdlavgearn

if (myid==0) print*,phi,dd

		CALL polex_compute_print(11)
		dd(0:2) = d012(1:3) *nmdlavgearn;
	ENDIF

! all 1
	IF (switch(3)==3 .OR. switch(3)==5) THEN

		dd = 0d0
		dd(1) = SUM(d012) *nmdlavgearn;

!		temp = 1d0 - 0.1d0*d012(2)/(d012(1)+d012(3))
!		dd = d012*temp*nmdlavgearn
!		dd(1) = d012(2)*1.1d0 *nmdlavgearn

if (myid==0) print*,phi,dd

		CALL polex_compute_print(12)
		dd(0:2) = d012(1:3) *nmdlavgearn;
	ENDIF

! all 2
	IF (switch(3)==4 .OR. switch(3)==5) THEN

		dd = 0d0
		dd(2) = SUM(d012) *nmdlavgearn;

!		temp = 1d0 - 0.1d0*d012(3)/SUM(d012(1:2))
!		dd = d012*temp*nmdlavgearn
!		dd(2) = d012(3)*1.1d0 *nmdlavgearn

if (myid==0) print*,phi,dd

		CALL polex_compute_print(13)
		dd(0:2) = d012(1:3) *nmdlavgearn;
	ENDIF


!--------------------------------------------------------------------------------------------
! 4. financial markets
!--------------------------------------------------------------------------------------------


! relax ig constraint
	IF (switch(4)==1 .OR. switch(4)==5) THEN

		skmin = -gg
		CALL gen_skgrid(skmin,skmax)

		CALL polex_compute_print(40)
	ENDIF

! relax lifecycle constraint
	IF (switch(4)==2 .OR. switch(4)==5) THEN

		gg = 2d0*gg
		bc_lc = -gg
		smin = bc_lc/(1d0+rr);
		CALL gen_sgrid(smin,smax)

		CALL polex_compute_print(41)
	ENDIF

! relax both
	IF (switch(4)==3 .OR. switch(4)==5) THEN

		gg = 0d0
		bc_lc = -gg
		smin = bc_lc/(1d0+rr);
		CALL gen_sgrid(smin,smax)

		CALL polex_compute_print(42)
	ENDIF

! eliminate lumpsum and tighten
	IF (switch(4)==4 .OR. switch(4)==5) THEN

		gg = 0d0
		bc_lc = -gg
		smin = bc_lc/(1d0+rr);
		CALL gen_sgrid(smin,smax)

		CALL polex_compute_print(43)
	ENDIF

!-----------------------------------------------------------
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!-----------------------------------------------------------
	IF (myid==0) CLOSE(100)
	CALL MPI_FINALIZE(ierr); STOP
!-----------------------------------------------------------

CONTAINS
!-----------------------------------------------------------
SUBROUTINE parallelize
IMPLICIT NONE
	INTEGER :: BLen(2)=(/6,111/)
	INTEGER :: Typs(2)=(/MPI_INTEGER,MPI_DOUBLE_PRECISION/)
	INTEGER(KIND=MPI_ADDRESS_KIND) :: LB,EXTENT
	INTEGER(KIND=MPI_ADDRESS_KIND) :: DSpl(2)

	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

	NIproc = NI/nprocs
	ALLOCATE( hhidproc(NIproc) )
	ALLOCATE( colfpproc(NIproc), colfkproc(NIproc), colfgproc(NIproc) )
	ALLOCATE( netincpproc(4:12,NIproc), netinckproc(4:12,NIproc), netincgproc(4:12,NIproc) )

	ALLOCATE( abi(NIproc,TT+2),epi(NIProc,TT+2,4:10) )
!	ALLOCATE( earnvec(NIproc,3:10) )

	nwp10 =		sspace10/nprocs; 
	nwp89 = 	sspace89/nprocs; 
	nwp678 = 	sspace678/nprocs; 
	nwp5 =		sspace5/nprocs
	nwp4 = 		sspace4/nprocs;

	ts10 = 	myid*nwp10; 	tf10 =	(myid+1)*nwp10-1
	ts89 = 	myid*nwp89; 	tf89 = 	(myid+1)*nwp89-1
	ts678 = myid*nwp678; 	tf678 = (myid+1)*nwp678-1
	ts5 =	myid*nwp5;		tf5 =	(myid+1)*nwp5-1
	ts4 = 	myid*nwp4; 		tf4 = 	(myid+1)*nwp4-1

! individuals	
	ALLOCATE(v10split(nwp10))
	ALLOCATE(w8v9split(nwp89),s910split(nwp89),n89split(nwp89))
	ALLOCATE(v678split(nwp678),s78split(nwp678),n67split(nwp678))
	ALLOCATE(v5split(nwp5),s6split(nwp5),n5split(nwp5))		
	ALLOCATE(v4split(nwp4),s5split(nwp4),n4split(nwp4))		
! kids
	ALLOCATE(l0split(nwp5),l12split(nwp678),c3split(nwp678),n3s4split(nwp89))

! MPI DERIVED TYPE
	CALL MPI_GET_ADDRESS(hhid(1)%col_p,DSpl(1),ierr)
	CALL MPI_GET_ADDRESS(hhid(1)%hh_p(1),DSpl(2),ierr)
	DSpl=DSpl-DSpl(1)
	CALL MPI_TYPE_CREATE_STRUCT(2,BLen,DSpl,Typs,MPI_AGENT,ierr)
	CALL MPI_TYPE_COMMIT(MPI_AGENT,ierr)

! check if allocated correctly
!	CALL MPI_TYPE_GET_EXTENT(MPI_AGENT,LB,EXTENT,IERR)
!	IF (myid==0) print*,lb,extent,dspl

RETURN
END SUBROUTINE
!-----------------------------------------------------------
SUBROUTINE initialize
IMPLICIT NONE

	INTEGER :: i

	IF (myid==0) THEN
		OPEN(11,FILE="output/guessreal.out")
		READ(11,*),params
		CLOSE(11)

		OPEN(22,FILE="output/guess.out")
		DO i=1,3;READ(22,*);ENDDO
		READ(22,*),plb
		READ(22,*),pub
		DO i=1,2;READ(22,*);ENDDO
!		READ(22,*),weights
		READ(22,*),outweights
		DO i=1,4;READ(22,*);ENDDO
		READ(22,*),hfrac,sfrac,hkfrac,skfrac
		DO i=1,2;READ(22,*);ENDDO
		READ(22,*),itmax,bign,switch
		CLOSE(22)

		start = MPI_wtime()
		OPEN(100,FILE="output/logfile.out",STATUS='replace')
	ENDIF
	CALL MPI_BCAST(params,nop,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(plb,nop,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(pub,nop,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!	CALL MPI_BCAST(weights,nop,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(outweights,nout,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(hfrac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(sfrac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(hkfrac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(skfrac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(itmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(bign,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(switch,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE

!-----------------------------------------------------------
! compute and print
SUBROUTINE polex_compute_print(polindx)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: polindx
	REAL(dp) :: calww

	calww = ww

	CALL shortrunpe(polindx)
	CALL longrunpe(polindx)
	CALL longrunge(polindx)

	rr = (1d0+r0)**per-1d0;

	gg = gp01(1) *nmdlavgearn;


	bc_lc =  -gg;
	smin = bc_lc/(1d0+rr)
	CALL gen_sgrid(smin,smax)

	skmin = 0d0
	CALL gen_skgrid(skmin,skmax)

	Windx = Windx0 *tfp
	hmax = nmdlavgearn/Windx*hfrac
	CALL gen_hgrid(hmax)
	hkmax = (hmax*hkfrac)
	CALL gen_hkgrid(hkmax)

	ww = calww
	ups = colshare/(1d0-colshare) *eprem
	ups = ww**ces *ups**(1d0-ces)
	ups = 1d0/(1d0+ups)

	perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0)
	perm_cons_denom = 1d0 + bta**(1d0/chi) *(1d0+(1d0-rtax)*rr)**(1d0/chi-1d0) *perm_cons_denom

	wskill(1) = (1d0-ups)**(1d0/(1d0-ces)) *ww**(ces/(ces-1d0))
	wskill(1) = ups**(1d0/(1d0-ces)) + wskill(1)
	wskill(1) = Windx *wskill(1)**((1d0-ces)/ces)
	wskill(2) = wskill(1) *ww

	DO i = 0,2
		lamb(i,:) = lamb0(i)/wskill**kamm(i)
	ENDDO
RETURN
END SUBROUTINE
!-----------------------------------------------------------
END PROGRAM MAIN
