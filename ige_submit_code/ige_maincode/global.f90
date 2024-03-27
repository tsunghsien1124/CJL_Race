!-----------------------------------------------------------------------------------------------------------
!!!!! WHEN DOING REVISION, CREATE AGE-DEPENDENT GRIDS?
!-----------------------------------------------------------------------------------------------------------

MODULE GLOBAL
USE MPI
IMPLICIT NONE
	PUBLIC
	INTEGER,PARAMETER :: dp = KIND(1d0)
	REAL(dp),PARAMETER :: goldrat = (1d0+sqrt(5d0))/2d0

! all in 2000 dollars
! PSID: avg hh earnings 38331.43 dollars // last period: 42578.63 // 36 years: 42971.49
! targets: average opportunity time cost, anchor: college enrollment
	REAL(dp),PARAMETER :: psidavgearn = 38674.2d0
	REAL(dp),PARAMETER :: nmdlavgearn = 1d0,zbeq=3d3/psidavgearn

	REAL(dp),PARAMETER :: cutb30 = 1d3/psidavgearn *4.07464954416234d0
	REAL(dp),PARAMETER :: cuta30 = 1.5d3/psidavgearn *4.07464954416234d0
	REAL(dp),PARAMETER :: timeb30 = 260d0/1400d0
	REAL(dp),PARAMETER :: timea30 = 520d0/1400d0


! SWITCH*4:
! 1: CALIBRATE?					0 (NO) 1 (NO,EQM)		2 (YES)
! 2: EXPERIMENTS:PARAMS:0-4?	0 (NO) 1 (YES FOR ALL)	2 (ADLT PARAMS)
!								3 (SHOCK PARAMS)		4 (ECE PARAMS)
! 3: EXPERIMENTS:GOVERN:0-4?	0 (NO) 1 (YES FOR ALL)	2 (ADLT POL)
!								3 (COLLEGE)				4 (ECE PARAMS)
!-----------------------------------------------------------------------------------------------------------
! PARAMETERS FOR COMPUTATION
!------------------------------------------------------------------------------------------------------------------
	INTEGER :: itmax,bign,switch(4)
	INTEGER :: iflt

	REAL(dp),PARAMETER :: infty=1d16, zero=1d-16 !HUGE(infty) EPSILON(zero)
	REAL(dp),PARAMETER :: lzero = LOG(zero)
	INTEGER,PARAMETER :: nop=13, maxfn=100, nom=11

	REAL(dp),PARAMETER :: rtol=5d-3, wtol=5d-2, atol=1d-1
	REAL(dp),PARAMETER :: ctol=1d-2, vtol=1d-1, tol=1d-8
	INTEGER,PARAMETER :: pn=15
	INTEGER,PARAMETER :: NI=120000, TT=200	! 10000 for each proc; TT generations

	INTEGER,PARAMETER :: agen=12,cn=2,an=3,en=2,hn=12,sn=24
	INTEGER,PARAMETER :: powh=2,pows=2

	REAL(dp),PARAMETER :: adjv = 1d0/goldrat ! .5d0 !

!	REAL(dp),PARAMETER :: beq_gdp = 2.65d0						! hendricks
!	REAL(dp),PARAMETER :: beq_nw = 0.31d0						! gale scholz 51.8 45
	REAL(dp),DIMENSION(nop) :: plb,pub,params,params_exp
	REAL(dp),DIMENSION(nom) :: moments,weights

	INTEGER,PARAMETER :: nout = 9
	REAL(dp),DIMENSION(nout) :: outweights

	REAL(dp) :: start,duration,finish,distance
	INTEGER :: piter,caliter
!------------------------------------------------------------------------------------------------------------------
! FIXED PARAMETERS: mostly from PSID
!------------------------------------------------------------------------------------------------------------------
	INTEGER,PARAMETER :: per=6
	REAL(dp),PARAMETER :: chi=2d0
	REAL(dp),PARAMETER :: aeq=1.7d0

	REAL(dp),PARAMETER :: alph=.33d0, delk=.067d0	!.376d0, delk=0.059d0 !, Dk=(1d0-delk)**per
	REAL(dp),PARAMETER :: ces = 1d0-1d0/1.441d0
	REAL(dp),PARAMETER :: r0=4d-2

	REAL(dp),PARAMETER :: nmdlavginc = nmdlavgearn* (1d0+ ((1d0+r0)**6-1d0)/((1d0+r0+delk)**6-1d0) *alph/(1d0-alph))

! tax/transfers
	REAL(dp),PARAMETER :: rtax = 0.15d0,taxgdp = 0.29d0
	REAL(dp),DIMENSION(3),PARAMETER :: taues = (/ 0.106d0,0.035d0,0.1111d0 /)
	!0d0,0d0,0d0,0d0/) !0d0,0d0,0d0/) ! 0d0 !2d0/15d0

	REAL(dp),DIMENSION(3),PARAMETER :: gp01 = (/ 0.02d0*nmdlavginc/nmdlavgearn,0.0687d0,0.3353d0/)
	! 0d0,0d0/) !0.02d0 !
! world bank vs PSID
!	REAL(dp),DIMENSION(3),PARAMETER :: d012 = (/0.1912d0/5d0,0.1912d0,0.22628d0/)/(1d0-taxgdp)
	REAL(dp),DIMENSION(3),PARAMETER :: d012 = (/908.3698D0,3493.836D0,3845.92D0/)/psidavgearn

! college costs: from EDUC STATS and CPS (hcaprisk)
	REAL(dp),PARAMETER :: ccost=7059.08333333333D0/psidavgearn
	REAL(dp),PARAMETER :: ctime(2,2) = (/0d0,1d0,	1d0/2d0,1d0/2d0/)
	REAL(dp) :: gg,p0,p1,kappa,dd(0:2)

! fixed parameters from CDS/PSID
	REAL(dp),PARAMETER :: kamm(0:2) = (/0.9033114d0,0.7108408d0,0.6771832d0/)
	REAL(dp),PARAMETER :: lamb0(0:2) = kamm**kamm * (1d0-kamm)**(1d0-kamm)
	REAL(dp) :: lamb(0:2,2)

!	REAL(dp),PARAMETER :: phi0 = 1d0
!	REAL(dp),PARAMETER :: phi(0:2) = (/1.0056035d0,-0.00376753d0,0.00147892d0/)
!	REAL(dp),PARAMETER :: phi(0:2) = (/0.18669905d0,-0.00402475d0,0.039025d0/)	!log(testscore)

	REAL(dp),PARAMETER :: phi(1:2) = 0.5d0

	REAL(dp),PARAMETER :: sig_e(2) = (/.122722d0,.2071547d0/)
!	REAL(dp),PARAMETER :: sig_e = .1770161d0	!SQRT((0.111d0**2)*6d0)

	REAL(dp) :: omega(2)
!	REAL(dp),PARAMETER :: omega(2) = (/0.6356687d0,0.4531415d0/)	! COBB-DOUGLAS
!	REAL(dp),PARAMETER :: omega(2) = (/0.67924753d0,0.42824605d0/)	! CES

! calibration targets from PSID
! slopes: age 48/24 for hs, 54/24 for col
! log e (10: age 60-65) earnings var
	REAL(dp),PARAMETER :: eprem = 1.519012d0	! this is exponentiate then sum
												! eprem7 is sum then exponentiate
	REAL(dp),PARAMETER :: colshare = 0.4798591d0

	! nlsy premia: 1.46786626627134d0
	! nlsy enrollment: 41.35666d0,6.85d0,0.32 // gap: 37.5888208d0 or 35.6494491d0
	! nlsy: 17.63729d0/12.50731d0,23.61649d0/15.3677d0 /) ! (/1.410158539,1.536761519/)

	REAL(dp),PARAMETER :: igtfr = 0.208d0	!+0.31d0	! intended transfers ! 208.0, 0.120, 0.310d0

!------------------------------------------------------------------------------------------------------------------
! CALIBRATED PARAMETERS / POLICY VARIABLES
!------------------------------------------------------------------------------------------------------------------
	REAL(dp) :: thet,qq			! altruism and implied utility pareto weights
	REAL(dp) :: rh_a			! ability and time investment in children
	REAL(dp) :: sig_a		! shock-corr with ability, ability shock sd
	REAL(dp) :: gamm(2) 		! ben-porath for hs-col ! (/0.832d0,0.871d0/) ! 0d0
	REAL(dp) :: zcols,zcolm		! z-col pref correlation
	REAL(dp) :: mloga
!	REAL(dp) :: nu				!=0.9d0

	REAL(dp) :: mu_a,sig_eta

! parameters calibrated in eqm: rr is fixed at r0. x1's are implied values.
! normalize avg earnings to 1. ww is w2/w2 (col/hs)
	REAL(dp) :: bta,ups
	REAL(dp) :: zeta
	REAL(dp) ::	rr, ww, Windx, Windx0
	REAL(dp) :: r1, ww1
	REAL(dp) :: ups1,eprem1,avgearn

	REAL(dp) :: wskill(2),perm_cons_denom

! taxes & transfers & SS for policy experiments
	REAL(dp) :: tau0,tau1,tauS
	REAL(dp) :: avgtax

! implied moments
	REAL(dp) :: tfp,gdp,cap,lab(2)
	REAL(dp) :: beq,earn12,cinv
	REAL(dp) :: enroll,prof(2),logvar(2)

! borrowing constraints
	REAL(dp) :: bc_lc

!------------------------------------------------------------------------------------------------------------------
! VF'S AND GRIDS
!------------------------------------------------------------------------------------------------------------------
	REAL(dp) :: hmin,hkmin,smin,skmin
	REAL(dp) :: hmax,hkmax,smax,skmax
	REAL(dp) :: hbins(hn+1),hkbins(hn+1)
	REAL(dp) :: sbins(sn+1),skbins(sn+1)

	REAL(dp) :: hstep,sstep,skstep,hkstep
	REAL(dp) :: hgrid(hn),hkgrid(hn),sgrid(sn),skgrid(sn)

	REAL(dp) :: agrid(an),astat(an),atrans(an,an)	!rouwenhorst
	REAL(dp) :: egrid(en,2)							!kennan

	REAL(dp) :: hfrac,sfrac,hkfrac,skfrac


! VFIT
! retirement
	REAL(dp),DIMENSION(sn,hn,cn) :: v10fun
! oldage, intervivos / after college
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn,cn) :: v9fun,s10fun,n9fun,s4fun
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn,cn) :: w8fun,s9fun,n8fun,n3fun
! midage, kid goes to school/college
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn) :: v6fun,v7fun,v8fun
! midage, kid in school
	REAL(dp),DIMENSION(sn,hn,hn,an,an,cn) :: s7fun,s8fun,n6fun,n7fun,l1fun,l2fun
	INTEGER,DIMENSION(sn,hn,hn,an,an,cn)  :: c3fun
! youngage, child-birth
	REAL(dp),DIMENSION(sn,hn,an,an,cn) :: v5fun,s6fun,n5fun,l0fun
! youngage, alone
	REAL(dp),DIMENSION(sn,hn,an,cn) :: v4fun,s5fun,n4fun


!------------------------------------------------------------------------------------------------------------------
! MPI
!------------------------------------------------------------------------------------------------------------------
	INTEGER :: myid,nprocs,ierr,ti
	INTEGER :: NIproc

	INTEGER :: nwp10,	ts10,	tf10
	INTEGER :: nwp89,	ts89,	tf89
	INTEGER :: nwp678,	ts678,	tf678
	INTEGER :: nwp5,	ts5,	tf5
	INTEGER :: nwp4,	ts4,	tf4

	INTEGER,PARAMETER :: sspace10	= sn*hn*cn
	INTEGER,PARAMETER :: sspace89	= sn*hn*cn *an *an*hn*cn
	INTEGER,PARAMETER :: sspace678	= sn*hn*cn *an *an*hn
	INTEGER,PARAMETER :: sspace5	= sn*hn*cn *an *an
	INTEGER,PARAMETER :: sspace4	= sn*hn*cn *an

! retirement, alone
	REAL(dp),ALLOCATABLE :: v10split(:)
! oldage, intervivos / after college
	REAL(dp),ALLOCATABLE :: w8v9split(:),s910split(:),n89split(:),n3s4split(:)
! midage, kid goes schoo1/college
	REAL(dp),ALLOCATABLE :: v678split(:),s78split(:),n67split(:),l12split(:)
	INTEGER,ALLOCATABLE :: c3split(:)
! youngage, kid born
	REAL(dp),ALLOCATABLE :: v5split(:),s6split(:),n5split(:),l0split(:)
! youngage, independent
	REAL(dp),ALLOCATABLE :: v4split(:),s5split(:),n4split(:)

!------------------------------------------------------------------------------------------------------------------
! MCMC: most of this crap i only need for experiments: FIX LATER
	INTEGER,ALLOCATABLE :: abi(:,:),epi(:,:,:)

! ALWAYS KEEP INTEGER-REAL ORDER FOR MPI DERIVED TYPE
! need up to grandchildren's exogenous shocks for short-run simulations
	TYPE :: AGENT
		INTEGER :: col_p,col_k,col_g						! college
		INTEGER :: ability_p,ability_k,ability_g			! ability shocks
		REAL(dp),DIMENSION(10) 		:: hh_p,hh_k,hh_g		! hcap states
		REAL(dp),DIMENSION(0:9) 	:: nn_p,nn_k,nn_g		! time choices
		REAL(dp),DIMENSION(3:10)	:: ee_p,ee_k,ee_g		! wealth states
		REAL(dp),DIMENSION(4:agen)	:: ss_p,ss_k,ss_g		! wealth states
!		REAL(dp) :: zy_p,zy_k,zy_g,zo_p,zo_k							! wlth states
!		REAL(dp) :: lftm_p,lftm_k,e58_p(4),e56_k(2)
!		REAL(dp) :: retw_p,w69_p(4),w57_k(3),coy_p(2),coy_k(2)			! ige states
!		REAL(dp) :: col_kip
	END TYPE AGENT
	TYPE(AGENT) :: hhid(NI)
	TYPE(AGENT),ALLOCATABLE :: hhidproc(:)
	INTEGER :: MPI_AGENT

	REAL(dp),DIMENSION(NI) :: lfkavge,lfkearn,lfkwlth,lfkutil
	REAL(dp),DIMENSION(NI) :: lfpavge,lfpearn,lfpwlth
	REAL(dp),DIMENSION(NI) :: lfgavge,lfgearn,lfgwlth

	REAL(dp),DIMENSION(NI) :: colfp,colfk,colfg
	REAL(dp),ALLOCATABLE :: colfpproc(:),colfkproc(:),colfgproc(:)

	REAL(dp),DIMENSION(4:12,NI) :: netincp,netinck,netincg
	REAL(dp),ALLOCATABLE :: netincpproc(:,:),netinckproc(:,:),netincgproc(:,:)

!	REAL(dp),ALLOCATABLE :: earnvec(:,:)



!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
! ALL TARGETS
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------

REAL(dp),DIMENSION(7),PARAMETER :: stdprof = (/ .4743609d0,.5059385d0,.5292786d0,.5638516d0,.586065d0,.6526362d0,.6857275d0 /)

REAL(dp),DIMENSION(7),PARAMETER :: avgprof = EXP( (/ -.5364066d0,-.2900843d0,-.1504221d0,-.0805986d0,-.0112411d0,-.0021376d0,-.1139525d0 /) + .0021376d0)

REAL(dp),DIMENSION(7),PARAMETER :: hscprof = EXP( (/ -.3044764d0,-.1199633d0,-.0441615d0,-.0250541d0,.0262443d0,-.0138407d0,-.1806566d0 /) + .0138407d0 )

REAL(dp),DIMENSION(7),PARAMETER :: colprof = EXP( (/ -.1662192d0,.1419125d0,.3454348d0,.4659746d0,.5533911d0,.6116829d0,.5548698d0 /) + .0138407d0 )	!	- .6116829d0 )

REAL(dp),DIMENSION(7),PARAMETER :: premprof = EXP( (/ -.1662192d0,.1419125d0,.3454348d0,.4659746d0,.5533911d0,.6116829d0,.5548698d0 /) - (/ -.3044764d0,-.1199633d0,-.0441615d0,-.0250541d0,.0262443d0,-.0138407d0,-.1806566d0 /) )

REAL(dp),DIMENSION(3),PARAMETER :: timeprof = (/.1573504d0,.1035712d0,.0825784d0/)

! ppige,educslope,nlsy educ pers,gs bequests
REAL(dp),DIMENSION(5),PARAMETER :: igcprof = (/0.341d0,0.675d0,0.634d0-colshare,igtfr,0.45d0/)
!0.356494491d0

REAL(dp),DIMENSION(6),PARAMETER :: testprof = (/ &
								4.982081d0,32.31565d0,49.53989d0, &
								2.94013d0,9.435945d0,11.66267d0 /)

REAL(dp),DIMENSION(7) :: mavgprof,mhscprof,mcolprof,mpremprof,mstdprof
REAL(dp) :: mtimeprof(3),migcprof(5)	!,mtestprof(6)


! bequests, ige, avg invest; hs-col earnings slope, enroll, gini
! ige: 0.344 is benchmark; 0.452 is 10-90 pctl
! educ slope: 0.675
	REAL(dp),PARAMETER :: eprem7 = EXP((LOG(premprof(1))+LOG(premprof(2))+LOG(premprof(3))+LOG(premprof(4))+LOG(premprof(5))+LOG(premprof(6))+LOG(premprof(7)))/7d0)
	REAL(dp),PARAMETER :: psidlogstd(2) = (/stdprof(1),stdprof(7)/)
	REAL(dp),PARAMETER :: slope(2) = (/ hscprof(6)/hscprof(1),colprof(6)/colprof(1) /)
	REAL(dp),DIMENSION(nom),PARAMETER :: targets = (/igtfr,0.675d0,0.1573504d0, &
											& eprem7,slope,psidlogstd, &
											& r0,colshare,1d0 /)*1d2
	! ig trans, ige, time, wgini=0.62 (hendricks)

	REAL(dp),PARAMETER :: mu_e(2) = LOG(avgprof(7))-LOG(avgprof(6)) + 0.0283492d0*(/-1D0,1D0/)/2D0

!	REAL(dp),PARAMETER :: mu_e(2) = (/ LOG(hscprof(7))-LOG(hscprof(6)), LOG(colprof(7))-LOG(colprof(6)) /)


SAVE
!-----------------------------------------------------------------------------------------------------------
END MODULE
