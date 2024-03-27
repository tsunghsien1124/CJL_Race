MODULE printresults
USE GLOBAL
IMPLICIT NONE

	REAL(dp),PARAMETER,PRIVATE :: ninfty = -1D4+1D-3	!9999999999999d0
	REAL(dp),PARAMETER,PRIVATE :: pinfty = 1D4-1D-3	!9999999999999d0

CONTAINS
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printguessreal
IMPLICIT NONE

	OPEN(11,FILE="output/guessreal.out",STATUS='replace')
	WRITE(11,'(20ES)'),params
	CLOSE(11)

RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printguess
IMPLICIT NONE
	CHARACTER(LEN=20) :: afmt,ffmt

	WRITE(afmt,'(A,I3,A)') '(',nop,'(A8))'
	WRITE(ffmt,'(A,I3,A)') '(',nop,'(F8.3))'

	OPEN(22,FILE="output/guess.out",STATUS='replace')

	WRITE(22,afmt),'thet','rho_a','zeta','omeg1','omeg2','bp','zcols','zcolm','mu_a','sig_a','bta','ww','tfp'
	WRITE(22,ffmt),params
	CALL printline(22)
	WRITE(22,ffmt),plb
	WRITE(22,ffmt),pub
	CALL printline(22)
	WRITE(22,afmt),'hsc','col','std','time','igc','eprem1','i.rate','enroll','avgearn'
	WRITE(22,ffmt),outweights
	CALL printline(22)
	WRITE(22,'(A15,F10.3)'),'DISTANCE:',distance
	CALL printline(22)
	WRITE(22,'(2A16)'),'hmax/smax*avge','hk/sk*max'
	WRITE(22,'(4F8.3)'),hfrac,sfrac,hkfrac,skfrac
	CALL printline(22)
	WRITE(22,'(20A8)'),'itmax','bign','switch'
	WRITE(22,'(20I8)'),itmax,bign,switch
	CALL printmoments(22)
!	CALL printprices(22)
	CALL printcheck(22)

	CLOSE(22)

RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printmoments(fileno)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: fileno

	CALL printline(fileno)
	WRITE(fileno,'(A20,10F10.3)'),'eprem/avge:',targets(4),targets(nom)
	WRITE(fileno,'(A20,10F10.3)'),'eprem/avge:',(/eprem1,avgearn/nmdlavgearn/)*1d2
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A20,10F10.3)'),'IR/enrl/ww/ups0:',targets(9:10),(/ww,ups/)*1d2
	WRITE(fileno,'(A20,10F10.3)'),'IR/enrl/ww/ups1:',(/r1,enroll,ww1,ups1/)*1d2
!	WRITE(fileno,'(A20,5F10.3)'),'/negretw:',1d2*(/beq,negretw/)

	CALL printline(fileno)
	WRITE(fileno,'(A20,10F10.3)'),'avg:',avgprof
	WRITE(fileno,'(A20,10F10.3)'),'mavg:',mavgprof
	WRITE(fileno,'(A20,10F10.3)'),'prem:',premprof
	WRITE(fileno,'(A20,10F10.3)'),'mprem:',mpremprof
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A20,10F10.3)'),'hsc:',hscprof
	WRITE(fileno,'(A20,10F10.3)'),'mhsc:',mhscprof
	WRITE(fileno,'(A20,10F10.3)'),'col:',colprof
	WRITE(fileno,'(A20,10F10.3)'),'mcol:',mcolprof
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A20,10F10.3)'),'logstd:',stdprof
	WRITE(fileno,'(A20,10F10.3)'),'mlogstd:',mstdprof
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A20,10F10.3)'),'time:',timeprof
	WRITE(fileno,'(A20,10F10.3)'),'mtime:',mtimeprof
!	WRITE(fileno,'(A20,10F10.3)'),'test:',testprof
!	WRITE(fileno,'(A20,10F10.3)'),'mtest:',mtestprof
	WRITE(fileno,'(A20,10F10.3)'),'igc4:',igcprof
	WRITE(fileno,'(A20,10F10.3)'),'migc4:',migcprof


RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printcheck(fileno)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: fileno
	CHARACTER(LEN=20) :: hfmt,sfmt

	WRITE(hfmt,'(A,I2,A)'),'(',hn+1,'(F6.3))'
	WRITE(sfmt,'(A,I2,A)'),'(',sn+1,'(F6.3))'

	CALL printline(fileno)
	WRITE(fileno,'(A15,F6.3)'),'hbins:',SUM(hbins)
	WRITE(fileno,*),'--------------'
	WRITE(fileno,hfmt),hbins
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A15,F6.3)'),'hkbins:',SUM(hkbins)
	WRITE(fileno,*),'--------------'
	WRITE(fileno,hfmt),hkbins
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A15,F6.3)'),'sbins:',SUM(sbins)
	WRITE(fileno,*),'--------------'
	WRITE(fileno,sfmt),sbins
	WRITE(fileno,*),'--------------'
	WRITE(fileno,'(A15,F6.3)'),'skbins:',SUM(skbins)
	WRITE(fileno,*),'--------------'
	WRITE(fileno,sfmt),skbins
RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printvalues
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: pfmt = '(100F10.4)'
	INTEGER :: ci,cii,ai,aii,hi,hii,si

! write out retirement valfunc
	OPEN(110,FILE="output/v10fun.out",STATUS='replace')
	DO hi=1,hn;DO ci=1,cn
		WRITE(110,pfmt),MIN(MAX(v10fun(:,hi,ci),ninfty),pinfty)
	ENDDO;ENDDO
	CLOSE(110)

! write out valfuncs with child
	OPEN(610,FILE="output/v6fun.out",STATUS='replace')
	OPEN(710,FILE="output/v7fun.out",STATUS='replace')
	OPEN(810,FILE="output/v8fun.out",STATUS='replace')
	OPEN(820,FILE="output/w8fun.out",STATUS='replace')
	OPEN(910,FILE="output/v9fun.out",STATUS='replace')
	DO hi=1,hn;DO hii=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(610,pfmt),MIN(MAX(v6fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(710,pfmt),MIN(MAX(v7fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(810,pfmt),MIN(MAX(v8fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		DO cii=1,cn
			WRITE(820,pfmt),MIN(MAX(w8fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
			WRITE(910,pfmt),MIN(MAX(v9fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
		ENDDO;
	ENDDO;ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(610);CLOSE(710);CLOSE(810);CLOSE(820);CLOSE(910)

! write out valfuncs with baby
	OPEN(510,FILE="output/v5fun.out",STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(510,pfmt),MIN(MAX(v5fun(:,hi,ai,aii,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(510)

! write out age 4 valfunc
	OPEN(410,FILE="output/v4fun.out" ,STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO ci=1,cn
		WRITE(410,pfmt),MIN(MAX(v4fun(:,hi,ai,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO
	CLOSE(410)

RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE printpolicy
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: pfmt = '(100F10.4)'
	CHARACTER(LEN=20),PARAMETER :: ifmt = '(100I10)'
	INTEGER,PARAMETER :: pinfti = INT(pinfty), ninfti = INT(ninfty)
	INTEGER :: ci,cii,ai,aii,hi,hii,si

! write out hc-time with child
	OPEN(610,FILE="output/n6fun.out",STATUS='replace')
	OPEN(710,FILE="output/n7fun.out",STATUS='replace')
	OPEN(810,FILE="output/n8fun.out",STATUS='replace')
	OPEN(910,FILE="output/n9fun.out",STATUS='replace')
	DO hi=1,hn;DO hii=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(610,pfmt),MIN(MAX(n6fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(710,pfmt),MIN(MAX(n7fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		DO cii=1,cn
			WRITE(810,pfmt),MIN(MAX(n8fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
			WRITE(910,pfmt),MIN(MAX(n9fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
		ENDDO;
	ENDDO;ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(610);CLOSE(710);CLOSE(810);CLOSE(910)

! write out hc-time with baby
	OPEN(510,FILE="output/n5fun.out",STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(510,pfmt),MIN(MAX(n5fun(:,hi,ai,aii,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(510)

! write out age 4 hc-time
	OPEN(410,FILE="output/n4fun.out" ,STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO ci=1,cn
		WRITE(410,pfmt),MIN(MAX(n4fun(:,hi,ai,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO
	CLOSE(410)
!---------------------------------
! write out savings with child
	OPEN(610,FILE="output/s7fun.out",STATUS='replace')
	OPEN(710,FILE="output/s8fun.out",STATUS='replace')
	OPEN(810,FILE="output/s9fun.out",STATUS='replace')
	OPEN(910,FILE="output/s10fun.out",STATUS='replace')
	DO hi=1,hn;DO hii=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(610,pfmt),MIN(MAX(s7fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(710,pfmt),MIN(MAX(s8fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		DO cii=1,cn
			WRITE(810,pfmt),MIN(MAX(s9fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
			WRITE(910,pfmt),MIN(MAX(s10fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
		ENDDO;
	ENDDO;ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(610);CLOSE(710);CLOSE(810);CLOSE(910)

! write out savings with baby
	OPEN(510,FILE="output/s6fun.out",STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(510,pfmt),MIN(MAX(n5fun(:,hi,ai,aii,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(510)

! write out age 4 savings
	OPEN(410,FILE="output/s5fun.out" ,STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO ci=1,cn
		WRITE(410,pfmt),MIN(MAX(n4fun(:,hi,ai,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO
	CLOSE(410)
!---------------------------------
! write out bequests, kid's hc-time and college choice
	OPEN(610,FILE="output/l1fun.out",STATUS='replace')
	OPEN(710,FILE="output/l2fun.out",STATUS='replace')
	OPEN(810,FILE="output/c3fun.out",STATUS='replace')
	OPEN(820,FILE="output/n3fun.out",STATUS='replace')
	OPEN(910,FILE="output/s4fun.out",STATUS='replace')
	DO hi=1,hn;DO hii=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(610,pfmt),MIN(MAX(l1fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(710,pfmt),MIN(MAX(l2fun(:,hi,hii,ai,aii,ci),ninfty),pinfty)
		WRITE(810,ifmt),MIN(MAX(c3fun(:,hi,hii,ai,aii,ci),ninfti),pinfti)
		DO cii=1,cn
			WRITE(820,pfmt),MIN(MAX(n3fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
			WRITE(910,pfmt),MIN(MAX(s4fun(:,hi,hii,ai,aii,ci,cii),ninfty),pinfty)
		ENDDO;
	ENDDO;ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(610);CLOSE(710);CLOSE(810);CLOSE(820);CLOSE(910)

! write out child-time
	OPEN(510,FILE="output/l0fun.out",STATUS='replace')
	DO hi=1,hn;DO ai=1,an;DO aii=1,an;DO ci=1,cn;
		WRITE(510,pfmt),MIN(MAX(l0fun(:,hi,ai,aii,ci),ninfty),pinfty)
	ENDDO;ENDDO;ENDDO;ENDDO
	CLOSE(510)

RETURN
END SUBROUTINE
!---------------------------------
SUBROUTINE printsimulation
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: pfmt = '(2I4,100F9.4)'
	INTEGER :: ind

	OPEN(1000,FILE="results/hhidgen1.out",STATUS='replace')
	OPEN(2000,FILE="results/hhidgen2.out",STATUS='replace')
	OPEN(3000,FILE="results/hhidgen3.out",STATUS='replace')
	DO ind=1,3
		WRITE(ind*1000,'(2A4,100A9)'),'abi','col',	&
						'hh_1','hh_2','hh_3','hh_4','hh_5','hh_6','hh_7','hh_8','hh_9','hh_10',	&
						'ee_3','ee_4','ee_5','ee_6','ee_7','ee_8','ee_9' ,'ee_10',	&
						'ss_4','ss_5','ss_6','ss_7','ss_8','ss_9','ss_10','ss_11','ss_12',	&
						'll_0','ll_1','ll_2'
	ENDDO
	DO ind=1,NI
		WRITE(1000,pfmt),hhid(ind)%ability_p,hhid(ind)%col_p,hhid(ind)%hh_p,hhid(ind)%ee_p,hhid(ind)%ss_p,hhid(ind)%nn_p(0:2)
		WRITE(2000,pfmt),hhid(ind)%ability_k,hhid(ind)%col_k,hhid(ind)%hh_k,hhid(ind)%ee_k,hhid(ind)%ss_k,hhid(ind)%nn_k(0:2)
		WRITE(3000,pfmt),hhid(ind)%ability_g,hhid(ind)%col_g,hhid(ind)%hh_g,hhid(ind)%ee_g,hhid(ind)%ss_g,hhid(ind)%nn_g(0:2)
	ENDDO
	DO ind=1,3;CLOSE(ind*1000);ENDDO

RETURN
END SUBROUTINE
!---------------------------------
SUBROUTINE printlifetime
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: pfmt = '(20F9.4)'
	INTEGER :: ind

	OPEN(1000,FILE="results/lifetime.out",STATUS='replace')
	WRITE(1000,'(10A9)'),'lfkavge','lfkearn','lfkwlth','lfkutil','lfpavge','lfpearn','lfpwlth','lfgavge','lfgearn','lfgwlth'
	DO ind=1,NI
		WRITE(1000,pfmt),lfkavge(ind),lfkearn(ind),lfkwlth(ind),lfkutil(ind),lfpavge(ind),lfpearn(ind),lfpwlth(ind),lfgavge(ind),lfgearn(ind),lfgwlth(ind)
	ENDDO
	CLOSE(1000)

RETURN
END SUBROUTINE
!---------------------------------
SUBROUTINE printnetinc
IMPLICIT NONE
	CHARACTER(LEN=20),PARAMETER :: pfmt = '(20F9.4)'
	INTEGER :: ind

	OPEN(1000,FILE="results/netincp.out",STATUS='replace')
	WRITE(1000,'(12A9)'),'inc4','inc5','inc6','inc7','inc8','inc9','inc10','inc11','inc12'
	DO ind=1,NI
		WRITE(1000,pfmt),netincp(:,ind)
	ENDDO
	CLOSE(1000)

	OPEN(1000,FILE="results/netinck.out",STATUS='replace')
	WRITE(1000,'(12A9)'),'inc4','inc5','inc6','inc7','inc8','inc9','inc10','inc11','inc12'
	DO ind=1,NI
		WRITE(1000,pfmt),netinck(:,ind)
	ENDDO
	CLOSE(1000)

	OPEN(1000,FILE="results/netincg.out",STATUS='replace')
	WRITE(1000,'(12A9)'),'inc4','inc5','inc6','inc7','inc8','inc9','inc10','inc11','inc12'
	DO ind=1,NI
		WRITE(1000,pfmt),netincg(:,ind)
	ENDDO
	CLOSE(1000)

RETURN
END SUBROUTINE
!---------------------------------
SUBROUTINE printline(fileno)
IMPLICIT NONE
	INTEGER,INTENT(IN) :: fileno
	WRITE(fileno,*),'--------------------------------------------------'
END SUBROUTINE
!---------------------------------
END MODULE
