! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see 
! the GNU General Public License.
      MODULE  Stand    
      USE DETDATA ; USE DISLIN
      CONTAINS   
! Common procedures for TOMB and CRCS
!--------------------------------------------------------
      SUBROUTINE cov(rx,ry,ilx,rcov) !  Correlation procedure     
      IMPLICIT NONE
      INTEGER,INTENT(IN)                :: ilx       ! Data count
      REAL(8),DIMENSION(ilx),INTENT(IN) :: rx,ry     ! Data fields
      REAL(8),INTENT(OUT)               :: rcov      ! Correlation coefficient
      REAL(8)                           :: mx,my,lx,sq
      lx=DBLE(ilx) ; mx=SUM(rx)/lx ; my=SUM(ry)/lx 
      sq=SQRT(SUM((rx-mx)**2)*SUM((ry-my)**2))
      IF (sq.GT.0.D0) THEN     ! Are series identical
        rcov=SUM((rx-mx)*(ry-my))/sq
      ELSE
        rcov=0.D0
      ENDIF
      RETURN
      END SUBROUTINE cov
!----------------------------------------------------------------------
      SUBROUTINE covmiss(rx,ry,rok,ilx,rcov,covok)  ! Correlation procedure     
      IMPLICIT NONE                                 ! Exclude absent values 
      INTEGER,INTENT(IN)                :: ilx      ! Series length
      REAL(8),INTENT(IN),DIMENSION(ilx) :: rx,ry    ! Input series
      LOGICAL,INTENT(IN),DIMENSION(ilx) :: rok      ! Rings OK or missing
      LOGICAL,INTENT(OUT)               :: covok    ! Procedure worked
      REAL(8),INTENT(OUT)               :: rcov     ! Correlation value 
      REAL(8),DIMENSION(ilx)            :: wk1,wk2  ! Working storage areas
      REAL(8)                           :: mx,my,lx,sq
      INTEGER                           :: i,j
      covok=.TRUE.
      IF (.NOT.ALL(rok(1:ilx))) THEN
        j=0
        DO i=1,ilx       ! Close gaps
          IF (rok(i)) THEN
            j=j+1 ; wk1(j)=rx(i) ; wk2(j)=ry(i)
          ENDIF
        ENDDO
        IF (j.GT.4) THEN
          lx=DBLE(j) ; mx=SUM(wk1(1:j))/lx ; my=SUM(wk2(1:j))/lx 
          sq=SQRT(SUM((wk1(1:j)-mx)**2)*SUM((wk2(1:j)-my)**2))
          IF (sq.GT.0.D0) THEN   ! Are series identical
            rcov=SUM((wk1(1:j)-mx)*(wk2(1:j)-my))/sq
          ELSE
            rcov=0.D0
          ENDIF
        ELSE
          rcov=2.D0 ; covok=.FALSE.
        ENDIF
      ELSE
        lx=DBLE(ilx) ; mx=SUM(rx)/lx ; my=SUM(ry)/lx 
        sq=SQRT(SUM((rx-mx)**2)*SUM((ry-my)**2))
        IF (sq.GT.0.D0) THEN   ! Are series identical
          rcov=SUM((rx-mx)*(ry-my))/sq
        ELSE
          rcov=0.D0
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE covmiss
 !------------------------------------------------------------------
      SUBROUTINE QCKSRT(N,ARR)   ! From ARSTAN
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: N
      REAL(8),DIMENSION(N),INTENT(INOUT) :: ARR
      INTEGER,PARAMETER  :: M=7,NSTACK=10000
      REAL(8),PARAMETER  :: FM=7875.D0,FA=211.D0
      REAL(8),PARAMETER  :: FC=1663.D0,FMI=1.2698413D-4
     INTEGER,DIMENSION(NSTACK) :: ISTACK
      INTEGER :: JSTACK,L,IR,J,I,IQ
      REAL(8) :: FX,A
      JSTACK=0 ; L=1 ; IR=N ; FX=0.D0
      D4: DO 
        IF(IR-L.LT.M)THEN
          DO J=L+1,IR
            A=ARR(J)
            I2: DO I=J-1,1,-1
              IF(ARR(I).LE.A) EXIT I2
              ARR(I+1)=ARR(I)
            ENDDO I2 
            IF (I.EQ.-1) I=0
            ARR(I+1)=A  
          ENDDO 
          IF(JSTACK.EQ.0) EXIT D4
          IR=ISTACK(JSTACK) ; L=ISTACK(JSTACK-1) ; JSTACK=JSTACK-2
        ELSE
          I=L ; J=IR ; FX=MOD(FX*FA+FC,FM)
          IQ=L+INT(DBLE(IR-L+1)*(FX*FMI))
          A=ARR(IQ)  ; ARR(IQ)=ARR(L)
          D2: DO  
            JDO: DO J=J,1,-1
              IF(A.GE.ARR(J)) EXIT JDO
            ENDDO JDO
            IF(J.LE.I)THEN
              ARR(I)=A ; EXIT D2
            ENDIF
            ARR(I)=ARR(J) ; I=I+1
            I1: DO I=I,N-1 ; IF(A.LE.ARR(I)) EXIT I1 ; ENDDO I1
            IF(J.LE.I)THEN
              ARR(J)=A ; I=J ; EXIT D2
            ENDIF
            ARR(J)=ARR(I) ; J=J-1 
          ENDDO D2
          JSTACK=JSTACK+2 
          IF(JSTACK.GT.NSTACK) &
            CALL out_err("QCKSRT - NSTACK must be larger.")
          IF(IR-I.GE.I-L)THEN
            ISTACK(JSTACK)=IR ; ISTACK(JSTACK-1)=I+1 ; IR=I-1
          ELSE
            ISTACK(JSTACK)=I-1 ; ISTACK(JSTACK-1)=L ; L=I+1
          ENDIF
        ENDIF
      ENDDO D4
      RETURN
      END SUBROUTINE QCKSRT
!-----------------------------------------------------------------
      SUBROUTINE MFS(X,N,XMED)  ! Computes the median of n observations
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: N     ! Number of observations 
      REAL(8),DIMENSION(N),INTENT(INOUT) :: X     ! Input/output series
      REAL(8),INTENT(OUT)                :: XMED  ! median of X
      INTEGER                            :: TEMP1,J,K
      CALL QCKSRT(N,X)       ! RANK THE OBSERVATIONS
      K=N ; J=(K/2)+1 ; TEMP1=N+1-J
      XMED=(X(J)+X(TEMP1))/2.D0                
      RETURN
      END SUBROUTINE MFS
!--------------------------------------------------------------------
!SBR BIWGT: BIWEIGHT ROBUST MEAN ESTIMATION USING THE NEWTON-RAPHSON ALGORITHM
!  PROGRAMMED BY: EDWARD R COOK * TREE-RING LABORATORY
!     DAT  - INPUT  - DATA VECTOR OF LENGTH N
!     WRK - INPUT  - WORK VECTOR OF LENGTH N
!      C  - INPUT  - TUNING CONSTANT USED FOR COMPUTING BIWEIGHT MEAN
!		     WHERE:  C USUALLY FALLS BETWEEN 6 AND 9, INCLUSIVE
!      N  - INPUT  - NUMBER OF OBSERVATIONS IN DATA VECTOR X
!    BIMN - OUTPUT - BIWEIGHT ROBUST MEAN
!    BISG - OUTPUT - BIWEIGHT STANDARD DEVIATION
!    BISE - OUTPUT - BIWEIGHT MEAN STANDARD ERROR
!     EFF - OUTPUT - EFFICIENCY OF BIWEIGHT MEAN RELATIVE TO GAUSSIAN MEAN
!    ITER - OUTPUT - NUMBER OF NEWTON-RAPHSON ITERATIONS
!    	SUBROUTINES CALLED - SORT
      SUBROUTINE BIWGT(DAT,C,N,BIMN,BISG)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: N
      REAL(8),DIMENSION(N),INTENT(INOUT) :: DAT
      REAL(8),INTENT(IN)                 :: C
      REAL(8),INTENT( OUT)               :: BIMN,BISG
      REAL(8),DIMENSION(N)               :: WRK
      INTEGER :: ITER,I,M,MM
      REAL(8) :: MED,MAD,RN,GMEAN,GVAR,SMED,SAV
      REAL(8) :: XSUM,YSUM,U,USQ,DUSQ,FUSQ,CMAD,ABSU,DEP
      REAL(8),PARAMETER :: TOL=0.001D0      !TOLERANCE CRITERION
      IF (N.EQ.1) THEN                      !IF N=1, JUST RETURN
        BIMN=DAT(1) ; BISG=0.D0 ; RETURN
      ENDIF
      RN=DBLE(N) ; GMEAN=SUM(DAT(1:N))/RN   !COMPUTE MEAN AND VARIANCE
      GVAR=SUM((DAT(1:N)-GMEAN)**2)/(RN-1.D0)  !GAUSSIAN VARIANCE (GVAR)
      CALL MFS(DAT,N,SMED)
      M=(N+1)/2 ; MM=N-M+1 ; MED=SMED ; BIMN=SMED !MEDIAN(MED) OF DATA VECTOR X
      SAV=SMED ; WRK(1:N)=ABS(DAT(1:N)-MED) !ABSOLUTE DEPARTURES FROM MED
      CALL QCKSRT(N,WRK)
      MAD=(WRK(M)+WRK(MM))/2.D0             !MEDIAN ABSOLUTE DEVIATION (MAD)
      IF (N.LE.5) THEN  !IF N .LE. 5, USE MEDIAN AND
        MAD=3.D0*MAD/2.D0 ; BISG=MAD ; RETURN         
      ENDIF   
      IF (MAD.LE.0.000001D0) THEN
        BIMN=MED ; BISG=MAD ; RETURN
      ENDIF
      bimn=smed
! COMPUTE BISQUARE MEAN USING THE MEDIAN AS THE INITIAL ESTIMATE
      CMAD=C*MAD                  ! WEIGHT MAD BY TUNING CONSTANT
      ITERDO: DO ITER=1,11        ! AND BEGIN ITERATIONS
        IF(ITER.GT.10)THEN
          CALL out_err("BIWGT - NO CONVERGENCE") ; STOP
        ENDIF
        XSUM=0.D0 ; YSUM=0.D0
        DO I=1,N
          U=(DAT(I)-BIMN)/CMAD ; ABSU=DABS(U)   !ROBUST DEPARTURES (U)
          IF(ABSU.LT.1.D0)THEN                       !IF ABS(U>=1) IGNORE, ELSE
            USQ=U*U ; DUSQ=1.D0-USQ ; FUSQ=1.D0-5.D0*USQ !COMPUTE AND SUM TERMS FOR
            XSUM=XSUM+U*DUSQ**2 ; YSUM=YSUM+DUSQ*FUSQ !THE BIWEIGHT INFLUENCE CURVE 
          ENDIF                                       !ITS FIRST-DERIVATIVE
        END DO          
        BIMN=BIMN+CMAD*(XSUM/YSUM)       !NEWTON-RAPHSON ESTIMATE
        IF(ABS(BIMN-SAV).GT.TOL)THEN     !TEST IT AGAINST TOL
          SAV=BIMN                       !IF NOT CLOSE ENOUGH, SAVE IT
        ELSE                             !AND ITERATE AGAIN
          EXIT ITERDO
        ENDIF
      ENDDO ITERDO                                    
      CMAD=DBLE(C*MAD) ; XSUM=0.D0 ;YSUM=0.D0  !COMPUTE THE BIWEIGHT ESTIMATE OF
      DO I=1,N                                   !SCALE AND STANDARD ERROR OF THE MEAN
        U=DBLE((DAT(I)-BIMN))/CMAD ; USQ=U*U     !ROBUST DEPARTURES (U)
        IF(DABS(U).LT.1.D0)THEN     !IF ABS(U>=1) IGNORE, ELSE
          DUSQ=1.D0-USQ ; FUSQ=1.D0-5.D0*USQ  !COMPUTE AND SUM TERMS
          DEP=DAT(I)-BIMN
          XSUM=XSUM+DEP**2*DUSQ**4  !NUMERATOR TERM OF SCALE ESTIMATE
          YSUM=YSUM+DUSQ*FUSQ       !DENOMINATOR TERM IN SCALE ESTIMATE
        ENDIF
      ENDDO
      BISG=DSQRT((RN*XSUM)/(YSUM*(YSUM-1.D0))) !BIWEIGHT STANDARD DEVIATION
      RETURN
      END SUBROUTINE BIWGT
!---------------------------------------------------------------------------------
! SBR CURVE: COMPUTE EXPONENTIAL CURVE COEFFICIENTS PROGRAMMED BY FRITTS, HUZAR,
! BOTTORF AND RAY. REVISED OCT 1968 (V BOCKMAN,1983) BY E R COOK, R L HOLMES 
! FEB 1984, FEB 1985. EXPONENTIAL CURVE EQUATION COMPUTED:
! Y = AH * EXP(-BH * T) + EH  (AH=EQ(1),BH=EQ(2),EH=EQ(5))
! EQ(6)= ITERATIONS REQUIRED TO CONVERGE, EQ(7)= SUM OF ARRAY Z
      SUBROUTINE curvet(N,DX1,DX2,EQ)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: N    ! Number of data values in z
      REAL(8),DIMENSION(1:N),INTENT(IN)  :: DX1  ! Data to which to fit curve
      REAL(8),DIMENSION(1:N),INTENT(OUT) :: DX2  ! Working array of length n
      REAL(8),DIMENSION(7),INTENT(OUT)   :: EQ   ! EXP. curve coefficients comPUTED
      CHARACTER(3)         :: MESS
      INTEGER              :: I,K,IT,MO,LO,MU      
      REAL(8)              :: BB,D,E,G,EH,AH,YF,YL,RN,DLT
      REAL(8),DIMENSION(3) :: BH,RSS
      LOGICAL              :: FAILED
      REAL(8),PARAMETER    :: DINC=0.0000005D0
      FAILED=.FALSE. ; DX2=0.D0
      IF(N .LT. 20)THEN
        WRITE(MESS,'(I3)') N
        CALL out_err("CURVE: Only "//MESS//" values, no exp fitted")
        FAILED = .TRUE.
      ELSE
        RN=DBLE(N) ; D=SUM(DBLE(DX1(1:N))) 
        YF=SUM(DBLE(DX1(1:10)))/10.D0
        YL=SUM(DBLE(DX1(N-9:N)))/10.D0 
        BH(2)=MAX((DLOG(YF)-DLOG(YL))/(RN-10.D0),1.D-8)
        DLT=0.0005D0 ; MO=0 ; MU=0
        ITDO: DO IT=1,34                   !  Iteration loop
          BH(1)=BH(2)-DLT ; BH(3)=BH(2)+DLT
          KDO: DO K=1,3
            IF(K.NE.2.OR.MU.LE.0) THEN 
              IF(ABS(BH(K))*RN.GT.112.8D0) THEN  ! presumes last number is greatest
                FAILED=.TRUE. ; EXIT ITDO
              ENDIF
              DX2(1:N)=DEXP(-BH(K)*(/(DBLE(I),I=1,N)/)) 
              BB=SUM(DBLE(DX2(1:N)))
              E=SUM(DBLE(DX2(1:N))**2)
              G=SUM(DBLE(DX1(1:N))*DBLE(DX2(1:N)))
              IF(BB.GT.1.D15.OR.E.GT.1.D15.OR.BB.LT.1.D-15.OR.E.LT.1.D-15) THEN
                FAILED=.TRUE. ; EXIT ITDO
              ENDIF
              EH=(D*E-BB*G)/(RN*E-BB**2) ; AH=(D-RN*EH)/BB
              RSS(K)=SUM((DBLE(DX1(1:N))-EH-AH*DBLE(DX2(1:N)))**2)
            ENDIF
          END DO KDO  
          LO=0 ; MU=1
          DO K=1,3,2
            IF(RSS(K).LT. RSS(2))THEN
              LO=-1 ; BH(2)=BH(K) ;  RSS(2)=RSS(K)
            ENDIF
          END DO  !30
          IF(LO.GE.0) MO=1
          IF(MO.LE.0)THEN
            DLT=DLT*2.D0
          ELSE
            DLT=DLT*0.5D0 ; IF(DLT.LE.DINC) EXIT ITDO
          ENDIF
        END DO ITDO !13
        IF (IT.GT.34) THEN
          CALL out_err("CURVE: Over 34 iterations, ")
          FAILED=.TRUE.
        ENDIF
      ENDIF
      IF (FAILED) THEN  ! Flag values if coefficients not computed
        DX2(1)=-999.0D0 ; EQ(2)=-999.D0
      ELSE
        IF(EH.GT.9999.D0.OR.AH.LE.1.D-15.OR.BH(2).LT.1.D-15) BH(2)=0.D0 
        EQ(1)=AH ; EQ(2)=BH(2) ; EQ(5)=EH
        EQ(6)=DBLE(IT) ; EQ(7)=DBLE(D)
        DX2(1:N)=AH*DEXP(BH(2)*(/(DBLE(-I),I=1,N)/))+EH ! Curve values
      END IF
      RETURN
      END SUBROUTINE curvet
!-----------------------------------------------------------------
      SUBROUTINE TREND(n,dx1,dx2,sl,yi,xi)  ! Best fit straight line
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n    ! Number of values in series
      REAL(8),DIMENSION(n),INTENT(IN)  :: dx1  ! Time series array
      REAL(8),DIMENSION(n),INTENT(OUT) :: dx2  ! Calculated trend line values 
      REAL(8),INTENT(OUT)              :: sl   ! Slope of trend line
      REAL(8),INTENT(OUT)              :: yi   ! Y-axis intercept
      REAL(8),INTENT(OUT)              :: xi   ! X-axis intercept
      REAL(8),DIMENSION(n)             :: wk
      REAL(8) :: xy,x,y,xx,an
      INTEGER :: i=1
      an=DBLE(n) ; wk(1:n)=(/(DBLE(i),i=1,n)/)
      xy=SUM(dx1(1:n)*wk(1:n)) ; x=SUM(wk(1:n)) ; y=SUM(dx1(1:n))
      xx=SUM(wk(1:n)**2) ; sl=(an*xy-x*y)/(an*xx-x**2)
      yi=y/an-sl*x/an ; dx2(1:n)=wk(1:n)*sl+yi
      IF (sl.NE.0.D0) xi=(-yi)/sl 
      RETURN
      END SUBROUTINE TREND
!----------------------------------------------------------
      SUBROUTINE TRENDmiss(n,dx1,dx2,okt,sl,yi,xi)  ! Best fit straight line
      IMPLICIT NONE   ! For series with missing values
      INTEGER,INTENT(IN)               :: n    ! Number of values in series
      LOGICAL,DIMENSION(n),INTENT(IN)  :: okt  ! Present or not
      REAL(8),DIMENSION(n),INTENT(IN)  :: dx1  ! Time series array
      REAL(8),DIMENSION(n),INTENT(OUT) :: dx2  ! Calculated trend line values 
      REAL(8),INTENT(OUT)              :: sl   ! Slope of trend line
      REAL(8),INTENT(OUT)              :: yi   ! Y-axis intercept
      REAL(8),INTENT(OUT)              :: xi   ! X-axis intercept
      REAL(8),DIMENSION(n)             :: wk
      REAL(8)    :: xy,x,y,xx,an
      INTEGER :: i=1
      an=DBLE(COUNT(okt(1:n))) ; wk(1:n)=(/(DBLE(i),i=1,n)/)
      xy=SUM(dx1(1:n)*wk(1:n),MASK=okt(1:n))
      x=SUM(wk(1:n),MASK=okt(1:n))
      y=SUM(dx1(1:n),MASK=okt(1:n))
      xx=SUM(wk(1:n)**2,MASK=okt(1:n))
      sl=(an*xy-x*y)/(an*xx-x**2)
      yi=y/an-sl*x/an ; dx2(1:n)=wk(1:n)*sl+yi
      IF (sl.NE.0.D0) xi=(-yi)/sl 
      RETURN
      END SUBROUTINE TRENDmiss
!----------------------------------------------------------
      SUBROUTINE GEXP(DX1,N,DX2,AH,BH,OK)
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: N       ! Number of years
      LOGICAL,DIMENSION(N),INTENT(IN)  :: OK      ! Valid ring
      REAL(8),DIMENSION(N),INTENT(IN)  :: DX1     ! Tree-ring series
      REAL(8),DIMENSION(N),INTENT(OUT) :: DX2     ! General exponential curve estimates
      REAL(8),INTENT(OUT)              :: AH,BH   ! Coefficients
      REAL(8),DIMENSION(N)             :: WKA,WKB ! Work vector
      INTEGER :: N1,J,I
      REAL(8) :: R,XA,AY,DN1,XT,YT,SXX,SYY,SXY,SGY,SGX,YXR
      AH=0.0D0 ; BH=0.0D0 ; N1=0
      DO I=1,N            ! Linearize the problem
        IF(DX1(I).GT.0.00001D0.AND.OK(I))THEN
          N1=N1+1 ; WKA(N1)=DBLE(I)
          WKB(N1)=DLOG(DX1(I))-DLOG(WKA(N1))
        ENDIF
      ENDDO               ! Now compute simple linear regression
      DN1=DBLE(N1) ; XA=SUM(WKA(1:N1))/DN1
      AY=SUM(WKB(1:N1))/DN1 ; SXX=0.D0 ; SYY=0.D0 ; SXY=0.D0
      DO J=1,N1
        XT=WKA(J)-XA ; YT=WKB(J)-AY
        SXX=SXX+XT**2 ; SYY=SYY+YT**2 ; SXY=SXY+XT*YT
      ENDDO
      SGY=DSQRT(SYY/(DN1-1.D0)) ; SGX=DSQRT(SXX/(DN1-1.D0))
      YXR=SGY/SGX ; R=SXY/DSQRT(SXX*SYY) ; BH=YXR*R ; AH=AY-BH*XA
      AH=DEXP(AH)   ! Unlinearize the model and calculate growth curve values
      DX2(1:N)=AH*(/(DBLE(I),I=1,N)/)*DEXP(BH*(/(DBLE(I),I=1,N)/))
      BH=-BH
      RETURN
      END SUBROUTINE GEXP
!-----------------------------------------------------------------------
      SUBROUTINE HUGHDI(DX1,NTOT,DX2,COEF,OK,POFF)
!     Y(I) = A * T(I)**CM * EXP(CK*T(I))   for I=1....NTOT
!     Hauptprogramm NCOEF muss 3 sein September 1989 /NOS/VE/ Braeker/Nogler
!     PC-Version 2.1/Juni 1992/P.Nogler Modified by E.R. Cook Oct 23, 1996
      IMPLICIT NONE
      INTEGER,INTENT(IN)                  :: NTOT  ! Number of years
      INTEGER,INTENT(IN)                  :: POFF  ! Pith offset years
      LOGICAL,DIMENSION(NTOT),INTENT(IN)  :: OK    ! Valid ring
      REAL(8),DIMENSION(NTOT),INTENT(IN)  :: DX1   ! Tree-ring series
      REAL(8),DIMENSION(NTOT),INTENT(OUT) :: DX2   ! Hugershoff curve estimates
      REAL(8),DIMENSION(7),INTENT(OUT)    :: COEF  ! 1:3 are hugershoff coefficients
      REAL(8),DIMENSION(NTOT)             :: WK    ! Time series as year numbers
      REAL(8)  :: A,CK,CM,SE,SF,CC,SUMT,SUMX2,SUMT2,SUMXT
      REAL(8)  :: SUMYT,SUMXY,SUMY,SUMX,XK,YK,XL,YL,RN1,RN2,RN3,RN
      REAL(8)  :: RZA,RZ2,RZM,RZ1,RZK,DXY
      INTEGER  :: NCOEF,NRED,K=1
!     wk(1:ntot)=(/(DBLE(k),k=1,ntot)/)
      wk(1:ntot)=(/(DBLE(k),k=POFF,POFF+ntot-1)/)
      CC=-999.0D0 ; NRED=0 ; SUMT=0.D0 ; SUMX=0.D0 ; SUMX2=0.D0
      SUMT2=0.D0 ; SUMXT=0.D0 ; SUMY=0.D0 ; SUMXY=0.D0
      SUMYT=0.D0
      DO K=1,NTOT
        XK=WK(K) ; YK=DX1(K)
        IF(OK(K).AND.(YK.GT.1.D-15).AND.(XK.GT.1.D-15)) THEN
          NRED=NRED+1 ; XL=DLOG(XK) ; YL=DLOG(YK) ; SUMT=SUMT+XK
          SUMT2=SUMT2+XK**2 ; SUMX=SUMX+XL ; SUMX2=SUMX2+XL**2
          SUMXT=SUMXT+XL*XK ; SUMY=SUMY+YL ; SUMXY=SUMXY+XL*YL
          SUMYT=SUMYT+YL*XK
        ENDIF
      END DO  
      RN1=(SUMX2*SUMT2-SUMXT**2) ; RN2=(SUMT*SUMXT-SUMX*SUMT2) 
      RN3=(SUMX*SUMXT-SUMT*SUMX2) ; RN=DBLE(NRED)*RN1+SUMX*RN2+SUMT*RN3
      RZA=SUMY*RN1+SUMXY*RN2+SUMYT*RN3 ; A=RZA/RN
      RZ2=(SUMX*SUMT-DBLE(NRED)*SUMXT)
      RZM=SUMY*RN2+SUMXY*(DBLE(NRED)*SUMT2-SUMT**2)+SUMYT*RZ2
      CM=RZM/RN ; RZ1=(SUMT*SUMX2-SUMX*SUMXT)
      RZ2=(DBLE(NRED)*SUMXT-SUMX*SUMT)
      RZK=SUMY*RZ1+SUMXY*RZ2+SUMYT*(SUMX**2-DBLE(NRED)*SUMX2)
      CK=RZK/RN ; NCOEF=3 ; COEF(1:3)=(/DEXP(A),CM,-CK/)
!     IF(COEF(2).LE.1.D-15)THEN
!       COEF(2)=-999.D0 ; DX2(1)=-999.D0 ; RETURN
!     ENDIF
      SE=0.D0 ; SF=0.D0
      DO K=1,NTOT 
        XK=DBLE(WK(K)) ;  YK=DX1(K)
        DXY=COEF(1)*XK**COEF(2)*DEXP(COEF(3)*XK)
        DX2(K)=DXY
        IF(OK(K))THEN 
          SE=SE+(YK-DXY) ; SF=SF+DXY 
        ENDIF
      ENDDO
      CC=SE/SF+1.D0
      IF(CC.NE.0.D0) THEN 
        COEF(1)=COEF(1)*CC ; DX2(1:NTOT)=DX2(1:NTOT)*CC 
      ENDIF
      RETURN 
      END SUBROUTINE HUGHDI
!-------------------------------------------------------------------- 
      SUBROUTINE splinec(n,rw,ssy,rwp)  ! Spline calculation
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n    ! Length of series 
      INTEGER,DIMENSION(n),INTENT(IN)  :: ssy  ! Spline stiffness each year
      REAL(8),DIMENSION(n),INTENT(IN)  :: rw   ! Series to fit 
      REAL(8),DIMENSION(n),INTENT(OUT) :: rwp  ! Spline curve
      INTEGER                          :: i
      REAL(8)                          :: arg 
      REAL(8),DIMENSION(-1:n)          :: a1,a2,a3,a4,p  ! Working arrays
      REAL(8),PARAMETER                :: pi = 3.1415926535897935D0 
      DO i=1,n-2                  ! Stiffness for each year
        arg=COS(2.D0*pi/DBLE(ssy(i)))
        p(i)=(2.D0*(arg-1.D0)**2)/(arg+2.D0)
      ENDDO
      a1=0.D0 ; a2=0.D0 ; a3=0.D0 ; a4=0.D0
      a1(3:n-2)=1.D0 ; a2(2:n-2)=-4.D0+p(2:n-2)
      a3(1:n-2)=6.D0+p(1:n-2)*4.D0
      a4(1:n-2)=rw(1:n-2)+rw(3:n)-2.D0*rw(2:n-1)
      DO i=1,n-2                     ! Solve matrix
        a1(i)=a1(i)*a3(i-2)
        a2(i)=(a2(i)-a1(i)*a2(i-1))*a3(i-1)
        a3(i)=1.D0/SQRT(a3(i)-a1(i)**2-a2(i)**2)
        a4(i)=(a4(i)-a4(i-2)*a1(i)-a4(i-1)*a2(i))*a3(i)
      ENDDO
      DO i=n-2,1,-1              
        a4(i)=(a4(i)-a4(i+1)*a2(i+1)-a4(i+2)*a1(i+2))*a3(i)
      END DO
      rwp(1:n)=rw(1:n)-a4(-1:n-2)+2.D0*a4(0:n-1)-a4(1:n)
      RETURN
      END SUBROUTINE splinec
!-----------------------------------------------------------------
      SUBROUTINE spline_miss(n,rw,ss,rwp,ok) ! Calculates fixed spline
      IMPLICIT NONE    ! Where data can be missing and closes gaps
      INTEGER,INTENT(IN)               :: n    ! Length of series 
      INTEGER,INTENT(IN)               :: ss   ! Spline stiffness
      REAL(8),DIMENSION(n),INTENT(IN)  :: rw   ! Series to fit 
      LOGICAL,DIMENSION(n),INTENT(IN)  :: ok   ! Series to fit 
      REAL(8),DIMENSION(n),INTENT(OUT) :: rwp  ! Spline curve
      REAL(8),DIMENSION(n)             :: wk   ! Working area
      INTEGER,DIMENSION(n)             :: ssy  ! Stiffness each year
      INTEGER :: i,j
      j=0
      DO i=1,n       ! Pack into gaps
        IF (ok(i)) THEN 
          j=j+1 ; wk(j)=rw(i)
        ENDIF
      ENDDO          ! Fit spline
      ssy(1:j-2)=ss ; CALL splinec(j,wk(1:j),ssy(1:j),rwp(1:j))      
      DO i=n,2,-1    ! Unpack across gaps
        IF (ok(i)) THEN
          rwp(i)=rwp(j) ; j=j-1
        ELSE
          rwp(i)=rwp(MIN(n,i+1))  ! Copy next value into gap
        ENDIF
      ENDDO          
      RETURN
      END SUBROUTINE spline_miss
!-----------------------------------------------------------------
      SUBROUTINE splinet(n,rw,ss,rwp) ! Calculates fixed spline
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n    ! Length of series 
      INTEGER,INTENT(IN)               :: ss   ! Spline stiffness
      REAL(8),DIMENSION(n),INTENT(IN)  :: rw   ! Series to fit 
      REAL(8),DIMENSION(n),INTENT(OUT) :: rwp  ! Spline curve
      INTEGER,DIMENSION(n)             :: ssy  ! Stiffness each year
      ssy(1:n-2)=ss ; CALL splinec(n,rw,ssy,rwp)      
      RETURN
      END SUBROUTINE splinet
!-----------------------------------------------------------------
      SUBROUTINE spline3(n,rws,cnt,ss,rwp,rise)  ! Age related spline, no rise
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n    ! Length of series 
      INTEGER,INTENT(IN)               :: ss   ! Age offset stiffness
      REAL(8),DIMENSION(n),INTENT(IN)  :: rws  ! Series to fit 
      INTEGER,DIMENSION(n),INTENT(IN)  :: cnt  ! Series counts 
      LOGICAL,INTENT(IN)               :: rise ! Tail can go up
      REAL(8),DIMENSION(n),INTENT(OUT) :: rwp  ! Spline curve
      INTEGER,DIMENSION(n)             :: ssy  ! Stiffness each year
      INTEGER,DIMENSION(1)             :: k
      INTEGER                          :: i=1,j,q
      JD: DO j=n,1,-1              ! 4 or more series only
        IF (cnt(j).GE.4) EXIT JD   ! at old end of curve
      ENDDO JD
      ssy(1:j)=(/(ss+i,i=1,j)/)    ! Spline stiffness
      CALL splinec(j,rws,ssy,rwp)  ! Fit spline
      IF (rise) THEN               ! Only basal area
        rwp(j+1:n)=rwp(j)          ! Linear extend
      ELSE                         ! No rise at end final third
        k=MINLOC(rwp(2*j/3:j)) ; q=k(1)+2*j/3-1 
        rwp(q+1:n)=rwp(q)          ! Linear extension
      ENDIF
      RETURN
      END SUBROUTINE spline3
!-----------------------------------------------------------------
      SUBROUTINE out_err(words) ! Output a message and wait for key press
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: words
      INTEGER :: i
      CALL SETCLR(black)
      CALL MESSAG(words//" - hit key to continue",310,1560)
      CALL CSRPOS(msx,msy,i)
      CALL SETCLR(white)
      CALL AREAF((/300,2100,2100,300/),(/1550,1550,1610,1610/),4)
      CALL SETCLR(black)
      RETURN
      END SUBROUTINE out_err
!---------------------------------------------------------------------
      FUNCTION io_err(words,fname) ! Display input output error message
      IMPLICIT NONE
      LOGICAL                   :: io_err
      CHARACTER(*),INTENT(IN)   :: words,fname
      CHARACTER(6)              :: numb
      INTEGER                   :: i
      io_err=FA
      IF (ios.NE.0) THEN
        CALL SETCLR(black)
        WRITE(numb,'(I6)') ios
        CALL MESSAG(words//TRIM(fname)//" fail err ="//numb,310,1560)
        CALL CSRPOS(msx,msy,i) ;  CALL SETCLR(white)
        CALL AREAF((/300,2100,2100,300/),(/1550,1550,1760,1760/),4)
        CALL SETCLR(black)
      ENDIF
      RETURN
      END FUNCTION io_err
!----------------------------------------------------------------------
      SUBROUTINE but_draw(i,label)   ! Draw a button on the screen
      IMPLICIT NONE
      INTEGER,INTENT(IN)      :: i        ! Button number
      CHARACTER(*),INTENT(IN) :: label    ! Button label
      INTEGER                 :: j
      IF (b(i)%on) THEN ; j=bhigh         ! Button highlight colour
      ELSE              ; j=silver        ! Do not highlight
      ENDIF
      IF (b(i)%ok) THEN                   ! If turned on
        CALL SETCLR(j)                    ! Background box
        CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                   (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
        CALL SETCLR(red)
        CALL SHDPAT(0)                    ! Outline pattern
        CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                   (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
        CALL SHDPAT(16)                   ! Filled pattern
        CALL SETCLR(black)
        IF (LEN_TRIM(label).LT.1) THEN
          CALL MESSAG(b(i)%lab,b(i)%x1+10,b(i)%y1+12)  ! Fixed label
        ELSE
          CALL MESSAG(label,b(i)%x1+10,b(i)%y1+12)     ! Variable label
        ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE but_draw
!--------------------------------------------------------
      SUBROUTINE ch_disp(maxcol,maxi,names,dhs,dhc,tok)  ! Displays list for selection
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)                      :: maxcol ! number of columns
      INTEGER,INTENT(IN)                      :: maxi   ! number of items in menu
      CHARACTER(*),DIMENSION(maxi),INTENT(IN) :: names  ! names of each entry
      INTEGER,INTENT(IN)                      :: dhs    ! number of selected item
      INTEGER,DIMENSION(maxi),INTENT(IN)      :: dhc    ! selected list
      INTEGER,INTENT(IN)                      :: tok    ! If trees
      INTEGER  :: i,j,m,imax
      INTEGER  :: ys,xs,ra,rb,rc 
      imax=MIN(maxcol*rw,maxi)          ! Limited by size of display area
      ys=(bbot-btop)/rw ; xs=(gright-cnx)/maxcol
      DO i=1,dhs                        ! For each selected item
        j=dhc(i)-sccc 
        IF (j.GE.0.AND.j.LT.imax) THEN  ! If selected name visible
          CALL SETCLR(yellow)           ! Highlight selected name
          ra=cnx+xs*(j/rw)              ! Column offset
          rb=btop+MOD(j,rw)*ys          ! row offset
          CALL AREAF((/ra,ra+xs,ra+xs,ra/),(/rb,rb,rb+ys,rb+ys/),4)
        ENDIF
      ENDDO
      ra=cnx ; rb=ra+xs
      DO i=1,(imax-1)/rw+1              ! For each column
        CALL SETCLR(green)
        j=MIN(rw,imax-rw*(i-1)) ; rc=btop-ys
        DO j=0,j                        ! Draw the lines
          rc=rc+ys ; CALL LINE(ra,rc,rb,rc)
        ENDDO
        CALL LINE(ra,btop,ra,rc)
        CALL LINE(rb,btop,rb,rc) ; CALL SETCLR(black) 
        rc=btop+12
        DO j=rw*(i-1),MIN(maxi-sccc+1,rw*i)-1  ! Write the names
          IF (tok.EQ.2) THEN                   ! Trees from sorted list
            CALL MESSAG(names(tre(j+sccc)),ra+5,rc)                
          ELSEIF (tok.EQ.3) THEN               ! RCS no + Trees 
            m=tre(j+sccc)
            WRITE(mess,'(A12,"R",I2," D",I2)') names(m),trr(m),trm(m)
            CALL MESSAG(mess(1:25),ra+5,rc)                
          ELSEIF (tok.EQ.1) THEN
            CALL MESSAG(names(j+sccc),ra+5,rc)                
          ENDIF
          rc=rc+ys     
        ENDDO
        ra=rb ; rb=ra+xs
      ENDDO
      RETURN 
      END SUBROUTINE ch_disp
!-----------------------------------------------------------------
      SUBROUTINE get_int(bno,numb)  ! Read an integer, optional sign
      IMPLICIT NONE
      INTEGER,INTENT(IN)     :: bno
      INTEGER,INTENT(INOUT)  :: numb
      INTEGER                :: i,j,m,kk,osign
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Input Integer")
      CALL SENDBF()  ! Update display
      CALL HEIGHT(30) ; kk=0 ; osign=1
      b(bno)%lab="" ; b(bno)%on=TR ; CALL but_draw(bno,"") ; m=0
      ID: DO i=1,6                     ! Sign and up to 5 digits
        CALL CSRPOS(msx,msy,j) ; b(bno)%lab(i+m:i+m)=CHAR(j)
        CALL MESSAG(b(bno)%lab(1:i+m),b(bno)%x1+10,b(bno)%y1+12)
        IF (i.EQ.1.AND.j.EQ.45) THEN   ! Minus sign
          osign=-1 ; m=1
          CALL CSRPOS(msx,msy,j) ; b(bno)%lab(i+m:i+m)=CHAR(j)
          CALL MESSAG(b(bno)%lab(1:2),b(bno)%x1+10,b(bno)%y1+12) 
        ENDIF
        IF (j.GT.47.AND.j.LT.58) THEN 
          kk=10*kk+(j-48)
        ELSE
          EXIT ID
        ENDIF
      ENDDO ID
      CALL HEIGHT(25) ; IF (j.EQ.13) numb=kk*osign
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Running")
      CALL SENDBF()  ! Update display
      RETURN
      END SUBROUTINE get_int
!---------------------------------------------------------------------
      SUBROUTINE get_real(bno,numb) ! Get user input decimal
      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: bno     ! Button box number
      REAL(8),INTENT(OUT) :: numb    ! Output value
      CHARACTER(60)       :: dum
      LOGICAL             :: dot
      INTEGER             :: i,j,osign
      REAL(8)             :: mult
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Input Real")
      CALL SENDBF()  ! Update display
      CALL HEIGHT(25) ; dot=TR ; numb=0.D0 ; mult=0.1D0
      osign=1 ; dum=b(bno)%lab ; b(bno)%lab=""
      b(bno)%on=TR ; CALL but_draw(bno,"")
      ID: DO i = 1,9                ! Sign and one point, up to 9 digits
        CALL CSRPOS(msx,msy,j) ; b(bno)%lab(i:i)=CHAR(j)
        CALL MESSAG(b(bno)%lab(1:i),b(bno)%x1+10,b(bno)%y1+12)
        IF (i.EQ.1.AND.j.EQ.45) THEN ! Minus sign
          osign=-1
        ELSEIF ((b(bno)%lab(i:i).EQ.".").AND.dot) THEN  ! Decimal point now
          dot=FA
        ELSEIF (j.GT.47.AND.j.LT.58) THEN
          IF (dot) THEN
            numb=10.D0*numb+DBLE(j-48)
          ELSE
            numb=numb+mult*DBLE(j-48) ; mult=mult*0.1D0
          ENDIF
        ELSE
          EXIT ID
        END IF
      ENDDO ID
      IF (j.EQ.13) THEN
        numb=numb*DBLE(osign)
      ELSE
        numb=-999999.D0
      ENDIF
      CALL HEIGHT(25)
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Running")
      CALL SENDBF()  ! Update display
      RETURN
      END SUBROUTINE get_real
!---------------------------------------------------------------------
      SUBROUTINE get_str(bno)      ! Gets a name, up to 12 chars
      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: bno   ! Button to use
      CHARACTER(60)       :: dum
      INTEGER             :: i,j,m
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Input Name")
      CALL SENDBF()  ! Update display
      dum=b(bno)%lab ; b(bno)%lab=" " ; CALL but_draw(bno,"") ; m=0
      ID: DO i=1,60      
        CALL CSRPOS(msx,msy,j) ; b(bno)%lab(i+m:i+m)=CHAR(j)
        CALL MESSAG(b(bno)%lab(1:i),b(bno)%x1+10,b(bno)%y1+12)
        IF (j.EQ.13) THEN
          b(bno)%lab(i:i) = " " ; EXIT ID
        ENDIF
      ENDDO ID   ! If less than 2 chars then not changed
      IF (i.LE.2) b(bno)%lab = dum 
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Running")
      CALL SENDBF()  ! Update display
      RETURN
      END SUBROUTINE get_str
!---------------------------------------------------------------------
      SUBROUTINE mwrite(i)        ! Write message to screen
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: i     ! Message number
      CALL SETCLR(black)
      CALL MESSAG(mtext(i),madd(i)%x1,madd(i)%y1+15)  
      RETURN 
      END SUBROUTINE mwrite
!----------------------------------------------------------
      SUBROUTINE plot_trees(cnt,trees)
      IMPLICIT NONE   ! Shaded tree count area
      INTEGER,INTENT(IN)                :: cnt    ! Data count
      INTEGER,DIMENSION(cnt),INTENT(IN) :: trees  ! Tree counts
      CHARACTER(5),DIMENSION(4)         :: labs
      REAL(8)                           :: yscal  ! Y scale
      REAL(8)                           :: xscal  ! X scale
      INTEGER,DIMENSION(0:cnt+1)        :: zz4,zz5
      INTEGER                           :: i,p,lab
      p=MAXVAL(trees(1:cnt))
      SELECT CASE (p)  ! Select suitable yscale
        CASE (   : 50)  ; lab=50
        CASE (51 :100)  ; lab=100
        CASE (101:160)  ; lab=160
        CASE (161:300)  ; lab=300
        CASE (301:500)  ; lab=500
        CASE (501:900)  ; lab=900
        CASE (901:1200) ; lab=1200
        CASE (1201:)    ; lab=(p/500+1)*500
      END SELECT
      CALL SETCLR(silver) ; CALL HEIGHT(18) ; CALL LINWID(1)
      yscal=DBLE(grb-grt)/DBLE(lab)
      xscal=DBLE(grr-grl)/DBLE(cnt)
      DO i=1,4 ; WRITE(labs(i),'(I5)') (lab*i)/4 ; ENDDO
      zz4(1:cnt)=grb-NINT(DBLE(trees)*yscal)
      zz4(0)=grb ; zz4(cnt+1)=grb
      zz5(1:cnt)=grl+NINT(xscal*(/(DBLE(i),i=1,cnt)/))
      zz5(0)=zz5(1) ; zz5(cnt+1)=zz5(cnt)
      CALL AREAF(zz5,zz4,cnt+2)
 ! Print yscale labels
      CALL SETCLR(black)
      i=5-LEN_TRIM(ADJUSTL(labs(4)))+1
      CALL MESSAG("Sample",grr+8,grt+((grb-grt)*38)/100)
      CALL MESSAG("Count",grr+8,grt+((grb-grt)*62)/100)
      CALL MESSAG(labs(4)(i:5),grr+8,grt)
      CALL MESSAG(labs(3)(i:5),grr+8,grt+((grb-grt)*25)/100)
      CALL MESSAG(labs(2)(i:5),grr+8,grt+((grb-grt)*50)/100)
      CALL MESSAG(labs(1)(i:5),grr+8,grt+((grb-grt)*75)/100)
      CALL HEIGHT(22) ; CALL LINWID(1)
      RETURN
      END SUBROUTINE plot_trees
!--------------------------------------------------------------------
      SUBROUTINE plot_trees3(cnt,trees,years,yax,off,scol)
      IMPLICIT NONE             ! Shaded tree count area
      INTEGER,INTENT(IN)                :: cnt    ! Data count
      INTEGER,INTENT(IN)                :: scol   ! Colour
      INTEGER,DIMENSION(cnt),INTENT(IN) :: trees  ! Tree counts
      REAL(8),DIMENSION(cnt),INTENT(IN) :: years  ! Year
      REAL(8),DIMENSION(0:cnt+1)        :: zz4    ! Year
      REAL(8),DIMENSION(0:cnt+1)        :: zz5    ! Counts
      REAL(8),INTENT(IN)                :: yax    ! Y axis size
      REAL(8),INTENT(IN)                :: off    ! offset
      CHARACTER(5),DIMENSION(4)         :: labs
      INTEGER                           :: i,p,lab
      REAL(8)                           :: yscal  ! Y scale
      p=MAXVAL(trees(1:cnt))
      SELECT CASE (p)  ! Select suitable yscale
        CASE (   : 50)  ; lab=50
        CASE (51 :100)  ; lab=100
        CASE (101:160)  ; lab=160
        CASE (161:300)  ; lab=300
        CASE (301:500)  ; lab=500
        CASE (501:900)  ; lab=900
        CASE (901:1200) ; lab=1200
        CASE (1201:)    ; lab=(p/500+1)*500
      END SELECT
      CALL SETCLR(scol) ; yscal=yax/DBLE(lab)
      DO i=1,4 ; WRITE(labs(i),'(I5)') (lab*i)/4 ; ENDDO
      zz4(1:cnt)=DBLE(trees)*yscal+off
      zz4(0)=off ; zz4(cnt+1)=off
      zz5(1:cnt)=years
      zz5(0)=zz5(1) ; zz5(cnt+1)=zz5(cnt)
      CALL RLAREA(zz5,zz4,cnt+2)
      ! Print yscale labels
      IF (scol.EQ.silver) CALL SETCLR(black) 
      i=5-LEN_TRIM(ADJUSTL(labs(4)))+1
      CALL MESSAG("Sample",grr+8,grb-((grb-grt)*62)/100)
      CALL MESSAG("Count",grr+8,grb-((grb-grt)*38)/100)
      CALL MESSAG(labs(4)(i:5),grr+8,grt)
      CALL MESSAG(labs(3)(i:5),grr+8,grt+((grb-grt)*25)/100)
      CALL MESSAG(labs(2)(i:5),grr+8,grt+((grb-grt)*50)/100)
      CALL MESSAG(labs(1)(i:5),grr+8,grt+((grb-grt)*75)/100)
      RETURN
      END SUBROUTINE plot_trees3
!--------------------------------------------------------------------
      SUBROUTINE plot_treesq(cnt,trees,nox)
      IMPLICIT NONE   ! Shaded tree count area
      INTEGER,INTENT(IN)                :: cnt    ! Data count
      INTEGER,DIMENSION(cnt),INTENT(IN) :: trees  ! Tree counts
      INTEGER,INTENT(IN)                :: nox    ! Maximum count
      CHARACTER(5),DIMENSION(4)         :: labs
      REAL(8)                           :: yscal  ! Y scale
      REAL(8)                           :: xscal  ! X scale
      INTEGER,DIMENSION(0:cnt+1)        :: zz4,zz5
      INTEGER                           :: i,p,lab
      p=MAX(MAXVAL(trees(1:cnt)),nox)
      SELECT CASE (p)  ! Select suitable yscale
        CASE (   : 50)  ; lab=50
        CASE (51 :100)  ; lab=100
        CASE (101:160)  ; lab=160
        CASE (161:300)  ; lab=300
        CASE (301:500)  ; lab=500
        CASE (501:900)  ; lab=900
        CASE (901:1200) ; lab=1200
        CASE (1201:)    ; lab=(p/500+1)*500
      END SELECT
      CALL SETCLR(silver) ; CALL HEIGHT(18) ; CALL LINWID(1)
      yscal=DBLE(grb-grt)/DBLE(lab)
      xscal=DBLE(grr-grl)/DBLE(cnt)
      DO i=1,4 ; WRITE(labs(i),'(I5)') (lab*i)/4 ; ENDDO
      zz4(1:cnt)=grb-NINT(DBLE(trees)*yscal)
      zz4(0)=grb ; zz4(cnt+1)=grb
      zz5(1:cnt)=grl+NINT(xscal*(/(DBLE(i),i=1,cnt)/))
      zz5(0)=zz5(1) ; zz5(cnt+1)=zz5(cnt)
      CALL AREAF(zz5,zz4,cnt+2)
 ! Print yscale labels
      CALL SETCLR(black)
      i=5-LEN_TRIM(ADJUSTL(labs(4)))+1
      CALL MESSAG("Sample",grr+8,grt+((grb-grt)*38)/100)
      CALL MESSAG("Count",grr+8,grt+((grb-grt)*62)/100)
      CALL MESSAG(labs(4)(i:5),grr+8,grt)
      CALL MESSAG(labs(3)(i:5),grr+8,grt+((grb-grt)*25)/100)
      CALL MESSAG(labs(2)(i:5),grr+8,grt+((grb-grt)*50)/100)
      CALL MESSAG(labs(1)(i:5),grr+8,grt+((grb-grt)*75)/100)
      CALL HEIGHT(22) ; CALL LINWID(1)
      RETURN
      END SUBROUTINE plot_treesq
!--------------------------------------------------------------------
      SUBROUTINE pair_sort(n,rdat,idat) ! Returns pointers in sequence
      IMPLICIT NONE                      
      INTEGER,INTENT(IN)                 :: n     ! Number of values
      REAL(8),DIMENSION(n),INTENT(INOUT) :: rdat  ! Data to sort
      INTEGER,DIMENSION(n),INTENT(OUT)   :: idat  ! Pointers in sorted order
      INTEGER,DIMENSION(1)               :: k  
      INTEGER  :: i=1,j,m
      REAL(8)  :: rr
      idat(1:n)=(/(i,i=1,n)/)      ! Numbers 1 to n
      DO i=n,2,-1                  
        k=MAXLOC(rdat(1:i)) ; j=k(1) 
        IF (j.LT.i) THEN
          rr=rdat(j) ; m=idat(j) ; rdat(j:i-1)=rdat(j+1:i)
          idat(j:i-1)=idat(j+1:i) ; rdat(i)=rr ; idat(i)=m
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE pair_sort
!-----------------------------------------------------------------
      SUBROUTINE growth_rate(gr) ! Growth rate relative to 1 RCS
      IMPLICIT NONE                 
      REAL(8),DIMENSION(nc),INTENT(OUT) :: gr   ! Growth rate
      INTEGER,DIMENSION(mxt)            :: cnt
      REAL(8),DIMENSION(mxt)            :: wk
      INTEGER                           :: i,p,q,r,u,v,wly
      wly=1 ; wk=0.D0 ; cnt=0
      DO i=1,nc  ! For each tree or core       
        p=ad(i) ; r=yr(i) ; q=p+r-1  ! Ring address
        IF (poo.EQ.1) THEN  ! Using pith offset
          u=fy(i)-pth(i)+1 ; v=u+yr(i)-1  ! Age address
          wk(1:u-1)=wk(1:u-1)+pthr(i)/DBLE(u-1) ! Include pith in RCS
        ELSE
          u=1 ; v=yr(i)                    ! Age address
        ENDIF
        cnt(1:u-1)=cnt(1:u-1)+1 ; wly=MAX(wly,v)
        WHERE (xok(p:q))  ! Valid measurement
          cnt(u:v)=cnt(u:v)+1 ; wk(u:v)=wk(u:v)+x(p:q) 
        END WHERE
      ENDDO                     ! Mean of each age (RCS curve)
      wk(1)=wk(1)/DBLE(MAX(cnt(1),1))
      DO i=2,wly                ! RCS curve diameter
        wk(i)=wk(i-1)+wk(i)/DBLE(MAX(cnt(i),1))
      ENDDO
      DO i=1,nc  ! Final Radius / RCS curve radius that year
        p=ad(i) ; q=p+yr(i)-1   ! Presumes missing values estimated 
        IF (poo.EQ.1) THEN      ! Use pith offset
          gr(i)=(pthr(i)+SUM(x(p:q)))/wk(ly(i)-pth(i)+1) 
        ELSE 
          gr(i)=(SUM(x(p:q)))/wk(yr(i)) 
        ENDIF
      ENDDO
      RETURN 
      END SUBROUTINE growth_rate
!------------------------------------------------------------------------
      SUBROUTINE tree_sort()  ! Sorts trees into pointer tre(1:nc) 
      IMPLICIT NONE                 
      INTEGER  :: i=1
      IF (tst.EQ.1) THEN
        tre(1:nc)=(/(i,i=1,nc)/) ; RETURN     ! No sort 
      ENDIF
      SELECT CASE (tst) 
      CASE (2)                                ! Sort by age
        IF (poo.EQ.1) THEN                    ! Use pith offset
          tso(1:nc)=DBLE(ly(1:nc)-pth(1:nc)+1) 
        ELSE                                  ! Do not use pith offset
          tso(1:nc)=DBLE(yr(1:nc)) 
        ENDIF
      CASE (3)                                ! Sort by diameter
        DO i=1,nc ; tso(i)=SUM(x(ad(i):ad(i+1)-1)) ; ENDDO
        IF (poo.EQ.1) tso(1:nc)=tso(1:nc)+pthr(1:nc) ! Use pith offset
      CASE (4) ; CALL growth_rate(tso(1:nc))  ! Sort by growth rate
      CASE (5)                                ! Sort by name
        DO i=1,nc
          tso(i)=DBLE(COUNT(LOGICAL(nam(1:nc).LT.nam(i),1)))
        ENDDO
      CASE (6)
        IF (poo.EQ.1) THEN                    ! Use pith offset
          tso(1:nc)=DBLE(pth(1:nc))           ! Sort by tree first year
        ELSE
          tso(1:nc)=DBLE(fy(1:nc))            ! Sort by tree first year
        ENDIF
      CASE (7) ; tso(1:nc)=DBLE(ly(1:nc))     ! Sort by tree last year
      END SELECT
      CALL pair_sort(nc,tso(1:nc),tre(1:nc))          
      RETURN
      END SUBROUTINE tree_sort
!-------------------------------------------------------------------
      SUBROUTINE write_pith(fnam)  ! Write pith offset file
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN) :: fnam
      INTEGER                 :: i,n
      OPEN(28,FILE=fnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",fnam)) RETURN
      DO n=1,nc
        i=tre(n)
        WRITE(28,FMT='(A10,I6,F8.1)',IOSTAT=ios) &
          nam(i)(1:10),pth(i),pthr(i)/10.D0 
        IF (io_err("Write",fnam)) RETURN
      ENDDO 
      CLOSE(28)
      RETURN
      END SUBROUTINE write_pith
!----------------------------------------------------------------
      SUBROUTINE write_raw(tnam,dat)             ! Sub writes RAW files
      IMPLICIT NONE   ! Output record name, year, Value(10)
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat   ! Data set of values
      CHARACTER(*),INTENT(IN)           :: tnam  ! Output file name
      INTEGER,DIMENSION(mxt)            :: dno   ! Integer core measures 
      INTEGER                           :: i,j,k,m,n,p,q,r,yrs
      OPEN(27,FILE=tnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",tnam)) RETURN
      DO n=1,nc                                  ! For each core
        i=tre(n) ; p=ad(i) ; r=yr(i) ; q=p+yr(i)-1
        dno(1:r)=NINT(dat(p:q)*1000.D0) ; dno(r+1)=-9999 ! End of core marker
        WHERE (dno(1:r).EQ.999) dno(1:r)=998
        WHERE (.NOT.xok(p:q)) dno(1:r)=-9990     ! Missing rings negative
        IF (fy(i).GE.0) THEN
          yrs=fy(i)-(fy(i)/10)*10                ! remnant
        ELSE
          yrs=fy(i)-((fy(i)+1)/10-1)*10          ! remnant
        ENDIF
        DO k=1-yrs,r+1,10                        ! For each 10 years 
          j=k+yrs ; m=MIN(k+9,r+1) 
          IF (j+fy(i)-1.LT.-999) THEN            ! Pre -999 name
            WRITE(27,IOSTAT=ios,FMT='(A7,I5,10I6)') &
               nam(i)(1:7),j+fy(i)-1,dno(j:m) 
          ELSE                                   ! Later than -1000 
            WRITE(27,IOSTAT=ios,FMT='(A8,I4,10I6)') &
               nam(i)(1:8),j+fy(i)-1,dno(j:m) 
          ENDIF
          yrs=0 ; IF (io_err("Write",tnam)) RETURN
        ENDDO                            
      ENDDO
      CLOSE(27)
      RETURN
      END SUBROUTINE write_raw
!-----------------------------------------------------------------
      SUBROUTINE write_heidel(tnam,dat) ! Write RAW heidelberg file
      IMPLICIT NONE   
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat   ! Data set of values
      CHARACTER(*),INTENT(IN)           :: tnam  ! Output file name
      INTEGER,DIMENSION(mxt)            :: wka
      CHARACTER(6)                      :: lab
      INTEGER                           :: i,k,n,p,q,r
      OPEN(47,FILE=tnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",tnam)) RETURN
      DO n=1,nc                                  ! For each core
        i=tre(n) 
        WRITE(47,'("HEADER:")')
        WRITE(47,'("KeyCode=",A12)') ADJUSTL(nam(i))
        WRITE(47,'("Unit=1/1000 mm")')
        WRITE(lab,'(I6)') ly(i)
        WRITE(47,'("DateEnd=",A6)') ADJUSTL(lab)
        WRITE(lab,'(I6)') yr(i)
        WRITE(47,'("Length=",A6)') ADJUSTL(lab)
        WRITE(47,'("DATA:Tree")')
        p=ad(i) ; r=yr(i) ; q=p+r-1
        WHERE (xok(p:q)) ! Presumes xok(-) is correct 
          wka(1:r)=NINT(dat(p:q)*1000.D0)
        ELSEWHERE  ! Missing values
          wka(1:r)=-9990
        END WHERE
        wka(yr(i)+1:yr(i)+10)=0 ; k=((yr(i)+9)/10)*10
        WRITE(47,'(10I6)',IOSTAT=ios) wka(1:k)
        IF (io_err("Write ",tnam)) STOP
      ENDDO
      CLOSE(47)
      RETURN
      END SUBROUTINE write_heidel
!-----------------------------------------------------------------
      SUBROUTINE write_compact(fnam,dat) ! Sub saves raw data (compact)
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN)           :: fnam           ! Raw data file name
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat            ! Data set of values
      CHARACTER(8),PARAMETER            :: frm="(19F4.0)" ! Data format
      INTEGER,DIMENSION(mxt)            :: dno            ! Integer core data 
      INTEGER                           :: i,n,p,q,r
      CHARACTER(3)                      :: prec  
      REAL(8)                           :: prc
      OPEN(24,FILE=fnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",fnam)) RETURN
      IF (ANY(dat(1:ad(nc)+yr(nc)-1).GT.9.8D0)) THEN
        prec=" -2" ; prc=100.D0
      ELSE
        prec=" -3" ; prc=1000.D0
      ENDIF
      DO n=1,nc                              ! For each core
        i=tre(n) ; p=ad(i) ; r=yr(i) ; q=p+r-1
        dno(1:r)=NINT(prc*dat(p:q))
        WHERE (.NOT.xok(p:q)) dno(1:r)=-999  ! Missing rings negative
        WRITE(24,IOSTAT=ios,FMT='(I8,"=N",I8,"=I ",A12,35X,A3,A8)') &
           r,fy(i),nam(i),prec,frm                ! Write header 
        WRITE(24,'(19I4)',IOSTAT=ios) dno(1:r)    ! Write data
        IF (io_err("write",nam(i))) STOP
      ENDDO
      CLOSE(24)
      RETURN
      END SUBROUTINE write_compact
!-----------------------------------------------------------------
      SUBROUTINE write_ind(n,zfy,dat,cnt,styp,fnam)  ! Sub writes index files
      IMPLICIT NONE   ! Output record - name, year, 10*(Value,Count),Typ
      INTEGER,INTENT(IN)              :: n     ! Number of data points
      INTEGER,INTENT(IN)              :: zfy   ! First year
      REAL(8),DIMENSION(n),INTENT(IN) :: dat   ! Data
      INTEGER,DIMENSION(n),INTENT(IN) :: cnt   ! Counts
      CHARACTER(3),INTENT(IN)         :: styp  ! File type
      CHARACTER(*),INTENT(IN)         :: fnam  ! Output file name
      TYPE pair ; INTEGER             :: dval,num ; END TYPE
      TYPE(pair),DIMENSION(0:n+19)    :: drec
      INTEGER                         :: j,p,q
      REAL(8)                         :: fscal ! Scaling factor
      IF (idb.EQ.2) THEN      ! Normalised indices (range -4 to +4)      
        fscal=100.D0
      ELSE                    ! Ratios indices (range 0 to +4) 
        fscal=1000.D0
      ENDIF
      IF (zfy.LT.0) THEN      ! Negative first year
        q=((zfy+1)/10-1)*10 ; p=zfy-q 
      ELSE                    ! Positive first year
        q=(zfy/10)*10 ; p=zfy-q
      ENDIF
      drec(0:n+19)%num=0      ! Null values
      drec(0:n+19)%dval=9990  
      WHERE (cnt(1:n).GE.1)
        drec(p:p+n-1)%dval=MIN(NINT(dat(1:n)*fscal),9999)
        drec(p:p+n-1)%num=MIN(cnt(1:n),999)
      END WHERE
      DO j=0,n+p,10          ! For each 10 years 
        IF (q.LT.-999) THEN   ! 5 char name
          WRITE(23,IOSTAT=ios,FMT='(A5,I5,10(I4,I3),"  ",A3)') &
            namc(1:5),q+p,drec(j:j+9),styp 
        ELSE                  ! 6 char name
          WRITE(23,IOSTAT=ios,FMT='(A6,I4,10(I4,I3),"  ",A3)') &
            namc(1:6),q+p,drec(j:j+9),styp 
        ENDIF
        IF (io_err("Write data",fnam)) STOP
        q=q+10 ; p=0
      ENDDO                                      
      RETURN
      END SUBROUTINE write_ind
!------------------------------------------------------------------------
      SUBROUTINE write_index(znam)  ! Sub writes crns.crn file
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN) :: znam  ! Output file name
      CHARACTER(3)            :: typ
      INTEGER                 :: r
      INQUIRE(FILE=znam,EXIST=fileok)  
      IF (append.AND.fileok) THEN
        OPEN(23,FILE=znam,IOSTAT=ios,STATUS="OLD",POSITION="APPEND")
        IF (io_err("Open App",znam)) STOP
      ELSE
        OPEN(23,FILE=znam,IOSTAT=ios,STATUS="REPLACE")
        IF (io_err("Open Replace",znam)) STOP
      ENDIF
      IF     (jrb.GE.2)   THEN ; typ=jrblab(jrb)(1:3)
      ELSEIF (idt.GE. 10) THEN ; typ=crtyp(8)
      ELSEIF (idt.LE.-10) THEN ; typ=crtyp(9) 
      ELSE                     ; typ=crtyp(idt)
      ENDIF
      CALL crn_head() ; r=cyr(cf) ; namc=wnam(cf)
      CALL write_ind(r,cfy(cf),crn(1:r,cf),num(1:r,cf), &
        typ,wnam(cf))
      CLOSE(23) 
      RETURN
      END SUBROUTINE write_index
!-----------------------------------------------------------------
      SUBROUTINE read_open(ref1,fnam)    ! Opens an existing file
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)      :: ref1
      CHARACTER(*),INTENT(IN) :: fnam
      INTEGER                 :: i
      INQUIRE(FILE=fnam,EXIST=fileok)  
      IF (.NOT.fileok) THEN
        i=LEN_TRIM(fnam) ; IF (fnam(i-2:i).NE."pth") &
        CALL out_err("'"//TRIM(fnam)//"' - Does not exist") ; RETURN
      ENDIF
      OPEN(ref1,FILE=fnam,IOSTAT=ios,STATUS="OLD")
      IF (io_err("Open",fnam)) STOP
      RETURN                                                         
      END SUBROUTINE read_open
!-----------------------------------------------------------------------------      
      SUBROUTINE mdet_miss(fr,to1)          ! Infill missing rings
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)       :: fr,to1    ! First/last tree in CRN
      REAL(8),DIMENSION(1:mxy) :: dat,drc   ! Chronology and RCS values
      INTEGER,DIMENSION(1:mxy) :: dn1,dn2   ! Chronology and RCS counts
      REAL(8)                  :: frc       ! Ratio tree mean / RCS mean
      INTEGER                  :: i,j,p,q,r,s,t,u,v,yfy,yly,yyr,rly
      yfy=MINVAL(fy(fr:to1)) ; yly=MAXVAL(ly(fr:to1))
      yyr=yly-yfy+1 ; rly=MAXVAL(ly(fr:to1)-pth(fr:to1)+1)
      dn1(1:yyr)=0 ; dat(1:yyr)=0.D0
      dn2(1:rly)=0 ; drc(1:rly)=0.D0 
      DO i=fr,to1         
        p=ad(i) ; r=yr(i) ; q=p+r-1 ; u=fy(i)-yfy+1 ; v=u+r-1 
        s=fy(i)-pth(i)+1 ; t=s+r-1 
        WHERE (xok(p:q))
          dat(u:v)=dat(u:v)+x(p:q) ; dn1(u:v)=dn1(u:v)+1 ! Ring values & counts
          drc(s:t)=drc(s:t)+x(p:q) ; dn2(s:t)=dn2(s:t)+1 ! RCS values & counts
        END WHERE
      ENDDO                                 ! Mean raw and RCS curve 
      WHERE (dn2(1:rly).GT.0) drc(1:rly)=drc(1:rly)/DBLE(dn2(1:rly))
      WHERE (dn1(1:yyr).GT.0) dat(1:yyr)=dat(1:yyr)/DBLE(dn1(1:yyr))
      DO i=fr,to1                           ! For each tree
        p=ad(i) ; r=yr(i) ; q=p+r-1
        IF (.NOT.ALL(xok(p:q))) THEN        ! Tree has missing ring(s)
          u=fy(i)-yfy+1 ; v=u+r-1 ; s=fy(i)-pth(i)+1 ; t=s+r-1 
          frc=SUM(x(p:q),MASK=xok(p:q))/SUM(drc(s:t),MASK=xok(p:q)) 
          DO j=p,q
            IF (xok(j)) THEN         ! No action
            ELSEIF (dn1(j-p+u).GE.1) THEN
              x(j)=dat(j-p+u)*frc   ! From other trees
            ELSE
              x(j)=drc(j-p+s)*frc   ! From RCS curve
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN 
      END SUBROUTINE mdet_miss
!------------------------------------------------------------------------
      SUBROUTINE read_pith(fnam,fr,to1)  ! Sub reads .pth file
      IMPLICIT NONE  
      CHARACTER(*),INTENT(IN) :: fnam
      CHARACTER(60)           :: tnam
      INTEGER,INTENT(IN)      :: fr,to1  ! From tree to tree
      CHARACTER(12)           :: snam
      INTEGER                 :: i,j
      j=LEN_TRIM(fnam) 
      JD: DO j=j,3,-1 ; IF (fnam(j:j).EQ.".") EXIT JD ; ENDDO JD
      tnam=fnam(1:j)//"pth"
      snam(9:12)="    "
      INQUIRE(FILE=tnam,EXIST=fileok)  
      IF (.NOT.fileok) THEN       ! Try upper case
        tnam=fnam(1:j)//"PTH"
        INQUIRE(FILE=tnam,EXIST=fileok)  
      ENDIF
      IF (.NOT.fileok) THEN       ! No pith file to read
        pth(fr:to1)=fy(fr:to1)-1
        pthr(fr:to1)=0.1D0        ! Use defaults
        RETURN
      ELSE
        CALL read_open(69,tnam)
        DO i=fr,to1
          READ(69,'(A10,I6,F8.0)',IOSTAT=ios) snam(1:10),pth(i),pthr(i) 
          IF (io_err("Read1",nam(i)//tnam)) STOP
          IF (snam(1:8).NE.nam(i)(1:8)) THEN ! Pith names not match
            IF (snam(1:7).NE.nam(i)(1:7).OR.nam(i)(8:8).NE."-") THEN
              CALL out_err("7Pith "//snam//" tree "//nam(i)//tnam)
              pth(i)=fy(i)-1 ; pthr(i)=0.1D0 
            ENDIF
          ENDIF
          IF (pth(i).GE.fy(i)) THEN    ! Pith year too recent
            CALL out_err( &
              "8Pith after first ring "//nam(i)//tnam)
            pth(i)=fy(i)-1 ; pthr(i)=0.1D0 
          ENDIF
          IF (pthr(i).LT.0.1D0) pthr(i)=0.1D0 ! No pith radius
        ENDDO
        CLOSE(69)
      ENDIF
      pthr(fr:to1)=pthr(fr:to1)*10.D0        ! Convert cm to mm
      RETURN
      END SUBROUTINE read_pith
!----------------------------------------------------------------
      SUBROUTINE read_raw(fnam,head)           ! Sub reads RAW file
      IMPLICIT NONE                            ! Format name, year, Value(10)
      CHARACTER(*),INTENT(IN) :: fnam          ! Raw data file name
      INTEGER,INTENT(IN)      :: head          ! Number of header records
      CHARACTER(8)            :: tnam          ! Core name
      INTEGER                 :: i,j,k,m,n
      INTEGER,DIMENSION(mxt)  :: dno           ! Integer core measures 
      CALL read_open(21,fnam)
      DO i=1,head                              ! Read header records
        READ(21,*,IOSTAT=ios)   
        IF (io_err("Read Headers ",fnam)) STOP
      ENDDO         
      nc=nc+1 ; ad(1)=1 ; nam(nc)=" "          ! May be appending trees
      READ(21,'(A8,I4,10I6)',IOSTAT=ios) nam(nc)(1:8),fy(nc),dno(1:10) 
      IF (io_err("Read First ",fnam)) STOP     ! Read 1st data record
      IF (nam(nc)(8:8).EQ."-") fy(nc)=fy(nc)*(-1)   ! Pre year -999
      ID: DO i=nc,mxs                          ! For each core
        IF (fy(i).GE.0) THEN
          m=1-fy(i)+(fy(i)/10+1)*10            ! Next data address +ve years
        ELSE
          m=1-fy(i)+((fy(i)+1)/10)*10          ! Next data address -ve years
        ENDIF
        JDO: DO j=1,mxt+21,10                  ! For each further decade
          READ(21,'(A8,I4,10I6)',IOSTAT=ios) tnam,k,dno(m:m+9)
          m=m+10
          IF (tnam.EQ."        ") THEN
            WRITE(mess,'(I5)') m-11+fy(i)
            CALL out_err("Null name " &
              //TRIM(fnam)//" "//nam(i)//mess(1:5))
            ios=-1       ! Null name so presume end of file
          ENDIF
          IF (ios.EQ.-1) EXIT JDO              ! End of file
          IF (io_err("Read2", &
            TRIM(fnam)//" "//TRIM(tnam)//nam(MAX(i-1,1)))) STOP
          IF (tnam.NE.nam(i)(1:8)) THEN        ! Probably a new core 
            IF (nam(i)(1:8).EQ.tnam(1:7)//"-") THEN
              nam(i)(8:8)=tnam(8:8)            ! Just shorter year
            ELSE
              EXIT JDO                         ! New core
            ENDIF
          ENDIF
          IF (nam(i)(8:8).EQ."-") k=k*(-1)
          IF (k.NE.fy(i)+m-11) EXIT JDO        ! Year changed, name did not???
        ENDDO JDO
        JD1: DO j=m-11,MAX(m-20,1),-1     ! Find core terminating value
          IF (dno(j).EQ.-9999.OR.dno(j).EQ.999) EXIT JD1
        ENDDO JD1
        j=MAX(j,1) 
        yr(i)=j-1 ; ly(i)=fy(i)+j-2 ; n=ad(i) ; ad(i+1)=n+yr(i)
        IF (n+j.GT.mxd-10) THEN
          dno(mxd-n+2)=dno(j) ; j=mxd-n+2
          yr(i)=j-1 ; ly(i)=fy(i)+j-2
          ad(i+1)=n+yr(i) ; ios=-1 ; WRITE(mess,'(I7)') mxd
          CALL out_err("Too many rings "//TRIM(fnam)//tnam//mess(1:7))
        ENDIF
        IF (dno(j).EQ.-9999) THEN
          x(n:n+j-2)=DBLE(dno(1:j-1))/1000.D0  ! Convert to mm 
        ELSEIF (dno(j).EQ.999) THEN
          x(n:n+j-2)=DBLE(dno(1:j-1))/100.D0   ! Convert to mm
        ELSE
          WRITE(mess,'(2I5)') k,ly(i)
          CALL out_err("Bad end " &
            //TRIM(fnam)//" "//tnam//mess(1:10)//nam(MAX(i-1,1)))
          STOP
        ENDIF
        IF (ios.EQ.-1) EXIT ID                 ! End of file
        dno(1:10)=dno(m-10:m-1) ; fy(i+1)=k ; nam(i+1)=tnam//"    "
        IF (tnam(8:8).EQ."-") fy(i+1)=k*(-1)   ! New core found
        IF (nam(i)(8:8).EQ."-") nam(i)(8:8)=" "
      ENDDO ID
      IF (i.GT.mxs) THEN
        i=mxs
        CALL out_err("Too many trees - Last tree read "//nam(i))
      ENDIF
      CLOSE(21) ; CALL read_pith(fnam,nc,i) ; nc=i  ! Read pith data
      cc=MIN(nc,1) 
      RETURN
      END SUBROUTINE read_raw
!----------------------------------------------------------------
      SUBROUTINE read_heidel(fnam)   
      IMPLICIT NONE
      CHARACTER(*)  :: fnam
      CHARACTER(32) :: head
      INTEGER       :: i,j,k,prec
      LOGICAL       :: bok
!     CHARACTER(4),DIMENSION(mxs) :: spec    ! Tree species
!     INTEGER,DIMENSION(mxs)      :: swood   ! sapwood
      CALL read_open(19,fnam)
      READ(19,'(A32)',IOSTAT=ios) head 
      IF (io_err("Readh1 ",fnam)) STOP
      pth(nc+1:mxs)=-1 ; ad(1)=1      ! May be appending trees
      ID: DO i=nc+1,mxs                   
        bok=FA ; prec=1000            ! True/value if length or unit read in
        JD: DO j=1,200
          SELECT CASE (head(1:4))

          CASE ("HEAD") 
          CASE ("Date")   ! Begin, End and Leng can be in any sequence
            IF (Head(1:10).EQ."DateBegin=") THEN  
              READ(head(11:16),'(I6)',IOSTAT=ios) fy(i) 
              WRITE(mess,'(I6)') ad(i)
              IF (io_err("Readh2 ",head//nam(i)//mess(1:6))) STOP
              IF (bok) THEN ; ly(i)=fy(i)+yr(i)-1
              ELSE ; ly(i)=-999999 ; ENDIF
            ELSEIF (Head(1:8).EQ."DateEnd=") THEN
              READ(head(9:14),'(I6)',IOSTAT=ios) ly(i) 
              IF (io_err("Readh3 ",head)) STOP
              IF (bok) THEN ; fy(i)=ly(i)-yr(i)+1
              ELSE ; fy(i)=-999999 ; ENDIF
            ENDIF
          CASE ("Leng") ; READ(head(8:11),'(I4)',IOSTAT=ios) yr(i)
            IF (io_err("Readh4 ",head)) STOP
            bok=TR ; IF (ly(i).EQ.-999999) ly(i)=fy(i)+yr(i)-1
            IF (fy(i).EQ.-999999) fy(i)=ly(i)-yr(i)+1
          CASE ("Unit")
            IF (Head(8:11).EQ."100 ") prec=100  
!         CASE ("Spec")	; spec(i)=head(9:12)
!         CASE ("Pith") ; READ(head(6:9),'(I4)',IOSTAT=ios) pth(i)
!           IF (io_err("Read ",head)) STOP
!         CASE ("SapW") ; READ(head(14:17),'(I4)',IOSTAT=ios) swood(i)
!           IF (io_err("Read ",head)) STOP
          CASE ("KeyC") ; nam(i)=head(9:20)
!         CASE ("Comm") ; WRITE(72,'(A12,"  Com ",A24)') nam(i),head(9:32)
          CASE ("DATA")
            IF (ad(i)+yr(i).GT.mxd-15) THEN
              yr(i)=mxd-ad(i)-10 ; ly(i)=fy(i)+yr(i)-1
              READ(19,'(10F6.0)',IOSTAT=ios) &
                (x(k),k=ad(i),ad(i)+yr(i)-1) 
              IF (io_err("Readh5 ",fnam)) STOP
              WRITE(mess,'(I6)') mxd
              CALL out_err ("Too many rings "//TRIM(fnam) &
                & //TRIM(nam(i))//mess(1:6))
              ios=-9999
            ELSE
              READ(19,'(10F6.0)',IOSTAT=ios) &
                (x(k),k=ad(i),ad(i)+yr(i)-1) 
              IF (io_err("Readh6 ",fnam)) STOP
            ENDIF
            EXIT JD
!         CASE DEFAULT ; WRITE(72,'(2I6,"  ",A12,"  ",A32)') i,yr(i),nam(i),head
          ENDSELECT
          READ(19,'(A32)',IOSTAT=ios) head 
          IF (io_err("Readh7 ",fnam)) STOP
        ENDDO JD
        j=ad(i) ; k=ad(i)+yr(i)-1 ; ad(i+1)=k+1
        x(j:k)=x(j:k)/DBLE(prec)
        IF (ios.EQ.-9999) EXIT ID      ! Ran out of space
        READ(19,'(A32)',IOSTAT=ios) head 
        IF (ios.EQ.-1) EXIT ID         ! End of File 
        IF (io_err("Readh8",fnam)) STOP
      ENDDO ID
      IF (i.GT.mxs) THEN
        i=mxs
        CALL out_err("Too many trees - Last tree read "//nam(i))
      ENDIF
      CLOSE(19) ; CALL read_pith(fnam,nc+1,i) ; nc=i ! Read pith data
      RETURN
      END SUBROUTINE read_heidel
!--------------------------------------------------------------
      SUBROUTINE read_compact(fnam) ! Reads compact data file
      IMPLICIT NONE                 
      CHARACTER(*),INTENT(IN) :: fnam         ! Raw data file name
      CHARACTER(8)            :: frm          ! Data format
      INTEGER                 :: i,j,p,q,prec
      ad(1)=1 ; CALL read_open(41,fnam)
      ID: DO i=nc+1,mxs                       ! For each core
        nam(i)=" "
        READ(41,IOSTAT=ios,FMT='(1X,I7,2X,I8,3X,A12,35X,I3,A8)') &
            yr(i),fy(i),nam(i),prec,frm 
        IF (nam(i)(1:8).EQ."        ") ios=-1 ! Null name so end of file
        IF (ios.EQ.-1) EXIT ID                ! End of File 
        IF (io_err("Readc1",fnam)) STOP
        IF (ad(i)+yr(i).GT.mxd-15) THEN
          yr(i)=mxd-ad(i)-10 ; ly(i)=fy(i)+yr(i)-1
          ad(i+1)=ad(i)+yr(i)
          WRITE(mess,'(I6)') mxd
          CALL out_err ("Too many rings "//TRIM(fnam) &
            & //TRIM(nam(i))//mess(1:6))
          READ(41,IOSTAT=ios,FMT=frm) (x(j),j=ad(i),ad(i+1)-1) 
          IF (io_err("Readc2",fnam)) STOP
          EXIT ID
        ELSE
          ly(i)=fy(i)+yr(i)-1 ; ad(i+1)=ad(i)+yr(i)
          READ(41,IOSTAT=ios,FMT=frm) (x(j),j=ad(i),ad(i+1)-1) 
          IF (io_err("Readc3",fnam)) STOP
        ENDIF
      ENDDO ID
      IF (i.GT.mxs) THEN
        i=mxs
        CALL out_err("Too many trees - Last tree read "//nam(i))
      ENDIF
      CLOSE(41) ; j=nc+1 ; p=ad(j) ; q=ad(i)-1
      IF (prec.EQ.0) THEN
        IF (frm(5:5).EQ."3") x(p:q)=x(p:q)*10.D0 ! Convert to mm 
        x(p:q)=x(p:q)/1000.D0  
      ELSE
        x(p:q)=x(p:q)*(10.D0**prec)
      ENDIF
      nc=i-1 
      CALL read_pith(fnam,j,nc)                ! Read pith data 
      RETURN                                                         
      END SUBROUTINE read_compact
!----------------------------------------------------------------
      SUBROUTINE read_rft(fnam)    ! Works out file type
      IMPLICIT NONE                    
      CHARACTER(*),INTENT(IN) :: fnam
      CHARACTER(85)           :: rec1
      INTEGER                 :: i,j
      INTEGER,DIMENSION(22)   :: wk          
      CALL read_open(37,fnam)
      READ(37,IOSTAT=ios,FMT='(A85)') rec1      
      IF (io_err("Read ",fnam)) STOP
      IF (rec1(1:7).EQ."HEADER:") THEN  ! Heidelberg format
        CLOSE(37) ; CALL read_heidel(fnam)
      ELSEIF (rec1(72:72)//rec1(79:79).EQ."()") THEN ! Compact format
        IF (rec1(76:76).EQ."5".OR.rec1(76:76).EQ."4" &
                              .OR.rec1(76:76).EQ."3") THEN    
          CLOSE(37) ; CALL read_compact(fnam)
        ELSE
          CALL out_err("Compact file unknown type "//fnam) 
          CLOSE(37) ; RETURN
        ENDIF
      ELSE
        ID: DO i=1,4                      ! Could be several headers
          j=LEN_TRIM(rec1)
          IF (j.GE.18.AND.j.LE.72) THEN   ! Range for Tuscon
            READ(rec1(12:12),'(I1)',IOSTAT=ios) wk(21)
            IF ((ios.EQ.0).AND.(j.EQ.6*(10-wk(21))+12)) THEN ! Length for Tuscon
              READ(rec1,'(8X,I4,10I6)',IOSTAT=ios) wk(1:11)    
              IF (ios.EQ.0) THEN          ! Data fits Tuscon file
                CLOSE(37) ; CALL read_raw(fnam,i-1) ; EXIT ID
              ENDIF
            ELSE
              READ(rec1(9:12),'(I4)',IOSTAT=ios) wk(22)  ! BC years?
              IF ((ios.EQ.0).AND.((wk(22).LT.0).OR.((wk(22).GE.1000) &
                .AND.(rec1(8:8).EQ."-"))).AND.(j.EQ.12+6*wk(21))) THEN 
                READ(rec1,'(8X,I4,10I6)',IOSTAT=ios) wk(1:11)    
                IF (ios.EQ.0) THEN          ! Data fits Tuscon file
                  CLOSE(37) ; CALL read_raw(fnam,i-1) ; EXIT ID
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          READ(37,IOSTAT=ios,FMT='(A85)') rec1 ! Could be header so read next line
          IF (ios.NE.0.OR.i.EQ.4) THEN
            CALL out_err("File type unknown "//fnam)
            CLOSE(37) ; STOP 
          ENDIF
        ENDDO ID
      ENDIF
      cc=1 ; i=ad(nc)+yr(nc)-1 ; xok(1:i)=NINT(1000.D0*x(1:i)).GE.0
      IF (.NOT.ALL(xok(1:i))) CALL mdet_miss(1,nc) ! Infill missing rings
      RETURN
      END SUBROUTINE read_rft
!------------------------------------------------------------------------
      SUBROUTINE crn_head()  ! Chronology header
      IMPLICIT NONE   
      CHARACTER(5) :: lab
      CHARACTER(6) :: fmat
      INTEGER      :: j
      WRITE(23,'(20X,A20)',IOSTAT=ios) wnam(cf)
      IF (io_err("Write head ","CRN file")) STOP
      IF     (idt.GE.10) THEN
        WRITE(lab,'(I5)') idt
        mess(1:14)=TRIM(ADJUSTL(lab))//"yr Spline"
      ELSEIF (idt.LE.-10) THEN
        WRITE(lab,'(I5)') idt
        mess(1:14)=TRIM(ADJUSTL(lab))//"% Spline "
      ELSE
        mess(1:14)=idtlab(idt)
      ENDIF
      j=LEN_TRIM(mess(1:14))+1
      mess(j:j+15)=", "//itnlab(itn)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//indlab(ind)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//krblab(krb)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//isblab(isb)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//poolab(poo)  
      j=LEN_TRIM(mess(1:j+15))+1
      IF (idt.EQ.-2) THEN
        IF     (rdt.EQ.6) THEN
          WRITE(lab,'(I5)') rdtno
          mess(j:j+15)=", "//TRIM(ADJUSTL(lab))//"yr Smo   "
        ELSEIF (rdt.EQ.7) THEN
          WRITE(lab,'(I5)') -rdtno
          mess(j:j+15)=", "//TRIM(ADJUSTL(lab))//"% Smo    "
        ELSE
          mess(j:j+15)=", "//rdtlab(rdt)
        ENDIF
        j=LEN_TRIM(mess(1:j+15))+1
      ENDIF
      mess(j:j+15)=", "//idblab(idb)  
      j=LEN_TRIM(mess(1:j+15))
      WRITE(fmat,'("(A",I3,")")') j

      WRITE(23,FMT=fmat,IOSTAT=ios) mess(1:j)
      IF (io_err("Write head2 ","CRN file")) STOP
      IF (idt.EQ.-1) THEN
        mess=TRIM(trclab(trc))//", "
        j=LEN_TRIM(trclab(trc))+3
      ELSE
        j=1
      ENDIF
      mess(j:j+13)=tstlab(tst)  
      j=LEN_TRIM(mess(1:j+13))+1
      mess(j:j+15)=", "//sfolab(sfo)  
      j=LEN_TRIM(mess(1:j+15))+1
      IF (sfo.EQ.2) THEN
        WRITE(lab,'(I5)') sfono
        mess(j:j+7)="("//TRIM(ADJUSTL(lab))//") "
        j=LEN_TRIM(mess(1:j+7))+1
      ENDIF
      IF (idt.EQ.-2) THEN      ! If RCS
        mess(j:j+15)=", "//srclab(src)  
        j=LEN_TRIM(mess(1:j+15))+1
        IF (src.EQ.2) THEN
          WRITE(lab,'(I5)') srcno
          mess(j:j+7)="("//TRIM(ADJUSTL(lab))//") "
          j=LEN_TRIM(mess(1:j+7))+1
        ENDIF
      ENDIF  
      mess(j:j+15)=", "//gtrlab(gtr)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//bfclab(bfc)  
      j=LEN_TRIM(mess(1:j+15))+1
      mess(j:j+15)=", "//jrblab(jrb)  
      j=LEN_TRIM(mess(1:j+15))
      WRITE(fmat,'("(A",I3,")")') j
      WRITE(23,FMT=fmat,IOSTAT=ios) mess(1:j)
      IF (io_err("Write head2"," CRN file")) STOP
      RETURN
      END SUBROUTINE crn_head
!------------------------------------------------------------------------
      SUBROUTINE save_tucson(fnam)  ! Sub writes Tuscon chronology
      IMPLICIT NONE   
      CHARACTER(2),DIMENSION(12),PARAMETER :: styp=(/"1 ","2 ", &
        "3 ","4 ","5 ","6 ","7 ","8 ","9 ","10","11","12"/)
      CHARACTER(60),INTENT(INOUT) :: fnam   ! Output CRN name
      CHARACTER(3)                :: typ
      INTEGER                     :: i,j,k,p,q,u
      k=LEN_TRIM(wnam(cf))
      JD: DO j=k,k-3,-1 ; IF (wnam(cf)(j:j).EQ.".") EXIT JD ; ENDDO JD
      IF (j.LT.k-4) j=k
      j=j+1 ; k=j+2
      IF (b(211)%on) THEN              ! Mean RCS
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="M"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            p=sfy(i) ; q=sly(i) ; u=q-p+1 
            CALL write_ind(u,p,mval(p:q,i),mcnt(p:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="RCM" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        p=sfy(mx) ; q=sly(mx) ; u=q-p+1 
        CALL write_ind(u,p,mval(p:q,mx),mcnt(p:q,mx),typ,fnam)
      ENDIF
      IF (b(212)%on) THEN              ! Mean RCS
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="S"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            p=sfy(i) ; q=sly(i) ; u=q-p+1 
            CALL write_ind(u,p,msmo(p:q,i),mcnt(p:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="RCS" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        p=sfy(mx) ; q=sly(mx) ; u=q-p+1 
        CALL write_ind(u,p,msmo(p:q,mx),mcnt(p:q,mx),typ,fnam)
      ENDIF
      IF (b(213)%on) THEN              ! RCS Standard Deviation
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="D"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            p=sfy(i) ; q=sly(i) ; u=q-p+1 
            CALL write_ind(u,p,mssd(p:q,i),mcnt(p:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="RCD" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        p=sfy(mx) ; q=sly(mx) ; u=q-p+1 
        CALL write_ind(u,p,mssd(p:q,mx),mcnt(p:q,mx),typ,fnam)
      ENDIF
      IF (b(214)%on) THEN              ! RCS Standard Error
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="E"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            p=sfy(i) ; q=sly(i) ; u=q-p+1 
            CALL write_ind(u,p,mserr(p:q,i),mcnt(p:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="RCE" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        p=sfy(mx) ; q=sly(mx) ; u=q-p+1 
        CALL write_ind(u,p,mserr(p:q,mx),mcnt(p:q,mx),typ,fnam)
      ENDIF
      q=xyr
      IF (b(215)%on) THEN              ! Standard Indices
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="C"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            CALL write_ind(q,xfy,xcrn(1:q,i),xnum(1:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ=jrblab(jrb)(1:3) ; wnam(cf)(j:k)=typ ; CALL crn_head()
        CALL write_ind(q,xfy,xcrn(1:q,mx),xnum(1:q,mx),typ,fnam)
      ENDIF
      IF (b(216)%on) THEN              ! CRN standard deviation
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="D"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            CALL write_ind(q,xfy,xcsd(1:q,i),xnum(1:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="CSD" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        CALL write_ind(q,xfy,xcsd(1:q,mx),xnum(1:q,mx),typ,fnam)
      ENDIF
      IF (b(217)%on) THEN               ! Smoothed CRN
        IF (b(218)%on) THEN
          DO i=1,srcno
            typ="S"//styp(i) ; wnam(cf)(j:k)=typ ; CALL crn_head()
            CALL write_ind(q,xfy,xcsm(1:q,i),xnum(1:q,i),typ,fnam)
          ENDDO
        ENDIF
        typ="CSM" ; wnam(cf)(j:k)=typ ; CALL crn_head()
        CALL write_ind(q,xfy,xcsm(1:q,mx),xnum(1:q,mx),typ,fnam)
      ENDIF
      RETURN
      END SUBROUTINE save_tucson
!-----------------------------------------------------------------
      SUBROUTINE save_column()   ! Sub column chronology
      IMPLICIT NONE   
      CHARACTER(2),DIMENSION(12),PARAMETER  :: styp= (/"1 ","2 ", &
        "3 ","4 ","5 ","6 ","7 ","8 ","9 ","10","11","12"/) 
      INTEGER            :: j,m,p,q
      CHARACTER(402)     :: slin      
      CHARACTER(6)       :: fmat
      IF (ANY(b(215:217)%on)) THEN   ! Chronology columns
        slin=" " ; p=1 ; q=6 ; slin(p:q)="  Year"  
        IF (b(218)%on) THEN          ! If multi 
          DO m=1,srcno               ! Counts
            p=q+1 ; q=p+4 ; slin(p:q)="  N"//styp(m)
          ENDDO
        ENDIF
        p=q+1 ; q=p+4 ; slin(p:q)="  num"
        IF (b(215)%on) THEN
          IF (b(218)%on) THEN          ! If multi 
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    C"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)="    "//jrblab(jrb)(1:3)
        ENDIF
        IF (b(216)%on) THEN
          IF (b(218)%on) THEN          ! If multi 
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    D"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)=" CRNdev"
        ENDIF
        IF (b(217)%on) THEN 
          IF (b(218)%on) THEN          ! If multi 
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    S"//styp(m)
            ENDDO
          ENDIF 
          p=q+1 ; q=p+6 ; slin(p:q)=" CRNsmo"
        ENDIF
        WRITE(fmat,'("(A",I3,")")') q
        WRITE(23,FMT=fmat) slin(1:q) 
        DO j=1,xyr
          slin=" " ; p=1 ; q=6 ; WRITE(slin(p:q),'(I6)') xfy-1+j   
          IF (b(218)%on) THEN          ! If multi 
            DO m=1,srcno
              p=q+1 ; q=p+4 ; WRITE(slin(p:q),'(I5)') xnum(j,m)
            ENDDO
          ENDIF 
          p=q+1 ; q=p+4 ; WRITE(slin(p:q),'(I5)') xnum(j,mx)
          IF (b(215)%on) THEN
            IF (b(218)%on) THEN          ! If multi 
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcrn(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcrn(j,mx)
          ENDIF
          IF (b(216)%on) THEN
            IF (b(218)%on) THEN          ! If multi 
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcsd(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcsd(j,mx)
          ENDIF

          IF (b(217)%on) THEN
            IF (b(218)%on) THEN          ! If multi 
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcsm(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') xcsm(j,mx)
          ENDIF
          WRITE(fmat,'("(A",I3,")")') q
          WRITE(23,FMT=fmat) slin(1:q) 
        ENDDO
      ENDIF
      IF (ANY(b(211:214)%on)) THEN      ! RCS columns
        WRITE(23,*)
        slin=" " ; p=1 ; q=6 ; slin(p:q)="   Age"  
        IF (b(218)%on) THEN             ! If multi 
          DO m=1,srcno                    ! Counts
            p=q+1 ; q=p+4 ; slin(p:q)="  N"//styp(m)
          ENDDO
        ENDIF
        p=q+1 ; q=p+4 ; slin(p:q)="  num"
        IF (b(211)%on) THEN
          IF (b(218)%on) THEN           ! If multi
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    M"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)=" RCmean"
        ENDIF
        IF (b(212)%on) THEN
          IF (b(218)%on) THEN           ! If multi   
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    S"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)="    RCS"
        ENDIF
        IF (b(213)%on) THEN
          IF (b(218)%on) THEN           ! If multi
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    D"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)=" RCSDev"
        ENDIF
        IF (b(214)%on) THEN
          IF (b(218)%on) THEN           ! If multi
            DO m=1,srcno
              p=q+1 ; q=p+6 ; slin(p:q)="    E"//styp(m)
            ENDDO
          ENDIF
          p=q+1 ; q=p+6 ; slin(p:q)=" RCSerr"
        ENDIF
        WRITE(fmat,'("(A",I3,")")') q
        WRITE(23,FMT=fmat) slin(1:q) 
        DO j=1,sly(mx)
          slin=" " ; p=1 ; q=6 ; WRITE(slin(p:q),'(I6)') j   
          IF (b(218)%on) THEN           ! If multi
            DO m=1,srcno
              p=q+1 ; q=p+4 ; WRITE(slin(p:q),'(I5)') mcnt(j,m)
            ENDDO
          ENDIF 
          p=q+1 ; q=p+4 ; WRITE(slin(p:q),'(I5)') mcnt(j,mx)
          IF (b(211)%on) THEN
            IF (b(218)%on) THEN           ! If multi
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mval(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mval(j,mx)
          ENDIF
          IF (b(212)%on) THEN
            IF (b(218)%on) THEN           ! If multi
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') msmo(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') msmo(j,mx)
          ENDIF
          IF (b(213)%on) THEN
            IF (b(218)%on) THEN           ! If multi
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mssd(j,m)
              ENDDO
            ENDIF  
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mssd(j,mx)
          ENDIF
          IF (b(214)%on) THEN
            IF (b(218)%on) THEN           ! If multi
              DO m=1,srcno
                p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mserr(j,m)
              ENDDO
            ENDIF
            p=q+1 ; q=p+6 ; WRITE(slin(p:q),'(F7.3)') mserr(j,mx)
          ENDIF
          WRITE(fmat,'("(A",I3,")")') q
          WRITE(23,FMT=fmat) slin(1:q) 
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE save_column
!-----------------------------------------------------------------
      SUBROUTINE save_crn()       ! Sub writes chronology files
      IMPLICIT NONE   
      CHARACTER(60)      :: fnam     ! Filename
      LOGICAL            :: openc
      INTEGER            :: i,j
      i=LEN_TRIM(dirc) ; j=LEN_TRIM(namc) 
      fnam=dirc(1:i)//"/"//namc(1:j)//cftext(cft)//" "
      j=LEN_TRIM(cnam(cf))
      ID: DO i=j,1,-1
        IF (cnam(cf)(i:i).EQ."/".OR.cnam(cf)(i:i).EQ.CHAR(92)) EXIT ID
      ENDDO ID
      j=MIN(j,i+20)
      wnam(cf)=cnam(cf)(i+1:j)   ! Short filename - no directory
      INQUIRE(FILE=fnam,EXIST=openc)  
      IF (append.AND.openc) THEN     ! Appending data to a file
        OPEN(23,FILE=fnam,IOSTAT=ios,STATUS="OLD",POSITION="APPEND")
      ELSE                           ! Open empty file
        OPEN(23,FILE=fnam,IOSTAT=ios,STATUS="REPLACE")
      ENDIF
      IF (io_err("Open",fnam)) STOP
      IF (cft.EQ.1) THEN             ! Tucson format 
        CALL save_tucson(fnam)
      ELSE                           ! Column format
        CALL crn_head() ; CALL save_column()
      ENDIF
      CLOSE(23) 
      RETURN
      END SUBROUTINE save_crn
!-----------------------------------------------------------------
      SUBROUTINE save_data()   ! Sub writes data files
      IMPLICIT NONE   
      CHARACTER(60) :: fnam    ! Filename
      INTEGER       :: i,j,k
      i=LEN_TRIM(dirr) ; j=LEN_TRIM(namr) 
      DO k=219,224             ! Save selected raw data
        IF (b(k)%on) THEN      ! Need to save data 
          fnam=dirr(1:i)//"/"//namr(1:j)//datext(k-218)
          SELECT CASE ((k-219)*3+rft)
          CASE ( 1) ; CALL write_raw(fnam,x)
          CASE ( 2) ; CALL write_heidel(fnam,x)
          CASE ( 3) ; CALL write_compact(fnam,x)
          CASE ( 4) ; CALL write_raw(fnam,tx)
          CASE ( 5) ; CALL write_heidel(fnam,tx)
          CASE ( 6) ; CALL write_compact(fnam,tx)
          CASE ( 7) ; CALL write_raw(fnam,fx)
          CASE ( 8) ; CALL write_heidel(fnam,fx)
          CASE ( 9) ; CALL write_compact(fnam,fx)
          CASE (10) ; CALL write_raw(fnam,cx)
          CASE (11) ; CALL write_heidel(fnam,cx)
          CASE (12) ; CALL write_compact(fnam,cx)
          CASE (13) ; CALL write_raw(fnam,dx)
          CASE (14) ; CALL write_heidel(fnam,dx)
          CASE (15) ; CALL write_compact(fnam,dx)
          CASE (16) ; CALL write_raw(fnam,ax)
          CASE (17) ; CALL write_heidel(fnam,ax)
          CASE (18) ; CALL write_compact(fnam,ax)
          END SELECT     
        ENDIF
      ENDDO
      IF (b(219)%on) THEN  ! Sub writes pith if raw data
        fnam=dirr(1:i)//"/"//namr(1:j)//".pth"
        CALL write_pith(fnam)  
      ENDIF
      RETURN
      END SUBROUTINE save_data
!-----------------------------------------------------------------
      SUBROUTINE mndf(x,nok,n,nseg,mn,df,iopt)
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n,iopt
      INTEGER,INTENT(OUT)              :: nseg
      LOGICAL,DIMENSION(n)             :: nok     ! Value present
      REAL(8),DIMENSION(n),INTENT(IN)  :: x
      REAL(8),DIMENSION(n),INTENT(OUT) :: mn,df
      INTEGER :: i
      REAL(8) :: xm,dif
      nseg=0
      DO i=1,n-1
        IF (nok(i).AND.nok(i+1)) THEN
          IF(iopt.EQ.1) xm=(x(i)+x(i+1))/2.D0
          IF(iopt.EQ.2) xm=x(i+1)
          DIF=ABS(x(i+1)-x(i))
          IF(xm.GT.0.000001D0.AND.dif.GT.0.000001D0)THEN
            nseg=nseg+1 ; mn(nseg)=DLOG(xm) ; df(nseg)=DLOG(dif)
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE mndf
!--------------------------------------------------------------------
      SUBROUTINE pearsn(X,Y,N,R,A,B)
      IMPLICIT NONE
      INTEGER,INTENT(IN)              :: N
      REAL(8),DIMENSION(N),INTENT(IN) :: X,Y
      REAL(8),INTENT(OUT)             :: R,A,B
      REAL(8)                         :: XA,AY,SXX,SYY,SXY,RN
      RN=DBLE(N) ; XA=SUM(X(1:N))/RN ; AY=SUM(Y(1:N))/RN
      SXX=SUM((X(1:N)-XA)**2) ; SYY=SUM((Y(1:N)-AY)**2)
      SXY=SUM((Y(1:N)-AY)*(X(1:N)-XA)) ; R=SXY/SQRT(SXX*SYY)
      B=R*SQRT(SYY/SXX) ; A=AY-B*XA
      RETURN
      END SUBROUTINE pearsn
!--------------------------------------------------------------------
      SUBROUTINE trnfrm() ! Power transformation (ARSTAN opt 4)
      IMPLICIT NONE       ! Converts tx to tx
      INTEGER :: j,m,n,q,nseg
      REAL(8) :: r,a,bb,p,rmin,rminp,range1,rangep
      REAL(8),DIMENSION(mxy) :: wka,wkb
      REAL(8),PARAMETER      :: epsilon=0.0001D0
      DO j=1,nc   
        m=ad(j) ; q=yr(j) ; n=m+q-1
        rmin=MINVAL(tx(m:n)) ; range1=MAXVAL(tx(m:n))-rmin
        CALL MNDF(tx(m:n),xok(m:n),q,nseg,wka(1:q),wkb(1:q),2)
        CALL PEARSN(wka(1:q),wkb(1:q),nseg,r,a,bb)
        p=1.D0-bb ; IF(p.LT.0.05D0) p=0.D0
        p=MIN(p,1.D0)
        IF(p.GT.epsilon)THEN
          tx(m:n)=MAX(tx(m:n),0.01D0)**p
        ELSE
          tx(m:n)=DLOG(MAX(tx(m:n),0.01D0))
        ENDIF
        rminp=MINVAL(tx(m:n)) ; rangep=(MAXVAL(tx(m:n))-rminp)
        tx(m:n)=(tx(m:n)-rminp)*range1/rangep+rmin
      ENDDO
      RETURN
      END SUBROUTINE trnfrm
!-----------------------------------------------------------------
      SUBROUTINE Rbar(dat,rb,sdev) ! RBar - mean and Sdev
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat
      REAL(8),INTENT(OUT)               :: rb,sdev
      INTEGER  :: i,j,p,q,r,u,v,s,t,nr
      REAL(8)  :: cor,sd,sm,lx,xx,yy
      sm=0.D0 ; sd=0.D0 ; nr=0
      DO i=1,nc-1 ; DO j=i+1,nc
        p=MAX(fy(i),fy(j)) ; q=MIN(ly(i),ly(j)) ; r=q-p+1
        IF(r.GT.20)THEN  ! >20 year overlap needed
          u=ad(i)+p-fy(i) ; v=ad(i)+q-fy(i)
          s=ad(j)+p-fy(j) ; t=ad(j)+q-fy(j)
          lx=DBLE(r) ; xx=SUM(dat(u:v))/lx
          yy=SUM(dat(s:t))/lx 
          cor=SUM((dat(u:v)-xx)*(dat(s:t)-yy))/ &
            (SQRT(SUM((dat(u:v)-xx)**2)*SUM((dat(s:t)-yy)**2)))
          sm=sm+cor ; sd=sd+cor**2 ; nr=nr+1
        ENDIF
      ENDDO ; ENDDO
      rb=sm/DBLE(nr) ; sdev=SQRT((sd-sm*rb)/DBLE(nr-1))  
      RETURN 
      END SUBROUTINE Rbar
!------------------------------------------------------------------------
      SUBROUTINE stabVar() ! Count weighted variance adjustment
      IMPLICIT NONE        ! Uses dx (for RBar), xcrn and xnum
      REAL(8),DIMENSION(mxy) :: rc      ! Effective counts
      REAL(8)                :: rbr,rn,xm1,sd1,xm2,sd2
      INTEGER                :: w 
      CALL Rbar(dx,rbr,rn) ; w=xyr
      IF(w.GT.1.AND.rbr.GE.epsi)THEN
        rn=DBLE(COUNT(cok(1:w,mx)))               ! Count
        xm1=SUM(xcrn(1:w,mx),MASK=cok(1:w,mx))/rn ! Mean
        sd1=SQRT(SUM((xcrn(1:w,mx)-xm1)**2,MASK=cok(1:w,mx))/(rn-1.D0)) ! SDev
        WHERE (cok(1:w,mx)) 
          rc(1:w)=DBLE(xnum(1:w,mx))/(1.D0+DBLE(xnum(1:w,mx)-1)*rbr)
          xcrn(1:w,mx)=(xcrn(1:w,mx)-xm1)*SQRT(rc(1:w))+xm1
        ELSEWHERE
          xcrn(1:w,mx)=1.D0    
        END WHERE
        xm2=SUM(xcrn(1:w,mx),MASK=cok(1:w,mx))/rn     ! Mean
        sd2=SQRT(SUM((xcrn(1:w,mx)-xm2)**2,MASK=cok(1:w,mx))/(rn-1.D0)) ! SDev
        xcrn(1:w,mx)=((xcrn(1:w,mx)-xm2)/sd2)*sd1+xm1 ! Reset mean and SDev
      ENDIF
      RETURN
      END SUBROUTINE stabVar
!----------------------------------------------------------------------
      SUBROUTINE stabVar2(sm) ! High frequency only variance correction
      IMPLICIT NONE           ! Count weighted variance adjustment
      INTEGER,INTENT(IN)     :: sm   ! Smoothing stiffness
      REAL(8),DIMENSION(mxd) :: zx
      REAL(8),DIMENSION(mxy) :: rc   ! Effective counts
      REAL(8)                :: rbr,rn,xm1,sd1,xm2,sd2
      INTEGER                :: i,p,q,r,w 
      DO i=1,nc     ! Need smoothed indices for RBar
        r=yr(i) ; p=ad(i) ; q=p+r-1
        CALL splinet(r,dx(p:q),sm,zx(p:q))
        zx(p:q)=dx(p:q)-zx(p:q)
      ENDDO
      CALL Rbar(zx,rbr,rn) ; w=xyr  ! Calculate RBar and smooth chronology
      CALL spline_miss(w,xcrn(1:w,mx),sm,zx(1:w),cok(1:w,mx))
      xcrn(1:w,mx)=xcrn(1:w,mx)-zx(1:w)   ! High pass chronology
      IF(w.GT.1.AND.rbr.GE.epsi)THEN
        rn=DBLE(COUNT(cok(1:w,mx)))               ! Count
        xm1=SUM(xcrn(1:w,mx),MASK=cok(1:w,mx))/rn ! Mean
        sd1=SQRT(SUM((xcrn(1:w,mx)-xm1)**2,MASK=cok(1:w,mx))/(rn-1.D0)) ! SDev
        WHERE (cok(1:w,mx)) 
          rc(1:w)=DBLE(xnum(1:w,mx))/(1.D0+DBLE(xnum(1:w,mx)-1)*rbr)
          xcrn(1:w,mx)=(xcrn(1:w,mx)-xm1)*SQRT(rc(1:w))+xm1
        ELSEWHERE
          xcrn(1:w,mx)=1.D0    
        END WHERE
        xm2=SUM(xcrn(1:w,mx),MASK=cok(1:w,mx))/rn     ! Mean
        sd2=SQRT(SUM((xcrn(1:w,mx)-xm2)**2,MASK=cok(1:w,mx))/(rn-1.D0)) ! SDev
        xcrn(1:w,mx)=((xcrn(1:w,mx)-xm2)/sd2)*sd1+xm1 ! Reset mean and SDev
      ENDIF
      xcrn(1:w,mx)=xcrn(1:w,mx)+zx(1:w)  ! Put back low frequency
      RETURN
      END SUBROUTINE stabVar2
!---------------------------------------------------------------------------
      SUBROUTINE albino(dat,nsok,n,phi)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: n
      REAL(8),DIMENSION(N),INTENT(INOUT) :: dat
      REAL(8),DIMENSION(12),INTENT(IN)   :: phi
      LOGICAL,DIMENSION(N),INTENT(IN)    :: nsok
      REAL(8),DIMENSION(N+12)            :: wrk
      REAL(8)                            :: xm
      INTEGER                            :: i,j,ip
      ip=NINT(phi(12))
      xm=SUM(dat,MASK=nsok)/DBLE(COUNT(nsok))
      WHERE (nsok)
        wrk(1+ip:n+ip)=dat(1:n)-xm
      ELSEWHERE
        wrk(1+ip:n+ip)=0.D0
      END WHERE
      DO i=ip,1,-1 ; wrk(i)=SUM(wrk(i+1:i+ip)*phi(1:ip)) ; ENDDO
      DO i=ip+1,ip+n
        dat(i-ip)=wrk(i)-SUM((/(phi(j)*wrk(i-j),j=1,ip)/))+xm
      ENDDO
      WHERE (.NOT.nsok) dat=1.D0
      RETURN
      END SUBROUTINE albino
!--------------------------------------------------------------
      SUBROUTINE commie(dat,nsok,n,phi)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: n
      REAL(8),DIMENSION(N),INTENT(INOUT) :: dat
      REAL(8),DIMENSION(12),INTENT(IN)   :: phi
      LOGICAL,DIMENSION(N),INTENT(IN)    :: nsok
      REAL(8),DIMENSION(N+12)            :: wrk
      REAL(8)                            :: xm
      INTEGER                            :: i,j,ip
      ip=NINT(phi(12))
      xm=SUM(dat,MASK=nsok)/DBLE(COUNT(nsok))
      WHERE (nsok)
        wrk(1+ip:n+ip)=dat(1:n)-xm
      ELSEWHERE
        wrk(1+ip:n+ip)=0.D0
      ENDWHERE
      DO i=ip,1,-1 ; wrk(i)=SUM(wrk(i+1:i+ip)*phi(1:ip)) ; ENDDO ! Backcast
      DO i=ip+1,ip+n
        IF(nsok(i-ip))THEN
          wrk(i)=wrk(i)+SUM((/(wrk(i-j)*phi(j),j=1,ip)/))
          dat(i-ip)=wrk(i)+xm
        ELSE
          dat(i-ip)=1.D0
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE commie
!----------------------------------------------------------------------- 
      SUBROUTINE memcof(DAT,NSok,N,M,PHI,IOPT)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                 :: N,M,IOPT
      REAL(8),DIMENSION(N),INTENT(INOUT) :: DAT
      LOGICAL,DIMENSION(N),INTENT(IN)    :: NSok
      REAL(8),DIMENSION(12),INTENT(OUT)  :: PHI
      REAL(8),DIMENSION(N)               :: WK1,WK2
      REAL(8),DIMENSION(12)              :: WKM
      REAL(8)    :: RN,RN1,XM,PMM,SST,AICM,AICMM,HT,RK,PNEUM,DENOM
      REAL(8)    :: AIC,PM
      INTEGER :: J,K,IP
      WKM=0.D0 ; PHI=0.D0 ; RN=DBLE(N) ; RN1=DBLE(COUNT(NSok))
      XM=SUM(DAT,MASK=NSok)/RN1
      WHERE(NSok) DAT=DAT-XM
      PM=SUM(DAT**2,MASK=NSok)/RN1
      PMM=PM ; SST=PMM ; AICM=RN1*DLOG(PM)+2.D0 ; AICMM=AICM
      WK1(1:N-1)=DAT(1:N-1) ; WK2(1:N-1)=DAT(2:N)
      KDO: DO K=1,M
        PNEUM=SUM(WK1(1:N-K)*WK2(1:N-K))
        DENOM=SUM(WK1(1:N-K)**2+WK2(1:N-K)**2)
        PHI(K)=2.D0*PNEUM/DENOM
        PM=PM*(1.D0-PHI(K)**2) ; RK=DBLE(K)
        HT=2.D0*(RK+1.D0)*(RK+2.D0)/(RN1-RK-2.D0)
        AIC=RN1*DLOG(PM)+2.D0*(RK+1.D0)+HT
        IF(K.NE.1) PHI(1:K-1)= &
          WKM(1:K-1)-(/(PHI(K)*WKM(J),J=K-1,1,-1)/)
        IF(IOPT.LT.1)THEN
          IF(K.GT.1.AND.AIC.GT.AICM)THEN
            IP=K-1 ; PM=1.D0-PM/SST ; AIC=AICM
            PHI(1:IP)=WKM(1:IP)
            IF(AICM.GT.AICMM)THEN
              IP=0 ; PM=0.D0 ; AIC=AICMM
            ENDIF
            PHI(11)=PM ; PHI(12)=DBLE(IP)
            WHERE(NSok) DAT=DAT+XM
            RETURN
          ENDIF
        ELSE
          IF(K.EQ.M) EXIT KDO
        ENDIF
        WKM(1:K)=PHI(1:K)
        DO J=1,N-K-1
          IF(NSok(J).AND.NSok(J+1))THEN
            WK1(J)=WK1(J)-WKM(K)*WK2(J)
            WK2(J)=WK2(J+1)-WKM(K)*WK1(J+1)
          ENDIF
        END DO
        PMM=PM ; AICM=AIC
      END DO KDO
      IP=M-1 ; PM=1.D0-PM/SST ; AIC=AICM
      PHI(1:IP)=WKM(1:IP) ; PHI(11)=PM ; PHI(12)=DBLE(IP)
      WHERE(NSok) DAT=DAT+XM
      RETURN
      END SUBROUTINE memcof
!------------------------------------------------------------
! By: Edward R Cook performs pooled autoregression analysis
! estimates common persistence structure of tree-ring series.
      SUBROUTINE pool(darps,ipr)
      IMPLICIT NONE
      REAL(8),DIMENSION(12),INTENT(OUT) :: DARPS  ! Pooled AR Coefficients
      INTEGER,INTENT(OUT)     :: ipr              ! Max AR lag years
      REAL(8)                 :: DP               ! Output as RSQ in original
      REAL(8),DIMENSION(12)   :: AVARPS,GR,AIC,WK
      INTEGER                 :: I,J,K,NN,FT,LT,p,q
      REAL(8)                 :: RNT,GSUM,VP,AICM
      LOGICAL                 :: TEST
      DO J=1,NC
        p=ad(j) ; q=p+yr(j)-1
        ax(p:q)=dx(p:q)-SUM(dx(p:q))/DBLE(yr(j)) ! Remove mean from series 
      ENDDO
      DO K=1,12         ! Calculate and pool lag product sums
        GSUM=0.D0
        DO I=1,NC 
          DO J=1,NC 
            FT=MAX(FY(I),FY(J)) ; LT=MIN(LY(I),LY(J))-K+1
            IF (LT-FT.GE.0) THEN       ! If overlap exists
              GSUM=GSUM+SUM(ax(ad(i)+FT-fy(i)    :ad(i)+LT-fy(i))* &
                            ax(ad(j)+FT-fy(j)+K-1:ad(j)+LT-fy(j)+K-1))  
            ENDIF       ! Calculate lag products
          ENDDO
        ENDDO 
        GR(K)=GSUM      ! Save lag product sums
      END DO
      AVARPS(1:11)=GR(1:11)/GR(1)
! Solve for pooled ar coefficients using levinson-durbin recursion
      VP=GR(1) ; RNT=DBLE(xyr) ; AIC(1)=RNT*DLOG(VP)+2.D0
      AICM=AIC(1) ; IPR=0 ; TEST=.FALSE. 
      DO NN=2,11
        VP=GR(1)-SUM(AVARPS(2:NN)*GR(2:NN))
        DP=GR(NN+1)-SUM((/(AVARPS(I)*GR(NN-I+2),I=2,NN)/))
        AVARPS(NN+1)=DP/VP
        AIC(NN)=RNT*DLOG(VP)+DBLE(2*NN)
        IF(.NOT.TEST) THEN
          IF(IPR.GT.0.AND.AIC(NN).GT.AICM)THEN
            TEST=.TRUE.
          ELSEIF(AIC(NN).LE.AICM)THEN
            IPR=NN-1 ; DARPS(1:IPR)=AVARPS(2:IPR+1) ; AICM=AIC(NN)
          ENDIF
        ENDIF      
        WK(2:NN)=AVARPS(2:NN)-(/(AVARPS(NN+1)*AVARPS(I),I=NN,2,-1)/)
        AVARPS(2:NN)=WK(2:NN)
      ENDDO
      DARPS(12)=DBLE(IPR)
      RETURN
      END SUBROUTINE pool
!---------------------------------------------------------------------
      SUBROUTINE armod(ipr)   
      IMPLICIT NONE
      INTEGER,INTENT(IN)    :: ipr        ! Max AR lag years
      REAL(8),DIMENSION(12) :: cof        ! AR coefficients
      INTEGER               :: j,p,q,r
      DO j=1,nc
        r=yr(j) ; p=ad(j) ; q=p+r-1
        CALL memcof(ax(p:q),xok(p:q),r,ipr+1,cof,1)
        CALL albino(ax(p:q),xok(p:q),r,cof)
        ax(p:q)=ax(p:q)-SUM(ax(p:q))/DBLE(r)
        ax(p:q)=MAX(ax(p:q)+1.D0,0.D0)    ! Sets Means to 1.0 ??
      ENDDO
      RETURN
      END SUBROUTINE armod
!--------------------------------------------------------------------
      SUBROUTINE AR_model()  ! Creates RES (or ARS) chronologies 
      IMPLICIT NONE           
      REAL(8),DIMENSION(12) :: darps  ! Pooled AR Coefficients
      REAL(8),DIMENSION(12) :: phi    ! Stored AR Coefficients
      INTEGER               :: ipr,m  ! Max AR lag years
      IF (jrb.EQ.1) THEN
        m=ad(nc)+yr(nc)-1 ; ax(1:m)=dx(1:m)
      ELSE      
        CALL pool(darps,ipr) ; CALL armod(ipr) 
        IF     (krb.EQ.2) THEN 
          CALL robust_mean(ax)          ! Robust overall chronology  
        ELSE                   
          CALL arith_mean(ax)           ! Arith overall chronology
        ENDIF
        WHERE (.NOT.cok(1:xyr,mx))  xcrn(1:xyr,mx)=1.D0
        CALL memcof(xcrn(1:xyr,mx),cok(1:xyr,mx),xyr,ipr+1,phi,1)
        IF (NINT(phi(12)).GT.0) THEN    ! IP.gt.0
          CALL albino(xcrn(1:xyr,mx),cok(1:xyr,mx),xyr,phi)
          IF (jrb.EQ.3) &                 ! ARS Chronology
            CALL commie(xcrn(1:xyr,mx),cok(1:xyr,mx),xyr,darps) 
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE AR_model
!------------------------------------------------------------------------
! Robert Renka, Oak Ridge Natl. lab.
! This subroutine uses an order n*log(n) quick sort to sort a real 
! array x into increasing order.  The algorithm is as follows.  ind is
! initialized to the ordered sequence of indices 1,...,n, and all interchanges
! are applied to ind.  x is divided into two portions by picking a central
! element t.  The first and last elements are compared with t, and
! interchanges are applied as necessary so that the three values are in
! ascending order.  Interchanges are then applied so that all elements
! greater than t are in the upper portion of the array and all elements
! less than t are in the lower portion.  The upper and lower indices of one
! of the portions are saved in local arrays, and the process is repeated
! iteratively on the other portion.  When a portion is completely sorted,
! the process begins again by retrieving the indices bounding another
! unsorted portion.
! The ordering on x is defined by y(i) = x(ind(i)).
! note -- iu and il must be dimensioned >= log(n) where log has base 2.
! (ok for n up to about a billon)
      SUBROUTINE quicksort(n,x,ind)
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n      ! Count
      REAL(8),DIMENSION(N),INTENT(IN)  :: x(n)   ! Input array
      INTEGER,DIMENSION(N),INTENT(OUT) :: ind(n) ! Sorted index
      INTEGER,DIMENSION(21)            :: iu,il
      INTEGER   :: m,i=1,j,k,l,ij,it,itt,indx
      REAL(8)   :: r,t
! iu,il =  temporary storage for the upper and lower
!            indices of portions of the array x
! m =      index for iu and il
! i,j =    lower and upper indices of a portion of x
! k,l =    indices in the range i,...,j
! ij =     randomly chosen index between i and j
! it,itt = temporary storage for interchanges in ind
! indx =   temporary index for x
! r =      pseudo random number for generating ij
! t =      central element of x
      IF (n.LE.0) RETURN
      ind(1:n)=(/(i,i=1,n)/)  ! Initialize
      m=1 ; i=1 ; j=n ; r=0.375D0
! TOP OF LOOP
20    IF (i.GE.j) GO TO 70
      IF (r.LE.0.5898437D0) THEN
        r=r+0.0390625D0
      ELSE
        r=r-0.21875D0
      ENDIF
! Initialize k
30    k=i
! Select a central element of x and save it in t
      ij=i+INT(r*DBLE(j-i)) ; it=ind(ij) ; t=x(it)
! If first element of the array is greater than t, interchange with t
      indx=ind(i)
      IF (x(indx).GT.t) THEN
        ind(ij)=indx ; ind(i)=it ; it=indx ; t=x(it)
      ENDIF
! Initialize l
      l = j
! If last element of the array is less than t, interchange with t
      indx=ind(j)
      IF (x(indx).GE.t) GOTO 50
      ind(ij)=indx ; ind(j)=it ; it=indx ; t=x(it)
! If first element of the array is greater than t, interchange with t
      indx = ind(i)
      IF (x(indx).LE.t) GOTO 50
      ind(ij)=indx ; ind(i)=it ; it=indx ; t=x(it)
      GO TO 50
! Interchange elements k and l
40    itt=ind(l) ; ind(l)=ind(k) ; ind(k)=itt
! Find an element in the upper part of the array which is not larger than t
50    l=l-1
      indx=ind(l)
      IF (x(indx).GT.t) GOTO 50
! Find an element in the lower part of the array whcih is not smaller than t
60    k=k+1 ; indx=ind(k)
      IF (x(indx).LT.t) GOTO 60
! If k <= l, interchange elements k and l
      IF (k.LE.l) GOTO 40
! Save upper and lower subscripts of the portion of array yet to be sorted
      IF (l-i.GT.j-k) THEN
        il(m)=i ; iu(m)=l ; i=k ; m=m+1
        GOTO 80
      ENDIF
      il(m)=k ; iu(m)=j ; j=l ; m=m+1
      GOTO 80
! Begin again on another unsorted portion of the array
70    m=m-1
      IF (m.EQ.0) RETURN
      i=il(m) ; j=iu(m)
80    IF (j-i.GE.11) GOTO 30
      IF (i==1) GOTO 20
      i=i-1
! Sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.
90    i=i+1
      IF (i.EQ.j) GOTO 70
      indx=ind(i+1) ; t=x(indx)
      it=indx ; indx=ind(i)
      IF (x(indx).LE.t) GOTO 90
      k=i
100   ind(k+1)=ind(k) ; k=k-1 ; indx=ind(k)
      IF (t.LT.x(indx)) GOTO 100
      ind(k+1)=it ; GOTO 90
      END SUBROUTINE quicksort
!-----------------------------------------------------------------
      SUBROUTINE RandNorm(n,dat)  ! Return sorted random normal array
      IMPLICIT NONE  ! See Wiki - Box–Muller transform
      INTEGER,INTENT(IN)                 :: n     ! Count of numbers 
      REAL(8),DIMENSION(mxd),INTENT(OUT) :: dat   ! Output numbers 
      REAL(8)                            :: x1,x2,w
      INTEGER                            :: j,m,p
      INTEGER,DIMENSION(12)              :: seed
      seed=12848 ; p=1 ; m=n+6
      CALL RANDOM_SEED(put=seed) ! Initialise random sequence
      ID: DO
        JD: DO j=1,1000
          IF (m.GE.n+5) THEN
            CALL RANDOM_NUMBER(dat(p:n+5)) ; m=p  ! Random Uniform
          ENDIF
          x1=2.D0*dat(m)-1.D0 ; x2=2.D0*dat(m+1)-1.D0
          w=x1*x1+x2*x2 ; m=m+2
          IF (w.LT.1.D0) EXIT JD
        ENDDO JD
        w=SQRT((-2.D0*LOG(w))/w) ; dat(p)=x1*w ; dat(p+1)=x2*w ; p=p+2
        IF (p.GT.n) EXIT ID
      ENDDO ID
      RETURN
      END SUBROUTINE RandNorm
!--------------------------------------------------------------
      SUBROUTINE norm_distrib() ! Convert dx to normal dist
      IMPLICIT NONE             ! Uses ax
      INTEGER,DIMENSION(mxd) :: nx,ny   
      INTEGER                :: i=1,m
      m=ad(nc)+yr(nc)-1
      CALL quicksort(m,dx,nx)
      CALL randnorm(m,ax)
      CALL quicksort(m,ax,ny)
      DO i=1,m ; dx(nx(i))=ax(ny(i)) ; ENDDO
      IF     (krb.EQ.2) THEN 
        CALL robust_means()     ! Robust mean individuals  
        CALL robust_mean(dx)    ! Robust overall chronology  
      ELSE                   
        CALL arith_means()      ! Arith mean individuals 
        CALL arith_mean(dx)     ! Arith overall chronology
      ENDIF
      RETURN
      END SUBROUTINE norm_distrib
!------------------------------------------------------------------
      SUBROUTINE rcs_sdev()  ! SDEv and Error of RCS curves
      IMPLICIT NONE
      INTEGER :: i,j,k,p,q,u,v,rf
      mssd=0.D0
      DO i=1,nc
        j=trr(i)
        IF (j.GE.1) THEN                       ! May not use this tree
          IF (poo.EQ.1) THEN                   ! Use pith offset
            p=fy(i)-pth(i)+1
          ELSE                                 ! Do not use pith offset
            p=1
          ENDIF
          q=p+yr(i)-1 ; u=ad(i) ; v=u+yr(i)-1
          WHERE (xok(u:v)) 
            mssd(p:q,j)=mssd(p:q,j)+fx(u:v)**2   ! Ring sum by age
            mssd(p:q,mx)=mssd(p:q,mx)+fx(u:v)**2 ! Single Ring sum by age
          END WHERE
        ENDIF
      ENDDO                       
      DO k=1,cur+1                             ! Each curve and single
        rf=k ; IF (rf.GT.cur) rf=mx                 
        p=sfy(rf) ; q=sly(rf)
        WHERE (mok(p:q,rf))                    ! Arithmetic mean
          mssd(p:q,rf)=SQRT(MAX(mssd(p:q,rf)-(mval(p:q,rf)**2)* &
            DBLE(mcnt(p:q,rf)),0.001D0)/DBLE(MAX(mcnt(p:q,rf)-1,1)))  ! SDev 
          mserr(p:q,rf)=mssd(p:q,rf)/SQRT(DBLE(mcnt(p:q,rf))) ! Std Error 
        ELSEWHERE
          mssd(p:q,rf)=0.D0 ; mserr(p:q,rf)=0.D0     
        END WHERE
      ENDDO
      RETURN
      END SUBROUTINE rcs_sdev
!------------------------------------------------------------------
      SUBROUTINE rcs_dsdev()  ! Create arithmetic mean RCS curves
      IMPLICIT NONE
      INTEGER :: i,j,k,m,p,q,rf
      REAL(8) :: rad
      mserr=0.D0
      DO i=1,nc
        j=trr(i)
        IF (j.GE.1) THEN                       ! May not use this tree
          IF (poo.EQ.1) THEN                   ! Use pith offset
            rad=pthr(i)
          ELSE                                 ! Do not use pith offset
            rad=1.D0
          ENDIF
          DO k=ad(i),ad(i)+yr(i)-1             ! For each ring
            m=NINT(rad) ; rad=rad+x(k)         ! Radius from rings
            IF (xok(k)) THEN         
              mserr(m,j)=mserr(m,j)+fx(k)**2   ! Square sum by age
              mserr(m,mx)=mserr(m,mx)+fx(k)**2 ! Single square sum by age
            ENDIF
          ENDDO
        ENDIF
      ENDDO                       
      DO k=1,cur+1                             ! Each curve and single
        rf=k ; IF (rf.GT.cur) rf=mx                 
        p=dfy(rf) ; q=dly(rf)
        WHERE (dok(p:q,rf))                    ! Arithmetic mean
          mssd(p:q,rf)=SQRT(MAX(mserr(p:q,rf)-(dval(p:q,rf)**2)* &
            DBLE(dcnt(p:q,rf)),0.001D0)/DBLE(MAX(dcnt(p:q,rf)-1,1)))  ! SDev 
          mserr(p:q,rf)=mssd(p:q,rf)/SQRT(DBLE(dcnt(p:q,rf))) ! Std Error 
        ELSEWHERE
          mssd(p:q,rf)=0.D0 ; mserr(p:q,rf)=0.D0     
        END WHERE
      ENDDO
      RETURN
      END SUBROUTINE rcs_dsdev
!------------------------------------------------------------------
      SUBROUTINE crn_sdev() ! Overall Chronologies SDevs
      IMPLICIT NONE    
      INTEGER                :: i,u,v,p,q
      REAL(8),DIMENSION(mxy) :: wka,wkb
      IF (src.GT.1.AND.bfc.GT.1) THEN  ! Two+ chronologies
      wka=0.D0 ; wkb=0.D0
      DO i=1,nc         
        u=fy(i)-xfy+1 ; v=u+yr(i)-1 ; p=ad(i) ; q=p+yr(i)-1
        WHERE (xok(p:q))
          wka(u:v)=wka(u:v)+dx(p:q)    ! Indices sum
          wkb(u:v)=wkb(u:v)+dx(p:q)**2 ! Indices sum squares
        END WHERE
      ENDDO                 
      WHERE (xnum(1:xyr,mx).GT.3)      ! Sdev of each year with 4+ values
        xcsd(1:xyr,mx)=SQRT((wkb(1:xyr)- &
          (wka(1:xyr)**2)/DBLE(xnum(1:xyr,mx))) &
          /DBLE(xnum(1:xyr,mx)-1))     ! SDev 
      ELSEWHERE
        xcsd(1:xyr,mx)=0.D0            ! Missing SDevs
      END WHERE                    
      ENDIF
      RETURN
      END SUBROUTINE crn_sdev
!------------------------------------------------------------------
      SUBROUTINE robust_means() ! Robust mean each chronology
      IMPLICIT NONE
      REAL(8),DIMENSION(mxy) :: wka
      INTEGER                :: i,j,k,n,p,q 
      DO k=1,cur                        ! For each chronology
        DO i=xfy+xfa(k)-1,xfy+xla(k)-1  ! For each year of chronology
          n=0 ; p=i-xfy+1
          DO j=1,nc                     ! For each tree
            IF(trm(j).EQ.k.AND.i.GE.fy(j).AND.i.LE.ly(j)) THEN
              q=ad(j)+i-fy(j)
              IF (xok(q)) THEN
                n=n+1 ; wka(n)=dx(q)
              ENDIF
            ENDIF
          ENDDO
          IF    (n.EQ.0)THEN   ! No rings
            xcrn(p,k)=1.D0 ; xcsd(p,k)=0.D0 
          ELSEIF(n.EQ.1)THEN   ! One ring
            xcrn(p,k)=wka(1) ; xcsd(p,k)=0.D0
          ELSE
            CALL biwgt(wka,7.5D0,n,xcrn(p,k),xcsd(p,k))
          ENDIF
        ENDDO
!       xcrn(1:xyr,k)=MAX(xcrn(1:xyr,k),epsi)  ! Minimum chronology value 
      ENDDO
      RETURN
      END SUBROUTINE robust_means
!-----------------------------------------------------------------
      SUBROUTINE arith_means() ! Arithmetic mean each chronology
      IMPLICIT NONE    
      REAL(8),DIMENSION(mxy) :: wka
      INTEGER                :: i,j,p,q,u,v
      xcrn=0.D0 ; xcsd=0.D0
      DO i=1,nc         
        j=trm(i) ; u=fy(i)-xfy+1 ; v=u+yr(i)-1 ; p=ad(i) ; q=p+yr(i)-1
        WHERE (xok(p:q))
          xcrn(u:v,j)=xcrn(u:v,j)+dx(p:q)      ! Individual sums
          xcsd(u:v,j)=xcsd(u:v,j)+dx(p:q)**2   ! Individual squares
        END WHERE
      ENDDO
      DO i=1,cur        
        wka(1:xyr)=xcrn(1:xyr,i)
        WHERE (cok(1:xyr,i)) &   ! Mean of each year with values
          xcrn(1:xyr,i)=wka(1:xyr)/DBLE(xnum(1:xyr,i))
        WHERE (xnum(1:xyr,i).GT.3)
          xcsd(1:xyr,i)=SQRT(MAX(xcsd(1:xyr,i)-xcrn(1:xyr,i)* &
            wka(1:xyr),0.001D0)/DBLE(xnum(1:xyr,i)-1)) ! SDev 
        ELSEWHERE
          xcsd(1:xyr,i)=0.D0
        END WHERE              
!       xcrn(1:xyr,i)=MAX(xcrn(1:xyr,i),epsi)  ! Minimum chronology value 
      ENDDO
      RETURN
      END SUBROUTINE arith_means
!----------------------------------------------------------------
      SUBROUTINE robust_mean(dat) ! Robust mean full chronology
      IMPLICIT NONE
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat ! Data set to use
      REAL(8),DIMENSION(mxy)            :: wka
      INTEGER :: i,j,n,p,q 
      DO i=1,xyr               ! For each year 
        n=0 ; p=xfy+i-1        ! Calendar date
        DO j=1,nc              ! For each tree
          IF(p.GE.fy(j).AND.p.LE.ly(j)) THEN
            q=ad(j)+p-fy(j)
            IF (xok(q)) THEN
              n=n+1 ; wka(n)=dat(q)
            ENDIF
          ENDIF
        ENDDO
        IF (n.EQ.0)THEN   ! No rings
          xcrn(i,mx)=1.D0 ; xcsd(i,mx)=0.D0 
        ELSEIF(n.EQ.1)THEN   ! One ring
          xcrn(i,mx)=wka(1) ; xcsd(i,mx)=0.D0
        ELSE
          CALL biwgt(wka,7.5D0,n,xcrn(i,mx),xcsd(i,mx))
        ENDIF
      ENDDO
!     xcrn(1:xyr,mx)=MAX(xcrn(1:xyr,mx),epsi)  ! Minimum chronology value 
      RETURN
      END SUBROUTINE robust_mean
!-----------------------------------------------------------------
      SUBROUTINE arith_mean(dat) ! Arithmetic mean overall chronology
      IMPLICIT NONE    
      REAL(8),DIMENSION(mxd),INTENT(IN) :: dat ! Data set to use
      REAL(8),DIMENSION(mxy)            :: wka
      INTEGER :: i,p,q,u,v
      xcrn(1:mxy,mx)=0.D0 ; xcsd(1:mxy,mx)=0.D0
      DO i=1,nc         
        u=fy(i)-xfy+1 ; v=u+yr(i)-1 ; p=ad(i) ; q=p+yr(i)-1
        WHERE (xok(p:q))
          xcrn(u:v,mx)=xcrn(u:v,mx)+dat(p:q)    ! Overall sum
          xcsd(u:v,mx)=xcsd(u:v,mx)+dat(p:q)**2 ! Overall squares
        END WHERE
      ENDDO
      wka(1:xyr)=xcrn(1:xyr,mx)
      WHERE (cok(1:xyr,mx)) &   ! Mean of each year with values
          xcrn(1:xyr,mx)=wka(1:xyr)/DBLE(xnum(1:xyr,mx))
      WHERE (xnum(1:xyr,mx).GT.3)
        xcsd(1:xyr,mx)=SQRT(MAX(xcsd(1:xyr,mx)-xcrn(1:xyr,mx)* &
          wka(1:xyr),0.001D0)/DBLE(xnum(1:xyr,mx)-1)) ! SDev 
      ELSEWHERE
        xcsd(1:xyr,mx)=0.D0
      END WHERE
!     xcrn(1:xyr,mx)=MAX(xcrn(1:xyr,mx),epsi)  ! Minimum chronology value 
      RETURN
      END SUBROUTINE arith_mean
!----------------------------------------------------------------
      SUBROUTINE basal()  ! Convert to basal area increments
      IMPLICIT NONE       ! Creates tx from x
      INTEGER :: i,j
      REAL(8) :: rad1,rad2,area1,area2
      DO i=1,nc                          ! For each tree
        IF (poo.EQ.1) THEN               ! Use pith offset
          rad1=pthr(i) ; area1=pi*rad1**2
        ELSE                             ! Do not use pith offset
          rad1=0.D0 ; area1=0.D0
        ENDIF
        DO j=ad(i),ad(i)+yr(i)-1         ! For each ring
          rad2=rad1+x(j) ; area2=pi*rad2**2
          tx(j)=area2-area1 ; rad1=rad2 ; area1=area2
        ENDDO                            ! tx = basal area increments
      ENDDO
      RETURN
      END SUBROUTINE basal
!------------------------------------------------------------------
      SUBROUTINE det_crnfy()  ! Range of CRN and counts
      IMPLICIT NONE    
      INTEGER :: i,p,q,u,v
      xfy=MINVAL(fy(1:nc)) ; xly=MAXVAL(ly(1:nc))
      xyr=xly-xfy+1 ; xnum(1:mxy,mx)=0 
      DO i=1,nc         
        u=fy(i)-xfy+1 ; v=u+yr(i)-1 ; p=ad(i) ; q=ad(i)+yr(i)-1
        WHERE (xok(p:q)) xnum(u:v,mx)=xnum(u:v,mx)+1  ! Ring counts
      ENDDO            
      cyr(cf)=xyr ; cfy(cf)=xfy ; cly(cf)=xly
      num(1:xyr,cf)=xnum(1:xyr,mx)
      cok(1:mxy,mx)=xnum(1:mxy,mx).GE.1
      RETURN
      END SUBROUTINE det_crnfy
!-----------------------------------------------------------------          
      SUBROUTINE multi_alloc() ! Sort trees, allocate to curves
      IMPLICIT NONE            ! Set up ring counts and logicals      
      INTEGER  :: i,j,k,m,p,q,u,v,s,t
      REAL(8)  :: rad
      CALL det_crnfy                   
      IF (tst.LT.1.OR.tst.GT.4) tst=4  ! Must be none, age, diam or growth 
      CALL tree_sort()                 ! Sort the trees
      IF (src.EQ.1) THEN               ! Single RCS curve
        cur=1 ; trr(1:nc)=1 ; trm(1:nc)=1 
      ELSEIF (src.EQ.2) THEN           ! Multiple RCS curves
        cur=srcno
        m=nc/cur+1                     ! Trees in each RCS curve
        DO i=1,cur
          DO j=(i-1)*m+1,MIN(i*m,nc)  
            trr(tre(j))=i   ! Trees to RCS
            trm(tre(j))=i   ! RCS to trees
          ENDDO
        ENDDO
      ENDIF
      mcnt=0 ; ccnt=0 ; xnum(1:mxy,1:mx-1)=0 
      DO i=1,nc             
        IF (poo.EQ.1) THEN             ! Use pith offset
          s=fy(i)-pth(i)+1
        ELSE
          s=1
        ENDIF
        p=fy(i)-xfy+1 ; q=p+yr(i)-1 ; u=ad(i) ; v=u+yr(i)-1
        t=s+yr(i)-1 
        WHERE (xok(u:v)) xnum(p:q,trm(i))=xnum(p:q,trm(i))+1  ! CRN curve counts
        IF (trr(i).GE.1) THEN
          WHERE (xok(u:v))
            mcnt(s:t,trr(i))=mcnt(s:t,trr(i))+1  ! RCS curve counts
            mcnt(s:t,mx)=mcnt(s:t,mx)+1  ! RCS total counts
          END WHERE
        ENDIF
      ENDDO
      mok=mcnt.GE.1 ; cok=xnum.GE.1
      DO i=1,cur
        xfa(i)=MINVAL(fy(1:nc),MASK=LOGICAL(trm(1:nc).EQ.i,1))-xfy+1
        xla(i)=MAXVAL(ly(1:nc),MASK=LOGICAL(trm(1:nc).EQ.i,1))-xfy+1  
        xaa(i)=xla(i)-xfa(i)+1
        P1: DO p=1,mxt ; IF (mok(p,i)) EXIT P1 ; ENDDO P1
        Q1: DO q=mxt,1,-1 ; IF (mok(q,i)) EXIT Q1 ; ENDDO Q1
        sfy(i)=p ; sly(i)=q
        WHERE (cok(1:xyr,i)) ccnt(1:xyr)=ccnt(1:xyr)+1
      ENDDO
      sfy(mx)=MINVAL(sfy(1:cur)) ; sly(mx)=MAXVAL(sly(1:cur))
      IF (trc.EQ.2.OR.trc.EQ.3) THEN   ! Use diameters
        dcnt=0
        DO i=1,nc
          IF (poo.EQ.1) THEN           ! Use pith offset
            rad=pthr(i)
          ELSE                         ! Do not use pith offset
            rad=1.D0
          ENDIF
          IF (trr(i).GE.1) THEN
            DO k=ad(i),ad(i)+yr(i)-1     ! For each ring
              m=NINT(rad) ; rad=rad+x(k) ! Radius from ring
              IF (xok(k)) THEN         
                dcnt(m,trr(i))=dcnt(m,trr(i))+1 ! dRCS counts
                dcnt(m,mx)=dcnt(m,mx)+1  ! Single dRCS counts
              ENDIF
            ENDDO
          ENDIF
        ENDDO                       
        dok=dcnt.GE.1                  ! Rings present logical
        DO i=1,cur
          P2: DO p=1,mxt ; IF (dok(p,i)) EXIT P2 ; ENDDO P2
          Q2: DO q=mxt,1,-1 ; IF (dok(q,i)) EXIT Q2 ; ENDDO Q2
          dfy(i)=p ; dly(i)=q
        ENDDO
        dfy(mx)=MINVAL(dfy(1:cur)) ; dly(mx)=MAXVAL(dly(1:cur))
      ENDIF
      RETURN 
      END SUBROUTINE multi_alloc
!--------------------------------------------------------------------
      SUBROUTINE rcs_dsmooth(rf) ! Smooth an RCS curve
      IMPLICIT NONE ! Fit tail to rcns >= 4 then extend linearly 
      INTEGER,INTENT(IN)     :: rf       ! RCS curve number
      INTEGER                :: j,p,q,r  !,m
      REAL(8),DIMENSION(7)   :: eq
      p=dfy(rf) ; q=dly(rf) 
      JD: DO j=q,p,-1 ; IF (dcnt(j,rf).GE.4) EXIT JD ; ENDDO JD 
      r=j-p+1                            ! Smooth >= 4 only
      SELECT CASE (rdt)                  ! Select on the method 
      CASE (1)                           ! Age dependent smoothing
        ! Presumes 11 year minimum stiffness
        CALL spline3(q-p+1,dval(p:q,rf),dcnt(p:q,rf),10, &
          dsmo(p:q,rf),FA)  
      CASE (2)                           ! Unsmoothed RCS curve
        dsmo(p:q,rf)=dval(p:q,rf) 
      CASE (3)
        CALL curvet(r,dval(p:j,rf),dsmo(p:j,rf),eq) ! NEG EXP Curve
        IF (eq(2).LT.epsi.OR.eq(5).LT.epsi) THEN
          CALL out_err("Negative Exp Fit Failed")
          dsmo(p:q,rf)=dval(p:q,rf)      ! No smoothing default
        ELSE
          dsmo(j+1:q,rf)=dsmo(j,rf)
        ENDIF
      CASE (4)                           ! Linear regression
        CALL trend(r,dval(p:j,rf),dsmo(p:j,rf),eq(3),eq(5),eq(4))
        IF((NINT(eq(4)).GE.-5.AND.NINT(eq(4)).LE.r+5)) THEN
          dsmo(p:q,rf)=dval(p:q,rf)      ! No smoothing default
          CALL out_err("Trend line Fit Failed")
        ELSE
          dsmo(j+1:q,rf)=dsmo(j,rf)
        ENDIF
      CASE (5)                           ! Hugershoff Growth Curve
        CALL hughdi(dval(p:j,rf),r,dsmo(p:j,rf),eq,dok(p:j,rf),p-1)  
        dsmo(j+1:q,rf)=dsmo(j,rf)
      CASE (6)                           ! Fixed Spline
        CALL splinet(r,dval(p:j,rf),rdtno,dsmo(p:j,rf))
        WHERE (dsmo(p:j,rf).LT.0.01D0) dsmo(p:j,rf)=0.01D0 ! Minimum value
        dsmo(j+1:q,rf)=dsmo(j,rf)
      CASE (7)                           ! %n Spline
        CALL splinet(r,dval(p:j,rf),-(rdtno*r)/100,dsmo(p:j,rf)) 
        WHERE (dsmo(p:j,rf).LT.0.01D0) dsmo(p:j,rf)=0.01D0 ! Minimum value
        dsmo(j+1:q,rf)=dsmo(j,rf)
      ENDSELECT
      dsmo(q+1:mxt,rf)=dsmo(q,rf)        ! Linearly extend 
      RETURN
      END SUBROUTINE rcs_dsmooth
!------------------------------------------------------------------
      SUBROUTINE rcs_diam()  ! Create arithmetic mean RCS curves
      IMPLICIT NONE
      INTEGER :: i,j,k,m,p,q,rf
      REAL(8) :: rad
      dval=0.D0 
      DO i=1,nc
        j=trr(i)
        IF (j.GE.1) THEN                  ! May not use this tree
          IF (poo.EQ.1) THEN              ! Use pith offset
            rad=pthr(i)
          ELSE                            ! Do not use pith offset
            rad=1.D0
          ENDIF
          DO k=ad(i),ad(i)+yr(i)-1        ! For each ring
            m=NINT(rad) ; rad=rad+x(k)    ! Radius from ring
            IF (xok(k)) THEN         
              dval(m,j)=dval(m,j)+fx(k)   ! Ring sum by diameter
              dval(m,mx)=dval(m,mx)+fx(k) ! Single Ring sum by diameter
            ENDIF
          ENDDO
        ENDIF
      ENDDO                       
      DO k=1,cur+1                        ! Each curve and single
        rf=k ; IF (rf.GT.cur) rf=mx                 
        p=dfy(rf) ; q=dly(rf)
        WHERE (dok(p:q,rf))  &            ! Arithmetic mean
          dval(p:q,rf)=dval(p:q,rf)/DBLE(dcnt(p:q,rf))
        IF (.NOT.ALL(dok(p:q,rf))) THEN
          MD: DO m=p,q ; IF (.NOT.dok(m,rf)) EXIT MD ; ENDDO MD
          dval(p:m-1,rf)=dval(m,rf)  ! Infill missing at start
          DO m=m+1,q
            IF (.NOT.dok(m,rf)) dval(m,rf)=dval(m-1,rf)  ! Infill other gaps
          ENDDO
        ENDIF 
        CALL rcs_dsmooth(rf)
      ENDDO
      RETURN
      END SUBROUTINE rcs_diam
!------------------------------------------------------------------
      SUBROUTINE rcs_asmooth(rf)  ! Smooth an RCS curve
      IMPLICIT NONE ! Fit tail to rcns >= 4 then extend linearly 
      INTEGER,INTENT(IN)   :: rf   ! RCS curve number
      INTEGER              :: j,p,q,r 
      REAL(8),DIMENSION(7) :: eq
      p=sfy(rf) ; q=sly(rf) 
      JD: DO j=q,p,-1 ; IF (mcnt(j,rf).GE.4) EXIT JD ; ENDDO JD 
      r=j-p+1                            ! Smooth >= 4 only
      SELECT CASE (rdt)                  ! Select on the method 
      CASE (1)                           ! Age dependent smoothing
        ! Presumes 11 year minimum stiffness
        CALL spline3(q-p+1,mval(p:q,rf),mcnt(p:q,rf),10, &
          msmo(p:q,rf),itn.EQ.3.OR.q.GT.1500)  ! Age dep smoothing 
      CASE (2)                           ! Unsmoothed RCS curve
        msmo(p:q,rf)=mval(p:q,rf) 
      CASE (3)
        CALL curvet(r,mval(p:j,rf),msmo(p:j,rf),eq) ! NEG EXP Curve
        IF (eq(2).LT.epsi) THEN    ! .OR.eq(5).LT.epsi
          CALL out_err("Negative Exp Fit Failed")
          msmo(p:q,rf)=mval(p:q,rf)      ! No smoothing default
        ELSE
          msmo(j+1:q,rf)=msmo(j,rf)
        ENDIF
      CASE (4)                           ! Linear regression
        CALL trend(r,mval(p:j,rf),msmo(p:j,rf),eq(3),eq(5),eq(4))
        IF((NINT(eq(4)).GE.-5.AND.NINT(eq(4)).LE.r+5)) THEN
          msmo(p:q,rf)=mval(p:q,rf)      ! No smoothing default
          CALL out_err("Trend line Fit Failed")
        ELSE
          msmo(j+1:q,rf)=msmo(j,rf)
        ENDIF
      CASE (5)                           ! Hugershoff Growth Curve
        CALL hughdi(mval(p:j,rf),r,msmo(p:j,rf),eq,mok(p:j,rf),p-1)  
        msmo(j+1:q,rf)=msmo(j,rf)
      CASE (6)                           ! Fixed Spline
        CALL splinet(r,mval(p:j,rf),rdtno,msmo(p:j,rf))
        msmo(j+1:q,rf)=msmo(j,rf)
      CASE (7)                           ! %n Spline
        CALL splinet(r,mval(p:j,rf),-(rdtno*r)/100,msmo(p:j,rf)) 
        msmo(j+1:q,rf)=msmo(j,rf)
      ENDSELECT
      msmo(q+1:mxt,rf)=msmo(q,rf)        ! Linearly extend 
      WHERE (msmo(p:q,rf).LT.0.02D0) msmo(p:q,rf)=0.02D0
      RETURN   ! Apply minumum value of smallest measurable ring
      END SUBROUTINE rcs_asmooth
!------------------------------------------------------------------
      SUBROUTINE rcs_age()  ! Create arithmetic mean RCS curves
      IMPLICIT NONE
      INTEGER :: i,j,k,m,p,q,u,v,rf
      mval=0.D0
      DO i=1,nc
        j=trr(i)
        IF (j.GE.1) THEN                       ! May not use this tree
          IF (poo.EQ.1) THEN                   ! Use pith offset
            p=fy(i)-pth(i)+1
          ELSE                                 ! Do not use pith offset
            p=1
          ENDIF
          q=p+yr(i)-1 ; u=ad(i) ; v=u+yr(i)-1
          WHERE (xok(u:v)) 
            mval(p:q,j)=mval(p:q,j)+fx(u:v)    ! Ring sum by age
            mval(p:q,mx)=mval(p:q,mx)+fx(u:v)  ! Single Ring sum by age
          END WHERE
        ENDIF
      ENDDO                       
      DO k=1,cur+1                             ! Each curve and single
        rf=k ; IF (rf.GT.cur) rf=mx                 
        p=sfy(rf) ; q=sly(rf)
        WHERE (mok(p:q,rf)) &                  ! Arithmetic mean
          mval(p:q,rf)=mval(p:q,rf)/DBLE(mcnt(p:q,rf))
        IF (.NOT.ALL(mok(p:q,rf))) THEN
          MD: DO m=p,q ; IF (mok(m,rf)) EXIT MD ; ENDDO MD
          mval(p:m-1,rf)=mval(m,rf)  ! Infill missing start with 1st value
          DO m=m+1,q                 ! Infill other with previous value
            IF (.NOT.mok(m,rf)) mval(m,rf)=mval(m-1,rf) 
          ENDDO
        ENDIF 
        CALL rcs_asmooth(rf)
      ENDDO
      RETURN
      END SUBROUTINE rcs_age
!------------------------------------------------------------------
      SUBROUTINE single_mean() ! Calculate scaling factors 
      IMPLICIT NONE ! to convert means to that of single RCS
      INTEGER :: i,j,p,q,u,v
      REAL(8) :: rad
      DO i=1,nc                        ! Set detrending curve values
        IF (poo.EQ.1) THEN             ! Use pith offset
          rad=pthr(i) ; p=fy(i)-pth(i)+1
        ELSE                           ! Do not use pith offset
          rad=1.D0 ; p=1               
        ENDIF
        u=ad(i) ; v=u+yr(i)-1
        IF (trc.EQ.1) THEN  ! Age or basal area based curve
          q=p+yr(i)-1 ; cx(u:v)=msmo(p:q,mx)
        ELSEIF (trc.EQ.2) THEN         ! Diameter based curve
          DO j=u,v                     ! Diameter in 1mm steps
            u=NINT(rad) ; cx(j)=dsmo(u,mx) ; rad=rad+x(j)
          ENDDO
        ELSEIF (trc.EQ.3) THEN         ! Mean of age and diameter curves
          p=p-ad(i)
          DO j=u,v                     ! Diameter in 1mm steps
            u=NINT(rad) ; rad=rad+x(j)
            cx(j)=(dsmo(u,mx)+msmo(p+j,mx))/2.D0 
          ENDDO
        ENDIF                          ! Mean value of tree indices
        bfdx(i)=SUM(tx(u:v)/cx(u:v),MASK=xok(u:v))/ &
                DBLE(COUNT(xok(u:v)))
      ENDDO
      RETURN
      END SUBROUTINE single_mean
!------------------------------------------------------------------
      SUBROUTINE rcs_det() ! Fit RCS curve to each tree 
      IMPLICIT NONE
      INTEGER :: i,j,k,p,q,u,v
      REAL(8) :: rad
      DO i=1,nc                        ! Set detrending curve values
        IF (poo.EQ.1) THEN             ! Use pith offset
          rad=pthr(i) ; p=fy(i)-pth(i)+1
        ELSE                           ! Do not use pith offset
          rad=1.0D0 ; p=1               
        ENDIF
        k=trm(i) ; u=ad(i) ; v=u+yr(i)-1
        IF (trc.EQ.1) THEN             ! Age or basal area based curve
          q=p+yr(i)-1 ; cx(u:v)=msmo(p:q,k)
        ELSEIF (trc.EQ.2) THEN         ! Diameter based curve
          DO j=u,v                     ! Diameter in 1mm steps
            u=NINT(rad) ; cx(j)=dsmo(u,k) ; rad=rad+x(j)
          ENDDO
        ELSEIF (trc.EQ.3) THEN         ! Mean of age and diameter curves
          p=p-ad(i)
          DO j=u,v                     ! Diameter in 1mm steps
            u=NINT(rad) ; rad=rad+x(j)
            cx(j)=(dsmo(u,k)+msmo(p+j,k))/2.D0 
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE rcs_det
!------------------------------------------------------------------
      SUBROUTINE tree_ind() ! Tree indices, dx from tx and cx
      IMPLICIT NONE    
      INTEGER :: i,p,q
      REAL(8) :: rr,mn,mn1,sd,sd1
      DO i=1,nc                  
        p=ad(i) ; q=p+yr(i)-1
        dx(p:q)=tx(p:q)/cx(p:q)     ! Divide to get indices
        IF (ind.EQ.2) THEN          ! or Subtract & set mean and Sdev to
          rr=DBLE(COUNT(xok(p:q)))  ! mean and Sdev of ratios
          mn1=SUM(dx(p:q),MASK=xok(p:q))/rr
          sd1=SQRT(SUM((dx(p:q)-mn1)**2,MASK=xok(p:q))/(rr-1.D0))
          dx(p:q)=tx(p:q)-cx(p:q)
          mn=SUM(dx(p:q),MASK=xok(p:q))/rr
          sd=SQRT(SUM((dx(p:q)-mn)**2,MASK=xok(p:q))/(rr-1.D0))
          dx(p:q)=((dx(p:q)-mn)/sd)*sd1+mn1
        ENDIF
      ENDDO 
      RETURN
      END SUBROUTINE tree_ind
!------------------------------------------------------------------
      SUBROUTINE rcs_single() ! Set index means to single RCS means 
      IMPLICIT NONE
      INTEGER :: i,u,v
      DO i=1,nc               ! For each tree
        u=ad(i) ; v=u+yr(i)-1
        dx(u:v)=dx(u:v)*bfdx(i)*DBLE(COUNT(xok(u:v)))/ &
          SUM(dx(u:v),MASK=xok(u:v))
      ENDDO
      RETURN
      END SUBROUTINE rcs_single
!------------------------------------------------------------------
      SUBROUTINE rcs_crn()    ! Create RCS chronologies
      IMPLICIT NONE    
      INTEGER               :: i,p,q
      REAL(8),DIMENSION(mx) :: ra
      IF     (krb.EQ.2) THEN       ! Individual chronologies
        CALL robust_means()        ! Robust mean individuals  
      ELSE                   
        CALL arith_means()         ! Arith mean individuals 
      ENDIF
      IF (src.EQ.1) THEN           ! Only one chronology
        xcrn(1:xyr,mx)=xcrn(1:xyr,1)
        xcsd(1:xyr,mx)=xcsd(1:xyr,1)
      ELSE                         ! More than one chronology
        SELECT CASE (bfc)
        CASE (1)                   ! Average of all trees
          IF     (krb.EQ.2) THEN 
            CALL robust_mean(dx)   ! Robust overall chronology  
          ELSE                   
            CALL arith_mean(dx)    ! Arith overall chronology
          ENDIF
        CASE (2)                   ! Average individual chronologies
          xcrn(1:xyr,mx)=0.D0      
          DO i=1,cur
            WHERE (cok(1:xyr,i)) &
              xcrn(1:xyr,mx)=xcrn(1:xyr,mx)+xcrn(1:xyr,i)
          ENDDO
          WHERE (cok(1:xyr,mx)) 
            xcrn(1:xyr,mx)=xcrn(1:xyr,mx)/DBLE(ccnt(1:xyr))
          ELSEWHERE
            xcrn(1:xyr,mx)=1.D0
          END WHERE 
        CASE (3)                   ! Means = 1.D0
          xcrn(1:xyr,mx)=0.D0
          DO i=1,cur         
            ra(i)=SUM(xcrn(1:xyr,i),MASK=cok(1:xyr,i))/ &
              DBLE(COUNT(cok(1:xyr,i)))
            xcrn(1:xyr,i)=xcrn(1:xyr,i)/ra(i) 
            WHERE (cok(1:xyr,i)) &
              xcrn(1:xyr,mx)=xcrn(1:xyr,mx)+xcrn(1:xyr,i)  
          ENDDO
          WHERE (cok(1:xyr,mx)) 
            xcrn(1:xyr,mx)=xcrn(1:xyr,mx)/DBLE(ccnt(1:xyr))  
          ELSEWHERE
            xcrn(1:xyr,mx)=1.D0
          END WHERE 
          DO i=1,nc
            p=ad(i) ; q=p+yr(i)-1 ; dx(p:q)=dx(p:q)*ra(trm(i))
          ENDDO
        ENDSELECT
      ENDIF
      RETURN
      END SUBROUTINE rcs_crn
!------------------------------------------------------------------
      SUBROUTINE det_sf_rcs()  ! Signal-free RCS detrend 
      IMPLICIT NONE           
      INTEGER :: m
      REAL(8) :: mn
      IF (trc.EQ.2.OR.trc.EQ.3) CALL rcs_diam()    ! RCS by diameter
      IF (trc.NE.2) CALL rcs_age()                 ! RCS by age 
      IF (gtr.EQ.2) CALL single_mean()             ! Rescaling factors
      CALL rcs_det()                               ! Fit RCS to trees
      CALL tree_ind()                              ! Tree indices
      IF (src.GT.1.AND.gtr.EQ.2) CALL rcs_single() ! Mean single
      m=ad(nc)+yr(nc)-1     ! Scale dx to produce CRN mean of 1.0
      mn=SUM(tx(1:m)/cx(1:m),mask=xok(1:m))/DBLE(COUNT(xok(1:m)))
      dx(1:m)=dx(1:m)/mn    ! Needed to prevent drift 
      CALL rcs_crn()                            ! Create chronology
      RETURN
      END SUBROUTINE det_sf_rcs
!------------------------------------------------------------------------
      SUBROUTINE rcs_detrend()     ! Performs CRCS detrend
      IMPLICIT NONE           
      LOGICAL,DIMENSION(mxy) :: qok        ! True if > 1 tree in year
      REAL(8),DIMENSION(mxy) :: dat 
      INTEGER                :: i,k,p,q,r,u,v,sfiter
      IF (journal) THEN  ! Save detrend method for each tree
        OPEN(23,FILE="Journal.prn",IOSTAT=ios,STATUS="OLD",POSITION="APPEND")
        IF (io_err("Open App","Journal.prn")) STOP
      ENDIF
      CALL multi_alloc()           ! Allocate trees to RCS curves
      qok(1:xyr)=xnum(1:xyr,mx).GT.1  ! > 1 tree
      sfiter=1 ; IF (sfo.EQ.2) sfiter=sfiter+sfono
      KDO: DO k=1,sfiter           ! Signal-free iterations
        IF (k.GT.1) THEN           ! Subsequent iterations
          DO i=1,nc                ! For each tree
            p=ad(i) ; q=p+yr(i)-1 ; u=fy(i)-xfy+1 ; v=u+yr(i)-1
            WHERE (xok(p:q).AND.qok(u:v).AND.xcrn(u:v,mx).GE.0.01D0) 
              fx(p:q)=tx(p:q)/xcrn(u:v,mx) ! Measurements divided by CRN if > 0.01
            ELSEWHERE 
              fx(p:q)=tx(p:q)      ! Measurements
            END WHERE 
          ENDDO
        ENDIF
        CALL det_sf_rcs()          ! Create chronology
        IF (journal) THEN  ! Save Signal-free convergence details
          WRITE(23,'("Det=",I3,"  Iter",I3,"  Diff ",F10.5)') idt,k, &
            SUM(ABS(xcrn(1:xyr,mx)-dat(1:xyr)))/DBLE(xyr) ! Mean difference
        ENDIF
        IF (k.NE.1.AND. &
          ALL(ABS(xcrn(1:xyr,mx)-dat(1:xyr)).LT.0.001D0)) EXIT KDO 
        dat(1:xyr)=xcrn(1:xyr,mx)  ! Store chronology
      ENDDO KDO
      IF (trc.EQ.2) THEN           ! RCS curve SDev by diameter
        CALL rcs_dsdev()              
      ELSE                         ! RCS curve SDev by age
        CALL rcs_sdev()         
      ENDIF
      DO i=1,cur ! srcno
        p=xfa(i) ; q=xla(i) ; r=xaa(i)
        CALL spline_miss(r,xcrn(p:q,i),CDsp,xcsm(p:q,i),cok(p:q,i))
      ENDDO
      CALL spline_miss(xyr,xcrn(1:xyr,mx),CDsp,xcsm(1:xyr,mx),cok(1:xyr,mx))
      IF (journal) CLOSE(23)  ! Close journal file
      RETURN
      END SUBROUTINE rcs_detrend
!------------------------------------------------------------------
      SUBROUTINE det_sf_cur()  ! Signal-free curve fitting 
      IMPLICIT NONE           
      REAL(8),DIMENSION(7) :: eq
      REAL(8)              :: sl,yi,xi,ah,bh
      INTEGER              :: i,j,m,p,q,r
      DO i=1,nc
        p=ad(i) ; r=yr(i) ; q=p+r-1 
        SELECT CASE (jdt(i))                  ! Select on the detrend method 
        CASE (:-10)                           ! Percentage of length spline
          m=INT(DBLE(-jdt(i))*DBLE(r)/100.D0) ! %n stiffness
          CALL splinet(r,fx(p:q),m,cx(p:q))   !  jdt(i)=9
        CASE (0) ;  cx(p:q)=1.D0              ! Mean RAW CRN
        CASE (1)                              ! Exp or any sloping line
          CALL curvet(r,fx(p:q),cx(p:q),eq)   ! Neg Exp Curve
          IF(eq(2).LE.epsi.OR.eq(5).LT.epsi)THEN
            CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi) 
            IF(cx(p).GE.epsi.AND.cx(q).GE.epsi)THEN
              jdt(i)=3
            ELSE                  ! Horizontal if line approaches zero
              jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)  
            ENDIF
          ENDIF
        CASE (2)    ! Exp or negative sloping line
          CALL curvet(r,fx(p:q),cx(p:q),eq)  ! Neg Exp Curve
          IF(eq(2).LE.epsi.OR.eq(5).LT.epsi)THEN
            CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi) 
            IF(sl.LE.epsi.AND.(cx(p).GE.epsi.AND.cx(q).GE.epsi))THEN
              jdt(i)=4
            ELSE          ! Horizontal if line approaches zero
              jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)  
            ENDIF
          ENDIF
        CASE (3)                            ! Any straight line
          CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi)  ! Best fit line
          IF(cx(p).LT.epsi.OR.cx(q).LT.epsi)THEN
            jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)
          ENDIF
        CASE (4)                            ! Negatively sloping line
          CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi) ! Best fit line
          IF(sl.GE.epsi.OR.cx(p).LT.epsi.OR.cx(q).LT.epsi)THEN
            jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)
          ENDIF
        CASE (5)                            ! Mean (horizontal) line
          cx(p:q)=SUM(fx(p:q))/DBLE(r)      ! Arithmetic mean
        CASE (6)    ! Hugershoff curve
          j=fy(i)-pth(i) ; IF (poo.EQ.2) j=1  ! Use/not pith offset
          CALL hughdi(fx(p:q),r,cx(p:q),eq,xok(p:q),j)  !HUGERSHOFF CURVE
          IF(eq(2).LE.epsi)THEN   ! Hug failed try Mod Neg Exp
            CALL curvet(r,fx(p:q),cx(p:q),eq)
            IF(eq(2).LE.epsi.OR.eq(5).LT.epsi)THEN
              CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi)
              IF(sl.LE.epsi.AND.(cx(p).GE.epsi.AND.cx(q).GE.epsi))THEN
                jdt(i)=4
              ELSE                ! Horizontal if line < zero
                jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)  
              ENDIF
            ELSE
              jdt(i)=2
            ENDIF
          ENDIF
        CASE (7)    ! General exponential curve
          CALL gexp(fx(p:q),r,cx(p:q),ah,bh,xok(p:q)) 
          IF(bh.LE.epsi)THEN
            jdt(i)=3 ; CALL trend(r,fx(p:q),cx(p:q),sl,yi,xi)
            IF(sl.GE.epsi.OR.cx(p).LT.epsi.OR.cx(q).LT.epsi) THEN
              jdt(i)=5 ; cx(p:q)=SUM(fx(p:q))/DBLE(r)
            ENDIF
          ENDIF
        CASE (9:)  !  ; CALL det_spline()     ! Fixed length spline
          CALL splinet(r,fx(p:q),jdt(i),cx(p:q))   ! jdt=8
        END SELECT
        cx(p:q)=MAX(cx(p:q),0.02D0)         ! Minimum curve value i.e 1 ring
!       cx(p:q)=MAX(tx(p:q)/4.D0,cx(p:q))   ! Max index value = 4.0 to prevent   
      ENDDO                                 ! excessive signal-free values
      CALL tree_ind()                       ! Tree indices
      IF     (krb.EQ.2) THEN 
        CALL robust_mean(dx)                ! Robust overall chronology  
      ELSE                   
        CALL arith_mean(dx)                 ! Arith overall chronology
      ENDIF
      RETURN
      END SUBROUTINE det_sf_cur
!------------------------------------------------------------------------
      SUBROUTINE cur_detrend() ! Curve-fitting detrend
      IMPLICIT NONE           
      INTEGER,PARAMETER             :: iter=40    ! Maximum signal free iterations 
      LOGICAL,DIMENSION(mxy)        :: qok        ! True if > 1 tree in year
      INTEGER,DIMENSION(mxy)        :: ssy        ! Variable spline stiffness
      REAL(8),DIMENSION(mxy)        :: dat,wk 
      REAL(8)                       :: rr,mn
      INTEGER,DIMENSION(mxs,0:iter) :: jnum
      REAL(8),DIMENSION(0:iter)     :: diff
      LOGICAL                       :: chok       ! Curve type changed
      INTEGER                       :: i,k,p,q,u,v,w,sfi
      IF (journal) THEN  ! Save detrend method for each tree
        OPEN(23,FILE="Journal.prn",IOSTAT=ios,STATUS="OLD",POSITION="APPEND")
        IF (io_err("Open App","Journal.prn")) STOP
      ENDIF
      CALL det_crnfy()                 ! Tree counts by year
      trr(1:nc)=1 ; trm(1:nc)=1        ! Only one chronology 
      CALL tree_sort()                 ! Sort the trees
      qok(1:xyr)=xnum(1:xyr,mx).GT.1   ! > 1 tree
      rr=DBLE(COUNT(qok(1:xyr)))       ! > 1 tree count
      dat(1:xyr)=1.D0 ; ssy=0                
      DO i=1,nc                                   
        w=yr(i) ; u=fy(i)-xfy+1 ; v=u+w-1 ; ssy(u:v)=ssy(u:v)+w
      ENDDO
      WHERE (xnum(1:xyr,mx).GE.1)  
        ssy(1:xyr)=ssy(1:xyr)/xnum(1:xyr,mx)  ! Mean tree age of CRN
      ELSE WHERE
        ssy(1:xyr)=MINVAL(yr(1:nc))  ! Min value for missing years
      END WHERE 
      sfi=0 ; IF (sfo.GT.1) sfi=MIN(sfono,iter)     
      jdt=idt ; diff=0.D0         ! Requested detrend 
      KD: DO k=0,sfi                           
        IF (k.GT.0) THEN                   ! Signal-free iterations
          DO i=1,nc                        ! For each tree
            p=ad(i) ; q=p+yr(i)-1 ; u=fy(i)-xfy+1 ; v=u+yr(i)-1
            WHERE (qok(u:v).AND.dat(u:v).GE.0.01D0) 
              fx(p:q)=tx(p:q)/dat(u:v)  ! Measurements divided by CRN if > 0.0
            ELSEWHERE 
              fx(p:q)=tx(p:q)           ! Measurements
            END WHERE 
            IF (COUNT(qok(u:v)).GT.1) THEN
              mn=SUM(tx(p:q),MASK=qok(u:v))/SUM(fx(p:q),MASK=qok(u:v))
              WHERE (qok(u:v)) fx(p:q)=fx(p:q)*mn  ! Reset mean
            ENDIF
          ENDDO
          WHERE (jdt(1:nc).NE.-66) jdt(1:nc)=idt ! Original curve option   
        ENDIF
        CALL det_sf_cur()                  ! Full chronology in crn
        jnum(1:nc,k)=jdt(1:nc)             ! Store detrend method used
        IF (k.GT.0) THEN
          chok=FA                          ! Presume no changes   
          DO i=1,nc          
            IF (jnum(i,k).NE.jnum(i,k-1)) THEN
              IF (k.GE.5) THEN             ! Changes in 5th+ iteration
                jdt(i)=-66 ; chok=TR
              ENDIF
            ENDIF
          ENDDO
          IF (chok) THEN
            CALL det_sf_cur()      ! Redo CRN if changes
            jnum(1:nc,k)=jdt(1:nc) ! Detrend method used
          ENDIF
        ENDIF
        dat(1:xyr)=xcrn(1:xyr,mx)-dat(1:xyr)  ! Difference this iteration
        CALL splinec(xyr,dat(1:xyr),ssy(1:xyr),wk(1:xyr))   
        dat(1:xyr)=dat(1:xyr)-wk(1:xyr)       ! High-pass difference 
        diff(k)=SUM(ABS(dat(1:xyr)),MASK=qok(1:xyr))/rr ! Mean difference
        IF (journal) THEN  ! Save signal-free convergence details
          WRITE(23,'("Det=",I3," Iter",I3," Diff",F10.5)') idt,k,diff(k)
        ENDIF 
        IF (diff(k).LT.0.0001D0) EXIT KD
        dat(1:xyr)=xcrn(1:xyr,mx)                  
      ENDDO KD
      IF (journal) THEN  ! Save detrend method for each tree
        WRITE(23,*) 
        DO i=1,nc ; WRITE(23,'(A8,I5)') nam(i)(1:8),jdt(i) ; ENDDO
        WRITE(23,*) ; CLOSE(23)
      ENDIF
      RETURN
      END SUBROUTINE cur_detrend
!------------------------------------------------------------------------
      SUBROUTINE detrend()     ! Performs detrend
      IMPLICIT NONE           
      INTEGER :: m,p,q,r
      m=ad(nc)+yr(nc)-1 
      tx(1:m)=x(1:m)               ! Transformed measures
      IF (itn.EQ.3) CALL basal()   ! Convert to basal area
      IF (itn.EQ.2) CALL trnfrm()  ! Power Transform
      fx(1:m)=tx(1:m)              ! Signal-free default
      IF (idt.EQ.-2) THEN
        CALL rcs_detrend()
        jdt(1:nc)=-2
      ELSE
        CALL cur_detrend()
      ENDIF
      CALL crn_sdev()                    ! Main chronology SDev
      IF (idb.EQ.2) CALL norm_distrib()  ! Normal distribution
      CALL AR_model()                    ! STD or ARS or RES chronology 
      IF (isb.EQ.2) CALL stabVar()       ! Stabilise variance 
      IF (isb.EQ.3) CALL stabVar2(30)    ! Stabilise high frequency variance 
      r=xyr ; crn(1:r,cf)=xcrn(1:r,mx)   ! Save final chronology
      okc(1:r,cf)=cok(1:r,mx)            ! Save logicals
      IF (idt.EQ.-2) THEN                ! Smooth chronologies
        CALL spline_miss(r,xcrn(1:r,mx),CDsp,xcsm(1:r,mx),cok(1:r,mx))
        MD: DO m=1,cur                   
          p=xfa(m) ; q=xla(m) ; r=xaa(m)
          CALL spline_miss(r,xcrn(p:q,m),CDsp,xcsm(p:q,m),cok(p:q,m))
          IF (src.EQ.1) EXIT MD
        ENDDO MD
      ENDIF
      RETURN
      END SUBROUTINE detrend
!------------------------------------------------------------------------
      SUBROUTINE par_valid(tnam)   ! Check the detrend parameters
      IMPLICIT NONE           
      CHARACTER(60),INTENT(IN) :: tnam
      IF (idt.LT.-999.OR.idt.GT.9999.OR.idt.EQ.-1.OR. &
        idt.EQ.8.OR.idt.EQ.9.OR.(idt.GE.-9.AND.idt.LE.-3)) THEN
        CALL out_err(TRIM(tnam)//" Invalid IDT, NEXP used")
        idt=2
      ENDIF
      IF (itn.LT.1.OR.itn.GT.3) THEN
        CALL out_err(TRIM(tnam)//" Invalid ITN, No transform used")
        itn=1
      ENDIF
      IF (idt.EQ.-2) THEN
        IF (rdt.LT.1.OR.rdt.GT.7) THEN
          CALL out_err(TRIM(tnam)//" Invalid RDT, Age Dependant used")
          rdt=1
        ENDIF
      ENDIF
      IF (ind.LT.1.OR.ind.GT.2) THEN
        CALL out_err(TRIM(tnam)//" Invalid IND, Ratio used")
        ind=1
      ENDIF
      IF (krb.LT.1.OR.krb.GT.2) THEN
        CALL out_err(TRIM(tnam)//" Invalid KRB, Arith mean used")
        krb=1
      ENDIF
      IF (isb.LT.1.OR.isb.GT.3) THEN
        CALL out_err(TRIM(tnam)//" Invalid ISB, Variance stabilised")
        isb=2
      ENDIF
      IF (sfo.LT.1.OR.sfo.GT.2) THEN
        CALL out_err(TRIM(tnam)//" Invalid SFO, Sig Free used")
        sfo=2
      ENDIF
      IF (poo.LT.1.OR.poo.GT.2) THEN
        CALL out_err(TRIM(tnam)//" Invalid POO, Pith offset used")
        poo=1
      ENDIF
      IF (idt.EQ.-2.AND.(src.LT.1.OR.src.GT.3)) THEN
        CALL out_err(TRIM(tnam)//" Invalid SRC, Single RCS used")
        src=1
      ENDIF
      IF (idt.EQ.-2.AND.trc.LT.1.OR.trc.GT.3) THEN
        CALL out_err(TRIM(tnam)//" Invalid TRC")
        trc=1
      ENDIF
      IF (gtr.LT.1.OR.gtr.GT.2.OR. &           ! Error as mean single needs 
        (gtr.EQ.2.AND.(idt.NE.-2.OR.src.EQ.1 & ! multi RCS, ratios
          .OR.itn.NE.2.OR.ind.EQ.1))) THEN     ! and power transform off
        gtr=1 ; CALL out_err(TRIM(tnam)//" Invalid GTR, No transform")
      ENDIF 
      IF (tst.LT.1.AND.tst.GT.7) THEN
        CALL out_err(TRIM(tnam)//" Invalid TST, Sort by growth rate")
        tst=4
      ENDIF
      IF (idt.NE.-2.OR.src.EQ.1) THEN
        bfc=1
      ELSEIF (bfc.LT.1.OR.bfc.GT.3) THEN
        CALL out_err(TRIM(tnam)//" Invalid BFC, Mean of trees used")
        bfc=1
      ENDIF
      IF (idb.LT.1.OR.idb.GT.2) THEN
        CALL out_err(TRIM(tnam)// &
          " Invalid IDB, Distribution unchanged")
        idb=1
      ENDIF
      IF (jrb.LT.1.OR.jrb.GT.3) THEN
        CALL out_err(TRIM(tnam)//" Invalid IDB, Use STD chronology")
        jrb=1
      ENDIF
      IF (srcno.LT.1.OR.srcno.GT.mx-1) sfono=1
      IF (sfono.LT.2.OR.sfono.GT.100) sfono=10
      IF     (rdt.EQ.6) THEN
        IF (rdtno.LT.5.OR.rdtno.GT.999) rdtno=60
      ELSEIF (rdt.EQ.7) THEN
        IF (rdtno.LT.-999.OR.rdtno.GT.-5) rdtno=-10
      ENDIF 
      RETURN
      END SUBROUTINE par_valid
!-------------------------------------------------------------------- 
      SUBROUTINE prog_help(hno,hm,bf,bl)
      IMPLICIT NONE   ! Displays help text from Help.prn
      INTEGER,INTENT(IN)               :: hno     ! Number of messages
      INTEGER,DIMENSION(60),INTENT(IN) :: hm      ! List of messages
      INTEGER,INTENT(IN)               :: bf,bl   ! First and last buttons
      INTEGER                          :: menu    ! Menu number
      INTEGER                          :: butno   ! Button number
      INTEGER                          :: messno  ! Message number
      INTEGER                          :: lines   ! Number of lines of text
      INTEGER                          :: i,j
      CHARACTER(80)                    :: htext
      INTEGER                          :: ra,rb,rc,rd ! Help box address
      LOGICAL                          :: okb
      ra=80 ; rb=1250 ; rc=1800 ; rd=2080 ; mous=0 ; okb=TR
      BDO: DO i=bf,bl     ! Is mouse in a button box
        IF (b(i)%ok.AND.msx.GE.b(i)%x1.AND.msx.LE.b(i)%x2 & 
           .AND.msy.GE.b(i)%y1.AND.msy.LE.b(i)%y2) THEN
           mous=i ; EXIT BDO
        ENDIF
      ENDDO BDO
      IF (mous.EQ.0) THEN ! Not button so is it a message
        okb=FA
        CDO: DO j=1,hno        
          i=hm(j)
          IF (msx.GE.madd(i)%x1.AND.msx.LE.madd(i)%x2 & 
            .AND.msy.GE.madd(i)%y1.AND.msy.LE.madd(i)%y2) THEN
            mous=i ; EXIT CDO
         ENDIF
        ENDDO CDO
      ENDIF
      IF (mous.NE.0) THEN
        OPEN(29,FILE="help.prn",IOSTAT=ios,STATUS="OLD")
        IF (io_err("Open","help.prn")) STOP
        READ(29,*,IOSTAT=ios)
        IF (io_err("Read Head","help.prn")) STOP
        READ(29,*,IOSTAT=ios)
        IF (io_err("Read Head","help.prn")) STOP ! Headers
        ID: DO i=1,500           ! Look for message or mouse number
          READ(29,'(I2,3I4,A80)',IOSTAT=ios) menu,butno,messno,lines,htext
          IF (io_err("Read","help.prn")) STOP
          IF ((okb.AND.mous.EQ.butno).OR. &
            ((.NOT.okb).AND.mous.EQ.messno)) EXIT ID  ! Message found
          IF (menu.EQ.0) EXIT ID           ! No message found
          DO j=1,lines
            READ(29,*,IOSTAT=ios)
            IF (io_err("Read Lines","help.prn")) STOP
          ENDDO
        ENDDO ID
        CALL SETCLR(silver)
        CALL AREAF((/ra,rb,rb,ra/),(/rc,rc,rd,rd/),4)
        CALL SHDPAT(0)        ! Outlined Area
        CALL SETCLR(red)              
        CALL AREAF((/ra,rb,rb,ra/),(/rc,rc,rd,rd/),4)
        CALL SHDPAT(16)        ! Filled Area
        CALL SETCLR(black) ; rc=rc+22
        CALL MESSAG(htext,ra+8,rc)  ! Fixed label
        CALL HEIGHT(20)               ! Character size
        DO j=1,lines
          rc=rc+50 ; READ(29,'(A80)') htext
          CALL MESSAG(htext,ra+8,rc)  ! Fixed label
        ENDDO
        CALL HEIGHT(25)               ! Character size
        CLOSE(29) 
      ENDIF
      RETURN
      END SUBROUTINE prog_help
!--------------------------------------------------------------
      END MODULE Stand                                        
