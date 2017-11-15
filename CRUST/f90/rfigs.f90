! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see 
! the GNU General Public License.
      MODULE  rfigs    
      USE crustutil
      IMPLICIT NONE
      CHARACTER(30),DIMENSION(fin) :: rfnam
      CONTAINS   
!--------------------------------------------------------
      SUBROUTINE rfig_val()  ! Figures displays
      IMPLICIT NONE                 
      rfnam(1:30)= &
        (/"(1) Exit Figures Menu         ", &
          "(2) Save this plot            ", &
          "(3) Create CRNs (crns.fil)    ", &
          "(4) Mean tree   (mean.fil)    ", &
          "(5) Raw Data Stats (ars.fil)  ", &
          "(6) Missing Ring Rep (ars.fil)", &
          "(7) Save Missing rings (cf)   ", &
          "(8)                           ", &
          "(9) Trees to CRN format (cf)  ", &
          "(10) CRN to cols (crns.crn)   ", &
          "(11)                          ", &
          "(12)                          ", &
          "(13) All trees fy (ars.fil)   ", &
          "(14) All trees name (ars.fil) ", &
          "(15) Trees by start (ars.fil) ", &
          "(16) Trees by end (ars.fil)   ", &
          "(17)                          ", &
          "(18) Pith Offset Report (cf)  ", &
          "(19)                          ", &
          "(20) ARSTAN pith (ars.fil)    ", &
          "(21)                          ", &
          "(22) Join All Sites (ars.fil) ", &
          "(23)                          ", &
          "(24)                          ", &
          "(25)                          ", &
          "(26) EPS Rep 2 RCS (ars.fil)  ", &
          "(27) EPS Rep 1 RCS (ars.fil)  ", &
          "(28)                          ", &
          "(29) MXD Adjust (ars.fil)     ", &
          "(30) TRW Adjust (ars.fil)     "/)
      rfnam(31:60)= &
        (/"(31)                          ", &
          "(32) Boot CRU (cf)            ", &
          "(33)                          ", &
          "(34)                          ", &
          "(35) Multi Site 1 RCS (ars)   ", &
          "(36) Multi Site 4 RCS (ars)   ", &
          "(37)                          ", &
          "(38)                          ", &
          "(39)                          ", &
          "(40)                          ", &
          "(41)                          ", &
          "(42)                          ", &
          "(43)                          ", &
          "(44)                          ", & 
          "(45)                          ", &
          "(46)                          ", &
          "(47)                          ", &
          "(48)                          ", &
          "(49)                          ", &
          "(50)                          ", &
          "(51)                          ", &
          "(52)                          ", &
          "(53)                          ", &
          "(54)                          ", &
          "(55)                          ", &
          "(56)                          ", &
          "(57)                          ", &
          "(58)                          ", &
          "(59)                          ", &
          "(60)                          "/)
      rfnam(61:90)= &
        (/"(61)                          ", &
          "(62)                          ", &
          "(63)                          ", &
          "(64)                          ", &
          "(65)                          ", &
          "(66)                          ", &
          "(67)                          ", &
          "(68)                          ", &
          "(69)                          ", &
          "(70)                          ", &
          "(71)                          ", &
          "(72)                          ", &
          "(73)                          ", & 
          "(74)                          ", &
          "(75)                          ", &
          "(76)                          ", &
          "(77)                          ", &
          "(78)                          ", &
          "(79)                          ", &
          "(80)                          ", &
          "(81)                          ", &
          "(82)                          ", &
          "(83)                          ", &
          "(84)                          ", &
          "(85)                          ", &
          "(86)                          ", &
          "(87)                          ", &
          "(88)                          ", &
          "(89)                          ", &
          "(90)                          "/)
      RETURN 
      END SUBROUTINE rfig_val
!----------------------------------------------------------
      SUBROUTINE rf_plot(ref1,plot)  
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: ref1   ! Chosen action
      INTEGER,INTENT(IN) :: plot   ! Plot number
      INTEGER            :: i
      CALL ERASE()
      SELECT CASE (ref1)
      CASE (2)    ! Save plot as .ps file
        figm="rfigs/rfig0"
        CALL open_ps(plot,11)
        SELECT CASE (plot)
        CASE (15:16) ; CALL KKtreed()
        CASE (26:27) ; CALL EPS_prep2d()
        CASE (29:30) ; CALL Adjust_MXDd()
        CASE (32)    ; CALL boot_crud()
        CASE (35) ; CALL mrcs_sited(1)  
        CASE (36) ; CALL mrcs_sited(2)  
        END SELECT
        CALL plot_psend()
        CASE ( 3) ; CALL crns_save()                ! Create CRNS
        CASE ( 4) ; CALL mean_tree()                ! Means of each tree
        CASE ( 5) ; CALL raw_stats()                ! Stats of raw files
        CASE ( 6) ; CALL miss_report()              ! Reports missing rings
        CASE ( 7) ; CALL save_missing()             ! Infilled missing values
        CASE ( 9) ; CALL tree_crn()                 ! Tree to CRN format
        CASE (10) ; CALL column_crn(TR)             ! col.ind crns to columns
        CASE (13) ; CALL list_corey()               ! By site cores sorted by fy/ly
        CASE (14) ; CALL list_coren()               ! All cores sorted by name
        CASE (15) ; CALL KKtree(1) ; CALL KKtreed() ! Tree distribution fy 
        CASE (16) ; CALL KKtree(2) ; CALL KKtreed() ! Tree distribution ly 
        CASE (18) ; CALL pith_year()                ! Pith Estimates from radius
        CASE (20) ; CALL ARSTAN_pith()              ! Save column format 
        CASE (22) ; CALL all_sites()                ! Join all raw files 
        CASE (26:27) ; CALL EPS_prep(ref1-25)
                       CALL EPS_prep2d()            ! RCS-EPS Polar Sites
        CASE (29:30) ; CALL Adjust_MXD(ref1-28) 
                       CALL Adjust_MXDd()           ! Adjust mean & SDev 
        CASE (32) ; CALL boot_cru() ; CALL boot_crud() ! Create .boot file
        CASE (35) ; CALL mrcs_site(1) ; CALL mrcs_sited(1)  
        CASE (36) ; CALL mrcs_site(2) ; CALL mrcs_sited(2)  
      END SELECT
      b(272)%ok=TR
      IDO: DO i=1,50
        CALL but_draw(272,"")  ! b(271)%on=TR ; 
        CALL mouse_click(6,272,272) ; IF (mous.EQ.272) EXIT IDO
      ENDDO IDO
      b(272)%ok=FA
      RETURN 
      END SUBROUTINE rf_plot
!------------------------------------------------------------------------
      SUBROUTINE pith_sort(ref1,gok)  ! Sorts for PO work
      IMPLICIT NONE                 
      INTEGER,INTENT(OUT),DIMENSION(mxs) :: ref1 ! List in order
      LOGICAL,INTENT(OUT),DIMENSION(mxs) :: gok  ! First core of tree
      CHARACTER(12),DIMENSION(mxs) :: snam       ! Adjusted Tree names
      INTEGER      :: i,j,k,m,n,p,q
      CHARACTER(7) :: znam
      DO i=1,nc           ! Convert all names to upper case
        snam(i)=nam(i)
        DO j=1,8                   
          p=ICHAR(snam(i)(j:j))
          IF (p.GE.97.AND.p.LE.122) snam(i)(j:j)=CHAR(p-32)
        ENDDO 
      ENDDO
      gok=TR ; n=nc+1 ; snam(n)="ZZZZZ       " 
      DO i=1,nc           ! Sort cores on name
        DO j=1,nc ; IF (gok(j).AND.snam(j)(1:8).LT.snam(n)(1:8)) n=j ; ENDDO
        ref1(i)=n ; gok(n)=FA ; n=nc+1
      ENDDO
      i=1
      QD: DO q=1,nc     ! Sort by first year of core as well
        n=ref1(i) ; gok(i)=TR ; m=LEN_TRIM(snam(n)(1:8))-1
        znam(1:m)=snam(n)(1:m) 
        JD: DO j=i+1,nc
          IF (LEN_TRIM(snam(ref1(j))(1:8)).NE.m+1) EXIT JD
          IF (znam(1:m).NE.snam(ref1(j))(1:m)) EXIT JD
        ENDDO JD        ! i:j cores with same name start
        j=j-1 ; p=i     ! p = core with earliest year
        DO k=i+1,j ; IF (fy(ref1(k)).LT.fy(ref1(p))) p=k ; ENDDO
        ref1(i)=ref1(p) ; ref1(p)=n ; i=j+1 ; IF (i.GT.nc) EXIT QD
      ENDDO QD
      RETURN 
      END SUBROUTINE pith_sort
!------------------------------------------------------------------------
      SUBROUTINE mean_tree1()  
      IMPLICIT NONE                      
      REAL(8),DIMENSION(mxy) :: wka
      INTEGER,DIMENSION(mxy) :: wki
      INTEGER                :: i,j,k,m,n,p,q,r,u,v,w
      INTEGER,DIMENSION(mxs) :: ref1,reas
      LOGICAL,DIMENSION(mxs) :: gok
      LOGICAL                :: ggok,olap
      reas=0 ; CALL pith_sort(ref1,gok) ; i=1
      DO w=1,nc    ! Look for pith discrepancy
        J3: DO j=i,nc-1 ; IF (gok(j+1)) EXIT J3 ; ENDDO J3
        IF (i.LT.j) THEN
          n=ref1(i) ; p=fy(n) ; q=NINT(10.D0*pthr(n))
          DO k=i+1,j                              
            IF (p.EQ.fy(ref1(k)).AND.q.NE.NINT(10.D0*pthr(ref1(k)))) &
              WRITE(35,'("Pith ??? ",A8,2X,A8)') nam(n),nam(ref1(k))
            IF (q.GT.NINT(10.D0*pthr(ref1(k)))) &
              WRITE(35,'("Pith problem ",A8,2X,A8)') &
                nam(n),nam(ref1(k))
          ENDDO
        ENDIF
        i=j+1
      ENDDO
      i=1
      W1: DO w=1,nc              ! For each tree
        J1: DO j=i,nc-1 ; IF (gok(j+1)) EXIT J1 ; ENDDO J1
        IF (i.LT.j) THEN         ! i:j Cores to average
          ggok=FA ; m=LEN_TRIM(nam(ref1(i)))
          IF (nam(ref1(i))(m:m).EQ."M") THEN      ! Already mean core 
            ggok=TR ; reas(i:j)=1
          ENDIF
          DO k=i+1,j
            IF (m.NE.LEN_TRIM(nam(ref1(k)))) THEN ! Name wrong length
              ggok=TR ; reas(i:j)=2
            ENDIF
            IF (nam(ref1(k))(m:m).EQ."M") THEN    ! Already mean core 
              ggok=TR ; reas(i:j)=1
            ENDIF
          ENDDO    
          gok(i+1:j)=ggok
          IF (.NOT.ggok) THEN  ! Only average cores that overlap
            DO k=i,j
              olap=FA
              MD: DO m=i,j
                IF (m.NE.k) olap=(fy(ref1(k)).LT.ly(ref1(m)).AND. &
                                  ly(ref1(k)).GT.fy(ref1(m)))
                IF (olap) EXIT MD
              ENDDO MD
              IF (.NOT.olap) reas(k)=3 
            ENDDO    
            p=i
            DO k=i,j
              IF (reas(k).EQ.3) THEN  ! No overlap core, swap with p
                q=ref1(p)  ; ref1(p)=ref1(k) ; ref1(k)=q
                q=reas(p) ; reas(p)=reas(k) ; reas(k)=q
                gok(k)=FA ; p=p+1 ; gok(p)=TR 
              ENDIF
            ENDDO
          ENDIF         
        ELSE
          reas(i)=4
        ENDIF
        i=j+1 ; IF (i.GT.nc) EXIT W1
      ENDDO W1 
      WRITE(35,'("Name     Start   End Years   POy  POcm  Message")') 
      DO i=1,nc
        j=ref1(i)
        SELECT CASE (reas(i))
          CASE (0) ; m=LEN_TRIM(nam(j)) ; mess=nam(j)(1:m-1)//"M"
          CASE (1) ; mess="One core is a Mean"
          CASE (2) ; mess="Name length misfit" 
          CASE (3) ; mess="Core - no overlap " 
          CASE (4) ; mess="One core this tree"
        END SELECT
        WRITE(35,'(A8,4I6,F6.1,2X,A18)') nam(j),fy(j), &
          ly(j),yr(j),fy(j)-pth(j),pthr(j)/10.D0,mess(1:18)
      ENDDO
      i=1 ; m=nc ; ad(m+1)=ad(nc)+yr(nc) 
      W2: DO w=1,nc              ! For each tree
        J2: DO j=i,nc-1 ; IF (gok(j+1)) EXIT J2 ; ENDDO J2
        m=m+1 ; n=ref1(i)
        IF (i.LT.j) THEN         ! i:j Cores to average
          wka=0.D0 ; wki=0 ; fy(m)=fy(n) ; ly(m)=ly(n) 
          DO k=i+1,j             ! Year range of mean tree
            fy(m)=MIN(fy(m),fy(ref1(k)))
            ly(m)=MAX(ly(m),ly(ref1(k)))
          ENDDO
          yr(m)=ly(m)-fy(m)+1 ; ad(m+1)=ad(m)+yr(m)
          r=LEN_TRIM(nam(n))-1
          nam(m)(1:r+1)=nam(n)(1:r)//"M"
          nam(m)(r+2:12)=" "
          DO k=i,j
            n=ref1(k) ; r=yr(n) ; p=fy(n)-fy(m)+1 ; q=p+r-1
            u=ad(n) ; v=u+r-1 
            IF (p.EQ.1) THEN     ! Pith from earliest core
              pth(m)=pth(n) ; pthr(m)=pthr(n)
            ENDIF
            WHERE (xok(u:v))
              wka(p:q)=wka(p:q)+x(u:v) ; wki(p:q)=wki(p:q)+1
            ENDWHERE
          ENDDO
          r=yr(m) ; p=ad(m) ; q=ad(m)+r-1 
          WHERE (wki(1:r).GT.0)
            x(p:q)=wka(1:r)/DBLE(wki(1:r)) ; xok(p:q)=TR
          ELSEWHERE
            x(p:q)=-9.99D0 ; xok(p:q)=FA
          END WHERE
        ELSE                     ! Save this single core
          n=ref1(i) ; fy(m)=fy(n) ; r=yr(n) ; yr(m)=r ; ly(m)=ly(n)
          pth(m)=pth(n) ; pthr(m)=pthr(n) ; ad(m+1)=ad(m)+r
          nam(m)=nam(n) ; u=ad(n) ; v=u+r-1 ; p=ad(m) ; q=p+r-1
          x(p:q)=x(u:v) ; xok(p:q)=xok(u:v)
        ENDIF
        i=j+1 ; IF (i.GT.nc) EXIT W2
      ENDDO W2 
      n=nc+1 ; r=m-n+1 ; u=ad(n) ; v=ad(m)+yr(m)-1 ; w=v-u+1 ; nc=r
      fy(1:r)=fy(n:m) ; ly(1:r)=ly(n:m)    ; yr(1:r)=yr(n:m)
      ad(1:r)=ad(n:m)-u+1 ; pth(1:r)=pth(n:m) ; pthr(1:r)=pthr(n:m)
      nam(1:r)=nam(n:m) ; x(1:w)=x(u:v) ; xok(1:w)=xok(u:v)
      WRITE(35,*) ; WRITE(35,'("MEAN FILE trees")') 
      WRITE(35,'("Name     Start   End Years   POy  POcm")') 
      DO j=1,nc
        WRITE(35,'(A8,4I6,F6.1)') nam(j),fy(j), &
        ly(j),yr(j),fy(j)-pth(j),pthr(j)/10.D0
      ENDDO
!      OPEN(74,FILE="zzz.prn",IOSTAT=ios,STATUS="REPLACE")
!      DO j=1,nc 
!        WRITE(74,'(A10,I6,L2,I6,2X,A10)') &
!          nam(j),j,gok(j),ref1(j),nam(ref1(j))
!      ENDDO
!      CLOSE(74) ; STOP
      RETURN 
      END SUBROUTINE mean_tree1
!----------------------------------------------------------------
      SUBROUTINE mean_tree() ! Replaces cores by mean tree values.
! Presumes last used digit in core names is core designation and replaces
! with M in output file. Adds M to output file name. Only averages cores
! that overlap and copies across other cores. Reports to mean.prn file. 
      IMPLICIT NONE                      
      INTEGER       :: i,j,k
      OPEN(64,FILE="mean.fil",IOSTAT=ios,STATUS='OLD')
      IF (io_err("Open","mean.fil")) RETURN
      OPEN(35,FILE="mean.prn",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open","mean.prn")) RETURN
      cf=1
      KD: DO k=1,mxs         ! For each file
        READ(64,'(A60)',IOSTAT=ios) cnam(cf) 
        IF (ios.EQ.-1) EXIT KD
        IF (io_err("Read","Mean Tree")) STOP
        i=LEN_TRIM(cnam(cf))
        WRITE(35,'(A60)') cnam(cf) 
        nc=0 ; CALL read_rft(cnam(cf))
        CALL mean_tree1()
        JD: DO j=i,i-6,-1   ! Look for filename extention 
          IF (cnam(cf)(j:j).EQ.".") EXIT JD
        ENDDO JD
        IF (j.LT.i-6) j=i+1
        cnam(cf)=cnam(cf)(1:j-1)//"m"//cnam(cf)(j:i)
        tre(1:nc)=(/(i,i=1,nc)/)
        CALL write_raw(cnam(cf),x)
        CALL write_pith(cnam(cf)(1:j+1)//"pth")
        WRITE(35,*)
      ENDDO KD
      CLOSE(64) ; CLOSE(35)
      RETURN 
      END SUBROUTINE mean_tree
!----------------------------------------------------------------
      SUBROUTINE raw_stats()  ! Some statistics of raw data files
      IMPLICIT NONE                      
      LOGICAL,DIMENSION(mxt) :: qok
      LOGICAL  :: covok   ! Correlation suceeded
      INTEGER  :: tall,rall,i,j,k,m,p,q,r,s,t,u,v,nr
      REAL(8)  :: mnc,mnr,sd,mraw
      OPEN(30,FILE="rawstats.prn",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open","rawstats.prn")) RETURN
      CALL det_default() ; cf=1 ; idt=30 ; tall=0 ; rall=0
      KDO: DO k=1,nfil          ! For each file in ars.fil
        nc=0 ; CALL read_rft(cnam(k)) ; CALL detrend()
        IF (MOD(k,40).EQ.1) WRITE(30,'("File name",52X, &
          & "Cores Start   End Years   Rings  Corr  RBar  MnRaw")')
        m=ad(nc)+yr(nc)-1 ; mnc=0.D0 ; mnr=0.D0 ; nr=0
        mraw=SUM(x(1:m),MASK=xok(1:m))/DBLE(COUNT(xok(1:m)))
        DO i=1,nc-1
          DO j=i+1,nc
            p=MAX(fy(i),fy(j)) ; q=MIN(ly(i),ly(j)) ; r=q-p+1
            u=ad(i)+p-fy(i) ; v=u+r-1
            s=ad(j)+p-fy(j) ; t=s+r-1
            qok(1:r)=xok(u:v).AND.xok(s:t)   ! Exclude missing rings
            IF(COUNT(qok(1:r)).GT.20) THEN
              CALL covmiss(dx(u:v),dx(s:t),qok(1:r),r,sd,covok)
              IF (covok) THEN
                mnr=mnr+sd ; nr=nr+1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF (nr.GE.1) THEN
          mnr=mnr/DBLE(nr)
        ELSE
          mnr=-2.D0   ! No valid correlations
        ENDIF
        nr=0 
        DO i=1,nc    ! Tree correlation with CRN
          r=yr(i) ; p=ad(i) ; q=p+r-1
          u=fy(i)-cfy(1)+1 ; v=u+r-1
          qok(1:r)=num(u:v,1).GT.1.AND.xok(p:q)  ! Some trees no overlap
          IF (COUNT(qok(1:r)).GT.20) THEN  ! Remove this tree
            crn(u:v,3)=crn(u:v,1)-dx(p:q)/DBLE(MAX(num(u:v,1),1))
            CALL covmiss(dx(p:q),crn(u:v,3),qok(1:r),r,sd,covok)
            IF (covok) THEN 
              mnc=mnc+sd ; nr=nr+1
            ENDIF
          ENDIF
        ENDDO
        IF (nr.GE.1) THEN
          mnc=mnc/DBLE(nr)   ! Mean correlation
        ELSE
          mnc=-2.D0   ! No valid correlations
        ENDIF
        WRITE(30,'(A60,4I6,I8,2F6.2,F7.3)') &
          cnam(k),nc,cfy(1),cly(1),cyr(1),m,mnc,mnr,mraw
        tall=tall+nc ; rall=rall+m
      ENDDO KDO
      WRITE(30,*)
      WRITE(30,'(16X,"Cores     Rings     Sites")') 
      WRITE(30,'("Totals     ",3I10)') tall,rall,k-1
      CLOSE(30)
      RETURN 
      END SUBROUTINE raw_stats
!----------------------------------------------------------------
      SUBROUTINE miss_report()  ! Missing or zero rings report
      IMPLICIT NONE                      
      INTEGER                :: i,j,k,m,n,p,q,r
      INTEGER,DIMENSION(mxt) :: wka
      OPEN(30,FILE="missing.prn",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open","missing.prn")) RETURN
      KDO: DO k=1,nfil          ! For each file
        nc=0 ; CALL read_rft(cnam(k)) 
        DO i=1,nc
          p=ad(i) ; r=yr(i) ; q=p+r-1 ; n=p
          wka(1:r)=NINT(x(p:q)*1000.D0)
          IF (ANY(wka(1:r).EQ.0).OR.ANY(xok(p:q).EQV.FA)) THEN
            JD: DO j=p,q
              N1: DO n=n,q
                IF (.NOT.xok(n).OR.wka(n-p+1).EQ.0) EXIT N1
              ENDDO N1
              M1: DO m=n+1,q
                IF (xok(m).AND.wka(m-p+1).GT.0) EXIT M1
              ENDDO M1
              IF (n.GT.q) EXIT JD
              IF (m.EQ.n+1.AND.wka(n-p+1).EQ.0) THEN
                WRITE(30,'(I4,"    zero ring  ",A8,"  ",A60)') &
                n-p+fy(i),nam(i),cnam(k)
              ELSEIF (m.EQ.n+1) THEN
                WRITE(30,'(I4," missing ring  ",A8,"  ",A60)') &
                  n-p+fy(i),nam(i),cnam(k)
              ELSEIF (wka(n-p+1).EQ.0) THEN
                WRITE(30,'(I4,"-",I4," are zero ",A8,"  ",A60)') &
                  n-p+fy(i),m-p+fy(i)-1,nam(i),cnam(k)
              ELSE
                WRITE(30,'(I4,"-",I4," missing  ",A8,"  ",A60)') &
                n-p+fy(i),m-p+fy(i)-1,nam(i),cnam(k)
              ENDIF
              n=m
            ENDDO JD
          ENDIF
       ENDDO
      ENDDO KDO
      CLOSE(30)
      RETURN 
      END SUBROUTINE miss_report
!----------------------------------------------------------------
      SUBROUTINE save_missing()  ! Save missing rings
      IMPLICIT NONE                      
      INTEGER :: m
      nc=0 ; CALL read_rft(cnam(cf)) ; m=ad(nc)+yr(nc)-1
      tst=1 ; CALL tree_sort()   ! Original sequence
      xok(1:m)=TR ; m=LEN_TRIM(cnam(cf))
      MD: DO m=m,3,-1 ; IF (cnam(cf)(m:m).EQ.".") EXIT MD ; ENDDO MD
      IF (cnam(cf)(m+1:m+2).EQ."fn") THEN
        CALL write_heidel(cnam(cf)(1:m)//"hei",x)
      ELSE
        CALL write_raw(cnam(cf)(1:m)//"mis",x)
      ENDIF
      RETURN 
      END SUBROUTINE save_missing
!----------------------------------------------------------------
      SUBROUTINE tree_crn() ! Save files in column format
      IMPLICIT NONE                 
      INTEGER  :: i,p,q,r
      OPEN(23,FILE="trees.crn",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open Replace","trees.crn")) STOP
      DO i=1,nc   ! Align tree as a CRN
        r=yr(i) ; cfy(cf)=fy(i) ; cly(cf)=ly(i) ; cyr(cf)=r
        p=ad(i) ; q=p+r-1
        crn(1:r,cf)=dx(p:q)
        WHERE (xok(p:q)) 
          num(1:r,cf)=1
        ELSEWHERE
          num(1:r,cf)=0
        END WHERE
        namc=nam(i)
        WRITE(23,'(20X,A8)',IOSTAT=ios) nam(i)
        IF (io_err("write ","trees.crn")) STOP
        CALL write_ind(r,cfy(cf),crn(1:r,cf),num(1:r,cf), &
          "TRE",nam(i))
      ENDDO
      CLOSE(23) 
      RETURN 
      END SUBROUTINE tree_crn
!------------------------------------------------------------------------
      SUBROUTINE column_crn(ref1) ! Save files in column format
      IMPLICIT NONE                 
      LOGICAL,INTENT(IN) :: ref1  ! Read data or not 
      INTEGER            :: i,j,k,p,q,r,u,v,w
      IF (ref1) THEN 
        icf=0 ; CALL read_index("crns.crn")          
        OPEN(23,FILE="Column.prn",IOSTAT=ios,STATUS="REPLACE")
      ENDIF
      p=MINVAL(cfy(1:icf)) ; q=MAXVAL(cly(1:icf)) ; r=q-p+1
      DO i=icf,1,-1    ! Align data and insert missing values
        k=i*2 ; j=k-1 ;    u=cfy(i)-p+1 ; w=cyr(i) ; v=u+w-1
        num(u:v,k)=num(1:w,i) ; num(u:v,j)=NINT(crn(1:w,i)*1000.0D0)
        num(1:u-1,k)=0 ; num(v+1:r,k)=0 ; num(1:u-1,j)=-9999
        num(v+1:r,j)=-9999 ; wnam(j)=wnam(i) ; wnam(k)=""
      ENDDO
      DO j=1,2*icf,60
        u=j ; v=MIN(j+59,icf*2)
        WRITE(23,'(I3,"  Year")') 0
        DO i=1,v-u+1,2
          WRITE(23,'(I3,2X,A20)') i+u-1,wnam(i+u-1)
        ENDDO
        WRITE(23,'(I6,30(I6," num"))') 0,(/(i,i=1,v-u+1,2)/)
        DO i=1,r
          WRITE(23,'(I6,30(I6,I4))') p-1+i,num(i,u:v)
        ENDDO
      ENDDO
      CLOSE(23)
      RETURN 
      END SUBROUTINE column_crn
!------------------------------------------------------------------------
      SUBROUTINE crns_save()  ! Creates and saves CRNs
      IMPLICIT NONE                 
      INTEGER       :: i,j,k,m
      CHARACTER(60) :: tnam
      CHARACTER(4)  :: typ
      INTEGER       :: xad,yad
      append=FA ; tnam=" " ; yad=230 ; xad=50
      CALL SETCLR(black) 
      OPEN(29,FILE="crns.fil",IOSTAT=ios,STATUS='OLD')
      IF (io_err("Open","crns.fil")) RETURN
      READ(29,*)                   ! File header
      IF (cft.EQ.2) THEN
        OPEN(23,FILE="Column.prn",IOSTAT=ios,STATUS="REPLACE")
        IF (io_err("Open","crns.fil")) RETURN
      ENDIF 
      ID: DO i=1,mxc               ! For each CRN               
        cf=i 
        READ(29,'(A60,A4)',IOSTAT=ios) tnam,typ
        IF (ios.EQ.-1) EXIT ID     ! End of file
        IF (io_err("Read","crns.fil")) RETURN
        READ(29,'(3X,15I4,4I6)',IOSTAT=ios) idt,itn,rdt,ind,krb,isb, &
          sfo,poo,src,trc,gtr,tst,bfc,idb,jrb,srcno,sfono,rdtno
        IF (io_err("Read",TRIM(tnam)//" crns.fil")) RETURN
        CALL par_valid(tnam)
        CALL MESSAG(tnam,xad,yad)
        CALL SENDBF()   ! Update screen 
        yad=yad+50 ; nc=0
        IF (yad.GT.1800) THEN ; xad=xad+500 ; yad=230 ; ENDIF
        CALL read_rft(tnam) ; CALL detrend()
        j=LEN_TRIM(tnam)
        MD: DO m=j,j-3,-1 ; IF (tnam(m:m).EQ.".") EXIT MD ; ENDDO MD
        IF (m.EQ.j-4) m=j
        KD: DO k=j,1,-1
          IF (tnam(k:k).EQ."/".OR.tnam(k:k).EQ.CHAR(92)) EXIT KD
        ENDDO KD
        wnam(cf)=tnam(k+1:m-1)//" "//typ
        IF (cft.EQ.1) THEN
          CALL write_index("crns.crn") ; append=TR
        ELSE
          CALL crn_head() 
        ENDIF
      ENDDO ID
      CLOSE(29)    
      IF (cft.EQ.2) THEN
        WRITE(23,*)
        icf=i-1 ; CALL column_crn(FA)
      ENDIF
      RETURN 
      END SUBROUTINE crns_save
!------------------------------------------------------------------------
      SUBROUTINE list_corey() 
      IMPLICIT NONE                 
      INTEGER,PARAMETER             :: maxc=20000
      INTEGER                       :: i,j,m,n,p,q
      CHARACTER(12),DIMENSION(maxc) :: gnam
      INTEGER,DIMENSION(maxc)       :: gnum,gfy,gly,gst,ref1
      REAL(8),DIMENSION(maxc)       :: gp,gsum
      LOGICAL,DIMENSION(maxc)       :: gok
      CHARACTER(100)                :: rec1
      m=0
      ID: DO i=1,nfil
        nc=0 ; CALL read_rft(cnam(i))
        DO j=1,nc
          IF (m.EQ.maxc) THEN
            WRITE(rec1,'("Max cores",I6," exceeded")') maxc
            CALL out_err(rec1) ; EXIT ID
          ENDIF
          m=m+1 ; gnam(m)=nam(j) ; gnum(m)=i ; gfy(m)=fy(j)
          gly(m)=ly(j) ; p=ad(j) ; q=p+yr(j)-1 ; gp(m)=pthr(j)/10.D0
          gsum(m)=SUM(x(p:q),MASK=xok(p:q))
          gst(m)=fy(j)*mxt+yr(j)
        ENDDO
      ENDDO ID
      OPEN(72,FILE="core_year.prn",IOSTAT=ios,STATUS='REPLACE')
      IF (io_err("Open","corey.prn")) STOP   
      WRITE(72,'("Files = ",I4,"  Trees = ",I6)') i-1,m
      WRITE(72,*) ; WRITE(72, &
        '("Name     Start   End  Pith Sum Rings  File Name")')
      WRITE(72,*) ; gok=TR ; gst(m+1)=10000*mxt ; n=m+1
      DO i=1,m                             ! Sort all the trees
        DO j=1,m
          IF (gok(j)) THEN
            IF (gst(j).LT.gst(n)) THEN     ! First key fy*mxt+ly 
              n=j
            ELSEIF (gst(j).EQ.gst(n)) THEN
              IF (gsum(j).LT.gsum(n)) n=j  ! Second key ring sum
            ENDIF
          ENDIF 
        ENDDO
        ref1(i)=n ; gok(n)=FA ; n=m+1
      ENDDO
      p=ref1(m)
      DO i=1,m                             ! Print all the trees
        n=ref1(i) 
        WRITE(rec1,'(A8,2I6,F6.1,F10.3,2X,A60)') &
          gnam(n),gfy(n),gly(n),gp(n),gsum(n),cnam(gnum(n)) 
        IF (gst(n).EQ.gst(p)) THEN ! Omit repeated fields
          IF (gfy(n).EQ.gfy(p))   rec1( 9:14)="      "
          IF (gly(n).EQ.gly(p))   rec1(15:20)="      "
          IF (gp(n).EQ.gp(p))     rec1(21:26)="      "
          IF (gsum(n).EQ.gsum(p)) rec1(27:36)="          "
        ENDIF
        p=n ; WRITE(72,'(A100)') rec1
      ENDDO
      CLOSE(72)
      RETURN 
      END SUBROUTINE list_corey
!------------------------------------------------------------------------
      SUBROUTINE list_coren() ! Cores from many files sorted by name
      IMPLICIT NONE                 
      INTEGER,PARAMETER             :: maxc=20000
      INTEGER                       :: i,j,k,m,n,p,q
      CHARACTER(12),DIMENSION(maxc) :: gnam
      INTEGER,DIMENSION(maxc)       :: gnum,gfy,gly,miss,ref1
      REAL(8),DIMENSION(maxc)       :: gp,gsum
      LOGICAL,DIMENSION(maxc)       :: gok
      CHARACTER(103)                :: rec1 
      m=0
      ID: DO i=1,nfil           ! Read all trees
        nc=0 ; CALL read_rft(cnam(i))
        DO j=1,nc
          IF (m.EQ.maxc) THEN
            WRITE(rec1,'("Max cores",I6," exceeded")') maxc
            CALL out_err(rec1) ; EXIT ID
          ENDIF
          m=m+1 ; gnam(m)=nam(j)
          DO k=1,8              ! Convert all names to upper case
            p=ICHAR(gnam(m)(k:k))
            IF (p.GE.97.AND.p.LE.122) gnam(m)(k:k)=CHAR(p-32)
          ENDDO 
          gnum(m)=i ; gfy(m)=fy(j) ; gly(m)=ly(j) ; p=ad(j)
          q=p+yr(j)-1 ; gp(m)=pthr(j)/10.D0
          gsum(m)=SUM(x(p:q),MASK=xok(p:q))/10.D0
          miss(m)=COUNT(LOGICAL(.NOT.xok(p:q),1))
        ENDDO
      ENDDO ID
      OPEN(72,FILE="core_names.prn",IOSTAT=ios,STATUS='REPLACE')
      IF (io_err("Open","coren.prn")) STOP   
      WRITE(72,'("Files = ",I4,"  Trees = ",I6)') i-1,m
      WRITE(72,*) ; WRITE(72, &
        '("Name     Start   End Miss  Pith Sum Rings  File Name")')
      WRITE(72,*) ; gok=TR ; gnam(m+1)="~~~~~~~~~~~~" ; n=m+1
      DO i=1,m                  ! Sort all the trees
        DO j=1,m ; IF (gok(j).AND.gnam(j).LT.gnam(n)) n=j ; ENDDO
        ref1(i)=n ; gok(n)=FA ; n=m+1
      ENDDO
      p=ref1(m)
      DO i=1,m                  ! Print all the trees
        n=ref1(i) 
        WRITE(rec1,'(A8,2I6,I5,F6.1,F10.3,2X,A60)') &
          gnam(n),gfy(n),gly(n),miss(n),gp(n),gsum(n),cnam(gnum(n)) 
        IF (miss(n).EQ.0) rec1(21:25)="     "
        IF (gnam(n).EQ.gnam(p)) THEN ! Omit repeated fields
          rec1( 1: 8)="        "
          IF (gfy(n).EQ.gfy(p))   rec1( 9:14)="      "
          IF (gly(n).EQ.gly(p))   rec1(15:20)="      "
          IF (miss(n).EQ.miss(p)) rec1(21:25)="     "
          IF (gp(n).EQ.gp(p))     rec1(26:31)="      "
          IF (gsum(n).EQ.gsum(p)) rec1(32:41)="          "
        ENDIF
        p=n ; WRITE(72,'(A103)') rec1
      ENDDO
      CLOSE(72)
      RETURN 
      END SUBROUTINE list_coren
!------------------------------------------------------------------------
      SUBROUTINE KKtreed() ! Display tree distribution
      IMPLICIT NONE                 
      INTEGER                :: i=1,j,p,q,r
      REAL(8)                :: rr
      nc=MIN(nc,500)
      p=cfy(cf) ; q=cly(cf) ; r=q-p+1
      grl=400 ; grr=2600 ; grt=140
      IF (nc.LT.199) THEN
         grb=1200
      ELSE
         grb=1800
      ENDIF
      CALL plot_trees(r,num(1:r,cf))  
      CALL NAME('Calendar Year','X')   ! Axis name
      CALL NAME('Each Tree','Y')       ! Axis name
      CALL tombox(p,q,0.D0,DBLE(nc))
      CALL SETCLR(grey) ; CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22) ; CALL LINWID(1)
      CALL SETCLR(blue) ; CALL MESSAG("Pith Estimates",400,grt-40)
      CALL SETCLR(red)  ; CALL MESSAG("Measured Rings",1100,grt-40)
      CALL LINWID(MAX((grb-grt)/(2*nc),1)) 
      DO i=1,nc
        j=tre(i) ; rr=DBLE(i)-0.5D0
        IF (pth(j).LT.fy(j)-1) THEN
          CALL SETCLR(blue)
          CALL RLINE(DBLE(pth(j)),rr,DBLE(fy(j)),rr)
          CALL SETCLR(red)
        ENDIF
        CALL RLINE(DBLE(fy(j)),rr,DBLE(ly(j)),rr)
      ENDDO
      CALL SETCLR(black) ; CALL LINWID(4)
      CALL MESSAG("Tree distribution through time "//cnam(cf),grl,grt-100)
      CALL ENDGRF()
      CALL LINWID(1)
      RETURN 
      END SUBROUTINE KKtreed
!-----------------------------------------------------------------
      SUBROUTINE KKtree(fl)  ! Plot tree distribution through time 
      IMPLICIT NONE
      INTEGER,INTENT(IN)     :: fl
      REAL(8),DIMENSION(mxs) :: wk    ! Data area
      IF (fl.EQ.1) THEN
        wk(1:nc)=DBLE(pth(1:nc))      ! Sort by pith year
      ELSE
        wk(1:nc)=DBLE(ly(1:nc))       ! Sort by last year
      ENDIF
      CALL pair_sort(nc,wk(1:nc),tre(1:nc))
      CALL det_crnfy()  
      RETURN
      END SUBROUTINE KKtree
!-------------------------------------------------------------------
      SUBROUTINE pith_year()  ! Plots PO and sums of relevant rings
      IMPLICIT NONE                 
      INTEGER                :: i,n,p,r
      INTEGER,DIMENSION(mxs) :: ref1
      LOGICAL,DIMENSION(mxs) :: gok
      REAL(8)                :: rr
      OPEN(52,FILE="Core_pith.prn",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open","Core_pith.prn")) STOP  
      nc=0 ; CALL read_rft(cnam(cf))
      CALL pith_sort(ref1,gok)
      WRITE(52,*) cnam(cf)
      WRITE(52,'("Name        Pith    Pith  Start    End      Pith", & 
               &  "     Missing   1st 50")')
      WRITE(52,'("            Year  Offset   Year   Year  Estimate", &
               &  "  Radius(cm)    rings")')
      WRITE(52,*)
      DO i=1,nc                 
        n=ref1(i) 
        IF (gok(i)) THEN  ! Earliest core sum pith years of rings
          r=n ; p=ad(n) ; rr=SUM(x(p:p-1+fy(n)-pth(n)))/10.D0
          WRITE(52,'(A8,I8,F8.1,2I7,F10.3,12X,F9.3)') nam(n),pth(n), &
            pthr(n)/10.D0,fy(n),ly(n),rr,SUM(x(p:p+49))/10.D0
          pthr(n)=pthr(n)/10.D0
        ELSE              ! Later core, sum rings of earliest core
          rr=SUM(x(p:p-1+fy(n)-fy(r)))/10.D0
          WRITE(52,'(A8,I8,F8.1,2I7,10X,F12.3)') &
            nam(n),pth(n),pthr(n)/10.D0,fy(n),ly(n),rr
          pth(n)=pth(r) ; pthr(n)=pthr(r)+rr
        ENDIF
      ENDDO
      OPEN(53,FILE="blank.pth",IOSTAT=ios,STATUS="REPLACE")
      DO i=1,nc 
        WRITE(53,'(A8,I8,F8.1)') nam(i),pth(i),pthr(i)
      ENDDO
      CLOSE(52) ; CLOSE(53)
      RETURN 
      END SUBROUTINE pith_year
!------------------------------------------------------------------------
      SUBROUTINE ARSTAN_pith()  ! Save column plot Arstan pith
      IMPLICIT NONE                      
      INTEGER  :: i,j,k
      CHARACTER(60) :: anam
      KDO: DO k=1,nfil          ! For each file
        j=LEN_TRIM(cnam(k))
        ID: DO i=j,j-3,-1 ; IF (cnam(k)(i:i).EQ.".") EXIT ID ; ENDDO ID
        IF (cnam(k)(i:i).EQ.".") THEN
          anam=cnam(k)(1:i)//"apth"
        ELSE
          anam=cnam(k)(1:j)//".apth"
        ENDIF
        nc=0 ; CALL read_rft(cnam(k))
        OPEN(69,FILE=anam,IOSTAT=ios,STATUS="REPLACE")
        IF (io_err("Open",anam)) RETURN
          DO i=1,nc
          WRITE(69,'(I4,4X,A8,I6)',IOSTAT=ios) &
            fy(i)-pth(i)-1,nam(i),fy(i)
          IF (io_err("Write",anam)) RETURN
        ENDDO
        CLOSE(69)
      ENDDO KDO
      RETURN 
      END SUBROUTINE ARSTAN_pith
!----------------------------------------------------------------
      SUBROUTINE all_sites() ! Create one raw file from RCS.fil 
      IMPLICIT NONE                      
      INTEGER :: i
      cf=1 ; nc=0
      DO i=1,nfil         ! For each file
        CALL read_rft(cnam(i))
      ENDDO  
      tre(1:nc)=(/(i,i=1,nc)/)  ! Not sorted
      CALL write_raw("allsites.raw",x)
      CALL write_pith("allsites.pth")
      RETURN 
      END SUBROUTINE all_sites
!----------------------------------------------------------------
      SUBROUTINE adjust_MXDd()  
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka
      REAL(8)                :: ra,rb,rc,rd
      INTEGER                :: i=1,p,q,r
      r=MAX(cyr(23),cyr(24))
      ra=MIN(MINVAL(crn(1:r,23),MASK=okc(1:r,23)), &
             MINVAL(crn(1:r,24),MASK=okc(1:r,24)), &
             MINVAL(crn(1:r,27),MASK=okc(1:r,23)), &
             MINVAL(crn(1:r,28),MASK=okc(1:r,24))) 
      rb=MAX(MAXVAL(crn(1:r,23),MASK=okc(1:r,23)), &
             MAXVAL(crn(1:r,24),MASK=okc(1:r,24)), &
             MAXVAL(crn(1:r,27),MASK=okc(1:r,23)), &
             MAXVAL(crn(1:r,28),MASK=okc(1:r,24))) 
      r=cyr(1)
      rc=MIN(MINVAL(crn(1:r,25),MASK=okc(1:r,25)), &
             MINVAL(crn(1:r,26),MASK=okc(1:r,26)), &
             MINVAL(crn(1:r,29),MASK=okc(1:r,25)), &
             MINVAL(crn(1:r,30),MASK=okc(1:r,26))) 
      rd=MAX(MAXVAL(crn(1:r,25),MASK=okc(1:r,25)), &
             MAXVAL(crn(1:r,26),MASK=okc(1:r,26)), &
             MAXVAL(crn(1:r,29),MASK=okc(1:r,25)), &
             MAXVAL(crn(1:r,30),MASK=okc(1:r,26))) 
      CALL NAME('Measures','Y') ! Axis name
      CALL NAME('','X') ! Axis name
      CALL LABELS('NONE','X')
      grl=200 ; grr=2400 ; grt=140 ; grb=540
      wka(1:1000)=(/(DBLE(i),i=1,1000)/)
      p=1 ; q=MINVAL(cly(23:24)) ; r=q-p+1
      CALL plot_trees(r,num(1:q,23))  
      CALL NAME('Mean MXD','Y') ! Axis name
      CALL tombox(1,q,ra,rb)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22) ; CALL SETCLR(blue)
      p=cfy(23) ; q=cly(23) ; r=cyr(23)
      CALL thickthin(r,wka(p:q),crn(p:q,23),num(p:q,23),3)
      CALL MESSAG(TRIM(wnam(1))//" (counts)",grl+800,grb-45)
      CALL SETCLR(red)
      p=cfy(24) ; q=cly(24) ; r=cyr(24)
      CALL thickthin(r,wka(p:q),crn(p:q,24),num(p:q,24),3)
      CALL MESSAG(wnam(2),grl+1500,grt+25)
      CALL SETCLR(black)
      CALL MESSAG("Original RCS Curves",grl+100,grb-45)
      CALL ENDGRF() 

      CALL NAME('Ring Age','X') ! Axis name
      CALL LABELS('FLOAT','X')
      grt=550 ; grb=950
      wka(1:1000)=(/(DBLE(i),i=1,1000)/)
      p=1 ; q=MINVAL(cly(23:24)) ; r=q-p+1
      CALL plot_trees(r,num(1:q,24))  
      CALL NAME('Mean MXD','Y') ! Axis name
      CALL tombox(1,q,ra,rb)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22) ; CALL SETCLR(blue)
      p=cfy(23) ; q=cly(23) ; r=cyr(23)
      CALL thickthin(r,wka(p:q),crn(p:q,27),num(p:q,23),3)
      CALL MESSAG(wnam(3),grl+800,grb-45)
      CALL SETCLR(red)
      p=cfy(24) ; q=cly(24) ; r=cyr(24)
      CALL thickthin(r,wka(p:q),crn(p:q,28),num(p:q,24),3)
      CALL MESSAG(TRIM(wnam(2))//" (counts)",grl+1500,grt+25)
      CALL SETCLR(black)
      CALL MESSAG("Adjusted RCS Curves",grl+100,grb-45)
      CALL ENDGRF()

      CALL NAME('MXD Indices','Y') ! Axis name
      CALL NAME('','X') ! Axis name
      CALL LABELS('NONE','X')
      grt=1060 ; grb=1460
      r=cyr(1) ; wka(1:r)=(/(DBLE(i),i=cfy(1),cly(1))/)
      CALL plot_trees(r,num(1:r,25))  
      CALL NAME('Index Value','Y') ! Axis name
      CALL tombox(cfy(1),cly(1),rc,rd)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22)
      CALL SETCLR(blue)
      p=cfy(25) ; q=cly(25) ; r=cyr(25)
      CALL thickthin(r,wka(p:q),crn(p:q,25),num(p:q,25),3)
      CALL MESSAG(TRIM(wnam(1))//" (counts)",grl+800,grt+25)
       CALL SETCLR(red)
      p=cfy(26) ; q=cly(26) ; r=cyr(26)
      CALL thickthin(r,wka(p:q),crn(p:q,26),num(p:q,26),3)
      CALL MESSAG(wnam(2),grl+1500,grt+25)
      CALL SETCLR(black)
      CALL MESSAG("Original Chronologies",grl+100,grt+25)
      CALL ENDGRF()

      CALL NAME('Calendar Year','X') ! Axis name
      CALL LABELS('FLOAT','X')
      grt=1470 ; grb=1870
      r=cyr(1) ; wka(1:r)=(/(DBLE(i),i=cfy(1),cly(1))/)
      CALL plot_trees(r,num(1:r,26))  
      CALL NAME('Index Value','Y') ! Axis name
      CALL tombox(cfy(1),cly(1),rc,rd)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22)
      CALL SETCLR(blue)
      p=cfy(25) ; q=cly(25) ; r=cyr(25)
      CALL thickthin(r,wka(p:q),crn(p:q,29),num(p:q,25),3)
      CALL MESSAG(wnam(3),grl+800,grt+25)
      CALL SETCLR(red)
      p=cfy(26) ; q=cly(26) ; r=cyr(26)
      CALL thickthin(r,wka(p:q),crn(p:q,30),num(p:q,26),3)
      CALL MESSAG(TRIM(wnam(2))//" (counts)",grl+1500,grt+25)
      CALL SETCLR(black)
      CALL MESSAG("Adjusted Chronologies",grl+100,grt+25)
      CALL ENDGRF()
      RETURN 
      END SUBROUTINE adjust_MXDd
!-----------------------------------------------------------------------
      SUBROUTINE adj_MXD(ref1,nc1,mna,sda)  ! Test signal free RCS
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)  :: ref1   ! MXD=1, TRW=2
      INTEGER,INTENT(IN)  :: nc1    ! Trees in first file
      REAL(8),INTENT(OUT) :: mna    ! Mean adjustment
      REAL(8),INTENT(OUT) :: sda    ! SDev adjustment
      REAL(8)  :: mn1,mn2,sd1,sd2,rr1,rr2,rb 
      INTEGER  :: i,r,u,v,s,t,nr1,nr2
      rr1=0.D0 ; rr2=0.D0 ; sd1=0.D0 ; sd2=0.D0
      nr1=0 ; nr2=0 ; r=cyr(cf)
      DO i=1,nc                     ! Use signal-free indices      
        u=ad(i) ; v=u+yr(i)-1
        s=fy(i)-cfy(cf)+1 ; t=s+yr(i)-1
        IF (i.LE.nc1) THEN
          rr1=rr1+SUM(dx(u:v)/crn(s:t,cf),MASK=xok(u:v))
          sd1=sd1+SUM((dx(u:v)/crn(s:t,cf))**2,MASK=xok(u:v))
          nr1=nr1+COUNT(xok(u:v))
        ELSE
          rr2=rr2+SUM(dx(u:v)/crn(s:t,cf),MASK=xok(u:v))
          sd2=sd2+SUM((dx(u:v)/crn(s:t,cf))**2,MASK=xok(u:v))
          nr2=nr2+COUNT(xok(u:v))
        ENDIF
      ENDDO
      mn1=rr1/DBLE(nr1)
      mn2=rr2/DBLE(nr2)
      sd1=SQRT((sd1-mn1*rr1)/DBLE(nr1-1))   ! SDev 
      sd2=SQRT((sd2-mn2*rr2)/DBLE(nr2-1))   ! SDev 
      WRITE(19,'(2I6,2F9.5)') nr1,nc1,mn1,sd1
      WRITE(19,'(2I6,2F9.5,2I7)') nr2,nc-nc1,mn2,sd2,cf
      r=ad(nc1)+yr(nc1)-1
      mna=1.D0+(mn2-1.D0)*DBLE(nr1+nr2)/DBLE(nr1)
      IF (ref1.EQ.1) THEN  ! Only adjust SDev for MXD 
        rb=SUM(x(1:r),MASK=xok(1:r))/DBLE(COUNT(xok(1:r)))
        x(1:r)=(x(1:r)-rb)*sd2/sd1+rb
        sda=sd2/sd1
      ELSE
        sda=1.D0           ! No adjustment TRW SDev
      ENDIF 
      WHERE (xok(1:r)) x(1:r)=x(1:r)*mna  ! Reset Mean
      RETURN 
      END SUBROUTINE adj_MXD
!-----------------------------------------------------------------------
      SUBROUTINE adjust_MXD(ref1)  ! Test signal free RCS
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: ref1   ! MXD=1, TRW=2
      REAL(8)            :: mna    ! Mean adjustment
      REAL(8)            :: sda    ! SDev adjustment
      REAL(8)            :: mn,sd  ! Accumulated values
      INTEGER            :: i,p,q,s,t,u,v,nc1
      OPEN(19,FILE="Adjust.col",IOSTAT=ios,STATUS="REPLACE")
      i=LEN_TRIM(cnam(1))
      cnam(3)=cnam(1)(1:i-4)//"adj"//cnam(1)(i-3:i)
      cnam(4)=cnam(3)(1:i)//"pth"
      i=LEN_TRIM(wnam(1))
      wnam(3)=wnam(1)(1:i-4)//"adj"//wnam(1)(i-3:i)
      WRITE(19,*) TRIM(cnam(1))
      WRITE(19,*) TRIM(cnam(2))
      WRITE(19,*) TRIM(cnam(3))
      CALL det_default() ; nc=0 ;cf=1
      CALL read_rft(cnam(1)) ; nc1=nc
      CALL read_rft(cnam(2)) 
      sfo=2 ; idt=-2 ; CALL detrend()
      cfy(23)=MINVAL(fy(1:nc1)-pth(1:nc1))        ! First by age
      cly(23)=MAXVAL(ly(1:nc1)-pth(1:nc1))
      cyr(23)=cly(23)-cfy(23)+1
      cfy(24)=MINVAL(fy(nc1+1:nc)-pth(nc1+1:nc))  ! Second by age
      cly(24)=MAXVAL(ly(nc1+1:nc)-pth(nc1+1:nc))
      cyr(24)=cly(24)-cfy(24)+1
      cfy(25)=MINVAL(fy(1:nc1))-cfy(cf)+1         ! First by year
      cly(25)=MAXVAL(ly(1:nc1))-cfy(cf)+1
      cyr(25)=cly(25)-cfy(25)+1
      cfy(26)=MINVAL(fy(nc1+1:nc))-cfy(cf)+1      ! Second by year
      cly(26)=MAXVAL(ly(nc1+1:nc))-cfy(cf)+1 
      cyr(26)=cly(26)-cfy(26)+1
      num(1:mxy,23:26)=0 ; crn(1:mxy,23:30)=0.D0
      DO i=1,nc          
        p=fy(i)-pth(i) ; q=p+yr(i)-1 ; u=ad(i) ; v=u+yr(i)-1
        s=fy(i)-cfy(cf)+1 ; t=s+yr(i)-1
        IF (i.LE.nc1) THEN
          WHERE (xok(u:v))
            crn(p:q,23)=crn(p:q,23)+fx(u:v)   ! By age
            num(p:q,23)=num(p:q,23)+1
            crn(s:t,25)=crn(s:t,25)+dx(u:v)
            num(s:t,25)=num(s:t,25)+1
          END WHERE 
        ELSE
          WHERE (xok(u:v))
            crn(p:q,24)=crn(p:q,24)+fx(u:v)   ! By year
            num(p:q,24)=num(p:q,24)+1
            crn(s:t,26)=crn(s:t,26)+dx(u:v)
            num(s:t,26)=num(s:t,26)+1
          END WHERE 
        ENDIF
      ENDDO
      okc(1:mxy,23:26)=num(1:mxy,23:26).GE.1
      DO i=23,26
        p=cfy(i) ; q=cly(i)
        WHERE(okc(p:q,i)) &
          crn(p:q,i)=crn(p:q,i)/DBLE(num(p:q,i))
      ENDDO
      WRITE(19,'(" Rings Cores  Mean SF     SDev   Iter")')
      CALL adj_MXD(ref1,nc1,mna,sda) ; mn=mna ; sd=sda
      cf=2 ; CALL detrend()
      CALL adj_MXD(ref1,nc1,mna,sda) ; mn=mn*mna ; sd=sd*sda
      cf=3 ; CALL detrend()
      CALL adj_MXD(ref1,nc1,mna,sda) ; mn=mn*mna ; sd=sd*sda
      WRITE(19,'("Mean adj",F9.4,"  SDev adj",F9.4)') mn,sd
      nc=nc1
      tre(1:nc)=(/(i,i=1,nc)/)
      CALL write_raw(cnam(3),x)   ! Save adjusted file    
      CALL write_pith(cnam(4))  
      nc=0 ; cf=4
      CALL read_rft(cnam(3))
      CALL read_rft(cnam(2))
      CALL detrend()
      DO i=1,nc          
        p=fy(i)-pth(i) ; q=p+yr(i)-1 ; u=ad(i) ; v=u+yr(i)-1
        s=fy(i)-cfy(cf)+1 ; t=s+yr(i)-1
        IF (i.LE.nc1) THEN
          WHERE (xok(u:v))
            crn(p:q,27)=crn(p:q,27)+fx(u:v)   ! By age
            crn(s:t,29)=crn(s:t,29)+dx(u:v)
          END WHERE 
        ELSE
          WHERE (xok(u:v))
            crn(p:q,28)=crn(p:q,28)+fx(u:v)   ! By year
            crn(s:t,30)=crn(s:t,30)+dx(u:v)
          END WHERE 
        ENDIF
      ENDDO
      DO i=27,30
        p=cfy(i-4) ; q=cly(i-4)
        WHERE(okc(p:q,i-4)) &
          crn(p:q,i)=crn(p:q,i)/DBLE(num(p:q,i-4))
      ENDDO
      WRITE(19,*)
      WRITE(19,'("Column  1 File 1 RCS counts")')
      WRITE(19,'("Column  2 File 2 RCS counts")')
      WRITE(19,'("Column  3 File 1 RCS values")')
      WRITE(19,'("Column  4 File 2 RCS values")')
      WRITE(19,'("Column  5 File 1 Adjusted values")')
      WRITE(19,'("Column  6 File 2 Adjusted values")')
      DO i=1,MAXVAL(cyr(23:24))
        WRITE(19,'(3I5,4F7.3)') i,num(i,23:24), &
          crn(i,23:24),crn(i,27:28)
      ENDDO
      WRITE(19,*)
      WRITE(19,'("Column  1 File 1 CRN counts")')
      WRITE(19,'("Column  1 File 2 CRN counts")')
      WRITE(19,'("Column  2 File 1 CRN values")')
      WRITE(19,'("Column  2 File 2 CRN values")')
      WRITE(19,'("Column  3 File 1 Adjusted CRN")')
      WRITE(19,'("Column  3 File 2 Adjusted CRN")')
      q=MINVAL(cfy(25:26)) ; p=MAXVAL(cly(25:26))-q+1 
      DO i=1,p
        WRITE(19,'(3I5,4F7.3)') i-1+q,num(i,25:26), &
          crn(i,25:26),crn(i,29:30)
      ENDDO
      CLOSE(19)
      RETURN 
      END SUBROUTINE adjust_MXD
!-------------------------------------------------------------------
      SUBROUTINE EPS_prep2d()  
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka
      INTEGER                :: i=1,p,q,r,u,v
      p=cfy(2) ; q=cly(2)-1 ; r=q-p+1
      u=cfy(1)+p-1 ; v=cfy(1)+q-1
      wka(1:r)=(/(DBLE(i),i=u,v)/) 
      grl=200 ; grr=2400 ; grt=140 ; grb=440
      CALL LABELS('NONE','X')
      CALL NAME('','X')              ! Axis name
      CALL NAME('RBar Values','Y')   ! Axis name
      CALL tombox(u,v,0.0D0,0.99D0)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22) ; CALL SETCLR(black)
      CALL line_miss(r,wka(1:r),crn(p:q,8),okc(p:q,8))
      WRITE(mess(1:21),'(I4,"-year Spline RBAR")') RBsp
      CALL MESSAG(Mess(1:21),grl+800,grt+30)
      CALL SETCLR(red)
      CALL line_miss(r,wka(1:r),crn(p:q,11),okc(p:q,8))
      CALL MESSAG("RCS RBar",grl+1300,grt+30)
      CALL SETCLR(black)
      CALL MESSAG(cnam(60),grl+300,grt-45)
      CALL MESSAG(cnam(61),grl+1200,grt-45)
      CALL ENDGRF() 
      grt=450 ; grb=1050
      CALL LABELS('FLOAT','X')
      CALL NAME('Calendar Year','X')               ! Axis name
      CALL NAME('EPS Values','Y')   ! Axis name
      CALL plot_trees(r,num(p:q,1))  
      CALL tombox(u,v,0.5D0,1.0D0)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(1) ; CALL HEIGHT(22)
      CALL SETCLR(black)
      CALL line_miss(r,wka(1:r),crn(p:q,16),okc(p:q,8))
      CALL MESSAG("Spline EPS",grl+1000,grb-45)
      CALL SETCLR(red) 
      CALL line_miss(r,wka(1:r),crn(p:q,17),okc(p:q,8))
      CALL MESSAG("RCS EPS",grl+1400,grb-45)
      CALL SETCLR(blue)
      CALL line_miss(r,wka(1:r),crn(p:q,19),okc(p:q,8))
      CALL MESSAG("RCS adj eff EPS",grl+1800,grb-45)
      CALL SETCLR(black) ; CALL ENDGRF() 
      RETURN 
      END SUBROUTINE EPS_prep2d
!-------------------------------------------------------------------
      SUBROUTINE EPS_prep(ref1)    ! Investigate RCS EPS 
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ref1
      INTEGER  :: i,cfs
      nc=0 ; cfs=cf                ! Store selection
      CALL read_rft(cnam(cf))      ! Read raw data
      i=LEN_TRIM(wnam(cf))
      wnam(61)=wnam(cf)      
      cnam(62)="EPS_"//wnam(cf)(1:i-4)//".prn"      
      OPEN(71,FILE=TRIM(cnam(62)),IOSTAT=ios,STATUS="REPLACE")
      CALL det_default()
      sfo=2 ; idt=-2
      IF (ref1.EQ.1) THEN
        src=2 ; srcno=2            ! Two RCS curves 
        cnam(61)="Two-Curve, Signal-Free RCS"      
      ELSE
        src=1 ; srcno=1            ! Single RCS curve 
        cnam(61)="One-Curve, Signal-Free RCS"      
      ENDIF
      cf=1 ; CALL EPS_prep1(RBsp)    ! Process data
      CLOSE(71) ; cf=cfs
      RETURN
      END SUBROUTINE EPS_prep
!--------------------------------------------------------------
      SUBROUTINE Boot_write(msnam,n,p,r,w) ! RCS Plots
      IMPLICIT NONE                 
      CHARACTER(40),INTENT(IN) :: msnam
      INTEGER,INTENT(IN)       :: n,p,r,w
      INTEGER                  :: j,k
      CHARACTER(6),PARAMETER   :: frm="(26I4)"
      WRITE(74,IOSTAT=ios,FMT= &
        '(I8,"=N",I8,"=I "," CRN",I4,42X,A6)') r,p,n,frm
      IF (io_err("write2",msnam)) STOP
      DO k=1,r,260   ! Some buffer limit encountered???
        WRITE(74,IOSTAT=ios,FMT=frm) (num(j,3),j=k,MIN(k+259,r)) 
        IF (io_err("write3",msnam)) STOP
      ENDDO
      WRITE(74,IOSTAT=ios,FMT= &
        '(I8,"=N",I8,"=I "," Num",I4,42X,A6)') r,p,n,frm
      IF (io_err("write4",msnam)) STOP
      DO k=1,r,260   ! Some buffer limit encountered???
        WRITE(74,IOSTAT=ios,FMT=frm) (num(j,1),j=k,MIN(k+259,r)) 
        IF (io_err("write5",msnam)) STOP
      ENDDO
      num(1:w,4)=NINT(1000.D0*crn(1:w,2)) ! RCS Curve
      WRITE(74,IOSTAT=ios,FMT= &
        '(I8,"=N",I8,"=I "," RCS",I4,42X,A6)') w,1,n,frm
      IF (io_err("write6",msnam)) STOP
      DO k=1,w,260   ! Some buffer limit encountered???
        WRITE(74,IOSTAT=ios,FMT=frm) (num(j,4),j=k,MIN(k+259,w)) 
        IF (io_err("write7",msnam)) STOP
      ENDDO
      WRITE(74,IOSTAT=ios,FMT= &
        '(I8,"=N",I8,"=I "," RCS",I4,42X,A6)') w,1,n,frm
      IF (io_err("write4",msnam)) STOP
      DO k=1,w,260   ! Some buffer limit encountered???
        WRITE(74,IOSTAT=ios,FMT=frm) (num(j,2),j=k,MIN(k+259,w)) 
        IF (io_err("write5",msnam)) STOP
      ENDDO
      RETURN 
      END SUBROUTINE Boot_write
!------------------------------------------------------------------------
      SUBROUTINE Boot_crud()  ! Plots Bootstrap mean and SD
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka 
      REAL(8),DIMENSION(4)   :: boxx,boxy 
      REAL(8)                :: ra
      INTEGER                :: i=1,p,q,r
      grl=200 ; grr=2400 ; grt=150 ; grb=650
      r=cyr(2) ; wka(1:r)=(/(DBLE(i),i=1,r)/) 
      p=cfy(2) ; ra=MAXVAL(crn(1:r,6))
      CALL NAME('Ring Age','X')     ! Axis name
      CALL NAME('Ring Width','Y')   ! Axis name
      CALL plot_trees(r,num(1:r,2))  
      CALL tombox(0,r+1,0.D0,ra)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(4) ; CALL HEIGHT(22)
      CALL SETCLR(cyan)
      DO i=p,r
         boxx=(/wka(i)-0.5D0,wka(i)+0.5D0,wka(i)+0.5D0,wka(i)-0.5D0/)
         boxy=(/crn(i,6),crn(i,6),crn(i,5),crn(i,5)/)
         CALL RLAREA(boxx,boxy,4)
      ENDDO
      CALL LINWID(1) ; CALL SETCLR(red)
      CALL CURVE(wka(p:r),crn(p:r,2),r-p+1)
      CALL SETCLR(black)
      CALL MESSAG(wnam(cf),grl+600,grt+30)
      CALL ENDGRF() ; CALL LINWID(1)
      grt=800 ; grb=1300
      r=cyr(1) ; wka(1:r)=(/(DBLE(i),i=cfy(1),cly(1))/) 
      p=1 ; q=r ; r=q-p+1     
      ra=MAXVAL(crn(1:r,4))
      CALL NAME('Calendar Year','X')
      CALL NAME('Index','Y')   ! Axis name
      CALL plot_trees(r,num(p:q,1))  
      CALL tombox(cfy(1),cly(1),-0.5D0,ra)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! Gridlines
      CALL LINWID(4) ; CALL HEIGHT(22)
      CALL SETCLR(cyan)
      DO i=p,q
         boxx=(/wka(i)-0.5D0,wka(i)+0.5D0,wka(i)+0.5D0,wka(i)-0.5D0/)
         boxy=(/crn(i,4),crn(i,4),crn(i,3),crn(i,3)/)
         CALL RLAREA(boxx,boxy,4)
      ENDDO
      CALL LINWID(1) ; CALL SETCLR(red)
      CALL CURVE(wka(p:q),crn(p:q,1),r)
      CALL SETCLR(black)
      CALL MESSAG(wnam(cf),grl+600,grt+30)
      CALL ENDGRF() ; CALL LINWID(1)
      RETURN 
      END SUBROUTINE Boot_crud
!-------------------------------------------------------------------
      SUBROUTINE Boot_cru() ! Bootstrap RCS 
      IMPLICIT NONE  ! Simple one-curve RCS (not signal-free)
      CHARACTER(60)           :: msnam
      INTEGER                 :: n,w1,w2,st
      st=cf ; n=LEN_TRIM(cnam(cf))
      msnam=cnam(cf)(1:n)//"bootc"
      OPEN(74,FILE=msnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",msnam)) RETURN
      nc=0 ; CALL read_rft(cnam(cf))        ! Read chronology
      cf=1 ; CALL det_crnfy() ; w1=cyr(1) 
      w2=MAXVAL(ly(1:nc)-pth(1:nc)+1)       ! Max RCS years
      cfy(2)=MINVAL(fy(1:nc)-pth(1:nc)+1)   ! Min RCS years
      cyr(2)=w2 ; cly(2)=w2 ; cf=st
      CALL Boot_cru1(w1,w2,msnam) 
      CLOSE(74)
      RETURN 
      END SUBROUTINE Boot_cru
!------------------------------------------------------------------------
      SUBROUTINE Boot_cru1(w1,w2,msnam) 
      IMPLICIT NONE  ! Bootstrap core selection with replacement
      INTEGER,INTENT(IN)           :: w1     ! Chronology length
      INTEGER,INTENT(IN)           :: w2     ! RCS curve length
      REAL(8),DIMENSION(1:w1,25,2) :: bcrn   ! High and low chronology values
      REAL(8),DIMENSION(1:w2,25,2) :: brcs   ! High and low RCS values 
      INTEGER,DIMENSION(1:w1)      :: bcno   ! Chronology count
      INTEGER,DIMENSION(1:w2)      :: brno   ! RCS count 
      REAL(8),DIMENSION(nc)        :: rrn
      INTEGER,DIMENSION(nc)        :: t1
      INTEGER,DIMENSION(12)        :: seed
      CHARACTER(60),INTENT(IN)     :: msnam
      REAL(8)                      :: rnc
      INTEGER                      :: i,j,n,p,q,r,s,t,u,v
      seed=12848 ; rnc=DBLE(nc)
      CALL RANDOM_SEED(put=seed) ! Initialise random sequence
      bcrn(1:w1,1:25,1)=-10.D0 ; bcrn(1:w1,1:25,2)=10.D0
      brcs(1:w2,1:25,1)=-10.D0 ; brcs(1:w2,1:25,2)=10.D0
      crn(1:w1,8)=0.D0 ; crn(1:w2,9)=0.D0
      num(1:w1,8)=0    ; num(1:w2,9)=0
      bcno=0 ; brno=0 
      DO n=1,1001        !  1000 chronologies with replacement
        IF (n.LT.1001) THEN 
          CALL RANDOM_NUMBER(rrn)     ! Random Uniform
          t1=INT(rrn*rnc)+1
        ELSE
          t1(1:nc)=(/(i,i=1,nc)/)  ! Finallly use all trees
        ENDIF 
        crn(1:w1,1)=0.D0 ; num(1:w1,1)=0 
        crn(1:w2,2)=0.D0 ; num(1:w2,2)=0 
        DO i=1,nc     ! Accumulate RCS curve values
          j=t1(i) ; r=yr(j) ; p=ad(j) ; q=p+r-1
          u=fy(j)-pth(j)+1 ; v=u+r-1
          WHERE (xok(p:q))
            crn(u:v,2)=crn(u:v,2)+x(p:q)
            num(u:v,2)=num(u:v,2)+1
          END WHERE
        ENDDO
        WHERE (num(1:w2,2).GE.1) & ! Mean RCS curve
          crn(1:w2,2)=crn(1:w2,2)/DBLE(num(1:w2,2))  
        PD: DO p=1,w2 ; IF (num(p,2).GE.1) EXIT PD ; ENDDO PD
        CALL spline3(w2-p+1,crn(p:w2,2),num(p:w2,2),4,crn(p:w2,2),FA)  
        DO i=1,nc       ! Accumulate CRN values
          j=t1(i) ; r=yr(j) ; p=ad(j) ; q=p+r-1   ! Ring address
          u=fy(j)-cfy(1)+1 ; v=u+r-1              ! CRN address 
          s=fy(j)-pth(j)+1 ; t=s+r-1              ! RCS address
          WHERE (xok(p:q))
            crn(u:v,1)=crn(u:v,1)+x(p:q)/crn(s:t,2)
            num(u:v,1)=num(u:v,1)+1
          END WHERE
        ENDDO
        WHERE (num(1:w1,1).GT.0)   ! Mean CRN
          crn(1:w1,1)=crn(1:w1,1)/DBLE(num(1:w1,1)) 
        ELSEWHERE
          crn(1:w1,1)=0.D0
        END WHERE
        num(1:w1,3)=NINT(1000.D0*crn(1:w1,1)) 
        CALL Boot_write(msnam,n,cfy(1),w1,w2)    ! Save bootstrap CRNs
        IF (n.LT.1001) THEN 
          WHERE (num(1:w1,1).GT.0)   ! Sum CRNs
            crn(1:w1,8)=crn(1:w1,8)+crn(1:w1,1)
            num(1:w1,8)=num(1:w1,8)+1
          END WHERE
          WHERE (num(1:w2,2).GT.0)  ! Sum RCS curves
            crn(1:w2,9)=crn(1:w2,9)+crn(1:w2,2)
            num(1:w2,9)=num(1:w2,9)+1
          END WHERE
          DO i=1,w1      ! CRN values in each year
            IF (num(i,1).GE.1) THEN
              bcno(i)=bcno(i)+1       ! Count of values
              IF (crn(i,1).GT.bcrn(i,1,1)) THEN  ! High 25
                J1: DO j=2,25
                  IF (crn(i,1).LT.bcrn(i,j,1)) EXIT J1
                ENDDO J1
                IF (j.GT.2) bcrn(i,1:j-2,1)=bcrn(i,2:j-1,1) 
                bcrn(i,j-1,1)=crn(i,1)
              ENDIF
              IF (crn(i,1).LT.bcrn(i,1,2)) THEN ! Low 25
                J2: DO j=2,25
                  IF (crn(i,1).GT.bcrn(i,j,2)) EXIT J2
                ENDDO J2
                IF (j.GT.2) bcrn(i,1:j-2,2)=bcrn(i,2:j-1,2) 
                bcrn(i,j-1,2)=crn(i,1)
              ENDIF
            ENDIF
          ENDDO
          DO i=1,w2      ! RCS values in each year
            IF (num(i,2).GE.1) THEN
              brno(i)=brno(i)+1     ! Count of values
              IF (crn(i,2).GT.brcs(i,1,1)) THEN  ! High 25
                J3: DO j=2,25
                  IF (crn(i,2).LT.brcs(i,j,1)) EXIT J3
                ENDDO J3
                IF (j.GT.2) brcs(i,1:j-2,1)=brcs(i,2:j-1,1) 
                brcs(i,j-1,1)=crn(i,2)
              ENDIF
              IF (crn(i,2).LT.brcs(i,1,2)) THEN ! Low 25
                J4: DO j=2,25
                  IF (crn(i,2).GT.brcs(i,j,2)) EXIT J4
                ENDDO J4
                IF (j.GT.2) brcs(i,1:j-2,2)=brcs(i,2:j-1,2) 
                brcs(i,j-1,2)=crn(i,2)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO i=1,w1
        IF (bcno(i).EQ.0) THEN
          crn(i,3:4)=crn(i,1)
        ELSE
          u=NINT(DBLE(bcno(i))/40.D0)
          u=MAX(u,1) ; u=MIN(25,u) ; u=26-u
          crn(i,3)=bcrn(i,u,2) ; crn(i,4)=bcrn(i,u,1)
        ENDIF
      ENDDO 
      DO i=1,w2
        IF (brno(i).EQ.0) THEN
          crn(i,5:6)=crn(i,2)
        ELSE
          u=NINT(DBLE(brno(i))/40.D0)
          u=MAX(u,1) ; u=MIN(25,u) ; u=26-u
          crn(i,5)=brcs(i,u,2) ; crn(i,6)=brcs(i,u,1)
        ENDIF
      ENDDO 
      WHERE (num(1:w1,8).GT.0) &  ! Mean CRNs
        crn(1:w1,8)=crn(1:w1,8)/DBLE(num(1:w1,8))
      WHERE (num(1:w2,9).GT.0) &  ! Mean RCS curves
        crn(1:w2,9)=crn(1:w2,9)/DBLE(num(1:w2,9))
      OPEN(77,FILE="Saved.bootc",IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open","Saved.bootc")) RETURN
      WRITE(77,'("Ring Age")')
      WRITE(77,'("Bootstrap RCS Counts")')
      WRITE(77,'("All Trees RCS Value")')
      WRITE(77,'("Mean Bootstrap RCS Value")')
      WRITE(77,'("2.5% Low RCS Value")')
      WRITE(77,'("2.5% High RCS Value")')
      DO i=1,w2       ! Plot RCS values
        WRITE (77,'(2I5,4F8.3)') &
          i,brno(i),crn(i,2),crn(i,9),crn(i,5:6)
      ENDDO
      WRITE(77,'("Calendar year")')
      WRITE(77,'("Bootstrap Chronology Counts")')
      WRITE(77,'("All Trees Chronology Value")')
      WRITE(77,'("Mean Bootstrap Chronology Value")')
      WRITE(77,'("2.5% Low Chronology Value")')
      WRITE(77,'("2.5% High Chronology Value")')
      DO i=1,w1    ! Plot CRN values
        WRITE (77,'(2I5,4F8.3)') &
          cfy(1)-1+i,bcno(i),crn(i,1),crn(i,8),crn(i,3:4)
      ENDDO
      CLOSE(77)
      RETURN 
      END SUBROUTINE Boot_cru1
!------------------------------------------------------------------------
      SUBROUTINE read_sites(sit)   ! Reads site groups
      IMPLICIT NONE
      INTEGER,DIMENSION(0:14) :: sit  ! Separate sites
      INTEGER                 :: i=1,j=1,k
      DO i=1,nfil
        k=LEN_TRIM(cnam(i))
        JD: DO j=k,1,-1
          IF (cnam(i)(j:j).EQ."/".OR.cnam(i)(j:j).EQ."\") EXIT JD
        ENDDO JD
        j=MAX(1,j+1)   ! 1st 3 chars of filename
        wnam(20+i)=cnam(i)(j:j+3)
      ENDDO 
      nc=0 ; sit(0)=0
      DO j=1,nfil
        CALL read_rft(cnam(j)) ; sit(j)=nc
      ENDDO  
      RETURN
      END SUBROUTINE read_sites
!-------------------------------------------------------------------
      SUBROUTINE mrcs_sited(ref1)  ! Plots Best Fit
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)     :: ref1
      REAL(8),DIMENSION(mxy) :: wka 
      INTEGER                :: i=1,p,q,r
      REAL(8)                :: ra,rb
      r=sly(mx) ; wka(1:r)=(/(DBLE(i),i=1,r)/) 
      grl=200 ; grr=2400 ; grt=200 ; grb=800
      CALL NAME('Mean Ring','Y') ! Axis name
      CALL NAME('Ring Age','X') ! Axis name
      CALL plot_trees(r,mcnt(1:r,mx))  
      IF (ref1.EQ.1) THEN
        ra=MAXVAL(crn(1:r,3+nfil:2+2*nfil))
      ELSE 
        ra=MAXVAL(mval(1:r,1:4))
      ENDIF
      CALL tombox(1,r,0.D0,ra)
      CALL SETCLR(grey) ; CALL GRID(1,1) ! GRIDLINES
      CALL LINWID(1) ; CALL HEIGHT(22)
      IF (ref1.EQ.1) THEN
        DO i=3+nfil,2+2*nfil
          IF (i-2-nfil.NE.7) THEN  ! Avoid silver
            CALL SETCLR(i-2-nfil) 
          ELSE
            CALL SETCLR(15) 
          ENDIF 
          CALL MESSAG(wnam(i+18),(i-2-nfil)*150+50,grt-45)
          p=cfy(i) ; q=MIN(cly(i),2500) ; r=q-p+1
          CALL thickthin(r,wka(p:q),crn(p:q,i),num(p:q,i),6)
        ENDDO
      ELSE
        DO i=1,4
          CALL SETCLR(i) 
          p=sfy(i) ; q=MIN(sly(i),2500) ; r=q-p+1
          CALL thickthin(r,wka(p:q),mval(p:q,i),mcnt(p:q,i),6)
        ENDDO
      ENDIF
      CALL SETCLR(black) 
      CALL MESSAG(cnam(20),grl+500,grt+30)
      CALL ENDGRF()  

      r=cyr(1) ; wka(1:r)=(/(DBLE(i),i=cfy(1),cly(1))/) 
      CALL NAME('Mean Index','Y')    ! Axis name
      CALL NAME('Calendar Year','X') ! Axis name
      grt=960 ; grb=1560
      ra=MINVAL(crn(1:r,3:2+nfil))
      rb=MAXVAL(crn(1:r,3:2+nfil))
      CALL plot_trees(r,num(1:r,1))  
      CALL tombox(cfy(1),cly(1),ra,rb)
      CALL SETCLR(grey) ; CALL GRID(1,1)       ! GRIDLINES
      CALL LINWID(1) ; CALL HEIGHT(22)
      DO i=3,2+nfil
        IF (i-2.NE.7) THEN ! Avoid silver
          CALL SETCLR(i-2) 
        ELSE
          CALL SETCLR(15) 
        ENDIF 
        p=cfy(i)-cfy(1)+1 ; q=cly(i) ; r=q-p+1
        CALL thickthin(r,wka(p:q),crn(p:q,i),num(p:q,i),6)
        CALL MESSAG(wnam(i+18)(1:4),(i-2)*150+50,grt-45)
      ENDDO
      CALL SETCLR(black) 
      CALL MESSAG(cnam(21),grl+150,grt+30)
      CALL ENDGRF() 
      RETURN 
      END SUBROUTINE mrcs_sited
!-------------------------------------------------------------------
      SUBROUTINE mrcs_site1(sit)  ! Investigate modern sample bias 
      IMPLICIT NONE
      INTEGER,DIMENSION(0:14),INTENT(IN) :: sit    ! Separate sites
      INTEGER :: i,j,p,q,r,s,t,u,v
      cf=1 ; CALL detrend() 
      r=cyr(1) ; crn(1:r,3:50)=0.D0 ; num(1:r,3:50)=0
      cfy(3:50)=mxy ; cly(3:50)=-3000
      DO j=3,2+nfil
        DO i=sit(j-3)+1,sit(j-2)         ! Each tree this site 
          p=ad(i) ; r=yr(i) ; q=p+r-1    ! Ring address
          u=fy(i)-cfy(1)+1 ; v=u+r-1     ! Chronology address
          s=fy(i)-pth(i) ; t=s+r-1       ! RCS address
          cfy(j)=MIN(cfy(j),u)           ! Calendar range
          cly(j)=MAX(cly(j),v)
          cfy(j+nfil)=MIN(cfy(j+nfil),s)   ! Age range
          cly(j+nfil)=MAX(cly(j+nfil),t)
          WHERE (xok(p:q)) 
            crn(u:v,j)=crn(u:v,j)+dx(p:q)         ! Sub chronologies     
            num(u:v,j)=num(u:v,j)+1       
            crn(s:t,j+nfil)=crn(s:t,j+nfil)+fx(p:q) ! Sub RCS Curves     
            num(s:t,j+nfil)=num(s:t,j+nfil)+1       
          END WHERE
        ENDDO 
      ENDDO 
      r=cyr(1) ; j=2+2*nfil ; cyr(3:j)=cly(3:j)-cfy(3:j)+1
      WHERE (num(1:r,3:j).GT.1) &  ! Mean chronologies
        crn(1:r,3:j)=crn(1:r,3:j)/DBLE(num(1:r,3:j))
      DO i=3,j     ! Smooth chronologies and RCS curves
        p=cfy(i) ; q=cly(i) ; r=q-p+1
        CALL splinet(r,crn(p:q,i),50,crn(p:q,i)) 
      ENDDO 
      OPEN(71,FILE="zzz.prn",IOSTAT=ios,STATUS="REPLACE")
      DO i=3,j
        WRITE(71,'(I4,3I6)') i,cfy(i),cly(i),cyr(i)
      ENDDO  
      WRITE(71,'(12X,32(A3,4X))') wnam(21:34)(1:3),wnam(21:34)(1:3)
      DO i=1,r
        WRITE(71,'(2I6,32F7.2)') i,i-1+cfy(1),crn(i,3:j) 
      ENDDO  
      CLOSE(71)
      RETURN
      END SUBROUTINE mrcs_site1
!-------------------------------------------------------------------
      SUBROUTINE mrcs_site(ref1)  ! Individual sites in a group
      IMPLICIT NONE
      INTEGER,INTENT(IN)      :: ref1
      INTEGER,DIMENSION(0:14) :: sit    ! Separate sites
      CHARACTER(12)           :: jnam
      INTEGER                 :: i,p,q,r
      IF (nfil.GT.14) THEN
        CALL out_err("Too many sites - 1st 14 used")
        nfil=14
      ENDIF
      IF (nfil.LT.2) THEN
        CALL out_err("Too few sites - 2 or more needed")
        STOP
      ENDIF
      CALL read_sites(sit)    ! Read raw data
      IF (ref1.EQ.2.AND.nc.LT.200) THEN
        CALL out_err("Too few trees - 200+ for four-curve RCS")
        STOP
      ENDIF
      CALL det_default()
      sfo=2 ; idt=-2 ; idb=2
      IF (ref1.EQ.1) THEN
        src=1                 ! Single RCS curve 
        cnam(20)="Mean Growth by Age for each site"
        cnam(21)="One-curve normal RCS, site means"
        jnam="fig35.col"
      ELSE
        src=2 ; srcno=4       ! 4 RCS curves
        cnam(20)="Mean Growth by Age for 4 Growth Rates"
        cnam(21)="Four-curve normal RCS, site means"
        jnam="fig36.col"
      ENDIF
      CALL mrcs_site1(sit)
      IF (ref1.EQ.2) THEN
        DO i=1,4   ! Smooth mean ring
          p=sfy(i) ; q=sly(i) ; r=q-p+1
          CALL splinet(r,mval(p:q,i),50,mval(p:q,i)) 
        ENDDO
      ENDIF

      OPEN(19,FILE=jnam,IOSTAT=ios,STATUS="REPLACE")
      r=sly(mx) 
      IF (ref1.EQ.1) THEN
        WRITE(19,'("RCS curve counts")') 
        WRITE(19,'("  Age",15(3X,A4))') wnam(21:20+nfil)(1:4)
        DO i=1,r
          WRITE(19,'(I5,15I7)') i,num(i,3+nfil:2+2*nfil)
        ENDDO
        WRITE(19,*)
        WRITE(19,'("RCS curve values")') 
        WRITE(19,'("  Age",15(3X,A4))') wnam(21:20+nfil)(1:4)
        DO i=1,r
          WRITE(19,'(I5,15F7.3)') i,crn(i,3+nfil:2+2*nfil)
        ENDDO
      ELSE
        WRITE(19,'(9X,"RCS curve counts",9X,"RCS curve values")') 
        WRITE(19,'()') 
        WRITE(19,'(5X,A45)') "  1st  2nd  3rd  4th       1st  2nd  3rd  4th"
        DO i=1,r
          WRITE(19,'(5I5,4F7.3)') i,mcnt(i,1:4),mval(i,1:4)
        ENDDO
      ENDIF
      r=cyr(1) 
      WRITE(19,*)
      WRITE(19,'("Chronology counts")') 
      WRITE(19,'("  Year",15(3X,A4))') wnam(21:20+nfil)(1:4)
      DO i=1,r
        WRITE(19,'(I6,15I7)') i-1+cfy(1),num(i,3:2+nfil)
      ENDDO
      WRITE(19,*)
      WRITE(19,'("Chronology values")') 
      WRITE(19,'("  Year",15(3X,A4))') wnam(21:20+nfil)(1:4)
      DO i=1,r
        WRITE(19,'(I6,15F7.3)') i-1+cfy(1),crn(i,3:2+nfil)
      ENDDO
      CLOSE(19)
      RETURN
      END SUBROUTINE mrcs_site
!-------------------------------------------------------------------
      END MODULE rfigs 
                          

