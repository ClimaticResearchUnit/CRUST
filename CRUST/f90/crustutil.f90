! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see 
! the GNU General Public License.
      MODULE CRUSTUTIL    ! Standardisation procs
      USE STAND ; USE detdata
      CONTAINS   
!--------------------------------------------------------
      SUBROUTINE TOM_help(menno) ! Displays help text from RCShelp.prn 
      IMPLICIT NONE
      INTEGER,INTENT(IN)    :: menno   ! Menu calling the help routine
      INTEGER               :: hno     ! Number of messages in this menu
      INTEGER,DIMENSION(60) :: hm      ! List of message numbers menu
      INTEGER               :: bf,bl   ! First and last buttons this menu
      INTEGER               :: i
      SELECT CASE (menno)
      CASE (1)      ! RCS menu
        hno=9 ; bf=1 ; bl=58 ; hm(1:9)=(/1,2,7,20,21,22,23,24,41/)
      CASE (2)      ! Tree menu
        hno=5 ; bf=55 ; bl=61 ; hm(1:5)=(/7,8,9,55,56/)  
      CASE (3)      ! Disp menu
        hno=44 ; bf=66 ; bl=193
        hm(1:5)=(/26,27,28,29,30/)
        hm(6:26)=(/(i,i=67,87)/)
        hm(27:38)=(/101,102,103,105,106,107,109,110,111,113,117,121/) 
        hm(39:50)=(/125,129,133,137,141,145,149,153,157,161,165,169/) 
        hm(51:58)=(/173,177,181,185,190,191,192,193/) 
      CASE (4)      ! Save menu
        hno=34 ; bf=199 ; bl=224
        hm(1:10)=(/(i,i=10,19)/) ; hm(11:21)=(/(i,i=199,209)/)
        hm(22:35)=(/(i,i=211,224)/)
      CASE (6)      ! Name menu
        hno=3 ; bf=264 ; bl=268 ; hm(1:3)=(/4,5,268/) 
      CASE (7)      ! Figures menu
        hno=2 ; bf=271 ; bl=273 ; hm(1:2)=(/271,273/)
      CASE (8)      ! Chron menu
        hno=6 ; bf=340 ; bl=359 ; hm(1:6)=(/(i,i=45,50)/)
      CASE (9)      ! Met menu
        hno=0 ; bf=365 ; bl=384 
      CASE (10)     ! Xcor menu
        hno=0 ; bf=391 ; bl=406 
      CASE (13)      ! Start menu
        hno=17 ; bf=451 ; bl=465 ; hm(1:17)=(/(i,i=303,319)/) 
!      CASE (8)      ! Name 3 menu
!        hno=2 ; bf=243 ; bl=257 ; hm(1:2)=(/254,255/)  
!      CASE (21)      ! Start menu
!        hno=16 ; bf=202 ; bl=214 ; hm(1:16)=(/(i,i=3,18)/)
!      CASE (22)      ! Tree menu
!        hno=5 ; bf=1 ; bl=42 ; hm(1:5)=(/54,55,56,57,58/)  
!      CASE (23)      ! Name menu
!        hno=3 ; bf=231 ; bl=236 ; hm(1:3)=(/97,98,99/) 
!      CASE (24)      ! Disp menu
!        hno=6 ; bf=60 ; bl=76 ; hm(1:6)=(/(i,i=45,50)/)
!      CASE (25)      ! Save menu
!        hno=36 ; bf=239 ; bl=264 ; hm(1:36)=(/(i,i=100,135)/)
      ENDSELECT
      CALL prog_help(hno,hm(1:hno),bf,bl)
      RETURN
      END SUBROUTINE TOM_help
!--------------------------------------------------------------
      SUBROUTINE det_default()   ! Default detrend parameters
      IMPLICIT NONE           
      IF ((idt.GE.-9.AND.idt.LE.-4).OR.idt.EQ.-1.OR. &
           idt.EQ.8.OR.idt.EQ.9) idt=-2
      itn=1 ; rdt=1 ; ind=1 ; krb=1 ; isb=1 ; sfo=2 ; poo=1
      src=1 ; trc=1 ; gtr=1 ; tst=4 ; bfc=1 ; idb=1 ; jrb=1
      srcno=1 ; sfono=10 ; rdtno=60 ; ignor=0
      cft=1 ; rft=1
      RETURN
      END SUBROUTINE det_default
!-------------------------------------------------------------------- 
      SUBROUTINE read_indsav(msnam) ! Reads saved index files (compact)
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN) :: msnam
      CHARACTER(12)           :: frm
      INTEGER                 :: i,j
      REAL(8),DIMENSION(mxy)  :: wka
      OPEN(24,FILE=msnam,IOSTAT=ios,STATUS="OLD")
      IF (io_err("Open",msnam)) RETURN
      ID: DO i=1,mxc                ! For each file
        READ(24,'(I4,2I6,2X,A20)',IOSTAT=ios) &
          icf,cfy(i),cly(i),wnam(i)
        IF (io_err("Read1",msnam)) STOP
        IF (icf.NE.i) EXIT ID
      ENDDO ID
      WRITE(frm,'("(I6,",I3,"F7.3)")')  icf
      DO i=MINVAL(cfy(1:icf)),MAXVAL(cly(1:icf))
        READ(24,IOSTAT=ios,FMT=frm) j,wka(1:icf) 
        IF (io_err("Read2",msnam)) STOP
        DO j=1,icf
          IF (cfy(j).LE.i.AND.cly(j).GE.i) crn(i-cfy(j)+1,j)=wka(j)
        ENDDO
      ENDDO
      CLOSE(24) ; chc(1)=1 ; chs=1 ; cyr(1:icf)=cly(1:icf)-cfy(1:icf)+1
      okc=(crn.NE.-9.99D0)
      RETURN
      END SUBROUTINE read_indsav
!-----------------------------------------------------------------
      SUBROUTINE write_indsav(msnam) ! Sub saves index files (compact)
      IMPLICIT NONE  
      CHARACTER(*),INTENT(IN) :: msnam
      INTEGER                 :: i,j,m,p,q
      CHARACTER(12)           :: frm
      REAL(8),DIMENSION(mxy)  :: wka
      OPEN(24,FILE=msnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",msnam)) RETURN
      p=100000 ; q=-100000
      DO i=1,chs                  ! For each file
        m=chc(i) ; p=MIN(cfy(m),p) ; q=MAX(cly(m),q)
        WRITE(24,'(I4,2I6,2X,A20)',IOSTAT=ios) i,cfy(m),cly(m),wnam(m)
        IF (io_err("write1",msnam)) STOP
      ENDDO
      WRITE(24,'(I4,2I6,2X,A20)',IOSTAT=ios) chs,p,q,"  "
      WRITE(frm,'("(I6,",I3,"F7.3)")')  chs
      DO i=p,q
        DO j=1,chs
          m=chc(j) 
          namc=wnam(m)
          IF (cfy(m).LE.i.AND.cly(m).GE.i) THEN
            wka(j)=crn(i-cfy(m)+1,m)
          ELSE 
            wka(j)=-9.999D0
          ENDIF
        ENDDO
        WRITE(24,IOSTAT=ios,FMT=frm) i,wka(1:chs) 
        IF (io_err("write2",msnam)) STOP
      ENDDO
      CLOSE(24)
      RETURN
      END SUBROUTINE write_indsav
!-----------------------------------------------------------------
      SUBROUTINE read_index(fnam) ! Sub read all the index files
      IMPLICIT NONE            
      CHARACTER(*),INTENT(IN) :: fnam
      INTEGER,DIMENSION(21)   :: wki
      INTEGER                 :: j
      CALL read_open(26,fnam)
      READ(26,IOSTAT=ios,FMT='(A85)') mess     ! Read ahead
      IF (io_err("Read3",fnam)) STOP
      ID: DO icf=icf+1,mxc
        IF (ios.EQ.-1) EXIT ID                 ! End of file
        wnam(icf)=mess(21:40)                  ! 1st header has file name 
        JD: DO j=1,3                           ! Optional header lines
          READ(26,'(A85)',IOSTAT=ios) mess     ! Read next 
          IF (ios.EQ.-1) EXIT ID
          IF (io_err("Read3",fnam)) STOP
          READ(mess,'(6X,I4,10(I4,I3))',IOSTAT=ios) wki(1:21) ! Is it valid data
          IF (ios.EQ.0) EXIT JD                ! Data line is valid
        ENDDO JD
        CALL read_ind()
      ENDDO ID
      IF (icf.GT.mxc) THEN
        WRITE(mess,'(I5)') mxc
        CALL out_err("Max CRN limit reached "//mess(1:5))
      ENDIF
      CLOSE(26) ; chs=1 ; chc(1)=1 ; icf=icf-1
      RETURN
      END SUBROUTINE read_index
!---------------------------------------------------------------------------
      SUBROUTINE read_ind()  ! Reads an index file to crn(icf)
      IMPLICIT NONE              
      TYPE pair ; INTEGER :: val,num ; END TYPE
      TYPE(pair),DIMENSION(0:10) :: rec1
      INTEGER                    :: j,p,q,r
      CHARACTER(1)               :: sign1
      READ(mess,IOSTAT=ios,FMT='(5X,A,I4,10(I4,I3))') &  ! First line
        sign1,j,rec1(0:9) 
      IF (io_err("Read1","MESS")) STOP
      IF (sign1.EQ."-") j=0-j                    ! Find first value
      BD: DO r=0,9 ; IF (rec1(r)%val.NE.9990) EXIT BD ; ENDDO BD
      IF (MOD(j,10).EQ.0) j=j+r
      cfy(icf)=j ; p=1 ; q=10-r 
      READ(26,IOSTAT=ios,FMT='(A85)') mess      ! Read ahead
      IF (io_err("Read2",wnam(icf))) STOP
      JD: DO j=j+q,j+mxy,10                     ! For each decade 
        crn(p:q,icf)=DBLE(rec1(r:9)%val)/1000.D0
        num(p:q,icf)=rec1(r:9)%num ; p=q+1 ; q=p+9 ; r=0
        READ(mess,IOSTAT=ios,FMT='(10X,10(I4,I3))') rec1(0:9) 
        IF (io_err("Read3",wnam(icf))) STOP
        READ(26,IOSTAT=ios,FMT='(A85)') mess 
        IF (ios.EQ.-1) EXIT JD                    ! End of file
        IF (io_err("Read4",wnam(icf))) STOP
        WRITE(mess(89:93),'(I5)') j+10    ! Next decade
        IF (mess(89:89).EQ."-") THEN
          IF (mess(6:10).NE.mess(89:93)) EXIT JD  ! Chronology ends
        ELSE
          IF (mess(7:10).NE.mess(90:93)) EXIT JD  ! Chronology ends
        ENDIF
      ENDDO JD                           
      RD: DO r=8,0,-1 ; IF(.NOT.rec1(r)%val.EQ.9990) EXIT RD ; ENDDO RD
      q=p+r ; crn(p:q,icf)=DBLE(rec1(0:r)%val)/1000.D0
      num(p:q,icf)=rec1(0:r)%num ; cly(icf)=cfy(icf)+q-1
      cyr(icf)=q 
      okc(1:q,icf)=num(1:q,icf).GE.1
      RETURN
      END SUBROUTINE read_ind
!-----------------------------------------------------------------
      SUBROUTINE mouse_click(menu,range1,range2) 
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: range1,range2,menu  ! Button range
      INTEGER            :: i,but
      LOGICAL            :: helpb        ! In help state
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Click Mouse")
      CALL SENDBF()  ! Update display
      but=0 ; helpb=FA
      ADO: DO
        IF (but.EQ.5) EXIT ADO             ! Left mouse click - action
        IF (but.EQ.6) CALL TOM_help(menu)  ! Right mouse click - help
        CALL CSRPOS(msx,msy,but)           ! Cursor XY 
      ENDDO ADO
      mous=0
      BDO: DO i=range1,range2        ! Is mouse in a button box
        IF ((b(i)%ok).AND.(msx.GE.b(i)%x1).AND.(msx.LE.b(i)%x2) & 
           .AND.(msy.GE.b(i)%y1).AND.(msy.LE.b(i)%y2)) THEN
           mous=i ; EXIT BDO
        ENDIF
      ENDDO BDO
      CALL mwrite(stabut) ; CALL but_draw(stabut,"Running")
      CALL SENDBF()  ! Update display
      RETURN
      END SUBROUTINE mouse_click
!-----------------------------------------------------------
      SUBROUTINE but_draw2(i)   ! Palette selection buttons
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i   ! Button number
      CALL SETCLR(lcol(lco))    ! Current colour
      IF (i.EQ.191) THEN        ! Line width
        CALL LINWID(lth) 
        CALL LINE(b(i)%x1,b(i)%y1+22,b(i)%x2,b(i)%y1+22)
        CALL LINWID(1)
      ELSEIF (i.EQ.192) THEN    ! Line style
        CALL LINTYP(lst)
        CALL LINE(b(i)%x1,b(i)%y1+22,b(i)%x2,b(i)%y1+22)
        CALL LINTYP(0)
      ELSEIF (i.EQ.193) THEN    ! Colour selection
        CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                   (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
      ENDIF
      CALL SETCLR(red)          ! Red outline
      CALL SHDPAT(0)            ! Shading pattern
      CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                 (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
      CALL SHDPAT(16)           ! Shading pattern
      CALL SETCLR(black)
      RETURN 
      END SUBROUTINE but_draw2
!--------------------------------------------------------
      SUBROUTINE but_draw3(i)   ! Palette selection buttons
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i   ! Button number
      CALL SETCLR(lins(i-86)%co) ! Current colour
      CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                 (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
      CALL SETCLR(red)           ! Red outline
      CALL SHDPAT(0)             ! Empty box
      CALL AREAF((/b(i)%x1,b(i)%x2,b(i)%x2,b(i)%x1/), &
                 (/b(i)%y1,b(i)%y1,b(i)%y2,b(i)%y2/),4)
      CALL SHDPAT(16)            ! Full box
      CALL SETCLR(black)
      IF (i.GE.101.AND.lins(i-86)%on) THEN               
        CALL MESSAG("On",b(i)%x1+5,b(i)%y1+12)
      ELSEIF (i.GE.101) THEN
        CALL MESSAG("Off",b(i)%x1+5,b(i)%y1+12)
      ENDIF
      RETURN 
      END SUBROUTINE but_draw3
!--------------------------------------------------------
      SUBROUTINE line_miss(n,xval,yval,zok)  ! Line with missing values
      IMPLICIT NONE
      INTEGER,INTENT(IN)                :: n     ! Data length
      REAL(8),INTENT(IN),DIMENSION(1:n) :: xval  ! Data x value
      REAL(8),INTENT(IN),DIMENSION(1:n) :: yval  ! Data y value
      LOGICAL,INTENT(IN),DIMENSION(1:n) :: zok   ! Data valid
      INTEGER                           :: i
      IF (ALL(zok)) THEN
        CALL CURVE(xval,yval,n)           ! Line
      ELSE
        DO i=2,n                                 ! Broken line
          IF (zok(i-1).AND.zok(i)) THEN
            CALL RLINE(xval(i-1),yval(i-1),xval(i),yval(i))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE line_miss
!----------------------------------------------------------------
      SUBROUTINE plot_line(n,xval,yval,zok,lin,off) ! Plots lines
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)              :: n      ! Number of years
      REAL(8),DIMENSION(n),INTENT(IN) :: xval   ! Year numbers
      REAL(8),DIMENSION(n),INTENT(IN) :: yval   ! Values
      LOGICAL,DIMENSION(n),INTENT(IN) :: zok    ! Missing values
      INTEGER,INTENT(IN)              :: lin    ! Line reference
      INTEGER,INTENT(IN)              :: off    ! % offset
      IF (lins(lin)%on) THEN
        CALL LINWID(lins(lin)%th)               ! Line width
        CALL LINTYP(lins(lin)%st)               ! Line style
        CALL SETCLR(lins(lin)%co)               ! Line colour
        CALL line_miss(n,xval,yval,zok)         ! Plot line
        IF (off.GE.0) &
          CALL MESSAG(lins(lin)%lab,grl+5+((grr-grl)*off)/100,grt-40)
        CALL LINWID(1)
        CALL LINTYP(0)
      ENDIF
      RETURN
      END SUBROUTINE plot_line
!-----------------------------------------------------------------
      SUBROUTINE xticks(pfy,ply,pft,pst)   ! Scale values
      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: pfy   ! Plot first year
      INTEGER,INTENT(IN)  :: ply   ! Plot last year
      INTEGER,INTENT(OUT) :: pft   ! First mark
      INTEGER,INTENT(OUT) :: pst   ! Space between marks
      SELECT CASE (ply-pfy)
        CASE (    :  10) ; pst=   1 
        CASE (  11:  50) ; pst=   5 
        CASE (  51: 100) ; pst=  10 
        CASE ( 101: 200) ; pst=  20 
        CASE ( 201: 500) ; pst=  50 
        CASE ( 501:1000) ; pst= 100 
        CASE (1001:2000) ; pst= 200 
        CASE (2001:5000) ; pst= 500 
        CASE (5001:    ) ; pst=1000 
      ENDSELECT
      IF (pfy.LT.0) THEN
        pft=(pfy/pst)*pst
      ELSE
        pft=(pfy/pst+1)*pst
      ENDIF
      RETURN 
      END SUBROUTINE xticks
!--------------------------------------------------------
      SUBROUTINE yticks(bot,top,tk,bk)   ! Scale values
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: bot  ! Plot base
      REAL(8),INTENT(IN)  :: top  ! Plot top
      REAL(8),INTENT(OUT) :: bk   ! First mark
      REAL(8),INTENT(OUT) :: tk   ! Space between marks
      SELECT CASE (NINT(100.D0*(top-bot)))
        CASE (     :   50)   ; bk= 0.1D0 
        CASE (   51:  100)   ; bk= 0.2D0 
        CASE (  101:  250)   ; bk= 0.5D0 
        CASE (  251:  600)   ; bk= 1.D0 
        CASE (  601: 1000)   ; bk= 2.D0 
        CASE ( 1001: 2000)   ; bk= 5.D0 
        CASE ( 2001: 5000)   ; bk=10.D0 
        CASE ( 5001:10000)   ; bk=20.D0 
        CASE (10001:20000)   ; bk=50.D0 
        CASE (20001:50000)   ; bk=100.D0 
        CASE (50001:100000)  ; bk=200.D0 
        CASE (100001:200000) ; bk= 500.D0 
        CASE (200001:500000) ; bk=1000.D0 
        CASE (500001:      ) ; bk=2000.D0 
      ENDSELECT
      IF (bot.LT.-epsi) THEN
        tk=DBLE(INT(bot/bk))*bk
      ELSE
        tk=DBLE(INT(bot/bk)+1)*bk
      ENDIF
      RETURN 
      END SUBROUTINE yticks
!--------------------------------------------------------
      SUBROUTINE tombox(gfy,gly,bot,top)        ! Plots display box
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: gfy,gly             ! X range
      REAL(8),INTENT(IN) :: bot,top             ! Y range
      INTEGER :: i,j
      REAL(8) :: tk,bk
      CALL AXSPOS(grl,grb)
      CALL AXSLEN(grr-grl,grb-grt)
      CALL HEIGHT(18) 
      CALL SETCLR(black)
      CALL xticks(gfy,gly,i,j)
      CALL yticks(bot,top,tk,bk)
      CALL GRAF(DBLE(gfy),DBLE(gly),DBLE(i),DBLE(j),bot,top,tk,bk)
      END SUBROUTINE tombox
!-------------------------------------------------------------------
      SUBROUTINE plot_core()      ! Plots individual core
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka   ! Working areas
      INTEGER                :: i=1,p,q,r,u,v,t
      REAL(8)                :: top
      t=tre(cc) ; p=ad(t) ; r=yr(t) ; q=p+r-1 
      IF (poo.EQ.1) THEN
        wka(1:r)=(/(DBLE(i),i=fy(t)-pth(t)+1,ly(t)-pth(t)+1)/)
      ELSE
        wka(1:r)=(/(DBLE(i),i=1,yr(t))/)
      ENDIF
      top=MAX(MAXVAL(tx(p:q),MASK=xok(p:q)), &
              MAXVAL(cx(p:q),MASK=xok(p:q)))
      grl=300 ; grr=2400 ; grt=350 ; grb=950
      CALL NAME('Ring Age','X')   ! Axis name
      CALL NAME('Ring mm','Y')    ! Axis name
      CALL tombox(0,NINT(wka(r)),0.D0,top)
      CALL SETCLR(grey)
      CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22)
      CALL LINWID(4)
      IF (psok) THEN
        CALL SETCLR(black) 
        CALL MESSAG(cnam(cf),grl,grt-70)
      ENDIF
      CALL plot_line(r,wka(1:r),cx(p:q),xok(p:q),12,10)
      IF (src.GE.2) THEN
        WRITE(mess,'("(RCS",I2,", DET",I2,")")') trr(t),trm(t)
        CALL MESSAG(mess(1:14),grl+5+((grr-grl)*30)/100,grt-40)
      ENDIF
      CALL plot_line(r,wka(1:r),tx(p:q),xok(p:q),11,60) 
      CALL LINWID(1)
      CALL ENDGRF()
      CALL SETCLR(black)

      u=fy(t)-xfy+1 ; v=u+r-1
      wka(1:r)=(/(DBLE(i),i=fy(t),ly(t))/)
      top=MAX(MAXVAL(dx(p:q),MASK=xok(p:q)), &
              MAXVAL(xcrn(u:v,mx),MASK=cok(u:v,mx)))
      grl=300 ; grr=2400 ; grt=1150 ; grb=1750
      CALL NAME('Calendar Year','X') ! Axis name
      CALL NAME('Indices','Y')       ! Axis name
      CALL tombox(NINT(wka(1)),NINT(wka(r)),0.D0,top)
      CALL SETCLR(grey)
      CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22)
      CALL LINWID(4)
      CALL plot_line(r,wka(1:r),xcrn(u:v,mx),cok(u:v,mx),14,50)
      CALL plot_line(r,wka(1:r),dx(p:q),xok(p:q),13,2) 
      CALL LINWID(1)
      CALL ENDGRF()
      CALL SETCLR(black) 
      RETURN
      END SUBROUTINE plot_core
!-------------------------------------------------------------------
      SUBROUTINE det_plot()      ! Plots curve-fitting core
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka   ! Working areas
      INTEGER                :: i=1,p,q,r,u,v,t
      REAL(8)                :: top,bot
      t=tre(cc) ; p=ad(t) ; r=yr(t) ; q=p+r-1 
      wka(1:r)=(/(DBLE(i),i=fy(t)-pth(t)+1,ly(t)-pth(t)+1)/)
      top=MAX(MAXVAL(tx(p:q),MASK=xok(p:q)), &
              MAXVAL(cx(p:q),MASK=xok(p:q)))
      grl=500 ; grr=2400 ; grt=350 ; grb=950
      CALL NAME('Ring Age','X')   ! Axis name
      CALL NAME('Ring mm','Y')    ! Axis name
      CALL tombox(NINT(wka(1)),NINT(wka(r)),0.D0,top)
      CALL SETCLR(grey) ; CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22) ; CALL LINWID(1) ; CALL SETCLR(black) 
      u=jdt(t) ; IF (u.GT.8) u=8
      IF (u.LT.-4) u=9
      CALL MESSAG(TRIM(cnam(cf))//", tree "//nam(t)//" "// &
        idtlab(u),grl,grt-90)
      CALL SETCLR(red)
      CALL MESSAG("Measurements",grl,grt-35)
      CALL SETCLR(blue)
      CALL MESSAG("Detrend Curve",grl+(grr-grl)/2,grt-35)
      CALL line_miss(r,wka(1:r),cx(p:q),xok(p:q))
      CALL SETCLR(red)
      CALL line_miss(r,wka(1:r),tx(p:q),xok(p:q)) 
      CALL SETCLR(black) ; CALL ENDGRF()

      u=fy(t)-xfy+1 ; v=u+r-1
      wka(1:r)=(/(DBLE(i),i=fy(t),ly(t))/)
      top=MAX(MAXVAL(dx(p:q),MASK=xok(p:q)), &
              MAXVAL(xcrn(u:v,mx),MASK=cok(u:v,mx)))
      bot=MIN(MINVAL(dx(p:q)),MINVAL(crn(u:v,cf)))  
      IF (bot.GE.-epsi) THEN
        bot=MAX(DBLE(INT((bot*10.D0)))/10.D0,0.001D0)
      ELSE
        bot=DBLE(INT((bot*10.D0)))/10.D0
      ENDIF
      grl=500 ; grr=2400 ; grt=1150 ; grb=1750
      CALL NAME('Calendar Year','X') ! Axis name
      CALL NAME('Indices','Y')       ! Axis name
      CALL tombox(NINT(wka(1)),NINT(wka(r)),bot,top)
      CALL SETCLR(grey) ; CALL GRID(1,1)   ! GRIDLINES
      CALL HEIGHT(22) ; CALL SETCLR(red) ; CALL LINWID(1)
      CALL MESSAG("Tree Indices",grl,grt-35)
      CALL SETCLR(blue)
      CALL MESSAG("Site Indices",grl+(grr-grl)/2,grt-35)
      CALL line_miss(r,wka(1:r),xcrn(u:v,mx),cok(u:v,mx))
      CALL SETCLR(red)
      CALL line_miss(r,wka(1:r),dx(p:q),xok(p:q)) 
      CALL SETCLR(black) ; CALL ENDGRF()
      RETURN
      END SUBROUTINE det_plot
!-------------------------------------------------------------------
      SUBROUTINE plot_deta()        ! Plots RCS curve display
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka    ! Working areas
      INTEGER                :: i=1,j,p,q,r,u,v
      REAL(8)                :: top,bot
      wka(1:sly(mx))=(/(DBLE(i),i=1,sly(mx))/)
      p=MAX(RDfy,sfy(mx)) ; q=MIN(RDly,sly(mx)) ; r=q-p+1
      IF (lins(4)%on) THEN
        top=MAXVAL(mval(p:q,mx)+mssd(p:q,mx))
        bot=MINVAL(mval(p:q,mx)-mssd(p:q,mx))
      ELSEIF (lins(5)%on) THEN
        top=MAXVAL(mval(p:q,mx)+mserr(p:q,mx))
        bot=MINVAL(mval(p:q,mx)-mserr(p:q,mx))
      ELSE
        top=MAXVAL(mval(p:q,mx))
        bot=MINVAL(mval(p:q,mx))
      ENDIF
      IF (src.GE.2) THEN
        DO i=1,srcno
          u=MAX(p,sfy(i)) ; v=MIN(q,sly(i))
          top=MAX(top,MAXVAL(msmo(u:v,i)))
          bot=MIN(bot,MINVAL(msmo(u:v,i)))
        ENDDO
      ENDIF
      bot=MAX(DBLE(INT(bot*10.D0))/10.D0,0.D0)
      grl=500 ; grr=2400 ; grt=250 ; grb=900
      IF (lins(6)%on) CALL plot_trees(r,mcnt(p:q,mx))  
      CALL NAME('Ring Age','X')      ! Axis name
      IF (itn.EQ.3) THEN             ! Basal Area
        CALL NAME('Ring Area','Y')   ! Axis name
      ELSE
        CALL NAME('Ring mm','Y')     ! Axis name
      ENDIF
      CALL TICKS(10,'X')      ! Number of ticks between labels
      CALL TICKS(5,'Y')       ! Number of ticks between labels
      CALL tombox(RDfy,RDly,bot,top)
      CALL SETCLR(grey) 
      CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22)
      IF (psok) THEN
        CALL SETCLR(black) 
        CALL MESSAG(cnam(cf),grl,grt-90)
      ENDIF
      IF (src.GE.2) THEN    ! Individual RCS Curves
        DO j=1,srcno
          i=11+j*4 ; u=MAX(p,sfy(j)) ; v=MIN(q,sly(j)) ; r=v-u+1
          CALL plot_line(r,wka(u:v),mval(u:v,j)+mssd(u:v,j), &
            mok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),mval(u:v,j)-mssd(u:v,j), &
            mok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),mval(u:v,j)+mserr(u:v,j), &
            mok(u:v,j),i+3,-10)
          CALL plot_line(r,wka(u:v),mval(u:v,j)-mserr(u:v,j), &
            mok(u:v,j),i+3,-10) 
          CALL plot_line(r,wka(u:v),mval(u:v,j),mok(u:v,j),i+1,-10)
          CALL plot_line(r,wka(u:v),msmo(u:v,j),mok(u:v,j),i,-10) 
        ENDDO
      ENDIF
      r=q-p+1
      CALL plot_line(r,wka(p:q),mval(p:q,mx)+mssd(p:q,mx), &
        mok(p:q,mx),4,-10)
      CALL plot_line(r,wka(p:q),mval(p:q,mx)-mssd(p:q,mx), &
        mok(p:q,mx),4,60)
      CALL plot_line(r,wka(p:q),mval(p:q,mx)+mserr(p:q,mx), &
        mok(p:q,mx),5,-10)
      CALL plot_line(r,wka(p:q),mval(p:q,mx)-mserr(p:q,mx), &
        mok(p:q,mx),5,80) 
      CALL plot_line(r,wka(p:q),mval(p:q,mx),mok(p:q,mx),2,30)
      CALL plot_line(r,wka(p:q),msmo(p:q,mx),mok(p:q,mx),3,2) 
      CALL ENDGRF()
      CALL SETCLR(black) ! Reset lines
      RETURN
      END SUBROUTINE plot_deta
!-------------------------------------------------------------------
      SUBROUTINE plot_detb()  ! Plots RCS Diam curve display
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka    ! Working areas
      INTEGER                :: i=1,j,p,q,r,u,v
      REAL(8)                :: top,bot
      IF (dfy(mx).LT.RDfy.OR.dly(mx).LE.dfy(mx)) RETURN
      wka(1:dly(mx))=(/(DBLE(i),i=1,dly(mx))/)
      p=MAX(RDfy,dfy(mx)) ; q=MIN(RDly,dly(mx)) ; r=q-p+1
      IF (lins(4)%on) THEN
        top=MAXVAL(dval(p:q,mx)+mssd(p:q,mx))
        bot=MINVAL(dval(p:q,mx)-mssd(p:q,mx))
      ELSEIF (lins(5)%on) THEN
        top=MAXVAL(dval(p:q,mx)+mserr(p:q,mx))
        bot=MINVAL(dval(p:q,mx)-mserr(p:q,mx))
      ELSE
        top=MAXVAL(dval(p:q,mx))
        bot=MINVAL(dval(p:q,mx))
      ENDIF
      IF (src.EQ.2) THEN
        DO i=1,srcno
          u=dfy(i) ; v=dly(i)
          top=MAX(top,MAXVAL(dsmo(u:v,i)))
          bot=MIN(bot,MINVAL(dsmo(u:v,i)))
        ENDDO
      ENDIF
      bot=MAX(DBLE(INT(bot*10.D0))/10.D0,0.D0)
      grl=500 ; grr=2400 ; grt=250 ; grb=900
      IF (lins(6)%on) CALL plot_trees(r,dcnt(p:q,mx)) ! Tree counts 
      CALL NAME('Radius mm','X')   ! Axis name
      CALL NAME('Ring mm','Y')     ! Axis name
      CALL tombox(RDfy,RDly,bot,top)
      CALL SETCLR(silver) 
      CALL SETCLR(grey) 
      CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22)
      IF (psok) THEN

        CALL SETCLR(black) 
        CALL MESSAG(cnam(cf),grl,grt-70)
      ENDIF
      IF (src.GT.1) THEN    ! Individual RCS Diam Curves
        DO j=1,srcno
          i=11+j*4 ; ; u=MAX(p,dfy(j)) ; v=MIN(q,dly(j)) ; r=v-u+1
          CALL plot_line(r,wka(u:v),dval(u:v,j)+mssd(u:v,j), &
            dok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),dval(u:v,j)-mssd(u:v,j), &
            dok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),dval(u:v,j)+mserr(u:v,j), &
            dok(u:v,j),i+3,-10)
          CALL plot_line(r,wka(u:v),dval(u:v,j)-mserr(u:v,j), &
            dok(u:v,j),i+3,-10) 
          CALL plot_line(r,wka(u:v),dval(u:v,j),dok(u:v,j),i+1,-10)
          CALL plot_line(r,wka(u:v),dsmo(u:v,j),dok(u:v,j),i,-10) 
        ENDDO
      ENDIF
      r=q-p+1
      CALL plot_line(r,wka(p:q),dval(p:q,mx)+mssd(p:q,mx), &
        dok(p:q,mx),4,-10)
      CALL plot_line(r,wka(p:q),dval(p:q,mx)-mssd(p:q,mx), &
        dok(p:q,mx),4,60)
      CALL plot_line(r,wka(p:q),dval(p:q,mx)+mserr(p:q,mx), &
        dok(p:q,mx),5,-10)
      CALL plot_line(r,wka(p:q),dval(p:q,mx)-mserr(p:q,mx), &
        dok(p:q,mx),5,80) 
      CALL plot_line(r,wka(p:q),dval(p:q,mx),dok(p:q,mx),2,30)
      CALL plot_line(r,wka(p:q),dsmo(p:q,mx),dok(p:q,mx),3,2) 
      CALL ENDGRF()
      CALL SETCLR(black)
      RETURN
      END SUBROUTINE plot_detb
!-------------------------------------------------------------------
      SUBROUTINE plot_detc()        ! Plots CRNs display
      IMPLICIT NONE                 
      REAL(8),DIMENSION(mxy) :: wka     ! Working areas
      INTEGER                :: i=1,j,p,q,r,u,v
      REAL(8)                :: top,bot
      wka(1:xyr)=(/(DBLE(i),i=xfy,xly)/)
      p=MAX(xfy,CDfy)-xfy+1 ; q=MIN(xly,CDly)-xfy+1 ; r=q-p+1
      top=MAXVAL(xcrn(p:q,mx),MASK=cok(p:q,mx))
      bot=MINVAL(xcrn(p:q,mx),MASK=cok(p:q,mx))
      IF (src.GT.1) THEN
        DO i=1,srcno
          top=MAX(top,MAXVAL(xcsm(p:q,i),MASK=cok(p:q,i)))
          bot=MIN(bot,MINVAL(xcsm(p:q,i),MASK=cok(p:q,i)))
        ENDDO
      ENDIF
      IF (bot.GE.-epsi) THEN
        bot=MAX(DBLE(INT(bot*10.D0))/10.D0,0.D0)
      ELSE
        bot=DBLE(INT((bot*10.D0)))/10.D0
      ENDIF
      grl=500 ; grr=2400 ; grt=1100 ; grb=1750
      IF (lins(9)%on) CALL plot_trees(r,xnum(p:q,mx))  
      CALL NAME('Calendar Years','X')  ! Axis name
      CALL NAME('Index','Y')           ! Axis name
      CALL tombox(CDfy,CDly,bot,top)
      CALL SETCLR(grey) 
      CALL GRID(1,1)        ! GRIDLINES
      CALL HEIGHT(22)  

      IF (src.GT.1) THEN    ! Individual CRNS
        DO j=1,srcno
          i=55+j*4 ; u=MAX(p,xfa(j)) ; v=MIN(q,xla(j)) ; r=v-u+1
          CALL plot_line(r,wka(u:v),xcrn(u:v,j)+xcsd(u:v,j), &
            cok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),xcrn(u:v,j)-xcsd(u:v,j), &
            cok(u:v,j),i+2,-10)
          CALL plot_line(r,wka(u:v),xcrn(u:v,j),cok(u:v,j),i+1,-10) 
          CALL plot_line(r,wka(u:v),xcsm(u:v,j),cok(u:v,j),i,-10) 
        ENDDO
      ENDIF
      r=q-p+1
      CALL plot_line(r,wka(p:q),xcrn(p:q,mx)+xcsd(p:q,mx), &
        cok(p:q,mx),8,-10)
      CALL plot_line(r,wka(p:q),xcrn(p:q,mx)-xcsd(p:q,mx), &
        cok(p:q,mx),8,70)
      CALL plot_line(r,wka(p:q),xcrn(p:q,mx),cok(p:q,mx),7,35) 
      CALL plot_line(r,wka(p:q),xcsm(p:q,mx),cok(p:q,mx),10,2) 
      CALL ENDGRF()
      CALL SETCLR(black)
      RETURN
      END SUBROUTINE plot_detc
!-------------------------------------------------------------------
      SUBROUTINE plot_det()        ! Plots RCS curve display
      IMPLICIT NONE                 
      IF (idt.EQ.-2) THEN     ! RCS
        IF (trc.EQ.2) THEN    ! Diameter
          CALL plot_detb()
        ELSE                  ! Age
          CALL plot_deta()
        ENDIF
        CALL plot_detc()      ! Chronologies
      ELSE
        CALL det_plot()       ! Curve fitting 
      ENDIF
      RETURN
      END SUBROUTINE plot_det
!-------------------------------------------------------------------
      SUBROUTINE plot_psbeg()     ! Saves current figure as an .ps file
      IMPLICIT NONE
      CALL DISINI()
      CALL SETVLT('VGA')          ! Select 16 colour table
      CALL PSFONT('Courier')
      CALL WINMOD('NONE')
      CALL HEIGHT(22)             ! Character height - plot coords
      CALL SHDPAT(16)             ! Fill areas
      CALL TICPOS('REVERS','XY')  ! Internal tick marks
      CALL HNAME(20)              ! Character height for Axis names
      CALL LABDIG(-1,'X')         ! Digits after decimal -1=integer
      CALL TICKS(10,'X')          ! Number of ticks between labels
      CALL TICKS(5,'Y')           ! Number of ticks between labels
      RETURN
      END SUBROUTINE plot_psbeg
!------------------------------------------------------------------------
      SUBROUTINE plot_psend()     ! Saves current figure as an .ps file
      IMPLICIT NONE
      CALL DISFIN()
      psok=FA                     ! Turn off .eps messages
      CALL METAFL('XWIN')
      CALL SETPAG('DA4L')
      CALL DISINI()
      CALL WINMOD('NONE')         ! DISINI drops out 
      CALL CSRMOD('READ','POS')   ! Control cursor reading
      CALL SETVLT('VGA')          ! Select 16 colour table
      CALL SIMPLX()
      CALL HEIGHT(22)             ! Character height - plot coords
      CALL SHDPAT(16)             ! Fill areas
      CALL TICPOS('REVERS','XY')  ! Internal tick marks
      CALL HNAME(18)              ! Character height for Axis names
      CALL LABDIG(-1,'X')         ! Digits after decimal -1=integer
      CALL TICKS(10,'X')          ! Number of ticks between labels
      CALL TICKS(5,'Y')           ! Number of ticks between labels
      RETURN
      END SUBROUTINE plot_psend
!------------------------------------------------------------------------
      SUBROUTINE plot_ps(fig)   ! Saves current figure as an .ps file
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: fig ! Figure to be saved
      CHARACTER(11) :: plnam
      psok=TR                   ! Turn on eps messages 
      CALL DISFIN()
      CALL METAFL('EPS')
      SELECT CASE (fig)
        CASE (1) ; plnam='detrend.eps'   ! Filename
        CASE (2) ; plnam='Trees.eps'  ! Filename
      ENDSELECT
      CALL SETFIL(TRIM(plnam))
      CALL plot_psbeg() 
      SELECT CASE (fig)
        CASE (1) ; CALL plot_det()
        CASE (2) ; CALL plot_core() 
      END SELECT
      CALL plot_psend
      RETURN
      END SUBROUTINE plot_ps
!------------------------------------------------------------------------
      SUBROUTINE read_list()   
      IMPLICIT NONE
      INTEGER :: i,j,k
      CALL read_open(17,lnam)
      ID: DO i=1,mxc                   
        READ(17,'(A60)',IOSTAT=ios) cnam(i) 
        IF (ios.EQ.-1) EXIT ID            ! End of file
        IF (io_err("Read",lnam)) STOP
        k=LEN_TRIM(cnam(i))
        IF (k.LT.3) THEN
          WRITE(mess(1:4),'(I4)') i
          CALL out_err("No"//mess(1:4)//" filename too short "//lnam)
          STOP
        ENDIF
        J1: DO j=k,1,-1
          IF (cnam(i)(j:j).EQ."/".OR.cnam(i)(j:j).EQ.CHAR(92)) EXIT J1
        ENDDO J1
        k=MIN(j+20,k)
        wnam(i)=cnam(i)(j+1:k)
      ENDDO ID
      CLOSE(17)
      IF (i.EQ.1) THEN
        CALL out_err("No data files in "//lnam)
      ELSE
        nfil=i-1 ; IF (cf.GT.nfil.OR.cf.LT.1) cf=1
        b(24)%lab=cnam(cf)
      ENDIF
      RETURN
      END SUBROUTINE read_list
!--------------------------------------------------------------
      SUBROUTINE disp_init()           ! Display Initialise
      IMPLICIT NONE                 
      WHERE (lins(2:102)%on)
        b(88:188)%lab="On"
      ELSEWHERE
        b(88:188)%lab="Off"
      END WHERE
      WHERE (lins(2:14)%on)
        b(68:80)%lab="On"
      ELSEWHERE
        b(68:80)%lab="Off"
      END WHERE
      WRITE(b(81)%lab,'(I4)') RDfy
      WRITE(b(82)%lab,'(I4)') RDly
      WRITE(b(83)%lab,'(I5)') CDfy
      WRITE(b(84)%lab,'(I4)') CDly
      WRITE(b(85)%lab,'(I4)') CDsp
      WRITE(b(86)%lab,'(I5)') TDfy
      WRITE(b(87)%lab,'(I5)') TDly
      b(190)%lab="Turn Buttons on/off" ; oocol=TR
      b(202)%lab=dirc ; b(203)%lab=namc
      b(207)%lab=dirr ; b(208)%lab=namr
      append=TR ; b(201)%lab="Append"
      lth=1 ; lst=0 ; lco=0  ! Default line style and colour
      RETURN 
      END SUBROUTINE disp_init          
!------------------------------------------------------------------------
      SUBROUTINE write_default()  ! Writes parameter default file
      IMPLICIT NONE
      CHARACTER(14),PARAMETER :: fnam = "RCSdefault.fil"  
      INTEGER                 :: i
      OPEN(29,FILE=fnam,STATUS='REPLACE',IOSTAT=IOS)
      IF (io_err("Open",fnam)) STOP   ! Open error?
      WRITE(29,'("DEFAULT PARAMETER VALUES FOR RCS PROGRAM")',IOSTAT=ios)
      IF (io_err("Write  Head",fnam)) STOP
      WRITE(29,'(I4," = Screen width")',IOSTAT=ios) screenw
      IF (io_err("Write  screenw",fnam)) STOP
      WRITE(29,'(I4," = Screen height")',IOSTAT=ios) screenh
      IF (io_err("Write  screenh",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (IDT)")',IOSTAT=ios) idt,idtlab(idt)
      IF (io_err("Write  idt",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (ITN)")',IOSTAT=ios) itn,itnlab(itn)
      IF (io_err("Write  itn",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (RDT)")',IOSTAT=ios) rdt,rdtlab(rdt)
      IF (io_err("Write  rdt",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (IND)")',IOSTAT=ios) ind,indlab(ind)
      IF (io_err("Write  ind",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (KRB)")',IOSTAT=ios) krb,krblab(krb)
      IF (io_err("Write  krb",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (ISB)")',IOSTAT=ios) isb,isblab(isb)
      IF (io_err("Write  isb",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (SFO)")',IOSTAT=ios) sfo,sfolab(sfo)
      IF (io_err("Write  sfo",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (POO)")',IOSTAT=ios) poo,poolab(poo)
      IF (io_err("Write  poo",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (SRC)")',IOSTAT=ios) src,srclab(src)
      IF (io_err("Write  src",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (TRC)")',IOSTAT=ios) trc,trclab(trc)
      IF (io_err("Write  trc",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (GTR)")',IOSTAT=ios) gtr,gtrlab(gtr)
      IF (io_err("Write  gtr",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (TST)")',IOSTAT=ios) tst,tstlab(tst)
      IF (io_err("Write  tst",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (BFC)")',IOSTAT=ios) bfc,bfclab(bfc)
      IF (io_err("Write  bfc",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (IDB)")',IOSTAT=ios) idb,idblab(idb)
      IF (io_err("Write  idb",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (JRB)")',IOSTAT=ios) jrb,jrblab(jrb)
      IF (io_err("Write  jrb",fnam)) STOP
      WRITE(29,'(I3,"  Number of RCS curves (SRCNO)")',IOSTAT=ios) srcno
      IF (io_err("Write  srcno",fnam)) STOP
      WRITE(29,'(I3,"  Max SF iterations    (SFONO)")',IOSTAT=ios) sfono
      IF (io_err("Write  sfono",fnam)) STOP
      WRITE(29,'(I3,"  Spline Stiffness     (RDTNO)")',IOSTAT=ios) rdtno
      IF (io_err("Write  rdtno",fnam)) STOP
      WRITE(29,'(I3,"  Ignore 1st n rings   (IGNOR)")',IOSTAT=ios) ignor
      IF (io_err("Write  ignor",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (RFT)")',IOSTAT=ios) rft,rftlab(rft)
      IF (io_err("Write  rft",fnam)) STOP
      WRITE(29,'(I3,2X,A14," (CFT)")',IOSTAT=ios) cft,cftlab(cft)
      IF (io_err("Write  cft",fnam)) STOP
      WRITE(29,'(I3,"  Display Spline Length")',IOSTAT=ios) CDsp
      IF (io_err("Write  CDsp",fnam)) STOP
      WRITE(29,'(I3,"  RBar for EPS segment length")',IOSTAT=ios) RBsp
      IF (io_err("Write  RBsp",fnam)) STOP
      WRITE(29,'("PLOTTING OPTIONS")',IOSTAT=ios)
      IF (io_err("Write  plot",fnam)) STOP
      DO i=1,102
        WRITE(29,'(L1,3I3,2X,A30)',IOSTAT=ios) lins(i)
        IF (io_err("Write  lines",fnam)) STOP
      ENDDO
      CLOSE(29)
      RETURN
      END SUBROUTINE write_default
!----------------------------------------------------------------
      SUBROUTINE read_default()  ! Reads parameter default file   
      IMPLICIT NONE
      CHARACTER(14),PARAMETER :: fnam = "RCSdefault.fil"  
      INTEGER                 :: i
      sccc=1
      CALL read_open(23,fnam)
      IF (.NOT.fileok) STOP
      DO i=1,3  
        READ(23,*,IOSTAT=ios)
        IF (io_err("Read Head",fnam)) STOP
      ENDDO
      READ(23,'(I4)',IOSTAT=ios) idt
      IF (io_err("Read itn",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) itn
      IF (io_err("Read itn",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) rdt
      IF (io_err("Read rdt",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) ind
      IF (io_err("Read ind",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) krb
      IF (io_err("Read krb",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) isb
      IF (io_err("Read isb",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) sfo
      IF (io_err("Read sfo",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) poo
      IF (io_err("Read poo",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) src
      IF (io_err("Read src",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) trc
      IF (io_err("Read trc",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) gtr
      IF (io_err("Read gtr",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) tst
      IF (io_err("Read tst",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) bfc
      IF (io_err("Read bfc",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) idb
      IF (io_err("Read idb",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) jrb
      IF (io_err("Read jrb",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) srcno
      IF (io_err("Read srcno",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) sfono 
      IF (io_err("Read sfono",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) rdtno
      IF (io_err("Read rdtsl",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) ignor
      IF (io_err("Read ignor",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) rft
      IF (io_err("Read rft",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) cft
      IF (io_err("Read cft",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) CDsp
      IF (io_err("Read CDsp",fnam)) STOP
      READ(23,'(I4)',IOSTAT=ios) RBsp
      IF (io_err("Read RBsp",fnam)) STOP
      READ(23,*,IOSTAT=ios)  ! Header line
      IF (io_err("Read Head",fnam)) STOP
      DO i=1,102
        READ(23,'(L1,3I3,2X,A30)',IOSTAT=ios) lins(i)
        IF (io_err("Read Line",fnam)) STOP
      ENDDO
      journal=lins(1)%on
      CLOSE(23)
      WRITE(b( 1)%lab,'(A14)') itnlab(itn)
      WRITE(b( 2)%lab,'(A14)') poolab(poo)
      WRITE(b( 3)%lab,'(A14)') indlab(ind)
      WRITE(b( 4)%lab,'(A14)') krblab(krb)
      WRITE(b( 5)%lab,'(A14)') sfolab(sfo)
      WRITE(b( 6)%lab,'(A14)') jrblab(jrb)
      WRITE(b( 7)%lab,'(A14)') idblab(idb)
      WRITE(b( 8)%lab,'(A14)') isblab(isb)
      WRITE(b( 9)%lab,'(A14)') rdtlab(rdt)
      WRITE(b(10)%lab,'(A14)') srclab(src)
      WRITE(b(11)%lab,'(A14)') trclab(trc)
      WRITE(b(12)%lab,'(A14)') gtrlab(gtr)
      WRITE(b(19)%lab,'(A14)') tstlab(tst)
      WRITE(b(209)%lab,'(A14)') tstlab(tst)
      WRITE(b(18)%lab,'(A14)') bfclab(bfc)
      WRITE(b(16)%lab,'(I4)') srcno
      b(16)%ok=(src.GE.2) 
      WRITE(b(20)%lab,'(I4)') sfono
      b(20)%ok=(sfo.EQ.2)               ! If signal free
      WRITE(b(21)%lab,'(I4)') ignor
      WRITE(b(13)%lab,'(I4)') rdtno
      b(13)%ok=(rdt.EQ.6.OR.rdt.EQ.7)   ! If RCS spline needed
      b(17)%ok=(src.NE.1.AND.gtr.GT.1)  ! If Multi and Ring Transform
      b(37)%ok=(src.NE.1)               ! If multiple RCS
      b(12)%ok=(src.NE.1)     ! If multiple RCS
      WRITE(b(205)%lab,'(A14)') rftlab(rft)
      WRITE(b(200)%lab,'(A14)') cftlab(cft)
      idtsl=100       ! Default spline
      idtsn=-66       ! Default % spline
      WRITE(b(51)%lab,'(I4)') idtsl
      WRITE(b(52)%lab,'(I4)') idtsn
      mtext(12)=cftext(cft)
      b(1:3)%ok=TR ; b(5)%ok=TR ; b(9)%ok=TR ; b(11)%ok=TR
      dirc="work" ; namc="crns" ; dirr="work" ; namr="rawdata"
      CALL disp_init()
      RETURN
      END SUBROUTINE read_default
!-------------------------------------------------------------------------
      SUBROUTINE write_srcs()  ! Sub writes SRCS file
      IMPLICIT NONE   
      CHARACTER(60)             :: fnam
      INTEGER,DIMENSION(mx*2+1) :: lin1
      INTEGER                   :: i=1,j,m,p
      fnam=cnam(cf) ; j=LEN_TRIM(fnam) 
      JD: DO j=j,3,-1 ; IF (fnam(j:j).EQ.".") EXIT JD ; ENDDO JD
      fnam=fnam(1:j)//"src"
      OPEN(73,FILE=fnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",fnam)) RETURN
      WRITE(73,'(I2,"  Number of RCS curves")') srcno
      WRITE(73,'(I2,"  Transform - BAI")') itn
      WRITE(73,'(I2,"  RCS smoothing Option")') rdt
      WRITE(73,'(I2,"  Ratios/Differences")') ind
      WRITE(73,'(I2,"  Signal Free on/off")') sfo
      WRITE(73,'(I2,"  Use/not Pith Offset")') poo
      WRITE(73,'(I2,"  RCS Curve Type")') trc
      WRITE(73,'("Tree        Build  Detrend")')
      WRITE(73,'("              RCS  Tree")')
      DO j=1,nc
        WRITE(73,FMT='(A12,2I4)',IOSTAT=ios) nam(j),trr(j),trm(j) 
        IF (io_err("Write3 ",TRIM(fnam)//" "//nam(j))) RETURN
      ENDDO 
      WRITE(73,'("Smoothed SRCS Curve Values and Counts")',IOSTAT=ios)
      m=srcno 
      WRITE(mess,'("Curve",11("  RCS ",I2))') (/(i,i=1,m)/)
      mess=mess(1:6+m*8)//"One RCS"
      WRITE(73,'(A122)') mess
      IF (trc.EQ.2) THEN    ! Diameter based RCS
        p=dly(mx)
        WRITE(73,'("Start",12I8)') dfy(1:m),dfy(mx)
        WRITE(73,'("End  ",12I8)') dly(1:m),dly(mx)
        DO j=1,p     ! For each year
          lin1(1)=j
          DO i=1,m
            lin1(i*2+1)=dcnt(j,i)
            lin1(i*2)=NINT(dsmo(j,i)*1000.D0)
          ENDDO
          lin1(m*2+3)=dcnt(j,mx)
          lin1(m*2+2)=NINT(dsmo(j,mx)*1000.D0)
          WRITE(73,'(I5,11(I6,I4))') lin1(1:m*2+3)
        ENDDO
      ELSE                  ! Age based or Basal area 
        p=sly(mx)
        WRITE(73,'("Start",12I8)') sfy(1:m),sfy(mx)
        WRITE(73,'("End  ",12I8)') sly(1:m),sly(mx)
        DO j=1,p     ! For each year
          lin1(1)=j
          DO i=1,m
            lin1(i*2+1)=mcnt(j,i) ; lin1(i*2)=NINT(msmo(j,i)*1000.D0)
          ENDDO
          lin1(m*2+3)=mcnt(j,mx) ; lin1(m*2+2)=NINT(msmo(j,mx)*1000.D0)
          WRITE(73,'(I5,12(I6,I4))') lin1(1:m*2+3)
        ENDDO
      ENDIF
      CLOSE(73)
      RETURN
      END SUBROUTINE write_srcs
!----------------------------------------------------------------
      SUBROUTINE read_srcs()  ! Read SRCS file
      IMPLICIT NONE   
      INTEGER,DIMENSION(mx*2+1) :: lin1
      CHARACTER(60)             :: fnam
      CHARACTER(12)             :: tnam
      INTEGER                   :: i,j
      fnam=cnam(cf) ; j=LEN_TRIM(fnam) 
      JD: DO j=j,3,-1 ; IF (fnam(j:j).EQ.".") EXIT JD ; ENDDO JD
      fnam=fnam(1:j)//"src"
      OPEN(71,FILE=fnam,IOSTAT=ios,STATUS="OLD")
      IF (io_err("Open",fnam)) STOP
      READ(71,'(I2)',IOSTAT=ios) i      ! srcno
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.2.OR.i.GE.mx) THEN
        WRITE(mess,'("Invalid RCS count ",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        srcno=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! itn
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.3) THEN
        WRITE(mess,'("ITN invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        itn=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! rdt
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.7) THEN
        WRITE(mess,'("RDT invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        rdt=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! ind
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.2) THEN
        WRITE(mess,'("IND invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        ind=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! sfo
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.2) THEN
        WRITE(mess,'("SFO invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        sfo=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! poo
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.2) THEN
        WRITE(mess,'("POO invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        poo=i
      ENDIF
      READ(71,'(I2)',IOSTAT=ios) i      ! trc
      IF (io_err("Read1",fnam)) RETURN
      IF (i.LT.1.OR.i.GT.3) THEN
        WRITE(mess,'("TRC invalid",I3)') i
        CALL out_err(mess) ; STOP
      ELSE
        trc=i
      ENDIF
      READ(71,*) ; READ(71,*)    ! Header lines
      DO j=1,nc     ! For each tree
        READ(71,FMT='(A12,2I4)',IOSTAT=ios) tnam,trr(j),trm(j) 
        IF (io_err("Read3",fnam)) RETURN
        IF (tnam.NE.nam(j)) THEN  ! SRCS names not match
          CALL out_err("SRCS "//TRIM(tnam)//" loaded tree "//nam(j))
          IF (src.EQ.3) src=2 
          RETURN
        ENDIF
        IF (trr(j).LT.0.OR.trr(j).GT.srcno.OR. &
            trm(j).LT.1.OR.trm(j).GT.srcno) THEN  ! Invalid RCS curve number
          CALL out_err("SRCS Invalid RCS no "//tnam)
          IF (src.EQ.3) src=2 
          RETURN
        ENDIF
      ENDDO 
      READ(71,*) ; READ(71,*)  ! Headers
      IF (trc.EQ.2) THEN       ! Diameter based RCS
        dfy=0 ; dly=0 ; dcnt=0 ; dsmo=0.D0 
        READ(71,'(5X,12(I6,4X))') dfy(1:srcno),dfy(mx)   ! First years
        READ(71,'(5X,12(I6,4X))') dly(1:srcno),dly(mx)   ! Last years
        DO j=1,dly(mx) ! For each year
          READ(71,'(I5,12(I6,I4))') lin1(1:srcno*2+2)
          IF (lin1(1).NE.j) THEN
            WRITE(mess,'("Years wrong ",2I5)') lin1(1),j
            CALL out_err(mess) ; STOP
          ENDIF
          DO i=1,srcno
            dcnt(j,i)=lin1(i*2+1) ; dsmo(j,i)=DBLE(lin1(i*2))/1000.D0
          ENDDO
          dcnt(j,mx)=lin1(srcno*2+3)
          dsmo(j,mx)=DBLE(lin1(srcno*2+2))/1000.D0
        ENDDO
        dok=dcnt.GE.1 ; dval=dsmo
      ELSE                  ! Age based or Basal area 
        sfy=0 ; sly=0 ; mcnt=0 ; msmo=0.D0 ; mssd=0.D0 ; mserr=0.D0
        READ(71,'(5X,12I8)') sfy(1:srcno),sfy(mx)   ! First years
        READ(71,'(5X,12I8)') sly(1:srcno),sly(mx)   ! Last years
        DO j=1,sly(mx) ! For each year
          READ(71,'(I5,12(I6,I4))') lin1(1:srcno*2+3)
          IF (lin1(1).NE.j) THEN
            WRITE(mess,'("Years wrong ",2I5)') lin1(1),j
            CALL out_err(mess) ; STOP
          ENDIF
          DO i=1,srcno
            mcnt(j,i)=lin1(i*2+1) ; msmo(j,i)=DBLE(lin1(i*2))/1000.D0
          ENDDO
          mcnt(j,mx)=lin1(srcno*2+3)
          msmo(j,mx)=DBLE(lin1(srcno*2+2))/1000.D0
        ENDDO
        mok=mcnt.GE.1 ; mval=msmo
      ENDIF
      CLOSE(71)
      RETURN
      END SUBROUTINE read_srcs
!----------------------------------------------------------------
      SUBROUTINE ignore_rings()  ! Removes 1st few rings from CRNs
      IMPLICIT NONE  ! Must retain half (or more) rings each tree               
      INTEGER :: i,p,q,r,u,v,ign,first,last
      nc=0 ; CALL read_rft(cnam(cf))
      IF (ignor.GT.0) THEN
!        OPEN(73,FILE="ignore.prn",IOSTAT=ios,STATUS="REPLACE")
        DO i=1,nc        ! For each tree
!          WRITE(73,'(A8,3I5,I7,3F8.3)') nam(i)(1:8),fy(i),ly(i),yr(i), &
!            ad(i),pthr(i),x(ad(i)),x(ad(i)+yr(i)-1)
          first=fy(i)-pth(i)+1              ! First age
          last=ly(i)-pth(i)+1-yr(i)/2       ! Last age can remove      
          IF (ignor.GE.first) THEN          ! Blank some rings
            ign=MIN(ignor,last)-first      
            p=ad(i) ; r=yr(i) ; q=p+r-1     ! Where tree was
            IF (i.GT.1) THEN
              u=ad(i-1)+yr(i-1)
            ELSE
              u=1
            ENDIF
            v=u+r-1-ign ! New address of rings
            fy(i)=fy(i)+ign                 ! New first year 
            yr(i)=yr(i)-ign                 ! New length
            pthr(i)=pthr(i)+SUM(x(p:p+ign)) ! Adjust pith radius
          ENDIF
          ad(i)=u
          xok(u:v)=xok(p+ign:q)             ! Adjust logicals  
          x(u:v)=x(p+ign:q)                 ! Adjust measurements
!          WRITE(73,'(A8,3I5,I7,3F8.3,I5)') nam(i)(1:8),fy(i),ly(i),yr(i), &
!            ad(i),pthr(i),x(ad(i)),x(ad(i)+yr(i)-1),ign
!          WRITE(73,*)
        ENDDO
        ad(nc+1)=ad(nc)+yr(nc)              ! Next ring address
      ENDIF
      ignok=FA  
!      CLOSE(73)
      RETURN 
      END SUBROUTINE ignore_rings          
!------------------------------------------------------------------------
      SUBROUTINE new_chron(loadok)  ! Builds new chronology 
      IMPLICIT NONE           
      LOGICAL,INTENT(IN) :: loadok
      LOGICAL :: dataok
      IF (loadok) THEN     ! May need to read raw data
        INQUIRE(FILE=cnam(cf),EXIST=dataok)  ! Load data
        IF (dataok) THEN
          nc=0 ; CALL read_rft(cnam(cf))
          IF (ignor.GT.0) CALL ignore_rings()
          srcno=MIN(MAX(nc/40,1),mx-1)  ! Suggested MRCS curves
          WRITE(b(16)%lab,'(I4)') srcno
!          tst=4 ; b(55)%lab=tstlab(tst) 
        ELSE
          CALL out_err(TRIM(cnam(cf))//" Raw data file does not exist")
          STOP
        ENDIF
      ENDIF
      IF (ignok) CALL ignore_rings()
      CALL detrend()
      IF (trc.EQ.2) THEN                ! Diameter curves
        RDfy=dfy(mx)   ; WRITE(b(81)%lab,'(I4)') RDfy
        RDly=dly(mx)   ; WRITE(b(82)%lab,'(I4)') RDly
      ELSE                              ! Age based curves
        RDfy=sfy(mx)   ; WRITE(b(81)%lab,'(I4)') RDfy
        RDly=sly(mx)   ; WRITE(b(82)%lab,'(I4)') RDly
      ENDIF
      CDfy=xfy         ; WRITE(b(83)%lab,'(I5)') CDfy
      CDly=xly         ; WRITE(b(84)%lab,'(I4)') CDly
      CDsp=40          ; WRITE(b(85)%lab,'(I4)') CDsp
      TDfy=fy(tre(cc)) ; WRITE(b(86)%lab,'(I5)') TDfy
      TDly=ly(tre(cc)) ; WRITE(b(87)%lab,'(I5)') TDly
      bred=FA          ! Red highlight turned off
      RETURN
      END SUBROUTINE new_chron
!------------------------------------------------------------------------
      SUBROUTINE start_init()    ! Initialise start menu
      IMPLICIT NONE                 
      INQUIRE(FILE="RCSdefault.fil",EXIST=fileok)  
      IF (.NOT.fileok) CALL write_default()
      CALL read_default()
      INQUIRE(FILE=lnam,EXIST=fileok)  
      IF (fileok) THEN
        CALL read_list() ; CALL new_chron(TR)   
      ELSE                 
        CALL out_err("Does not exist "//lnam) ; STOP
      ENDIF
      RETURN
      END SUBROUTINE start_init
!------------------------------------------------------------------
      SUBROUTINE thickthin(n,xval,yval,cnt,cut)  
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)              :: n      ! Number of years
      REAL(8),DIMENSION(n),INTENT(IN) :: xval   ! Year numbers
      REAL(8),DIMENSION(n),INTENT(IN) :: yval   ! Values
      INTEGER,DIMENSION(n),INTENT(IN) :: cnt
      INTEGER,INTENT(IN)              :: cut    ! Thin count
      LOGICAL,DIMENSION(n)            :: nok
      CALL LINWID(4) ; nok=cnt.GT.cut 
      CALL line_miss(n,xval,yval,nok)
      CALL LINWID(1) ; nok=cnt.GT.0 
      CALL line_miss(n,xval,yval,nok)
      RETURN 
      END SUBROUTINE thickthin
!-------------------------------------------------------------------
      SUBROUTINE pgcirc(xcent,ycent,rx,ry,fill)  ! Ex PGCIRC
      IMPLICIT NONE
      REAL(8),INTENT(IN)       :: xcent,ycent ! Centre of circle
      REAL(8),INTENT(IN)       :: rx,ry       ! Radius each direction 
      LOGICAL,INTENT(IN)       :: fill        ! Filled or outline
      INTEGER,PARAMETER        :: pts=72      ! Points around circle
      INTEGER                  :: i
      REAL(8)                  :: angle
      REAL(8),DIMENSION(0:pts) :: x,y
      DO i=0,pts                                                            
         angle=DBLE(i)/11.46D0                                       
         x(i)=xcent+rx*COS(angle)                                       
         y(i)=ycent+ry*SIN(angle)                                       
      ENDDO                                                                  
      IF (fill) THEN
        CALL RLAREA(x,y,72)
      ELSE
        CALL SHDPAT(0)             ! Shading pattern
        CALL RLAREA(x,y,72)                                                    
        CALL SHDPAT(16)            ! Shading pattern
      ENDIF                                              
      END SUBROUTINE pgcirc
!-----------------------------------------------------------------------        
      SUBROUTINE plot_fill(n,bfy,yh,yl) ! Plot block of data
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)              :: n      ! Number of years
      INTEGER,INTENT(IN)              :: bfy    ! First X value
      REAL(8),DIMENSION(n),INTENT(IN) :: yh,yl  ! Yhigh and Ylow 
      REAL(8),DIMENSION(n)            :: wkb
      INTEGER                         :: i=1
      wkb(1:n)=(/(DBLE(i),i=bfy,bfy+n-1)/)
      CALL SHDCRV(wkb(1:n),yh,n,wkb(1:n),yl,n)
      RETURN 
      END SUBROUTINE plot_fill
!------------------------------------------------------------------------
      SUBROUTINE open_ps(plot,chno)   ! Open postscript file
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: plot ! Plot number
      INTEGER,INTENT(IN) :: chno ! Character number
      CALL DISFIN()
      CALL METAFL('EPS')
      IF     (plot.GT.9) THEN
        WRITE(figm(chno:chno+5),'(I2,".eps")') plot
      ELSEIF (plot.GT.0) THEN
        WRITE(figm(chno+1:chno+5),'(I1,".eps")') plot
      ENDIF 
      CALL SETFIL(TRIM(figm))
      CALL plot_psbeg() 
      RETURN 
      END SUBROUTINE open_ps
!---------------------------------------------------------------
      SUBROUTINE norm1(n,dat,dsd,mn,sd)         ! Normalise data    
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n     ! Years
      REAL(8),DIMENSION(n),INTENT(IN)  :: dat   ! Data
      REAL(8),DIMENSION(n),INTENT(OUT) :: dsd   ! Normalised data
      REAL(8),INTENT(OUT)              :: mn,sd ! Mean and S.Dev
      mn=SUM(dat)/DBLE(n)
      sd=SQRT(MAX(SUM((dat-mn)**2),0.001D0)/DBLE(MAX(n-1,1)))
      dsd=(dat-mn)/sd
      RETURN
      END SUBROUTINE norm1
!----------------------------------------------------------------------
      SUBROUTINE norm2(n,dsum,dsq,cnt,dmn,dsd) ! Normalise series of data    
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n     ! Years
      REAL(8),DIMENSION(n),INTENT(IN)  :: dsum  ! Sums of data
      REAL(8),DIMENSION(n),INTENT(IN)  :: dsq   ! Sums of squares
      INTEGER,DIMENSION(n),INTENT(IN)  :: cnt   ! Data counts
      REAL(8),DIMENSION(n),INTENT(OUT) :: dmn   ! Mean of data
      REAL(8),DIMENSION(n),INTENT(OUT) :: dsd   ! Standard deviations
      REAL(8),DIMENSION(n)             :: wk,mn
      wk=DBLE(MAX(cnt,1)) ; mn=dsum/wk ; wk=DBLE(MAX(cnt-1,1))
      dsd=SQRT(MAX(dsq-mn*dsum,0.001D0)/wk) ; dmn=mn 
      RETURN
      END SUBROUTINE norm2
!----------------------------------------------------------------------
      SUBROUTINE ea_bfmsort(fs)  ! Create logicals for fast/slow,
      IMPLICIT NONE               ! young/old and big/small
      LOGICAL,INTENT(OUT),DIMENSION(mxs,6) :: fs 
      INTEGER  :: i
      fs=FA
      tst=4 ; CALL tree_sort()                 ! Fast and slow
      DO i=1,nc/2 ; fs(tre(nc+1-i),1)=TR ; fs(tre(i),2)=TR ; ENDDO 
      tst=2 ; CALL tree_sort()                 ! Old and young
      DO i=1,nc/2 ; fs(tre(nc+1-i),4)=TR ; fs(tre(i),3)=TR ; ENDDO 
      tst=3 ; CALL tree_sort()                 ! Big and small
      DO i=1,nc/2 ; fs(tre(nc+1-i),6)=TR ; fs(tre(i),5)=TR ; ENDDO 
      RETURN
      END SUBROUTINE ea_bfmsort
!------------------------------------------------------------------------
      SUBROUTINE read_met()   ! Program reads and stores met data
      IMPLICIT NONE
      INTEGER       :: j,p,r,s
      CHARACTER(60) :: fnam
      CALL read_open(27,"clim.fil") ; met=0.D0 ; mcf=0
      JD: DO j=1,maxm/12
        READ(27,'(A60)',IOSTAT=ios) metn(j) 
        IF (ios.EQ.-1) EXIT JD            ! End of file
        IF (io_err("Read","clim.fil")) STOP
        r=(j-1)*12+1 ; s=r+11 ; fnam=TRIM(metn(j))
        CALL read_metf(fnam) ; p=myr(r)
      ENDDO JD
      CLOSE(27) ; mfil=MIN(j-1,maxm/12) ; fm=MINVAL(mfy(1:mcf))
      lm=MAXVAL(mly(1:mcf)) ; chs=1 ; chc(1)=1
      RETURN
      END SUBROUTINE read_met
!--------------------------------------------------------------
      SUBROUTINE read_metf(fnam)       ! Program reads and stores met data
      IMPLICIT NONE   ! Missing values = -999
      INTEGER      :: i,j,jj,k,m,p,q,r,s
      CHARACTER(*) :: fnam
      CALL read_open(56,fnam)
      READ(56,'(56X,2I4)',iostat=ios) p,q ! Read header first/last year
      IF (io_err("Read ",fnam)) STOP
      r=mcf+1 ; s=r+11 ; mcf=s
      mfy(r:s)=p ; mly(r:s)=q ; p=q-p+1
      myr(r:s)=p ; m=LEN(TRIM(fnam))
      JD: DO jj=m,1,-1 ; IF (fnam(jj:jj).EQ.".") EXIT JD ; ENDDO JD
      KD: DO k=m,1,-1
        IF (fnam(k:k).EQ."/".OR.fnam(k:k).EQ.CHAR(92)) EXIT KD
      ENDDO KD
      k=MAX(k+1,1) ; IF (jj.LE.1) jj=m
      DO i=1,12 ; mnam(r-1+i)=fnam(k:jj)//mth(i) ; ENDDO
      DO i=1,p            ! Read each year
        READ(56,'(4X,12F5.0)',iostat=ios) met(i,r:s)
        IF (io_err("Read ",fnam)) STOP
      ENDDO
      CLOSE(56) ; p=myr(r) ; okm(1:p,r:s)=FA
      DO i=r,s
        J1: DO j=1,p-1   ! Remove leading -999
          IF (NINT(met(j,i)).NE.-999) EXIT J1
        ENDDO J1
        IF (j.GT.1) THEN  
          mfy(i)=mfy(i)-1+j ; myr(i)=myr(i)-j+1
          met(1:myr(i),i)=met(j:p,i) 
        ENDIF
        J2: DO j=p,2,-1  ! Remove trailing -999
          IF (NINT(met(j,i)).NE.-999) EXIT J2
        ENDDO J2
        IF (j.LT.p) THEN   
          mly(i)=mly(i)-p+j ; myr(i)=myr(i)-p+j
        ENDIF
        j=myr(i) ; okm(1:j,i)=(NINT(met(1:j,i)).NE.-999) 
      ENDDO
      IF (fnam(jj-1:jj-1).EQ."T") &
        WHERE (okm(1:p,r:s)) met(1:p,r:s)=met(1:p,r:s)*0.1D0 
      RETURN
      END SUBROUTINE read_metf
!--------------------------------------------------------------
      SUBROUTINE bootadjust() ! Bootstrap smoothing adjustment
      IMPLICIT NONE                      
      INTEGER                :: i,m,r
      REAL(8),DIMENSION(mxy) :: wka
      REAL(8)                :: aut
      r=cyr(cf) ; m=MAX(r/40,20)                  ! Spline stiffness m
      CALL cov(crn(1:r-1,cf),crn(2:r,cf),r-1,aut) ! Autocorrelation aut
      DO i=cf,cf+2                  ! Smooth CRN & limits with spline
        CALL splinet(r,crn(1:r,i),m,crn(1:r,i))
      ENDDO
      wka(1:r)=SQRT((1.D0+aut)/((1.D0-aut)*DBLE(m))) ! Annual adjustment
      DO i=cf+1,cf+2       ! Adjust lower and upper limits
        crn(1:r,i)=crn(1:r,cf)+(crn(1:r,i)-crn(1:r,cf))*wka(1:r) 
      ENDDO
      RETURN 
      END SUBROUTINE bootadjust
!----------------------------------------------------------------
      SUBROUTINE plot_boot(col1,bfy)  ! Creates output displays
      IMPLICIT NONE                 
      INTEGER,INTENT(IN)     :: col1  ! Colour
      INTEGER,INTENT(IN)     :: bfy   ! First year to plot
      INTEGER                :: i,j,k,m,n
      REAL(8),DIMENSION(mxy) :: wka,wkb
      CALL SETCLR(col1) ; n=125 ; k=bfy ; m=k-cfy(cf)+1
      DO i=1,(cyr(cf)-m)/n+1    
        DO j=1,n
          wka(j)=crn(m-1+j,cf+1) ; wka(j+n)=crn(n+m-j,cf+2)
          wkb(j)=DBLE(k-1+j) ; wkb(j+n)=DBLE(n+k-j)
        ENDDO
        CALL RLAREA(wkb(1:n+n),wka(1:n+n),n+n)
        k=k+n-1 ; m=m+n-1 ; n=MIN(n,cly(cf)-k+1)
      ENDDO  
      RETURN 
      END SUBROUTINE plot_boot
!------------------------------------------------------------------------
      SUBROUTINE boots()
! Tree indices in array dx(1:mxd), nc=tree count, cf=current chronology
! fy(1:nc), ly(1:nc), yr(1:nc), ad(1:nc)=Address of first ring each tree
! mxy=max chronology years, mxs=max number of trees. cfy(cf), cly(cf), 
! cyr(cf), crn(1:mxy,cf)=chronology, crn(1:mxy,cf+1)=lower confidence 
! limits, crn(1:mxy,cf+2)=upper confidence limits
      IMPLICIT NONE
      INTEGER,PARAMETER        :: nsims=1000 ! Generated bootstrap samples
      REAL(8),DIMENSION(mxs)   :: boot       ! Tree indices for year
      REAL(8),DIMENSION(nsims) :: simval     ! Samples - mean values
      REAL(8)                  :: dlow,dhi   ! Bootstrap confidence limits
      INTEGER                  :: mvar       ! Count of tree indices for year
      INTEGER                  :: i,j
      INTEGER,DIMENSION(12)    :: seed
      REAL(8)                  :: acc,xave
      seed=12848 ; CALL random_seed(put=seed) ! Initialise random sequence
      DO i=cfy(cf),cly(cf)                ! For each chronology year
!        IF (MOD(i,20).EQ.0) WRITE(19,*) i
        mvar=0
        DO j=1,nc                         ! For each tree
          IF (fy(j).LE.i.AND.ly(j).GE.i) THEN ! Tree has ring for year
            mvar=mvar+1 ; boot(mvar)=dx(i-fy(j)+ad(j))
          ENDIF
        ENDDO
        j=i-cfy(cf)+1
        IF (mvar.GE.6) THEN      ! Limited count of 6
! generate the boostrap psuedo-statistics and return all of them in msim
          CALL bootstrap(boot,simval,mvar,nsims)
! Calculate the acceleration constants for the data set using the jacknife 
! estimation of frangos and schucany comp. stat. & data anal. 9 (1990) 271-281
          xave=SUM(boot(1:mvar))/DBLE(mvar)
          boot(1:mvar)=boot(1:mvar)-xave
          acc=SUM(boot(1:mvar)**3)/(6.D0*(SUM(boot(1:mvar)**2)**1.5D0))
          CALL monte (xave,nsims,acc,simval,dlow,dhi)
! BCa or accelated confidence limits are calculated as they are most general
! and produce asymmetric limits eg those for the non-zero correlation coefficient
          crn(j,cf+1)=dlow       ! Lower limit
          crn(j,cf+2)=dhi        ! Upper limit
        ELSE                     ! Default to CRN value
          crn(j,cf+1)=crn(j,cf)  
          crn(j,cf+2)=crn(j,cf)  
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE boots
!--------------------------------------------------------------
      SUBROUTINE bootstrap(boot,simval,mvar,nsims)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                   :: mvar,nsims
      REAL(8),DIMENSION(mvar),INTENT(IN)   :: boot
      REAL(8),DIMENSION(nsims),INTENT(OUT) :: simval ! Means of mvar samples 
      INTEGER :: i,j
      REAL(8) :: ur,vr
      simval=0.D0 ; vr=DBLE(mvar)
      DO j=1,nsims
        DO i=1,mvar             ! Sum of mvar random samples
          CALL RANDOM_NUMBER(ur)
          simval(j)=simval(j)+boot(MIN(INT(ur*DBLE(mvar))+1,mvar))
        ENDDO
        simval(j)=simval(j)/vr  ! Arithmetic mean of samples
      ENDDO
      RETURN
      END SUBROUTINE bootstrap
!--------------------------------------------------------------
      SUBROUTINE monte(xave,nsims,acc,simval,dlow,dhi)
      IMPLICIT NONE
! As 214 appl. stat. (1985) vol. 34, no.3. Sets up monte carlo confidence
! intervals uses function ppnd - 111 and alnorn as 66. User supplies the
! estimate xave of a parameter, and monte obtains the confidence intervals
      INTEGER,INTENT(IN)                     :: nsims   
      REAL(8),INTENT(IN)                     :: xave,acc
      REAL(8),INTENT(OUT)                    :: dlow,dhi
      REAL(8),DIMENSION(nsims),INTENT(INOUT) :: simval
      INTEGER                                :: i,j,limita1,limita2,m    
      INTEGER,DIMENSION(1)                   :: k
! Calculate bias adjustment, so that sorting limits known?
      CALL biasad(xave,nsims,acc,simval,limita1,limita2,m)
      IF (limita1.LE.0.OR.limita1.GT.nsims) THEN
        dlow=simval(1)
      ELSE
        DO i=1,limita1  ! Select the nth lowest value       
          k=MINLOC(simval(1:nsims)) ; j=k(1) 
          dlow=simval(j) ; simval(j)=xave
        ENDDO
      ENDIF
      IF (limita2.LE.0.OR.limita2.GT.nsims) THEN
        dhi=simval(nsims)
      ELSE
        DO i=nsims,limita2,-1  ! Select the nth highest value       
          k=MAXLOC(simval(1:nsims)) ; j=k(1) 
          dhi=simval(j) ; simval(j)=xave
        ENDDO
      ENDIF
!      WRITE(19,'(4F10.3,4I6)') xave,acc,dlow,dhi,limita1,limita2,m
      RETURN
      END SUBROUTINE monte
!-------------------------------------------------------------------
      SUBROUTINE biasad (xave,nsims,acc,simval,limita1,limita2,m)
      IMPLICIT NONE
! Finds adjustment required to implement bias-corrected percentile
! method uses functions alnorm and ppnd algorithms as 66 and as 111
      INTEGER,INTENT(IN)                  :: nsims
      REAL(8),INTENT(IN)                  :: xave,acc
      REAL(8),DIMENSION(nsims),INTENT(IN) :: simval
      INTEGER,INTENT(OUT)                 :: limita1,limita2
      REAL(8),PARAMETER :: half=0.5D0,two=2.D0,hun=100.D0,conf=95.D0
      INTEGER     :: j,k,m  ! Confidence limit 0.0>conf<100.0
      REAL(8)     :: z,zed,fn1,xzed1,xzed2
      z=-ppnd(half*(hun-conf)/hun)
      j=0 ; k=0
! m=(count values < xave) + (half count values = xave), rounded 
      DO m=1,nsims
        IF (simval(m).LT.xave) THEN
          j=j+1 ; k=k+1
        ELSEIF (simval(m).EQ.xave) THEN
          k=k+1
        ENDIF
      ENDDO  
      m=(j+k+1)/2
      IF (m.LE.0.OR.m.GE.nsims) THEN
        CALL out_err("Biasad - Sim values below/above average ??")
        STOP
      ENDIF
      zed=two*ppnd(DBLE(m)/DBLE(nsims))
      fn1=DBLE(nsims+1)
! Calculate limits from the acceleration constant
      xzed1=zed/two+(zed/two-z)/(1.D0-acc*(zed/two-z))
      xzed2=zed/two+(zed/two+z)/(1.D0-acc*(zed/two+z))
      limita1=NINT(fn1*alnorm(xzed1,.false.))
      limita2=NINT(fn1*alnorm(xzed2,.false.))
      RETURN
      END SUBROUTINE biasad
!--------------------------------------------------------------------
      REAL(8) FUNCTION ppnd(p)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: p
! Algorithm as 111  appl. stat. (1977), vol. 26, no.1. Produces normal
! deviate corresponding to lower tail area of p real version for 
! eps = 2**(-31). Functions  abs,alog and sqrt are used.

      REAL(8),PARAMETER :: zero=1.D-15,split=0.42D0,half=0.5D0,one=1.0D0
      REAL(8),PARAMETER :: a0=   2.50662823884D0
      REAL(8),PARAMETER :: a1= -18.61500062529D0
      REAL(8),PARAMETER :: a2=  41.39119773534D0
      REAL(8),PARAMETER :: a3= -25.44106049637D0
      REAL(8),PARAMETER :: b1=  -8.47351093090D0
      REAL(8),PARAMETER :: b2=  23.08336743743D0
      REAL(8),PARAMETER :: b3= -21.06224101826D0
      REAL(8),PARAMETER :: b4=   3.13082909833D0
! Hash sum ab 143.70383558076
      REAL(8),PARAMETER :: c0=  -2.78718931138D0
      REAL(8),PARAMETER :: c1=  -2.29796479134D0
      REAL(8),PARAMETER :: c2=   4.85914127135D0
      REAL(8),PARAMETER :: c3=   2.32121276858D0
      REAL(8),PARAMETER :: d1=   3.54388924762D0
      REAL(8),PARAMETER :: d2=   1.63706781897D0
! Hash sum cd 17.43746520924
      REAL(8)           :: q,r,pl
      pl=DBLE(p) ; q=pl-half
      IF (DABS(q).LE.split) THEN
        r=q*q
        ppnd=q*(((a3*r+a2)*r+a1)*r+a0)/  &
            ((((b4*r+b3)*r+b2)*r+b1)*r+one)
        RETURN
      ENDIF
      r=pl
      IF (q.GT.zero) r=one-pl
      IF (r.GT.zero) THEN
        r=DSQRT(-DLOG(r))
        ppnd=(((c3*r+c2)*r+c1)*r+c0)/((d2*r+d1)*r+one)
        IF (q.LT.1.D-15) ppnd=-ppnd
        RETURN
      ENDIF
      CALL out_err("Ppnd - r less than zero ??")
      STOP
      END FUNCTION ppnd
!-------------------------------------------------------------------
      REAL(8) FUNCTION alnorm(x,upper)
      IMPLICIT NONE
! Algorithm as 66  appl. stat. (1973) vol. 22, no.3. Evaluates the 
! tail area of the standardised normal curve from x to infinity if
! upper is .true. or from minus infinity to x if upper is .false.
      REAL(8),INTENT(IN) :: x
      REAL(8)            :: z,y,al
      LOGICAL,INTENT(IN) :: upper
      LOGICAL            :: up
! ltone and utzero must be set to suit the particular computer 
      REAL(8),PARAMETER  :: zero=1.D-15, half=0.5D0, one=1.D0, con=1.28D0
      REAL(8),PARAMETER  :: ltone=7.D0, utzero=18.66D0
      up=upper ; z=DBLE(x)
      IF (z.LT.zero) THEN ; up=.NOT.up ; z=-z ; ENDIF
      IF ((z.LE.ltone).OR.(up.AND.z.LE.utzero)) THEN
        y=half*z*z
        IF (z.LE.con) THEN
          al=half-z*(.398942280444D0-.399903438504D0*y/ &
                 (y+5.75885480458D0-29.8213557808D0 /       &
                 (y+2.62433121679D0+48.6959930692D0 /       &
                 (y+5.92885724438D0))))
        ELSE
          al=.398942280385D0*DEXP(-y)/              &
                 (z+3.805D-8+1.00000615302D0/      &
                 (z+3.98064794D-4+1.98615381364D0/ &
                 (z-.151679116635D0+5.29330324926D0/ &
                 (z+4.8385912808D0-15.1508972451D0/  &
                 (z+.742380924027D0+30.789933034D0/  &
                 (z+3.99019417011D0))))))
        ENDIF
      ELSE
        al=zero
      ENDIF
      IF (.not.up) al=one-al
      alnorm=al
      RETURN
      END FUNCTION alnorm
!-------------------------------------------------------------
      SUBROUTINE randgen(ref1) ! Create random trees
      IMPLICIT NONE                
      INTEGER,INTENT(IN)    :: ref1
      INTEGER               :: i,j,p,q,r
      REAL(8)               :: rn,inc,rn1,inc1
      CHARACTER(11)         :: fnam
      INTEGER,DIMENSION(12) :: seed
      seed=11 ; CALL random_seed(put=seed) ! Initialise random sequence
      ad(1)=1 ; cly(cf)=2000
      IF (MOD(ref1-1,32).LE.15) THEN     ! Modern chronology
        nc=101 ; cfy(cf)=1700 ; cyr(cf)=301 ; wnam(1)(3:3)="M"
      ELSE                              ! Sub-fossil chronology
        nc=301 ; cfy(cf)=1300 ; cyr(cf)=701 ; wnam(1)(3:3)="S"
      ENDIF
      DO i=1,nc                         ! Tree names
        IF (i.LE.99) THEN
          WRITE(nam(i),'("RAND",I2,"  ")') i
        ELSE
          WRITE(nam(i),'("RAND",I3," ")') i
        ENDIF
        fy(i)=cfy(cf)+(i-1)*2           ! Trees every 2 years
      ENDDO 
      IF (MOD(ref1-1,16).LE.7) THEN      ! Same length
        yr(1:nc)=301 ; fy(nc-100:nc)=fy(nc-100) ; wnam(1)(4:4)="S"
      ELSE                              ! Variable length
        yr(nc-100:nc)=2001-fy(nc-100:nc) ; wnam(1)(4:4)="V"

        DO i=1,nc-101
          CALL RANDOM_NUMBER(rn) ; yr(i)=INT(rn*201.D0)+100
        ENDDO
      ENDIF  
      ly(1:nc)=fy(1:nc)+yr(1:nc)-1
      DO i=1,nc ; ad(i+1)=ad(i)+yr(i) ; ENDDO
      pth(1:nc)=fy(1:nc)-1 ; pthr(1:nc)=0.1D0
      r=cyr(cf) ; crn(1:r,cf)=1.D0
      IF (MOD(ref1-1,128).LE.63) THEN
        SELECT CASE (MOD(ref1-1,4))       ! Signal at end
        CASE (0) ; wnam(1)(5:10)="EUstep" ! Step up
!         IF (ref.EQ.45) THEN       ! Sine wave option
!           DO i=1,r
!             crn(i,cf)=1.0+0.1*SIND(DBLE(i)*9.0)
!           ENDDO
!         ELSE
          p=1951-cfy(cf) ; crn(p:r,cf)=1.2D0              
!         ENDIF
        CASE (1) ; wnam(1)(5:10)="EUslop"   ! Slope up
          p=1901-cfy(cf) 
          crn(p:r,cf)=(/(1.D0+DBLE(i-p+1)/500.D0,i=p,r)/) 
        CASE (2) ; wnam(1)(5:10)="EDstep"   ! Step down
          p=1951-cfy(cf) ; crn(p:r,cf)=0.8D0              
        CASE (3) ; wnam(1)(5:10)="EDslop"   ! Slope down
          p=1901-cfy(cf) 
          crn(p:r,cf)=(/(1.D0-DBLE(i-p+1)/500.D0,i=p,r)/) 
        ENDSELECT                    
      ELSE
        SELECT CASE (MOD(ref1-1,4))         ! Signal at start
        CASE (0) ; wnam(1)(5:10)="SUstep"   ! Step up
          crn(1:50,cf)=1.2D0                     
        CASE (1) ; wnam(1)(5:10)="SUslop"   ! Slope up
          crn(1:100,cf)=(/(1.0_8+DBLE(101-i)/500.0_8,i=1,100)/) 
        CASE (2) ; wnam(1)(5:10)="SDstep"   ! Step down
          crn(1:50,cf)=0.8D0                   
        CASE (3) ; wnam(1)(5:10)="SDslop"   ! Slope down
          crn(1:100,cf)=(/(1.D0-DBLE(101-i)/500.D0,i=1,100)/) 
        ENDSELECT                    
      ENDIF
      IF (MOD(ref1-1,64).LT.31) THEN         ! Linear
        CALL RANDOM_NUMBER(rn1) ; CALL RANDOM_NUMBER(inc1)
        DO i=1,nc
          rn=(rn1-0.5D0)/2.5D0+1.2D0       !  1.0 to 1.4 start value
          inc=(inc1-0.5D0)/2.5D0+0.5D0     !  0.4 to 0.8 end value
          p=ad(i) ; q=ad(i)+yr(i)-1 ; r=fy(i)-cfy(cf)+1-p
          inc=(rn-inc)/DBLE(yr(i))    
          DO j=p,q                          ! Measurement = line x signal
            x(j)=rn*crn(j+r,cf) ; rn=rn-inc
          ENDDO
        ENDDO
        wnam(1)(2:2)="L"
      ELSE                                  ! Exponential
        DO i=1,nc
          r=yr(i) ; p=ad(i)-1 ; q=fy(i)-cfy(cf) 
          CALL RANDOM_NUMBER(rn) ; CALL RANDOM_NUMBER(inc)
          inc=(inc-0.5D0)/2.5D0+1.2D0     ! 1.0 to 1.4 start value
          rn=DBLE(r+1)+rn*DBLE(r/2)       ! Tree lives longer
          DO j=1,r
            x(p+j)=(inc*EXP(-DBLE(j)/rn)+0.1D0)*crn(q+j,cf)
          ENDDO
        ENDDO
        wnam(1)(2:2)="E"
      ENDIF
      IF (MOD(ref1-1,8).LE.3) THEN           ! No Noise
        wnam(1)(1:1)="N"
      ELSE                                  ! Noise +/- 20%
        DO i=1,ad(nc+1)-1
          CALL RANDOM_NUMBER(rn) ; x(i)=x(i)*((rn-0.5D0)/2.5D0+1.D0)
        ENDDO
        wnam(1)(1:1)="Y"   
      ENDIF
      CALL det_crnfy() ; xok=TR ; wnam(1)(11:14)=".raw"
      WRITE(fnam,'("Rand",I3,".raw")') ref1
      CALL write_raw(fnam,x) 
      RETURN 
      END SUBROUTINE randgen
!------------------------------------------------------------------
      SUBROUTINE read_metsav(msnam) ! Reads saved met files (compact)
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN) :: msnam
      CHARACTER(12)           :: frm
      INTEGER                 :: i,j  
      REAL(8),DIMENSION(mmx)  :: wka
      OPEN(24,FILE=msnam,IOSTAT=ios,STATUS="OLD")
      IF (io_err("Open",msnam)) RETURN
      ID: DO i=1,maxm                ! For each file
        READ(24,'(I4,2I6,2X,A20)',IOSTAT=ios) &
          mcf,mfy(i),mly(i),mnam(i)
        IF (io_err("Read1",msnam)) STOP
        IF (mcf.NE.i) EXIT ID
      ENDDO ID
      WRITE(frm,'("(I6,",I3,"F7.1)")')  mcf
        fm=MINVAL(mfy(1:mcf)) ; lm=MAXVAL(mly(1:mcf))
      DO i=fm,lm


        READ(24,IOSTAT=ios,FMT=frm) j,wka(1:mcf) 
        IF (io_err("Read2",msnam)) STOP
        DO j=1,mcf
          IF (mfy(j).LE.i.AND.mly(j).GE.i) met(i-mfy(j)+1,j)=wka(j)
        ENDDO
      ENDDO
      CLOSE(24) ; chc(1)=1 ; chs=1 
      okm=(met.NE.-999.D0) 
      RETURN
      END SUBROUTINE read_metsav
!-----------------------------------------------------------------
      SUBROUTINE write_crnsav(msnam) ! Sub saves index files (Tucson)
      IMPLICIT NONE  
      CHARACTER(*),INTENT(IN) :: msnam
      INTEGER                 :: i,m,r
      OPEN(23,FILE=msnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",msnam)) RETURN
      DO i=1,chs                  ! For each file
        m=chc(i) ; r=cyr(m)
        WRITE(23,'(20X,A20)',IOSTAT=ios) wnam(cf)
        CALL write_ind(r,cfy(m),crn(1:r,m),num(1:r,m),"CRN",wnam(m))
      ENDDO
      CLOSE(23)
      RETURN
      END SUBROUTINE write_crnsav
!-----------------------------------------------------------------
      SUBROUTINE write_metsav(msnam) ! Sub saves met files (column)
      IMPLICIT NONE   
      CHARACTER(*),INTENT(IN) :: msnam
      INTEGER                 :: i,j,m,p,q
      CHARACTER(12)           :: frm
      REAL(8),DIMENSION(mmx)  :: wka
      OPEN(24,FILE=msnam,IOSTAT=ios,STATUS="REPLACE")
      IF (io_err("Open",msnam)) RETURN
      p=100000 ; q=-100000
      DO i=1,chs                  ! For each file
        m=chc(i) ; p=MIN(mfy(m),p) ; q=MAX(mly(m),q)
        WRITE(24,'(I4,2I6,2X,A20)',IOSTAT=ios) i,mfy(m),mly(m),mnam(m)
        IF (io_err("write1",msnam)) STOP
      ENDDO
      WRITE(24,'(I4,2I6,2X,A20)',IOSTAT=ios) chs,p,q,"  "
      WRITE(frm,'("(I6,",I3,"F7.1)")')  chs
      DO i=p,q         ! First year to last year
        DO j=1,chs     ! For each series
          m=chc(j) 
          IF (mfy(m).LE.i.AND.mly(m).GE.i) THEN
            IF (okm(i-mfy(m)+1,m)) THEN  ! In range and valid
              wka(j)=met(i-mfy(m)+1,m)
            ELSE 
              wka(j)=-999.0D0
            ENDIF
          ELSE
            wka(j)=-999.0D0
          ENDIF
        ENDDO
        WRITE(24,IOSTAT=ios,FMT=frm) i,wka(1:chs) 
        IF (io_err("write2",msnam)) STOP
      ENDDO
      CLOSE(24)
      RETURN
      END SUBROUTINE write_metsav
!-----------------------------------------------------------------
      SUBROUTINE TRENDW(n,dx1,dx2,ns,sl,yi,xi) ! Weighted best fit line
      IMPLICIT NONE
      INTEGER,INTENT(IN)               :: n    ! NUMBER OF VALUES IN SERIES
      REAL(8),DIMENSION(N),INTENT(IN)  :: dx1  ! TIME SERIES ARRAY
      REAL(8),DIMENSION(N),INTENT(OUT) :: dx2  ! CALCULATED TREND LINE VALUES 
      INTEGER,DIMENSION(N),INTENT(IN)  :: ns   ! Series counts
      REAL(8),INTENT(OUT)              :: sl   ! SLOPE OF TREND LINE
      REAL(8),INTENT(OUT)              :: yi   ! Y-AXIS INTERCEPT
      REAL(8),INTENT(OUT)              :: xi   ! X-AXIS INTERCEPT
      REAL(8),DIMENSION(n)             :: wk
      REAL(8) :: xy,x,y,xx,an
      INTEGER :: i=1
      an=DBLE(SUM(ns(1:n))) ; wk(1:n)=(/(DBLE(i),i=1,n)/)
      xy=SUM(dx1(1:n)*wk(1:n)*DBLE(ns(1:n)))
      x=SUM(wk(1:n)*DBLE(ns(1:n)))
      y=SUM(dx1(1:n)*DBLE(ns(1:n)))
      xx=SUM((wk(1:n)**2)*DBLE(ns(1:n)))
      sl=(an*xy-x*y)/(an*xx-x**2)
      yi=y/an-sl*x/an ; dx2(1:n)=wk(1:n)*sl+yi
      IF (sl.NE.0.D0) xi=(-yi)/sl 
      RETURN
      END SUBROUTINE TRENDW
!----------------------------------------------------------
      SUBROUTINE det_slope(r,dat,cnt,slop) ! Count weighted slope
      IMPLICIT NONE   
      INTEGER,INTENT(IN)              :: r         ! Length of series
      REAL(8),DIMENSION(r),INTENT(IN) :: dat        ! Values of series
      INTEGER,DIMENSION(r),INTENT(IN) :: cnt        ! Count each year
      REAL(8),INTENT(OUT)             :: slop       ! Output slope
      INTEGER                         :: i=1
      REAL(8)                         :: z1        ! Distance from centre
      REAL(8),DIMENSION(r)            :: k1,k2 
      k1(1:r)=(/(DBLE(i),i=1,r)/)                  ! Distance 
      k2(1:r)=DBLE(cnt(1:r))                       ! Counts
      z1=SUM(k2(1:r)*k1(1:r))/DBLE(SUM(cnt(1:r)))  ! Count weighted centre
      k1(1:r)=z1-k1(1:r) ; k2(1:r)=k2(1:r)*k1(1:r) ! Value * counts * dist
      slop=-SUM(dat(1:r)*k2(1:r))/SUM(k2(1:r)*k1(1:r)) ! Slope
      RETURN 
      END SUBROUTINE det_slope
!----------------------------------------------------------------
      SUBROUTINE crn_slope(rfy,rly,slop)  ! Reset slope over given period
      IMPLICIT NONE                      
      INTEGER,INTENT(IN)     :: rfy,rly   ! Period of slope
      REAL(8),INTENT(IN)     :: slop      ! Slope change 
      INTEGER                :: i=1,p,q,r
      REAL(8),DIMENSION(mxy) :: k1,k2     
      REAL(8)                :: z1,m1      ! Count weighted centre, mean
      r=cyr(cf) ; p=MAX(rfy-cfy(cf)+1,1) ; q=MIN(rly-cfy(cf)+1,r)  
      m1=DBLE(SUM(num(1:r,cf)))    ; k1(1:r)=(/(DBLE(i),i=1,r)/)  
      k2(1:r)=DBLE(num(1:r,cf)) ; z1=SUM(k2(p:q)*k1(p:q))/m1 
      crn(1:r,cf)=(1.0_8-slop*(z1-k1(1:r)))*crn(1:r,cf)
      m1=SUM(k2(1:r)*crn(1:r,cf))/m1      ! Mean
      crn(1:r,cf)=crn(1:r,cf)/m1          ! Reset weighted mean to 1.0
      RETURN 
      END SUBROUTINE crn_slope
!------------------------------------------------------------------
      SUBROUTINE EPS_prep1(lseg)         ! Investigate RCS EPS 
      IMPLICIT NONE
      INTEGER,INTENT(IN)     :: lseg     ! Segment length  
      REAL(8),DIMENSION(mxd) :: qx       ! Raw data storage
      INTEGER :: i,j,k,p,q,r,s,t,u,v,nr1,nr2
      LOGICAL :: corok
      REAL(8) :: segr,mnk,mnr,mnr2,corr,corr2,ef1,ef2
      CALL detrend() ; r=cyr(1)          ! 1 = Chronology
      crn(1:r,2)=xcsd(1:r,mx)            ! 2 = SDev
      okc(1:r,2)=num(1:r,1).GT.3
      IF (idb.EQ.1) THEN   ! SDev scaled if not Normal
        WHERE (okc(1:r,2).AND.crn(1:r,1).GT.0.D0)  ! 3 = Scaled SDev
          crn(1:r,3)=crn(1:r,2)/crn(1:r,1) 
        ELSE WHERE
          crn(1:r,3)=0.D0
        END WHERE
      ELSE
        crn(1:r,3)=crn(1:r,2)
      ENDIF 
      i=ad(nc)+yr(nc)-1 ; qx(1:i)=dx(1:i) ! Save RCS Indices
      cf=4 ; idt=lseg ; CALL detrend()    ! 4 = Spline Chronology
      cf=1 ; crn(1:r,5)=xcsd(1:r,mx)      ! 5 = Spline SDev
      IF (idb.EQ.1) THEN   ! SDev scaled if not Normal
        WHERE (okc(1:r,2).AND.crn(1:r,4).GT.0.D0)  ! 6 = Spine Scaled SDev
          crn(1:r,6)=crn(1:r,5)/crn(1:r,4) 
        ELSE WHERE
          crn(1:r,6)=0.D0
        END WHERE
      ELSE
        crn(1:r,6)=crn(1:r,5)
      ENDIF 
      num(1:r,8)=0 ; num(1:r,11)=0
      J2: DO j=r,1,-1 ; IF (okc(j,2)) EXIT J2 ; ENDDO J2
      I2: DO i=1,r    ; IF (okc(i,2)) EXIT I2 ; ENDDO I2
      cfy(2)=i+lseg/2                  ! Range >3 start 
      cly(2)=j-(lseg+1)/2+1            ! Range >3 end 
      cyr(2)=cly(2)-cfy(2)+1
      WRITE(71,'(A20)') wnam(61)
      WRITE(71,'(A)') TRIM(cnam(61))
      WRITE(71,*)   
      WRITE(71,'("Chronology",2I5,"  >4 Range",2I5,"RBar Range  fy ",2I5)') &
        cfy(1),cly(1),i+cfy(1)-1,j+cfy(1)-1,cfy(2)+cfy(1)-1,cly(2)+cfy(1)-1
      WRITE(71,'("Scnt    = Count for Rbar of Spline")')
      WRITE(71,'("Rcnt    = Count for Rbar of RCS")')
      WRITE(71,'("sp RBar = Rbar of",I4,"-year Spline")') lseg
      WRITE(71,'("RC RBar = Rbar of RCS")')
      WRITE(71,'("sp EPS  = EPS for",I4,"-year Spline")') lseg
      WRITE(71,'("RC EPS  = EPS for RCS")')
      WRITE(71,'("efEPS   = Effective EPS for RCS")')
      WRITE(71,'("aefEPS  = Adjusted effective EPS for RCS")')
      WRITE(71,*)   
      WRITE(71,'(2A34)') "  Year Count  Scnt  Rcnt sp RBar R", &
        "C RBar sp EPS RC EPS  efECS aefEPS"
      DO k=cfy(2),cly(2)              ! Rbar 50% overlap
        IF (okc(k,2)) THEN            ! Need enough samples for SDev calculation
          mnr=0.D0 ; mnr2=0.D0 ; nr1=0 ; nr2=0
          p=k-lseg/2+cfy(1) ; q=p+lseg-1
          DO i=1,nc-1           ! First tree
            IF (fy(i).LE.p.AND.ly(i).GE.q) THEN     ! Tree in segment
              DO j=i+1,nc       ! Second tree
                IF (fy(j).LE.p.AND.ly(j).GE.q) THEN ! Tree in segment
                  u=ad(i)+p-fy(i) ; v=ad(i)+q-fy(i)
                  s=ad(j)+p-fy(j) ; t=ad(j)+q-fy(j)
                  CALL covmiss(dx(u:v),dx(s:t),xok(u:v).AND.xok(s:t), &
                    lseg,corr,corok)
                  IF (corok) THEN
                    mnr=mnr+corr ; nr1=nr1+1        ! Sum Spline correlations
                    CALL covmiss(qx(u:v),qx(s:t),xok(u:v).AND.xok(s:t), &
                      lseg,corr2,corok)
                    IF (corok) THEN
                      mnr2=mnr2+corr2 ; nr2=nr2+1   ! Sum RCS correlations
                    ENDIF
                  ENDIF 
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          p=k-lseg/2 ; q=p+lseg-1
          num(k, 8)=nr1                     ! Rbar count spline
          crn(k, 8)=mnr/DBLE(MAX(nr1,1))    ! Rbar for spline
          num(k,11)=nr2                     ! Rbar count RCS
          crn(k,11)=mnr2/DBLE(MAX(nr2,1))   ! Rbar for RCS
          segr=DBLE(COUNT(okc(p:q,2)))      ! Number of SDevs in window
          crn(k, 9)=SUM(crn(p:q,5))/segr    ! Mean SDev spline 
          crn(k,12)=SUM(crn(p:q,2))/segr    ! Mean SDev RCS
          crn(k,10)=SUM(crn(p:q,6))/segr    ! Mean Adj SDev spline 
          crn(k,13)=SUM(crn(p:q,3))/segr    ! Mean Adj SDev RCS
          crn(k,14)=crn(k,12)/crn(k,9)      ! Ratio SDev
          crn(k,15)=crn(k,13)/crn(k,10)     ! Ratio Adj SDev 
          mnk=DBLE(num(k,1))
          ef1=mnk/(crn(k,12)/crn(k,9))**2   ! eff count
          ef2=mnk/(crn(k,13)/crn(k,10))**2  ! adj eff count
          IF (crn(k,8).GT.0.01D0) THEN      ! Positive correlation
            crn(k,16)=crn(k,8)*mnk/(1.D0+(mnk-1.D0)*crn(k,8)) ! Spline EPS
          ELSE
            crn(k,16)=0.D0
          ENDIF 
          IF (crn(k,11).GT.0.01D0) THEN     ! Positive correlation
            crn(k,17)=crn(k,11)*mnk/(1.D0+(mnk-1.D0)*crn(k,11))  ! RCS EPS
            crn(k,18)=crn(k,11)*ef1/(1.D0+(ef1-1.D0)*crn(k,11))  ! RCS EPS
            crn(k,19)=crn(k,11)*ef2/(1.D0+(ef2-1.D0)*crn(k,11))  ! RCS EPS
          ELSE
            crn(k,17:19)=0.D0
          ENDIF 
          WRITE(71,'(4I6,2F8.2,4F7.2)') k+cfy(1)-1,num(k,1), &
            num(k,8),num(k,11),crn(k,8),crn(k,11),crn(k,16:19)
        ENDIF
      ENDDO     
      WRITE(71,'(" ")')   
      WRITE(71,'("CRN    = RCS chronology indices")')
      WRITE(71,'("Rdev   = RCS standard deviation")')
      WRITE(71,'("ScRdev = RCS scaled standard deviation")')
      WRITE(71,'("SCRN   = Spline chronology indices")')
      WRITE(71,'("Sdev   = Spline standard deviation")')
      WRITE(71,'("ScSdev = Spline scaled standard deviation")')
      WRITE(71,'("Year Count     CRN    Rdev  ScRdev    SCRN    SDev  ScSdev")') 
      DO i=1,cyr(1)
        WRITE(71,'(2I5,6F8.4)')  cfy(1)-1+i,num(i,2),crn(i,1:6) 
      ENDDO
      r=cyr(1) ; okc(1:r,8)=num(1:r,11).GE.1.AND.num(1:r,8).GE.1
      RETURN
      END SUBROUTINE EPS_prep1
!--------------------------------------------------------------
      END MODULE CRUSTUTIL 
