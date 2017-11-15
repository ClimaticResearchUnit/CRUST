! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see 
! the GNU General Public License.
      MODULE crustprocs    ! Arstan - with graphics and no stats
      USE rfigs
      CONTAINS   
!-----------------------------------------------------------------
      SUBROUTINE crust_setup()    ! Finds PC details and sets defaults
      IMPLICIT NONE
      INTEGER :: i,ra,rb,rc,rd,re
! Next madd free ????  
      madd(35)=pos2(cnx+100,1800,cnx+1200,1850)
      madd(36)=pos2(cnx+1200,1800,cnx+2600,1850)
      madd(37)=pos2(cnx+1200,1850,cnx+2600,1900)
      madd(38)=pos2(cnx+1200,1900,cnx+2600,1950)
      mtext(35)="Copyright (C) 2013 Thomas Melvin - Program CRUST"
      mtext(36)="This program comes with ABSOLUTELY NO WARRANTY. This is"
      mtext(37)="free software, and you are welcome to redistribute it "
      mtext(38)="under certain conditions (GNU General Public License)."
! RCS main menu
      ra=2800 ; rb=2885 ; rc=2560 ; rd=100
      madd(1)=pos2(rc,rd,rb,rd+45) ; rd=rd+60   
      mtext(1)="Chronology Options" 
      b( 1)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Transform
      b( 2)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Pith offset
      b( 3)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Tree Index calc
      b( 4)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Mean chronology
      b( 5)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Signal free
      b( 6)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! AR - STD, RES or ARS
      b( 7)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Distrib adjust
      b( 8)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Variance stabilise
      b( 9)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Detrend RCS curve
      b(10)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Single RCS
      b(11)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! RCS type
      b(12)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! Ring transform
      b(13)=posn(ra+5,b( 9)%y1,rb,b( 9)%y2,"",FA,TR)  ! Detrend Spline sl
      b(14:15)%ok=FA  ! b(14:15)=spare
      b(16)=posn(ra+5,b(10)%y1,rb,b(10)%y2,"",FA,TR)  ! Number of RCS curves
      b(17)%ok=FA     ! b(17) Not used
      b(18)=posn(rc,rd,ra,rd+45,"",TR,TR) ; rd=rd+60  ! CRNs mean option
      b(19)=posn(2060,95,2300,140,"No Sort",TR,TR)    ! Sort option
      madd(20)=pos2(1780,95,2060,140)   
      mtext(20)="Sort Trees by"
      b(20)=posn(ra+5,b( 5)%y1,rb,b( 5)%y2,"",FA,TR)  ! Max SF iterations      
      b(21)=posn(1980,155,2080,200,"  0",TR,TR)       ! Ignore 1st 'n'
      madd(21)=pos2(1780,155,1980,200)   
      mtext(21)="Ignore first" 
      madd(22)=pos2(2100,155,2200,200)   
      mtext(22)="rings" 
      b(22)%ok=FA     ! b(22) Not Used
      b(23)=posn(2700,30,2820,75,"Quit",TR,TR)   ! Quit option
      madd(23)=pos2(2450,30,2700,75)
      mtext(23)="Exit Program"
      b(24)=posn(1300,30,2200,75," ?",TR,TR)     ! Chronology name
      madd(24)=pos2(1100,30,1300,75)
      mtext(24)="Data File"
      ra=60 ; rb=320 ; rd=335 
      madd(2)=pos2(ra,rd,rb,rd+45)   
      mtext(2)="Menu Options" 
      madd(3)=pos2(rb,rd,rb+440,rd+45) ; rd=rd+80
      mtext(3)=" - does not exist?"
      b(25:31)%ok=FA  ! b(25:31)=spare
      b(32)=posn(ra,rd,rb,rd+45,"Calc CRN",TR,TR)    ; rd=rd+80  ! Calculate CRN
      b(33)=posn(ra,rd,rb,rd+45,"Save Data",TR,TR)   ; rd=rd+80  ! Save Data
      b(34)=posn(ra,rd,rb,rd+45,"Save Plot",TR,TR)   ; rd=rd+80  ! Save Plot
      b(35)=posn(ra,rd,rb,rd+45,"Display Opt",TR,TR) ; rd=rd+80  ! Display Opt
      b(36)=posn(ra,rd,rb,rd+45,"Save Params",TR,TR) ; rd=rd+80  ! Save Parameters
      b(37)=posn(ra,rd,rb,rd+45,"Select RCS",TR,TR)  ; rd=rd+80  ! Select RCS curves
      b(38)=posn(ra,rd,rb,rd+45,"Tree Menu",TR,TR)   ; rd=rd+80  ! Tree Menu
      b(39:40)%ok=FA     ! b(39:40) spare
!  Curve fitting options
      ra=2830 ; rb=2960 ; rc=2550 ; rd=1000
      madd(41)=pos2(rc,rd,rb,rd+45)   
      mtext(41)="Curve-fitting Options" 
      DO i=0,9
        rd=rd+60 
        b(i+41)=posn(rc,rd,ra,rd+45,"",TR,FA)
        b(i+41)%lab=idtlab(i)
      ENDDO
      b(51)=posn(ra+5,b(49)%y1,rb,b(49)%y2,"",FA,TR)  ! Spline n
      b(52)=posn(ra+5,b(50)%y1,rb,b(50)%y2,"",FA,TR)  ! Spline %n
      b(53:54)%ok=FA   ! b(53:54) spare      
! Tree Menu     
      b(55)=posn(2180,120,2420,165,"No Sort",TR,TR) ! Sort option
      madd(55)=pos2(1900,120,2180,165)   
      mtext(55)="Sort Trees by"
      madd(56)=pos2(1600,30,2200,75)   ! Chronology name
      mtext(56)="Chronology"
      madd(9)=pos2(1250,30,1600,75)
      mtext(9)="Current Chronology"
      b(56)=posn(960,60,1160,105,"Previous",TR,TR) ! Previous tree name
      b(57)=posn(960,120,1080,165,"Next",TR,TR) ! Next tree name
      b(58)=posn(540,90,880,135,"",TR,TR)     ! Core name
      madd(7)=pos2(180,90,540,135)
      mtext(7)="Current Tree Name"
      b(59)=posn(2800,30,2920,75,"Exit",TR,TR) ! Quit option
      madd(8)=pos2(2500,30,2800,75)
      mtext(8)="Leave Tree Menu"
      b(60)=posn(1440,90,1760,135,"Save Plot",TR,TR)   ! Save Plot
      b(61)=posn(1440,150,1760,195,"Display Opt",TR,TR) ! Display Opt
      b(62:65)%ok=FA   ! b(62:65) spare      
! Display Menu     
      b(66)=posn(2360,220,2840,265,"Save as Default",TR,TR) ! Save Parameters
      b(67)=posn(2760,30,2880,75,"Exit",TR,TR) ! Quit option
      madd(67)=pos2(2360,30,2760,75)
      mtext(67)="Leave Display Menu"
      ra=80 ; rb=520 ; rc=640 ; rd=200
      madd(27)=pos2(ra+120,rd,rc,rd+45) ; rd=rd+60  
      mtext(27)="RCS Curve Display Options" 
      b(68)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Mean RCS curve
      madd(68)=pos2(ra,rd,rb,rd+45)   
      mtext(68)="Mean RCS curve" ; rd=rd+60 
      b(69)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Smooth RCS curve
      madd(69)=pos2(ra,rd,rb,rd+45)   
      mtext(69)="Smooth RCS curve" ; rd=rd+60 
      b(70)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! RCS standard deviation
      madd(70)=pos2(ra,rd,rb,rd+45)   
      mtext(70)="RCS Stand. Dev." ; rd=rd+60 
      b(71)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! RCS standard error
      madd(71)=pos2(ra,rd,rb,rd+45)   
      mtext(71)="RCS Stan. Error" ; rd=rd+60 
      b(72)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! RCS counts
      madd(72)=pos2(ra,rd,rb,rd+45)   
      mtext(72)="RCS ring count" ; rd=rd+60 
      b(81)=posn(rb,rd,rc,rd+45,"",TR,TR)  ! RCS first display year
      madd(81)=pos2(ra,rd,rb,rd+45)   
      mtext(81)="RCS First Year" ; rd=rd+60 
      b(82)=posn(rb,rd,rc,rd+45,"",TR,TR)  ! RCS last display year
      madd(82)=pos2(ra,rd,rb,rd+45)   
      mtext(82)="RCS Last Year" ; rd=rd+120 
      madd(26)=pos2(ra+120,rd,rc,rd+45) ; rd=rd+60
      mtext(26)="Chronology Display Options" 
      b(73)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Chronology
      madd(73)=pos2(ra,rd,rb,rd+45)   
      mtext(73)="Chronology" ; rd=rd+60 
      b(74)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! CRN standard deviation
      madd(74)=pos2(ra,rd,rb,rd+45)   
      mtext(74)="CRN Stand. Dev." ; rd=rd+60 
      b(75)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! CRN counts
      madd(75)=pos2(ra,rd,rb,rd+45)   
      mtext(75)="CRN ring count" ; rd=rd+60 
      b(76)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! CRN smooth
      madd(76)=pos2(ra,rd,rb,rd+45)   
      mtext(76)="Smooth CRN" ; rd=rd+60 
      b(83)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! CRN first display year
      madd(83)=pos2(ra,rd,rb,rd+45)   
      mtext(83)="CRN First Year" ; rd=rd+60  
      b(84)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! CRN last display year
      madd(84)=pos2(ra,rd,rb,rd+45)   
      mtext(84)="CRN Last Year" ; rd=rd+60 
      b(85)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Spline stiffness
      madd(85)=pos2(ra,rd,rb,rd+45)   
      mtext(85)="Spline Stiffness" ; rd=rd+120 
      madd(28)=pos2(ra+120,rd,rc,rd+45) ; rd=rd+60  
      mtext(28)="Tree Display Options" 
      b(77)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Tree raw measures
      madd(77)=pos2(ra,rd,rb,rd+45)   
      mtext(77)="Tree Measurments" ; rd=rd+60 
      b(78)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Tree detrending Curve
      madd(78)=pos2(ra,rd,rb,rd+45)   
      mtext(78)="Detrending Curve" ; rd=rd+60 
      b(79)=posn(rb,rd,rc,rd+45,"On",TR,TR)  ! Tree indices
      madd(79)=pos2(ra,rd,rb,rd+45)   
      mtext(79)="Tree Indices" ; rd=rd+60 
      b(80)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Tree chronology
      madd(80)=pos2(ra,rd,rb,rd+45)   
      mtext(80)="Tree chronology" ; rd=rd+60  
      b(86)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Tree first display year
      madd(86)=pos2(ra,rd,rb,rd+45)   
      mtext(86)="Tree First Year" ; rd=rd+60 
      b(87)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Tree last display year
      madd(87)=pos2(ra,rd,rb,rd+45)   
      mtext(87)="Tree Last Year" ; rd=rd+60 
      rb=710 ; rc=810 ; rd=200
      madd(29)=pos2(rb,rd,rc,rd+45) ; rd=rd+60
      mtext(29)="Apply" 
      b(88) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Mean RCS curve
      b(89) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Smooth RCS curve
      b(90) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! RCS standard deviation
      b(91) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! RCS standard error
      b(92) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+300 ! RCS counts
      b(93) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Chronology
      b(94) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! CRN standard deviation
      b(95) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! CRN counts
      b(96) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+360 ! CRN smoothed
      b(97) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Tree raw measures
      b(98) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Tree detrending Curve
      b(99) =posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Tree indices
      b(100)=posn(rb,rd,rc,rd+45,"",TR,TR) ; rd=rd+60  ! Tree chronology
      ra=1100 ; rb=1300 ; rc=1500 ; re=1700 ; rd=200 
      madd(102)=pos2(rb,rd,re,rd+45)   
      mtext(102)="Individual RCS Curves" ; rd=rd+60 
      madd(106)=pos2(ra,rd,rc,rd+45)   
      mtext(106)="Smooth    Mean"
      madd(110)=pos2(rc,rd,re+120,rd+45)   
      mtext(110)="S.Dev      S.Err" ; rd=rd+60 
      DO i=101,141,4
        b(i  )=posn(ra,rd,ra+100,rd+45,"",TR,TR)   ! Smooth CRN
        b(i+1)=posn(rb,rd,rb+100,rd+45,"",TR,TR)   ! Mean CRN
        b(i+2)=posn(rc,rd,rc+100,rd+45,"",TR,TR)   ! Standard deviation
        b(i+3)=posn(re,rd,re+100,rd+45,"",TR,TR)   ! Standard error
        madd(i)=pos2(re+160,rd,re+320,rd+45)   
        rd=rd+60
      ENDDO
      rd=980
      madd(103)=pos2(rb,rd,re,rd+45)   
      mtext(103)="Individual Chronologies" ; rd=rd+60 
      madd(107)=pos2(ra,rd,rc,rd+45)   
      mtext(107)="Smooth    Mean"
      madd(111)=pos2(rc,rd,re+120,rd+45)   
      mtext(111)="S.Dev      S.Err" ; rd=rd+60 
      DO i=145,185,4
        b(i  )=posn(ra,rd,ra+100,rd+45,"",TR,TR)   ! Smooth RCS
        b(i+1)=posn(rb,rd,rb+100,rd+45,"",TR,TR)   ! Mean RCS
        b(i+2)=posn(rc,rd,rc+100,rd+45,"",TR,TR)   ! Standard deviation
        b(i+3)=posn(re,rd,re+100,rd+45,"",TR,TR)   ! Standard error
        madd(i)=pos2(re+160,rd,re+320,rd+45)   
        rd=rd+60
      ENDDO
      mtext(101)="1st RCS"  ; mtext(105)="2nd RCS"  
      mtext(109)="3rd RCS"  ; mtext(113)="4th RCS"  
      mtext(117)="5th RCS"  ; mtext(121)="6th RCS"  
      mtext(125)="7th RCS"  ; mtext(129)="8th RCS"  
      mtext(133)="9th RCS"  ; mtext(137)="10th RCS"  
      mtext(141)="11th RCS" ; mtext(145)="1st CRN"  
      mtext(149)="2nd CRN"  ; mtext(153)="3rd CRN"  
      mtext(157)="4th CRN"  ; mtext(161)="5th CRN"  
      mtext(165)="6th CRN"  ; mtext(169)="7th CRN"  
      mtext(173)="8th CRN"  ; mtext(177)="9th CRN"  
      mtext(181)="10th CRN" ; mtext(185)="11th CRN"  
      ra=2280 ; rb=2620 ; rc=2900 ; rd=1200
      madd(190)=pos2(ra,rd,rc,rd+45)  ! on.off colour option 
      mtext(190)="Turn on/off or change line style" ; rd=rd+60 
      b(190)=posn(ra+60,rd,rb+100,rd+45,"",TR,TR)  ; rd=rd+120
      madd(30)=pos2(ra+60,rd,rc,rd+45)   
      mtext(30)="Line description to apply" ; rd=rd+60 
      b(191)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Line thickness
      madd(191)=pos2(ra,rd,rb,rd+45)   
      mtext(191)="Line Thickness" ; rd=rd+60 
      b(192)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Line Style
      madd(192)=pos2(ra,rd,rb,rd+45)   
      mtext(192)="Line Style" ; rd=rd+60 
      b(193)=posn(rb,rd,rc,rd+45,"",TR,TR)    ! Line Colour
      madd(193)=pos2(ra,rd,rb,rd+45)   
      mtext(193)="Line Colour" ; rd=rd+60 
      b(194:198)%ok=FA   ! b(194:198) spare      
! Save Menu     
      b(199)=posn(2600,30,2740,75,"Exit",TR,TR) ! Quit option
      madd(199)=pos2(2250,30,2600,75) 
      mtext(199)="Leave Save Menu"
      ra=880 ; rb=1360 ; rc=1840 ; rd=530
      b(200)=posn(rb-80,rd,rb+240,rd+45,"",TR,TR) ! Output CRN format
      madd(200)=pos2(ra,rd,ra+360,rd+45) 
      mtext(200)="Chronology Format"
      b(201)=posn(rc+240,rd,rc+400,rd+45,"Append",TR,TR) ! Append option
      madd(201)=pos2(rc-160,rd,rc+200,rd+45)
      mtext(201)="New File or Append" 
      b(202)=posn(ra,rd+270,ra+400,rd+315,"work",TR,TR) ! CRN Directory
      madd(202)=pos2(ra+80,rd+210,ra+480,rd+255) 
      mtext(202)="Directory"
      b(203)=posn(rb,rd+270,rb+400,rd+315,"crns",TR,TR) ! CRN File Name
      madd(203)=pos2(rb+80,rd+210,rb+480,rd+255) 
      mtext(203)="File Name"
      b(204)=posn(rc+400,rd+270,rc+600,rd+315,"Save CRN",TR,TR) ! Save CRN
      madd(204)=pos2(rc+340,rd+210,rc+660,rd+255)
      mtext(204)="Save Chronology"
      madd(11)=pos2(rc,rd+210,rc+240,rd+255) 
      mtext(11)="Extension"
      madd(12)=pos2(rc+40,rd+270,rc+240,rd+315) 
      mtext(12)=".crn"
      rd=1310
      b(205)=posn(rb-80,rd-60,rb+200,rd-15,"",TR,TR)  ! Output RAW format
      madd(205)=pos2(ra,rd-60,ra+360,rd-15) 
      mtext(205)="RAW Data Format"
      b(206)=posn(rc+400,rd+120,rc+600,rd+165,"Save RAW",TR,TR) ! Save RAW
      madd(206)=pos2(rc+340,rd+60,rc+660,rd+105)
      mtext(206)="Save RAW Data"
      b(207)=posn(ra,rd+120,ra+400,rd+165,"work",TR,TR) ! RAW Directory
      madd(207)=pos2(ra+80,rd+60,ra+480,rd+105) 
      mtext(207)="Directory"
      b(208)=posn(rb,rd+120,rb+400,rd+165,"rawdata",TR,TR) ! RAW File Name
      madd(208)=pos2(rb+80,rd+60,rb+480,rd-105) 
      mtext(208)="File Name"
      b(209)=posn(rc+680,rd-60,rc+920,rd-15,"No Sort",TR,TR) ! Sort option
      madd(209)=pos2(rc+400,rd-60,rc+680,rd-15)   
      mtext(209)="Sort Trees by"
      b(210)%ok=FA   ! b(210) spare      
      madd(13)=pos2(rc,rd-60,rc+60,rd-15) 
      mtext(13)="Extension"
      madd(14)=pos2(rc+40,rd,rc+240,rd+45) 
      mtext(14)=datext(1)
      madd(15)=pos2(rc+40,rd+60,rc+240,rd+105) 
      mtext(15)=datext(2)
      madd(16)=pos2(rc+40,rd+120,rc+240,rd+165) 
      mtext(16)=datext(3)
      madd(17)=pos2(rc+40,rd+180,rc+240,rd+225) 
      mtext(17)=datext(4)
      madd(18)=pos2(rc+40,rd+240,rc+240,rd+285) 
      mtext(18)=datext(5)
      madd(19)=pos2(rc+40,rd+300,rc+240,rd+345) 
      mtext(19)=datext(6)
      ra=160 ; rb=600 ; rc=720 ; rd=240
      madd(10)=pos2(ra+40,rd,rc,rd+45)  
      mtext(10)="Series to Save"   ; rd=rd+60 
      b(211)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Mean RCS curve
      madd(211)=pos2(ra,rd,rb,rd+45)   
      mtext(211)="Mean RCS curve" ; rd=rd+60 
      b(212)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Smooth RCS curve
      madd(212)=pos2(ra,rd,rb,rd+45)   
      mtext(212)="Smooth RCS curve" ; rd=rd+60 
      b(213)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! RCS standard deviation
      madd(213)=pos2(ra,rd,rb,rd+45)   
      mtext(213)="RCS Stand. Dev." ; rd=rd+60 
      b(214)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! RCS standard error
      madd(214)=pos2(ra,rd,rb,rd+45)   
      mtext(214)="RCS Stan. Error" ; rd=rd+60 
      b(215)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Chronology
      madd(215)=pos2(ra,rd,rb,rd+45)   
      mtext(215)="Chronology" ; rd=rd+60 
      b(216)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! CRN standard deviation
      madd(216)=pos2(ra,rd,rb,rd+45)   
      mtext(216)="CRN Stand. Dev." ; rd=rd+60 
      b(217)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! CRN smooth
      madd(217)=pos2(ra,rd,rb,rd+45)   
      mtext(217)="Smooth Chronology" ; rd=rd+60 
      b(218)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Multiple series
      madd(218)=pos2(ra,rd,rb,rd+45)   
      mtext(218)="Save Multiple Series" ; rd=rd+180
      b(219)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Tree raw measures
      madd(219)=pos2(ra,rd,rb,rd+45)   
      mtext(219)="Tree Measurments" ; rd=rd+60 
      b(220)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Transformed measurements
      madd(220)=pos2(ra,rd,rb,rd+45)   
      mtext(220)="Tree Transformed" ; rd=rd+60 
      b(221)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Sig-Free measurements
      madd(221)=pos2(ra,rd,rb,rd+45)   
      mtext(221)="Sig-Free Measures" ; rd=rd+60 
      b(222)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Tree detrending Curve
      madd(222)=pos2(ra,rd,rb,rd+45)   
      mtext(222)="Detrending Curve" ; rd=rd+60 
      b(223)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Tree indices
      madd(223)=pos2(ra,rd,rb,rd+45)   
      mtext(223)="Tree Indices" ; rd=rd+60 
      b(224)=posn(rb,rd,rc,rd+45,"Off",TR,FA)  ! Tree AR Model indices
      madd(224)=pos2(ra,rd,rb,rd+45)   
      mtext(224)="Tree AR Indices" ; rd=rd+60 
      b(225:234)%ok=FA   ! spare      
! Select Menu     
      b(235)=posn(2360,160,2600,205,"No Sort",TR,TR) ! Sort option
      madd(235)=pos2(2080,160,2360,205)   
      mtext(235)="Sort Trees by"
      madd(236)=pos2(1800,30,2320,75)   ! Chronology name
      mtext(236)="Chronology"
      madd(239)=pos2(1400,30,1800,75)
      mtext(239)="Current Chronology"
      b(236)=posn(960,70,1160,115,"Previous",TR,TR) ! Previous tree name
      b(237)=posn(960,130,1080,175,"Next",TR,TR) ! Next tree name
      b(238)=posn(540,90,880,135,"",TR,TR)     ! Core name
      madd(237)=pos2(80,90,540,135)
      mtext(237)="Current Tree Name"
      b(239)=posn(2840,30,2960,75,"Exit",TR,TR) ! Quit option
      madd(238)=pos2(2480,30,2840,75)
      mtext(238)="Leave Tree Menu"
      ra=2740 ; rb=2900 
      b(240)=posn(1640,100,1960,145,"Save SRCS",TR,TR) ! Save RCS file
      b(241)=posn(1640,160,1960,205,"Load SRCS",TR,TR) ! Load RCS file
      b(242)=posn(2360,100,2600,145,"Calc CRN",TR,TR)  ! Recalculate
      b(243)=posn(ra,390,rb,435,"DET  1",TR,FA)  ! RCS Curves
      b(244)=posn(ra,450,rb,495,"DET  2",TR,FA)  ! RCS Curves
      b(245)=posn(ra,510,rb,555,"DET  3",TR,FA)  ! RCS Curves
      b(246)=posn(ra,570,rb,615,"DET  4",TR,FA)  ! RCS Curves
      b(247)=posn(ra,630,rb,675,"DET  5",TR,FA)  ! RCS Curves
      b(248)=posn(ra,690,rb,735,"DET  6",TR,FA)  ! RCS Curves
      b(249)=posn(ra,750,rb,795,"DET  7",TR,FA)  ! RCS Curves
      b(250)=posn(ra,810,rb,855,"DET  8",TR,FA)  ! RCS Curves
      b(251)=posn(ra,870,rb,915,"DET  9",TR,FA)  ! RCS Curves
      b(252)=posn(ra,930,rb,975,"DET 10",TR,FA)  ! RCS Curves
      b(253)=posn(ra,990,rb,1035,"DET 11",TR,FA)  ! RCS Curves
      b(254)=posn(ra,240,rb,285,"as DET",TR,TR)  ! RCS DET/not
      madd(249)=pos2(ra-200,240,ra,285)
      mtext(249)="RCS curve"
      ra=2600 ; rb=2940
      madd(250)=pos2(ra,1280,rb,1325)
      mtext(250)="Selected Options"
      madd(243)=pos2(ra,1340,rb,1385)
      madd(244)=pos2(ra,1385,rb,1430)
      madd(245)=pos2(ra,1430,rb,1475)
      madd(246)=pos2(ra,1475,rb,1520)
      madd(247)=pos2(ra,1520,rb,1565)
      madd(248)=pos2(ra,1565,rb,1610)
      madd(257)=pos2(cnx,btop-45,cnx+120,btop)
      mtext(257)="R = RCS curve, D = Detrend Curve"
      b(255)=posn(cnx,btop,gright,bbot,"",TR,FA)  ! Choice table area
      b(256)=posn(2080,190,2200,235,"Exit",TR,TR) ! Exit option
      madd(255)=pos2(1600,190,2080,235)
      mtext(255)="Leave RCS select"
      madd(254)=pos2(ra+50,330,rb,375)
      mtext(254)="Selected DET"
      b(257)=posn(1000,btop-60,1320,btop-15,"Prev Page",FA,TR) 
      b(258)=posn(1400,btop-60,1720,btop-15,"Next Page",FA,TR)   
      b(259:270)%ok=FA   ! spare      
! Name menu
      b(264)=posn(cnx,btop,gright,bbot,"",TR,FA)      ! Choice table area
      madd(4)=pos2(cnx,btop-50,cnx+480,btop)
      mtext(4)="Choose a File/Tree"
      b(265)=posn(2080,30,2200,75,"Exit",TR,TR) ! Exit option
      madd(5)=pos2(1700,30,2080,75)
      mtext(5)="Leave file select"
      b(266)=posn(1000,btop-80,1320,btop-15,"Prev Page",FA,TR) 
      b(267)=posn(1400,btop-80,1720,btop-15,"Next Page",FA,TR)   
      b(268)=posn(150,30,400,75," ",TR,FA)     ! Status button
      madd(268)=pos2(10,30,150,75)
      mtext(268)="Status"
      b(269:270)%ok=FA   ! b(49:50) spare      
! Figures Menu
      b(271)=posn(cnx,btop,gright,bbot,"",TR,TR)    ! Choice table area
      b(272)=posn(2700,30,2900,75,"Continue",TR,TR) ! Next action
      madd(271)=pos2(cnx,btop-45,cnx+480,btop)
      mtext(271)="Choose a Plot/Action"
      b(273)=posn(1400,30,2350,75," ?",TR,TR)     ! Chronology name
      madd(273)=pos2(1000,30,1400,75)
      mtext(273)="Current chronology"
      ! b(274+) spare      
! Start Menu
      ra=600 ; rb=1000 ; rc=90 ; rd=135
      b(451)=posn(2600,160,2800,205,"Journal Off",TR,TR)  ! Journal on/off
      b(452)=posn(ra,rc+30,rb,rd+30," ?",TR,TR)       ! List name
      madd(303)=pos2(ra-440,rc+30,ra,rd+30)
      mtext(303)="File List Filename"
      madd(306)=pos2(rb,rc+30,rb+440,rd+30)
      mtext(306)=" - does not exist?"
      b(453)=posn(ra,rc+90,rb+480,rd+90," ?",TR,TR) ! chronology name
      madd(304)=pos2(ra-440,rc+90,ra,rd+90)
      mtext(304)="Current Data File"
      madd(305)=pos2(rb+480,rc+90,rb+920,rd+90)
      mtext(305)=" - does not exist?"
      b(454)=posn(2600,30,2800,75,"Quit",TR,TR)   ! Quit option
      madd(307)=pos2(2120,30,2600,75)
      mtext(307)="Leave the Program"
      ra=2600 ; rb=2800 ; rc=210 ; rd=255
      madd(319)=pos2(ra-320,rc+30,ra,rd+30)
      mtext(319)="Choose a Plot/Action"
      b(455)=posn(ra,rc+120,rb,rd+120,"Detrend",TR,TR)  ! Curve detrend
      madd(308)=pos2(ra-420,rc+120,ra,rd+120)
      mtext(308)="Curve-fit Detrend"
      rc=rc+180 ; rd=rd+180
      b(456)=posn(ra,rc,rb,rd,"Detrend",TR,TR) ! RCS Detrend
      madd(309)=pos2(ra-420,rc,ra,rd)
      mtext(309)="RCS Detrend"
      rc=rc+60 ; rd=rd+60
      b(457)=posn(ra,rc,rb,rd,"",FA,FA)  ! Was chrons
      b(458)=posn(ra,rc,rb,rd,"Figures",TR,TR) ! Figures
      madd(311)=pos2(ra-420,rc,ra,rd)
      mtext(311)="Reports/Figures menu"
      rc=rc+60 ; rd=rd+60
! Plot and Alps grid box data
      b(468)=posn(640,775,1080,820,"Write Grid data",TR,TR)
      b(469)=posn(640,890,1000,935," ",TR,TR)
      madd(33)=pos2(240,890,640,935) ; mtext(33)="Latitude"
      b(470)=posn(640,1025,1000,1070," ",TR,TR) 
      madd(34)=pos2(240,1025,640,1070) ; mtext(34)="Longitude"
      b(471)=posn(200,180,460,225,"Temperature",TR,TR)  
      b(472)=posn(200,240,460,285,"Precipitation",TR,FA)  
      b(473)=posn(200,300,460,345,"Cloudiness",TR,FA)  
      b(474)=posn(2700,30,2900,75,"Quit",TR,TR) ! Quit option
! Odds and Ends
      bhigh=yellow   ! Default button highlight colour
      bred=FA        ! Red highlight off
      cur=8 ; dirc="work" ; dirr="work"
      namc="crns" ; namr="rawdata" 
      cdsp=30 ; idtsl=100 ; idtsn=-66
      RETURN
      END SUBROUTINE crust_setup
!--------------------------------------------------------
      SUBROUTINE save_menu()            ! Save files 
      IMPLICIT NONE                 
      INTEGER :: j,k
      b(209)%lab=tstlab(tst)
      b(211:214)%ok=idt.EQ.-2           ! Only if RCS
      b(218)%ok=idt.EQ.-2.AND.src.GT.1  ! Only multi RCS
      IF (.NOT.b(218)%ok) b(218)%on=FA
      b(221)%ok=sfo.EQ.2                ! Only if Sig-free
      b(224)%ok=jrb.GT.1                ! Only if AR modelling
      ADO: DO j=1,10000                 ! Menu attempts
        CALL erase()                    ! Clear screen 
        DO k=199,224                   
          IF (b(k)%ok) THEN
            CALL but_draw(k,"") ; CALL mwrite(k)
          ENDIF
        ENDDO  
        DO k=10,13 ; CALL mwrite(k) ; ENDDO
        DO k=14,19 ; IF (b(k+205)%ok) CALL mwrite(k) ; ENDDO
        CALL mouse_click(4,199,224)     ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                        ! Mouse missed "on" buttons
        CASE (199) ; EXIT ADO           ! Leave display menu
        CASE (200)                      ! Output CRN format
          IF (cft.EQ.2) THEN ; cft=1 ; ELSE ; cft=cft+1 ; ENDIF
          b(mous)%lab=cftlab(cft) ; mtext(12)=cftext(cft)  ! File extention
        CASE (201)
          IF (append) THEN              ! Append or not
            append=FA ; b(mous)%lab="New File"
          ELSE
            append=TR ; b(mous)%lab="Append"
          ENDIF
        CASE (202) ; CALL get_str(mous) ! Current CRN directory
          IF (LEN_TRIM(b(mous)%lab).GT.2) THEN
            dirc=b(mous)%lab(1:20)
          ELSE                          ! Default directory
            dirc="work" ; b(mous)%lab=dirc 
          ENDIF
        CASE (203) ; CALL get_str(mous) ! Current CRN filename
          IF (LEN_TRIM(b(mous)%lab).GT.2) THEN
            namc=b(mous)%lab(1:20)
          ELSE                          ! Default filename
            namc="crns" ; b(mous)%lab=namc 
          ENDIF
        CASE (204) ; CALL save_crn() ! Write CRN file
        CASE (205)                      ! Output Raw format
          IF (rft.EQ.3) THEN ; rft=1 ; ELSE ; rft=rft+1 ; ENDIF
          b(mous)%lab=rftlab(rft)
        CASE (206) ; CALL save_data() ! Write RAW file
        CASE (207) ; CALL get_str(mous) ! Current RAW directoty
          IF (LEN_TRIM(b(mous)%lab).GT.2) THEN
            dirr=b(mous)%lab(1:20)
          ELSE                          ! Default directory
            dirr="work" ; b(mous)%lab=dirr  
          ENDIF
        CASE (208) ; CALL get_str(mous) ! Current RAW filename
          IF (LEN_TRIM(b(mous)%lab).GT.2) THEN
            namr=b(mous)%lab(1:20)
          ELSE                          ! Default filename 
            namr="rawdata" ; b(mous)%lab=namr
          ENDIF
        CASE (209)                      ! Tree sort option 
          IF (tst.EQ.7) THEN ; tst=1 ; ELSE ; tst=tst+1 ; ENDIF
          b(mous)%lab=tstlab(tst) ; CALL tree_sort()
        CASE (210)                      ! Not used
        CASE (211:224)                  ! Data selection
          IF (b(mous)%on) THEN
            b(mous)%on=FA ; b(mous)%lab="Off"
          ELSE
            b(mous)%on=TR ; b(mous)%lab="On"
          ENDIF
        END SELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE save_menu          
!------------------------------------------------------------------------
      SUBROUTINE disp_menu()               ! Display Options
      IMPLICIT NONE                 
      INTEGER :: i,j,m,r
      b(66)%ok=idt.EQ.-2
      b(190)%lab="Turn Buttons on/off" ; oocol=TR
      lth=1 ; lst=1 ; lco=1
      ADO: DO j=1,10000                    ! Menu attempts
        CALL erase()                       ! Clear screen 
        CALL but_draw(66,"")               ! Save parameters
        DO i=67,87   ; CALL but_draw(i,"") ; CALL mwrite(i) ; ENDDO
        DO i=88,100  ; CALL but_draw3(i) ; ENDDO
        IF (src.GT.1) THEN                 ! Optional Multi RCS
          CALL mwrite(102) ; CALL mwrite(103)
          CALL mwrite(106) ; CALL mwrite(107)
          CALL mwrite(110) ; CALL mwrite(111)
          DO i=101,100+srcno*4
            CALL but_draw3(i) ; CALL but_draw3(i+44)
            IF (MOD(i,4).EQ.1) THEN
              CALL mwrite(i) ; CALL mwrite(i+44)
            ENDIF
          ENDDO
        ENDIF
        CALL but_draw(190,"") ; CALL mwrite(190)
        IF (.NOT.oocol) THEN
          DO i=191,193 ; CALL but_draw2(i) ; CALL mwrite(i) ; ENDDO
          CALL mwrite(30)
        ENDIF
        DO i=26,29 ; CALL mwrite(i) ; ENDDO
        CALL mouse_click(3,66,193)         ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                           ! Mouse missed "on" buttons
        CASE (66) ; CALL write_default()   ! Save parameter values
        CASE (67) ; EXIT ADO               ! Leave display menu
        CASE (68:80) ; m=mous-66 ; lins(m)%on=.NOT.lins(m)%on
          IF(lins(m)%on)THEN ; b(mous)%lab="On"
          ELSE ; b(mous)%lab="Off" ; ENDIF
        CASE (81) ; CALL get_int(mous,i)   ! RCS first year
          IF (i.GE.0.AND.i.LE.RDly-10) RDfy=i
          WRITE(b(mous)%lab,'(I4)') RDfy
        CASE (82) ; CALL get_int(mous,i)   ! RCS last year
          IF (i.GE.RDfy+10.AND.i.LE.sly(1)) RDly=i
          WRITE(b(mous)%lab,'(I4)') RDly
        CASE (83) ; CALL get_int(mous,i)   ! CRN first year
          IF (i.GE.xfy.AND.i.LE.CDly-10) CDfy=i
          WRITE(b(mous)%lab,'(I5)') CDfy
        CASE (84) ; CALL get_int(mous,i)   ! CRN last year
          IF (i.GE.CDfy+10.AND.i.LE.xly) CDly=i
          WRITE(b(mous)%lab,'(I4)') CDly
        CASE (85) ; CALL get_int(mous,i)   ! CRN spline stiffness
          IF (i.GE.5.AND.i.LE.9999) CDsp=i
          r=xyr ; WRITE(b(mous)%lab,'(I4)') CDsp
          CALL splinet(r,xcrn(1:r,mx),CDsp,xcsm(1:r,mx))
        CASE (86) ; CALL get_int(mous,i)   ! Tree first year
          IF (i.GE.fy(tre(cc)).AND.i.LE.TDly-10) TDfy=i
          WRITE(b(mous)%lab,'(I5)') TDfy
        CASE (87) ; CALL get_int(mous,i)   ! Tree last year
          IF (i.GE.TDfy+10.AND.i.LE.ly(tre(cc))) TDly=i
          WRITE(b(mous)%lab,'(I5)') TDly
        CASE (88:188) ; m=mous-86 
          IF (oocol) THEN
            IF(lins(m)%on) THEN
              lins(m)%on=FA ; b(mous)%lab="Off"
            ELSE
              lins(m)%on=TR ; b(mous)%lab="On"
            ENDIF
          ELSE
            lins(m)%th=lth ; lins(m)%st=lst ; lins(m)%co=lcol(lco) 
          ENDIF
        CASE (190)                         ! On/off or colour option
          IF (oocol) THEN
            b(mous)%lab="Change Line Style"   ; oocol=FA
          ELSE
            b(mous)%lab="Turn Buttons on/off" ; oocol=TR
          ENDIF
        CASE (191)                          ! Line thickness
          IF (lth.EQ.8)  THEN ; lth=1 ; ELSE ; lth=lth+1 ; ENDIF
        CASE (192)                          ! Line style
          IF (lst.EQ.6)  THEN ; lst=0 ; ELSE ; lst=lst+1 ; ENDIF
        CASE (193)                          ! Line colour
          IF (lco.EQ.15) THEN ; lco=0 ; ELSE ; lco=lco+1 ; ENDIF
        END SELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE disp_menu          
!------------------------------------------------------------------------
      SUBROUTINE name_menu(ref1)  ! Get new data/tree
      IMPLICIT NONE                 
      INTEGER,INTENT(IN) :: ref1
      INTEGER :: j,k,p,q
      b(266:267)%ok=(ref1.EQ.2.AND.nc.GT.rw*ncol)
      IF (ref1.EQ.2) sccc=1+ &    ! Current tree in view
        MIN(MAX(0,cc/rw-2),MAX(nc/rw+1-ncol,0))*rw
      ADO: DO j=1,10000           ! Menu attempts
        CALL erase()              ! Clear screen 
        IF (ref1.EQ.1) THEN
          chc(1)=cf ; chs=1 ; sccc=1
          CALL ch_disp(ccol,nfil,cnam(1:nfil),chs,chc,1)
        ELSE
          chc(1)=cc ; chs=1 
          CALL ch_disp(ncol,nc,nam(1:nc),chs,chc,2)
        ENDIF
        CALL mwrite(4) ; CALL mwrite(5)
        CALL but_draw(265,"")      ! Exit button
        IF (sccc.GT.1) CALL but_draw(266,"")          ! Previous page
        IF (sccc.LT.nc-rw*ncol) CALL but_draw(267,"") ! Next page
        CALL mouse_click(6,264,267) ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                  ! Mouse missed "on" buttons
        CASE (264)                 
          IF (ref1.EQ.1) THEN     ! Choosing a new name
            k=rw*INT(DBLE(ccol)*DBLE(msx-cnx)/DBLE(gright-cnx))+ &
               INT(DBLE(rw*(b(mous)%y1-msy))/DBLE(b(mous)%y1-b(mous)%y2))+1
            IF (k.GE.1.AND.k.LE.nfil) THEN  ! Valid CRN
              cf=k ; q=LEN_TRIM(cnam(cf)) ; p=MAX(1,q-19)
              b(24)%lab=cnam(cf)(p:q)
            ENDIF
          ELSEIF (ref1.EQ.2) THEN ! Choosing a new tree
            k=rw*INT(DBLE(ncol)*DBLE(msx-cnx)/DBLE(gright-cnx))+ &
              INT(DBLE(rw*(b(mous)%y1-msy))/DBLE(b(mous)%y1-b(mous)%y2))+1
            k=k+sccc-1
            IF (k.GE.1.AND.k.LE.nc) cc=k 
          ENDIF
        CASE (265) ; EXIT ADO      ! Leave selection process
        CASE (266) ; sccc=MAX(1,sccc-rw*(ncol-1)) ! Previous page
        CASE (267) ; sccc=MIN(sccc+rw*(ncol-1),nc-rw*ncol+1) ! Next page
        END SELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE name_menu          
!------------------------------------------------------------------------
      SUBROUTINE tree_menu()    ! Display individual trees
      IMPLICIT NONE                 
      INTEGER  :: j,k
      cc=1 ; sccc=1 ; mtext(56)=cnam(cf) ; CALL tree_sort()
      ADO: DO j=1,10000           ! Menu attempts
        b(56)%ok=cc.GT.1 ; b(57)%ok=cc.LT.nc
        CALL erase()              ! Clear screen 
        CALL plot_core()    
        CALL mwrite(9) ; CALL mwrite(56)  ! Chronology name
        DO k=55,61                ! Draw button boxes
          IF (b(k)%ok) THEN
            SELECT CASE (k)
              CASE (55) ; CALL but_draw(k,tstlab(tst)) ; CALL mwrite(55)
              CASE (56:57) ; CALL but_draw(k,"") 
              CASE (58) ; CALL but_draw(k,nam(tre(cc))) ; CALL mwrite(7)
              CASE (59) ; CALL but_draw(k,"") ; CALL mwrite(8)
              CASE (60:61) ; CALL but_draw(k,"") 
           ENDSELECT
          ENDIF
        ENDDO
        CALL mouse_click(2,55,61) ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                  ! Mouse missed "on" buttons
        CASE (55)                 ! Tree sort option 
          IF (tst.EQ.7) THEN ; tst=1 ; ELSE ; tst=tst+1 ; ENDIF
          b(mous)%lab=tstlab(tst) ; CALL tree_sort()
        CASE (56) ; IF (cc.GT.1) cc=cc-1  ! Previous Tree
        CASE (57) ; IF (cc.LT.nc) cc=cc+1 ! Next Tree
        CASE (58) ; CALL name_menu(2)     ! Wants select a tree
        CASE (59) ; EXIT ADO              ! Wants to exit menu
        CASE (60) ; CALL plot_ps(2)       ! Save a .ps file
        CASE (61) ; CALL disp_menu()      ! Display menu
        END SELECT
        TDfy=fy(tre(cc)) ; WRITE(b(86)%lab,'(I5)') TDfy
        TDly=ly(tre(cc)) ; WRITE(b(87)%lab,'(I5)') TDly
      ENDDO ADO
      RETURN 
      END SUBROUTINE tree_menu          
!------------------------------------------------------------------------
      SUBROUTINE name_menu3()     ! Allocate RCS curve
      IMPLICIT NONE                 
      INTEGER :: i,j,k,m,dd
      b(257:258)%ok=(nc.GT.rw*ncol)
      sccc=1+MIN(MAX(0,cc/rw-2),MAX(nc/rw+1-ncol,0))*rw
      dd=trm(tre(cc))               ! Selected curve
      ADO: DO j=1,10000             ! Menu attempts
        b(243:253)%on=FA ; b(dd+242)%on=TR
        CALL erase() ; chc(1)=cc ; chs=1 
        CALL ch_disp(ccol,nc,nam(1:nc),chs,chc,3)
        DO i=243,250 ; CALL mwrite(i) ; ENDDO
        DO i=254,257 ; CALL mwrite(i) ; ENDDO
        DO i=243,253 ; CALL but_draw(i,"") ; ENDDO
        CALL but_draw(254,trrlab(trrok)) 
        CALL but_draw(256,"")       ! Exit button
        IF (sccc.GT.1) CALL but_draw(257,"")          ! Previous page
        IF (sccc.LT.nc-rw*ncol) CALL but_draw(258,"") ! Next page
        CALL mouse_click(5,243,258) ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                    ! Mouse missed "on" buttons
        CASE (243:253) ; dd=mous-242  ! Pick DET curve
        CASE (254)                  ! Pick RCS as DET/not
          IF (trrok.EQ.1) THEN ; trrok=2 ; ELSE ; trrok=1 ; ENDIF
        CASE (255)                  ! Allocating RCS curve
          k=rw*INT(DBLE(ccol)*DBLE(msx-cnx)/DBLE(gright-cnx))+ &
            INT(DBLE(rw*(b(mous)%y1-msy))/DBLE(b(mous)%y1-b(mous)%y2))+1
          k=k+sccc-1
          IF (k.GE.1.AND.k.LE.nc) THEN
            cc=k ; m=tre(cc) ; trm(m)=dd
            IF (trrok.EQ.2) THEN ; trr(m)=dd ; ELSE ; trr(m)=0 ; ENDIF
          ENDIF
        CASE (256) ; EXIT ADO       ! Leave selection process
        CASE (257) ; sccc=MAX(1,sccc-rw*(ncol-1)) ! Previous page
        CASE (258) ; sccc=MIN(sccc+rw*(ncol-1),nc-rw*ncol+1) ! Next page
        ENDSELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE name_menu3          
!------------------------------------------------------------------------
      SUBROUTINE srcs_menu()  ! Select trees for RCS curves
      IMPLICIT NONE                 
      INTEGER  :: j,k,m
      IF (src.NE.3) trr(1:nc)=trm(1:nc)
      cc=1 ; sccc=1 ; mtext(236)=cnam(cf) ! Open or set current as default
      DO j=1,mx-1 ; b(242+j)%ok=j.LE.srcno ; ENDDO ! RCS curves
      mtext(243)=itnlab(itn) ; mtext(244)=rdtlab(rdt)
      mtext(245)=indlab(ind) ; mtext(246)=sfolab(sfo)
      mtext(247)=poolab(poo) ; mtext(248)=trclab(trc)
      ADO: DO j=1,10000             ! Menu attempts
        m=tre(cc)
        IF (trr(m).EQ.0) THEN ; trrok=1 ; ELSE ; trrok=2 ; ENDIF  
        b(243:253)%on=FA ; b(242+trm(m))%on=TR
        b(236)%ok=cc.GT.1 ; b(237)%ok=cc.LT.nc
        CALL erase() ; CALL plot_core()    
        CALL mwrite(239) ; CALL mwrite(236) ; CALL mwrite(254) 
        DO k=243,250 ; CALL mwrite(k) ; ENDDO
        DO k=235,254                ! Draw button boxes
          IF (b(k)%ok) THEN
            SELECT CASE (k)
            CASE (235)
              CALL but_draw(k,tstlab(tst)) ; CALL mwrite(235)
            CASE (236:237) ; CALL but_draw(k,"") 
            CASE (238)
              CALL but_draw(k,nam(m)) ; CALL mwrite(237)
            CASE (239) ; CALL but_draw(k,"") ; CALL mwrite(238)
            CASE (240:253) ; CALL but_draw(k,"") 
            CASE (254) ; CALL but_draw(k,trrlab(trrok)) 
           ENDSELECT
          ENDIF
        ENDDO
        CALL mouse_click(7,235,254) ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                    ! Mouse missed "on" buttons
        CASE (235)                  ! Tree sort option 
          IF (tst.EQ.7) THEN ; tst=1 ; ELSE ; tst=tst+1 ; ENDIF
          b(mous)%lab=tstlab(tst) ; CALL tree_sort()
        CASE (236) ; IF (cc.GT.1) cc=cc-1  ! Previous Tree
          m=tre(cc) 
        CASE (237) ; IF (cc.LT.nc) cc=cc+1 ! Next Tree
          m=tre(cc)
        CASE (238) ; CALL name_menu3()     ! Wants select a tree
        CASE (239) ; EXIT ADO              ! Wants to exit menu
        CASE (240) ; CALL write_srcs()     ! Save SRCS file
        CASE (241) ; CALL read_srcs()      ! Load SRCS file
        CASE (242) ; CALL new_chron(FA)    ! Calculate chronology
        CASE (243:253) ; trm(m)=mous-242   ! Pick DET curve
          IF (trrok.EQ.2) THEN ; trr(m)=mous-242
          ELSE ; trr(m)=0 ; ENDIF
        CASE (254)                         ! Pick RCS as DET/not
          IF (trrok.EQ.2) THEN ; trrok=1 ; trr(m)=0
          ELSE ; trrok=2  ; trr(m)=trm(m) ; ENDIF
        ENDSELECT
        TDfy=fy(m) ; WRITE(b(86)%lab,'(I5)') TDfy
        TDly=ly(m) ; WRITE(b(87)%lab,'(I5)') TDly
      ENDDO ADO
      srcok=FA
      RETURN 
      END SUBROUTINE srcs_menu          
!------------------------------------------------------------------------
      SUBROUTINE detrend_menu() ! Curve fitting menu 
      IMPLICIT NONE                 
      INTEGER       :: i,j,k
      idt=2 ; gtr=1
      b(51:52)%ok=FA 
      b(41:50)%on=FA ; b(idt+41)%on=TR
      b(35:36)%ok=FA  ! Display options and save pars only for RCS
      CALL det_default()
      CALL new_chron(FA) 
      ADO: DO j=1, 10000    ! Menu attempts
        b(56)%ok=cc.GT.1 ; b(57)%ok=cc.LT.nc
        CALL erase()        ! Clear screen 
        CALL plot_det()
        CALL mwrite(1) ; CALL mwrite(2)
        DO k=1,8 ; CALL but_draw(k,"") ; ENDDO  ! Chronology Options
        CALL but_draw(19,"")  ! Sort option
        CALL but_draw(20,"")  ! Max SF iterations
        CALL but_draw(21,"")  ! Ignore 1st 'n' option
        DO k=20,24 ; CALL mwrite(k) ; ENDDO  ! Chronology Options
        CALL but_draw(23,"")
        CALL but_draw(24,cnam(cf))
        IF (bred) bhigh=red 
        CALL but_draw(32,"")
        bhigh=yellow
        CALL but_draw(33,"")  ! Chronology Options
        CALL but_draw(34,"")
        DO k=41,52 ; CALL but_draw(k,"") ; ENDDO 
        CALL mwrite(41)
        DO k=56,58                ! Draw button boxes
          IF (b(k)%ok) THEN
            SELECT CASE (k)
              CASE (56:57) ; CALL but_draw(k,"") 
              CASE (58) ; CALL but_draw(k,nam(tre(cc))) ; CALL mwrite(7)
           ENDSELECT
          ENDIF
        ENDDO
        CALL mouse_click(1,1,58)  ! Which button will be pressed
        IF ((mous.GE.1.AND.mous.LE.21).OR.  &
            (mous.GE.41.AND.mous.LE.52)) bred=TR ! bred=recalculate needed
        SELECT CASE (mous)      
        CASE (0)               ! Mouse missed "on" buttons
        CASE (1)               ! Transform
          IF (itn.EQ.1) THEN
            itn=2 ; ind=2      ! Transform and residuals
          ELSEIF (itn.EQ.2) THEN
            itn=3 ; ind=1      ! Basal Area Increment
          ELSEIF (itn.EQ.3) THEN
            itn=1 ; ind=1      ! No transform and ratios
          ENDIF
          b(3)%lab=indlab(ind) ;  b(mous)%lab=itnlab(itn)
        CASE (2)               ! Pith Offset / not
          IF (poo.EQ.2) THEN ; poo=1 ; ELSE ; poo=poo+1 ; ENDIF
          b(mous)%lab=poolab(poo)
        CASE (3)               ! Tree index calc
          IF (ind.EQ.2) THEN ; ind=1 ; ELSE ; ind=ind+1 ; ENDIF
          b(mous)%lab=indlab(ind)
        CASE (4)               ! Mean chronology
          IF (krb.EQ.2) THEN ; krb=1 ; ELSE ; krb=krb+1 ; ENDIF
          b(mous)%lab=krblab(krb)
        CASE (5)               ! Signal free / not
          IF (sfo.EQ.2) THEN ; sfo=1 ; ELSE ; sfo=sfo+1 ; ENDIF
          b(mous)%lab=sfolab(sfo) ; b(20)%ok=(sfo.EQ.2) 
        CASE (6)               ! Chron AR Model 
          IF (jrb.EQ.3) THEN ; jrb=1 ; ELSE ; jrb=jrb+1 ; ENDIF
          b(mous)%lab=jrblab(jrb)
        CASE (7)               ! Distribution adjust
          IF (idb.EQ.2) THEN ; idb=1 ; ELSE ; idb=idb+1 ; ENDIF
          b(mous)%lab=idblab(idb) 
        CASE (8)               ! CRN variance stabilize
          IF (isb.EQ.2) THEN ; isb=1 ; ELSE ; isb=isb+1 ; ENDIF
          b(mous)%lab=isblab(isb)
        CASE (19)                 ! Tree sort option 
          IF (tst.EQ.7) THEN ; tst=1 ; ELSE ; tst=tst+1 ; ENDIF
          b(mous)%lab=tstlab(tst) ; CALL tree_sort()
        CASE (20) ; CALL get_int(mous,i)    ! Max SF iterations
          IF (i.GE.2) sfono=i
          WRITE(b(mous)%lab,'(I4)') sfono
        CASE (21) ; CALL get_int(mous,i)    ! Ignore 1st 'n' rngs 
          IF (i.GE.0) THEN
            ignor=i ; ignok=TR
          ENDIF 
          WRITE(b(mous)%lab,'(I3)') ignor
        CASE (23) ; EXIT ADO                ! Quit program
        CASE (24) ; CALL name_menu(1)
          CALL new_chron(TR) ! New chronology name
        CASE (32) ; CALL new_chron(FA)    ! Calculate chronology
        CASE (33) ; CALL save_menu()      ! Save data or CRN
        CASE (34) ; CALL plot_ps(1)       ! Save current plot
!       CASE (35) ; CALL disp_menu()      ! Display menu
        CASE (41:48) ; b(41:50)%on=FA ; b(mous)%on=TR ; idt=mous-41  
          b(51:52)%ok=FA 
        CASE (49) ;    b(41:50)%on=FA ; b(mous)%on=TR ; idt=idtsl
          b(51)%ok=TR ; b(52)%ok=FA 
        CASE (50) ;    b(41:50)%on=FA ; b(mous)%on=TR ; idt=idtsn
          b(52)%ok=TR ; b(51)%ok=FA 
        CASE (51) ; CALL get_int(mous,i)   ! Detrend spline N
          IF (i.GT.9.OR.i.LT.9999) idtsl=i
          idt=idtsl ; WRITE(b(mous)%lab,'(I5)') idtsl
        CASE (52) ; CALL get_int(mous,i)   ! Detrend spline %N
          IF (i.LT.-9.OR.i.GT.-999) idtsn=i
          idt=idtsn ; WRITE(b(mous)%lab,'(I5)') idtsn
        CASE (56) ; IF (cc.GT.1) cc=cc-1  ! Previous Tree
        CASE (57) ; IF (cc.LT.nc) cc=cc+1 ! Next Tree
        CASE (58) ; CALL name_menu(2)     ! Wants select a tree
        END SELECT
      ENDDO ADO
      b(35:36)%ok=TR  ! Save pars back on
      RETURN 
      END SUBROUTINE detrend_menu          
!--------------------------------------------------------
      SUBROUTINE RCS_menu() ! Main RCS menu 
      IMPLICIT NONE                 
      CHARACTER(60) :: fnam
      LOGICAL       :: openc
      INTEGER       :: i,j,k
      CALL start_init()
      ADO: DO j=1, 10000    ! Menu attempts
        CALL erase()        ! Clear screen 
        CALL plot_det()
        CALL mwrite(1) ; CALL mwrite(2)
        DO k=1,21 ; CALL but_draw(k,"") ; ENDDO ! Chronology Options
        DO k=20,24 ; CALL mwrite(k) ; ENDDO     ! Chronology Options
        CALL but_draw(23,"")
        CALL but_draw(24,cnam(cf))
        IF (bred) bhigh=red 
        CALL but_draw(32,"")
        bhigh=yellow
        DO k=33,38 ; CALL but_draw(k,"") ; ENDDO ! Menu Options
        CALL mouse_click(1,1,38) ! Which button will be pressed
        IF (mous.GE.1.AND.mous.LE.21) bred=TR  ! bred=recalculate needed
        SELECT CASE (mous)      
        CASE (0)     ! Mouse missed "on" buttons
        CASE (1)                 ! Transform
          IF (itn.EQ.3) THEN
            itn=1 ; ind=1        ! No transform and ratios
          ELSEIF (itn.EQ.1) THEN
            itn=2 ; ind=2        ! Transform and residuals
            IF (gtr.NE.1) THEN
              gtr=1 ; b(12)%lab=gtrlab(gtr) ! Ring transform off         
            ENDIF
          ELSEIF (itn.EQ.2) THEN
            itn=3 ; ind=1        ! Basal Area Increment 
          ENDIF
          b(3)%lab=indlab(ind) ;  b(mous)%lab=itnlab(itn)
        CASE (2)               ! Pith Offset / not
          IF (poo.EQ.2) THEN ; poo=1 ; ELSE ; poo=poo+1 ; ENDIF
          b(mous)%lab=poolab(poo)
        CASE (3)               ! Tree index calc
          IF (ind.EQ.2) THEN ; ind=1 ; ELSE ; ind=ind+1 ; ENDIF
          IF (ind.EQ.2.AND.gtr.NE.1) gtr=1           
          b(mous)%lab=indlab(ind)
        CASE (4)               ! Mean chronology
          IF (krb.EQ.2) THEN ; krb=1 ; ELSE ; krb=krb+1 ; ENDIF
          b(mous)%lab=krblab(krb)
        CASE (5)               ! Signal free / not
          IF (sfo.EQ.2) THEN ; sfo=1 ; ELSE ; sfo=sfo+1 ; ENDIF
          b(mous)%lab=sfolab(sfo) ; b(20)%ok=(sfo.EQ.2) 
        CASE (6)               ! Chron AR Model 
          IF (jrb.EQ.3) THEN ; jrb=1 ; ELSE ; jrb=jrb+1 ; ENDIF
          b(mous)%lab=jrblab(jrb)
        CASE (7)               ! Distribution adjust
          IF (idb.EQ.2) THEN ; idb=1 ; ELSE ; idb=idb+1 ; ENDIF
          b(mous)%lab=idblab(idb) 
        CASE (8)               ! CRN variance stabilize
          IF (isb.EQ.2) THEN ; isb=1 ; ELSE ; isb=isb+1 ; ENDIF
          b(mous)%lab=isblab(isb)
        CASE (9)               ! Detrend RCS curve
          IF (rdt.EQ.7) THEN ; rdt=1 ; ELSE ; rdt=rdt+1 ; ENDIF
          b(mous)%lab=rdtlab(rdt)
          IF     (rdt.EQ.6) THEN
            IF (rdtno.LT.5.OR.rdtno.GT.999) rdtno=60
          ELSEIF (rdt.EQ.7) THEN
            IF (rdtno.LT.-999.OR.rdtno.GT.-5) rdtno=-10
          ENDIF 
          b(13)%ok=(rdt.GE.6) ; WRITE(b(13)%lab,'(I4)') rdtno
        CASE (10)              ! Single RCS / not
          IF     (src.EQ.3) THEN
            src=1 ; b(12)%ok=FA ; b(1:3)%ok=TR ; b(37)%ok=FA
            b(5)%ok=TR ; b(9)%ok=TR
            b(11)%ok=TR ; b(20)%ok=TR ; b(16:18)%ok=FA 
          ELSEIF (src.EQ.2) THEN
            fnam=cnam(cf) ; i=LEN_TRIM(fnam) 
            ID: DO i=i,3,-1 ; IF (fnam(i:i).EQ.".") EXIT ID ; ENDDO ID
            fnam=fnam(1:i)//"src"
            INQUIRE(FILE=fnam,EXIST=openc)  
            IF (openc) THEN   ! OK to use selected RCS
              CALL read_srcs()
              b(1:3)%ok=FA ; b(5)%ok=FA ; b(9)%ok=FA ; b(11)%ok=FA
              b(16)%ok=FA ; b(20)%ok=FA ; b(37)%ok=FA ; src=3
            ELSE
!             CALL out_err(TRIM(fnam)//" does not exist")
              src=1 ; b(12)%ok=FA ; b(16:18)%ok=FA ; b(37)%ok=FA
            ENDIF
          ELSEIF (src.EQ.1) THEN
            src=2 ; b(12)%ok=TR ; b(37)%ok=TR 
            b(16)%ok=TR ; b(18)%ok=TR 
          ENDIF
          b(mous)%lab=srclab(src)
        CASE (11)              ! RCS type (age, diam, basal etc)
          IF (trc.EQ.3) THEN ; trc=1 ; ELSE ; trc=trc+1 ; ENDIF
          b(mous)%lab=trclab(trc)
        CASE (12)              ! Ring transform
          IF (gtr.EQ.2) THEN ; gtr=1 ; ELSE ; gtr=2 ; ENDIF
          IF (src.EQ.1.AND.gtr.NE.1) gtr=1  ! Need Multi on
          IF ((itn.NE.1.OR.ind.NE.1).AND.gtr.NE.1) gtr=1  ! Power trans/ratios
          b(mous)%lab=gtrlab(gtr) 
        CASE (13) ; CALL get_int(mous,i)    ! Spline stiffness
          IF ((rdt.EQ.6.AND.i.GE.5.AND.i.LE.999) .OR. &
              (rdt.EQ.7.AND.i.LE.-5.AND.i.GE.-999)) THEN
            rdtno=i ; WRITE(b(mous)%lab,'(I4)') rdtno
          ENDIF
!       CASE (14) ! Not used 
!       CASE (15) ! Not used 
        CASE (16) ; CALL get_int(mous,i)    ! Multi RCS curves
          IF (i.GE.2.AND.i.LE.mx-1) srcno=MIN(i,MAX(nc/40,1))
          WRITE(b(mous)%lab,'(I4)') srcno
!       CASE (17) ! Not Used
        CASE (18)                           ! BFM crn options
          IF (bfc.EQ.3) THEN ; bfc=1 ; ELSE ; bfc=bfc+1 ; ENDIF
          b(mous)%lab=bfclab(bfc)
        CASE (19)                           ! Sort options
          IF (tst.EQ.7) THEN ; tst=1 ; ELSE ; tst=tst+1 ; ENDIF
          b(mous)%lab=tstlab(tst)
        CASE (20) ; CALL get_int(mous,i)    ! Max SF iterations
          IF (i.GE.2) sfono=i
          WRITE(b(mous)%lab,'(I4)') sfono
        CASE (21) ; CALL get_int(mous,i)    ! Ignore 1st 'n' rngs 
          IF (i.GE.0) THEN
            ignor=i ; ignok=TR
          ENDIF 
          WRITE(b(mous)%lab,'(I3)') ignor
!       CASE (22) ; Not Used
        CASE (23) ; EXIT ADO              ! Quit program
        CASE (24) ; CALL name_menu(1)
          CALL read_default()
          CALL new_chron(TR)              ! New chronology name
!       CASE (25) ; Not Used
        CASE (32) ; CALL new_chron(FA)    ! Calculate chronology
        CASE (33) ; CALL save_menu()      ! Save data or CRN
        CASE (34) ; CALL plot_ps(1)       ! Save current plot
        CASE (35) ; CALL disp_menu()      ! Display menu
        CASE (36) ; CALL write_default()  ! Save parameter values
        CASE (37) ; CALL srcs_menu()      ! Select RCS curves
        CASE (38) ; CALL tree_menu()      ! Look at trees
        END SELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE RCS_menu          
!--------------------------------------------------------
      SUBROUTINE tfigs_menu(ref1)  ! Creates output displays
      USE rfigs  
      IMPLICIT NONE
      INTEGER :: ref1    ! Menu selection
      INTEGER :: i,j,k,plold
      sccc=1 ; k=3 ; chs=1 ; b(272)%ok=FA ; plold=k
      JD: DO j=1,1000   ! Menu attempts  
        CALL ERASE()
        CALL mwrite(271) ; chc(1)=k 
        SELECT CASE (ref1) 
          CASE (1) ; CALL ch_disp(3,fin,rfnam(1:fin),chs,chc,1)
        ENDSELECT
        CALL but_draw(273,cnam(cf)) ; CALL mwrite(273)
        IF (b(272)%ok) CALL but_draw(272,"")
        CALL mouse_click(7,271,273)
        SELECT CASE (mous)
        CASE (271) 
          i=rw*INT(DBLE(3*(msx-cnx))/DBLE(gright-cnx))+ &
            INT(DBLE(rw*(b(mous)%y1-msy))/DBLE(b(mous)%y1-b(mous)%y2))+1
          IF (i.EQ.1) THEN
            EXIT JD
          ELSEIF (i.GT.1.AND.i.LE.fin) THEN
            plold=k ; k=i
            SELECT CASE (ref1) 
              CASE (1) ; CALL rf_plot(k,plold)
           ENDSELECT
          ENDIF
        CASE (272) ; EXIT JD
        CASE (273) ; CALL name_menu(1)
          CALL read_default()
          CALL new_chron(TR)    ! New chronology name
        ENDSELECT
      ENDDO JD
      RETURN 
      END SUBROUTINE tfigs_menu
!------------------------------------------------------------------------
      SUBROUTINE start_initt()    ! Initialise start menu
      USE rfigs
      IMPLICIT NONE                 
      CALL rfig_val()
      INQUIRE(FILE="RCSdefault.fil",EXIST=fileok)  
      IF (.NOT.fileok) THEN
        WRITE(*,'("File RCSdefault.fil missiing")')
        STOP
      ELSE
        OPEN(23,FILE="RCSdefault.fil",IOSTAT=ios,STATUS="OLD")
        IF (ios.NE.0) THEN
          WRITE(*,'("Error",I7," open RCSdefault.fil")') ios
          STOP 
        ENDIF
        READ(23,*,IOSTAT=ios)   ! Header line
        IF (ios.NE.0) THEN
          WRITE(*,'("Error",I7," read header RCSdefault.fil")') ios
          STOP 
        ENDIF
        READ(23,'(I4)',IOSTAT=ios) screenw
        IF (ios.NE.0) THEN
          WRITE(*,'("Error",I7," read width RCSdefault.fil")') ios
          STOP 
        ENDIF
        READ(23,'(I4)',IOSTAT=ios) screenh
        IF (ios.NE.0) THEN
          WRITE(*,'("Error",I7," read height RCSdefault.fil")') ios
          STOP 
        ENDIF
        CLOSE(23) 
      ENDIF
      CALL METAFL('XWIN')
      CALL WINDOW(0,0,screenw,screenh)  
      CALL SETPAG('DA4L')
      CALL SCRMOD('REVERS')       ! White background
      CALL DISINI()
      CALL WINMOD('NONE')         ! DISINI drops out 
      CALL CSRMOD('READ','POS')   ! Control cursor reading
      CALL SETVLT('VGA')          ! Select 16 colour table
      CALL SIMPLX()
      CALL HEIGHT(25)             ! Character height - plot coords
      CALL SHDPAT(16)             ! Fill areas
      CALL TICPOS('REVERS','XY')  ! Internal tick marks
      CALL HNAME(20)              ! Character height for Axis names
      CALL LABDIG(-1,'X')         ! Digits after decimal -1=integer
      CALL TICKS(10,'X')          ! Number of ticks between labels
      CALL TICKS(5,'Y')           ! Number of ticks between labels
      CALL read_default() 
      INQUIRE(FILE=lnam,EXIST=lfileok)  
      IF (lfileok) THEN
        CALL read_list()    
        INQUIRE(FILE=cnam(cf),EXIST=chronok)  
        IF (chronok) THEN
          nc=0 ; CALL read_rft(cnam(cf))
          CALL detrend()
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE start_initt
!------------------------------------------------------------------
      SUBROUTINE start_menu()           ! Start page menu 
      IMPLICIT NONE                 
      INTEGER  :: j,k
      lnam="ars.fil"
      CALL start_initt()
      IF (journal) THEN     
        b(451)%lab="Journal On" ; b(347)%lab=" On" 
      ELSE
        b(451)%lab="Journal Off" ; b(347)%lab="Off" 
      ENDIF
      ADO: DO j=1, 10000    ! Menu attempts
        CALL ERASE() ; CALL but_draw(452,lnam) ; CALL mwrite(303)     
        CALL but_draw(451,"")    
        IF (.NOT.lfileok) CALL mwrite(306)
        CALL but_draw(453,cnam(cf)) ; CALL mwrite(304)
        IF (.NOT.chronok) CALL mwrite(305)
        CALL but_draw(454,"") ; CALL mwrite(307) ; CALL mwrite(319)
        CALL but_draw(455,"") ; CALL mwrite(308)
        CALL but_draw(456,"") ; CALL mwrite(309)
        CALL but_draw(458,"") ; CALL mwrite(311)
        CALL mouse_click(13,451,458) ! Which button will be pressed
        SELECT CASE (mous)      
        CASE (0)                     ! Mouse missed "on" buttons
        CASE (451)
          IF (journal) THEN          ! Turn the journal on/off
            journal=FA ; b(mous)%lab="Journal Off"  
          ELSE
            journal=TR ; b(mous)%lab="Journal On"  
          ENDIF
        CASE (452) ; CALL get_str(mous)  ! Change list file name
          IF (LEN(TRIM(b(mous)%lab)).GT.2) THEN
            lnam=b(mous)%lab(1:20) ; CALL start_initt()
          ENDIF
        CASE (453)                   ! Select a chronology
          CALL name_menu(1) ; IF (sccc.GT.nfil-rw*5) sccc=1
          INQUIRE(FILE=cnam(cf),EXIST=chronok)  
        CASE (454) ; EXIT ADO        ! Wants to quit program
        CASE (455)              
          CALL detrend_menu()        ! Curve fitting options       
          CALL read_default()        ! Reset for RCS
          CALL new_chron(TR)
        CASE (456) ; CALL RCS_menu()             
        CASE (458)
          CALL tfigs_menu(1) ; CALL detrend()
        END SELECT
      ENDDO ADO
      RETURN 
      END SUBROUTINE start_menu          
!------------------------------------------------------------------------
     END MODULE crustprocs 
