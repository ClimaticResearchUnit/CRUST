! Copyright (C) 2013, Thomas M. Melvin and Keith R. Briffa, see 
! the GNU General Public License.
      MODULE DETDATA     ! Common data for CRU RCS
      IMPLICIT NONE
      INTEGER,PARAMETER   :: numbut =     480  ! Number of button boxes
      INTEGER,PARAMETER   :: mxc    =     120  ! Max number of chronology names
      INTEGER,PARAMETER   :: mxs    =    6000  ! Max number of trees
      INTEGER,PARAMETER   :: mxd    = mxs*250  ! Max rings
      INTEGER,PARAMETER   :: mxt    =    4700  ! Max years in a tree
      INTEGER,PARAMETER   :: mxy    =   10000  ! Max years in chronology
      INTEGER,PARAMETER   :: maxm   =     360  ! Max met series 
      INTEGER,PARAMETER   :: mmx    =     360  ! Max met years 
      INTEGER,PARAMETER   :: ccol   =       3  ! Max crn columns
      INTEGER,PARAMETER   :: ncol   =       7  ! Max tree columns
      INTEGER,PARAMETER   :: rw     =      30  ! Number of rows in selection box
      INTEGER,PARAMETER   :: stabut =     268  ! Status button number
      INTEGER,PARAMETER   :: fin =         90  ! Figures in menu choice
      
      INTEGER                      :: chs      ! Number of items selected
      INTEGER,DIMENSION(mxs)       :: chc      ! Numbers of items selected

      REAL(8),DIMENSION(mxd)       :: x        ! Raw data storage
      REAL(8),DIMENSION(mxd)       :: tx       ! Transformed raw
      REAL(8),DIMENSION(mxd)       :: fx       ! Transformed raw
      REAL(8),DIMENSION(mxd)       :: cx       ! Detrend curves
      REAL(8),DIMENSION(mxd)       :: dx       ! Tree indices
      REAL(8),DIMENSION(mxd)       :: ax       ! AR model indices
      LOGICAL,DIMENSION(mxd)       :: xok      ! False for missing rings

      INTEGER,DIMENSION(mxy)       :: iwk1     ! Integer array
! For each tree	  
      INTEGER,DIMENSION(mxs)       :: ad       ! Address first ring each tree
      INTEGER,DIMENSION(mxs)       :: fy       ! First year of each tree
      INTEGER,DIMENSION(mxs)       :: ly       ! Last year of each tree
      INTEGER,DIMENSION(mxs)       :: yr       ! Length of each tree
      INTEGER,DIMENSION(mxs)       :: pth      ! First pith year each tree
      REAL(8),DIMENSION(mxs)       :: pthr     ! Pith radius for each tree
      REAL(8),DIMENSION(mxs)       :: tso      ! Sorted tree values
      INTEGER,DIMENSION(mxs)       :: tre      ! Sorted tree order
      INTEGER,DIMENSION(mxs)       :: trr      ! RCS built using this tree
      INTEGER,DIMENSION(mxs)       :: trm      ! RCS to detrend this tree
      REAL(8),DIMENSION(mxs)       :: bfdx     ! Mean index factors
      INTEGER,DIMENSION(mxs)       :: jdt      ! Core detrend params
      INTEGER,DIMENSION(mxs)       :: ralt     ! Relative altitude each tree
      LOGICAL,DIMENSION(mxs)       :: raok     ! Relative altitude OK
      CHARACTER(12),DIMENSION(mxs) :: nam      ! Tree names
! For each chronology 
      REAL(8),DIMENSION(mxy,mxc)   :: crn      ! Chronology indices
      INTEGER,DIMENSION(mxy,mxc)   :: num      ! Chronology tree count
      LOGICAL,DIMENSION(mxy,mxc)   :: okc      ! Chronology index present
      INTEGER,DIMENSION(mxc)       :: cfy,cly,cyr ! CRN first,last and years
      CHARACTER(20),DIMENSION(mxc) :: wnam     ! short names of data files
      CHARACTER(60),DIMENSION(mxc) :: cnam     ! Names of raw data files
      INTEGER,DIMENSION(15,mxc)    :: pars     ! Parameter settings
      CHARACTER(16)                :: figm     ! .ps file name
! Multiple RCS fields
      INTEGER,PARAMETER            :: mx=12    ! Max number of curves + 1
      INTEGER                      :: cur      ! Number of SARCS curves
! Age based RCS curves
      INTEGER,DIMENSION(mx)        :: sfy,sly  ! First/last year by age
      REAL(8),DIMENSION(mxt,mx)    :: msmo     ! Smoothed RCS curves
      REAL(8),DIMENSION(mxt,mx)    :: mssd     ! Standard deviation
      REAL(8),DIMENSION(mxt,mx)    :: mserr    ! Standard error
      REAL(8),DIMENSION(mxt,mx)    :: mval     ! Mean RCS values
      INTEGER,DIMENSION(mxt,mx)    :: mcnt     ! Counts of rings
      LOGICAL,DIMENSION(mxt,mx)    :: mok      ! Rings present
! Diameter based RCS curves
      INTEGER,DIMENSION(mx)        :: dfy,dly  ! First/last year by diam
      REAL(8),DIMENSION(mxt,mx)    :: dval     ! Mean RCS diameter
      REAL(8),DIMENSION(mxt,mx)    :: dsmo     ! Smoothed RCS diameter
      INTEGER,DIMENSION(mxt,mx)    :: dcnt     ! Counts of rings
      LOGICAL,DIMENSION(mxt,mx)    :: dok      ! Rings present
! Chronologies each RCS curve
      INTEGER                      :: xfy,xly,xyr ! First,last and years
      INTEGER,DIMENSION(mx)        :: xfa,xla,xaa ! First,last and years
      REAL(8),DIMENSION(mxy,mx)    :: xcrn     ! Current chronology
      REAL(8),DIMENSION(mxy,mx)    :: xcsd     ! Chronology Stan Dev
      REAL(8),DIMENSION(mxy,mx)    :: xcsm     ! Chronology Smoothed
      INTEGER,DIMENSION(mxy,mx)    :: xnum     ! Current counts
      LOGICAL,DIMENSION(mxy,mx)    :: cok      ! False for missing values
      INTEGER,DIMENSION(mxy)       :: ccnt     ! Chronology counts

      CHARACTER(60)   :: lnam        ! List of names filename
      CHARACTER(60)   :: namc        ! Output CRN name 
      CHARACTER(60)   :: namr        ! Output RAW name 
      CHARACTER(60)   :: dirc        ! CRN directory name 
      CHARACTER(60)   :: dirr        ! RAW directory name 
      LOGICAL         :: Append      ! Add to existing file or not 
      CHARACTER(160)  :: mess        ! General message field
      INTEGER         :: ios         ! Input/Output error 
       
      INTEGER  :: nfil         ! Number raw data files
      INTEGER  :: icf          ! Number of raw files  
      INTEGER  :: cf           ! Current raw data file  
      INTEGER  :: nc           ! Number of trees
      INTEGER  :: cc           ! Currently open tree number
      INTEGER  :: sccc         ! First scroll tree number
      LOGICAL  :: srcok        ! Selected RCS files loaded
      LOGICAL  :: fileok       ! File list file exists
      LOGICAL  :: journal      ! If on write to journal
      INTEGER  :: mous         ! Mouse parameters
      INTEGER  :: msx,msy      ! button selection x,y address
      INTEGER  :: ecfy,ecly    ! Chronologies earliest and latest years   
      INTEGER  :: mcf          ! Last entry Index and Met file
      LOGICAL  :: lfileok      ! File list file exists
      LOGICAL  :: chronok      ! Current chronology file exists

      INTEGER  :: IDT = -2     ! Default value
      INTEGER  :: ITN          ! Transform + BAI
      INTEGER  :: RDT          ! Detrend RCS curve
      INTEGER  :: IND          ! Tree index calc
      INTEGER  :: KRB          ! Mean chronology
      INTEGER  :: ISB          ! CRN variance stabilize
      INTEGER  :: SFO          ! Signal free / not
      INTEGER  :: JRB          ! Chron AR Model 
      INTEGER  :: POO          ! Pith Offset / not
      INTEGER  :: SRC          ! Single RCS / multiple / Selected
      INTEGER  :: TRC          ! RCS type (age, diam, age-diam)
      INTEGER  :: GTR          ! Ring transform
      INTEGER  :: TST          ! Tree sort option
      INTEGER  :: BFC          ! BFM chronology options
      INTEGER  :: IDB          ! Tree index manipulation
      INTEGER  :: RFT          ! Output Raw format
      INTEGER  :: CFT          ! Output CRN format

      INTEGER  :: IDTSL        ! Curve-fit spline length
      INTEGER  :: IDTSN        ! Curve-fit % spline length
      INTEGER  :: RDTNO        ! Spline length
      INTEGER  :: SRCNO        ! Number Multiple RCS curves
      INTEGER  :: SFONO        ! Maximum signal-free iterations
      INTEGER  :: TRROK        ! RCS use DET or not
      INTEGER  :: IGNOR        ! Ignore 1st 'n' rings
      LOGICAL  :: IGNOK        ! Need to call ignore

      INTEGER  :: LTH          ! Line thickness
      INTEGER  :: LST          ! Line style 
      INTEGER  :: LCO          ! Line colour
       
      CHARACTER(14),DIMENSION(3),PARAMETER :: itnlab = &
        (/"No Transform  ","Adaptive Power","Basal Area    "/)
      CHARACTER(14),DIMENSION(7),PARAMETER :: rdtlab = &
        (/"Age Dep Smo   ","No Smooth     ",&
          "Neg Exp Smo   ","Trend Smooth  ",&
          "Huger Smo     ","Splin 50% Var ","Splin 50% N   "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: indlab = &
        (/"Ratios        ","Residuals     "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: krblab = &
        (/"Arith Mean    ","Robust Mean   "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: isblab = &
        (/"Var.Stab Off  ","RBAR Stab On  ","HF Stabilise  "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: sfolab = &
        (/"Sig Free Off  ","Sig Free On   "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: jrblab = &
        (/"STD Chron     ","RES Chron     ","ARS Chron     "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: poolab = &
        (/"PO on         ","PO off        "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: srclab = &
        (/"Single RCS    ","Multi RCS     ","Select RCS    "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: trclab = &
        (/"Age RCS       ","Diam RCS      ","Age+Diam RCS  "/) 
      CHARACTER(14),DIMENSION(2),PARAMETER :: gtrlab = &
        (/"No Adjust     ","Mean Single   "/)
      CHARACTER(14),DIMENSION(7),PARAMETER :: tstlab = & ! Tree sort option
        (/"Not Sorted    ","Tree Age      ","Tree Size     ", &
          "Growth Rate   ","Tree Name     ","Pith Year     ", &
          "Last Year     "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: bfclab = &
        (/"Mean trees    ","Arith Mn CRN  ","CRN Mn=1.0    "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: idblab = &
        (/"Ratios CRN    ","Normal CRN    "/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: rftlab = &
        (/"Tucson Data   ","Heidelberg    ","Compact Data  "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: cftlab = &
        (/"Tucson CRN    ","Column Format "/)
      CHARACTER(4),DIMENSION(6),PARAMETER  :: datext = &
        (/".rww",".trn",".sfm",".crv",".ind",".res"/)
      CHARACTER(4),DIMENSION(2),PARAMETER  :: cftext=(/".crn",".col"/)
      CHARACTER(14),DIMENSION(3),PARAMETER :: lthlab = &
        (/"Thin Line     ","Medium Line   ","Thick Line    "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: lstlab = &
        (/"Solid Line    ","Dashed Line   "/)
      CHARACTER(14),DIMENSION(2),PARAMETER :: trrlab = &
        (/"No RCS        ","As DET        "/)
      CHARACTER(3),DIMENSION(-2:9),PARAMETER :: crtyp = &
        (/"RCS","-  ","RAW","EXP","EXN","ANY", &
          "NEG","HOR","HUG","GEX","SPL","%SP"/)
      CHARACTER(14),DIMENSION(-2:9),PARAMETER :: idtlab = &
       (/"RCS method    ","Not used      ","No Detrend    ", &
         "Exp/Any Line  ","Exp/Neg Line  ","Sloping line  ", &
         "Neg Slope Line","Horiz Line    ","Hugershoff    ", &
         "General Exp   ","Splin 50% Var ","Splin 50% -N  "/)

      INTEGER,PARAMETER  :: white   = 0
      INTEGER,PARAMETER  :: navy    = 1
      INTEGER,PARAMETER  :: green   = 2
      INTEGER,PARAMETER  :: teal    = 3
      INTEGER,PARAMETER  :: maroon  = 4
      INTEGER,PARAMETER  :: purple  = 5
      INTEGER,PARAMETER  :: brown   = 6
      INTEGER,PARAMETER  :: silver  = 7
      INTEGER,PARAMETER  :: grey    = 8
      INTEGER,PARAMETER  :: blue    = 9
      INTEGER,PARAMETER  :: lime    = 10
      INTEGER,PARAMETER  :: cyan    = 11
      INTEGER,PARAMETER  :: red     = 12
      INTEGER,PARAMETER  :: pink    = 13
      INTEGER,PARAMETER  :: yellow  = 14
      INTEGER,PARAMETER  :: black   = 15

      INTEGER,DIMENSION(8),PARAMETER  :: col = &  
            (/yellow,cyan,grey,green,maroon,red,black,blue/)
      INTEGER,DIMENSION(16),PARAMETER :: lcol = &
        (/15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,0/)  ! Line colours
 
      LOGICAL,PARAMETER :: TR=.TRUE. , FA=.FALSE. ! True and false
      REAL(8),PARAMETER :: epsi = 0.000000001D0   ! Small number
      REAL(8),PARAMETER :: pi = 3.1415926535897935D0 

      TYPE posn                       ! Button box parameters   
        INTEGER       :: x1,y1,x2,y2  ! button corners
        CHARACTER(60) :: lab          ! button label
        LOGICAL       :: ok           ! display status
        LOGICAL       :: on           ! on/off status
      END TYPE
      TYPE pos2                       ! Message box    
        INTEGER       :: x1,y1,x2,y2  ! Message corners
      END TYPE
      TYPE(posn),DIMENSION(numbut) :: b      ! Mouse buttons
      TYPE(pos2),DIMENSION(320)    :: madd   ! Message addresses
      CHARACTER(60),DIMENSION(320) :: mtext  ! Message text
     
      INTEGER            :: maxxs, maxys   ! size of PC screen - pixels
      INTEGER,PARAMETER  :: sizex = 2000   ! size of graphics screen
      INTEGER,PARAMETER  :: sizey = 3000
      INTEGER            :: screenh,screenw  ! Pixels to use for screen
      INTEGER,PARAMETER  :: cnx   = 200    ! bottom left corner
      INTEGER            :: screen         ! device number for graphics screen
      INTEGER            :: gright = 2400  ! Lower graph = G2
      INTEGER,PARAMETER  :: btop = 200     ! button table button top/bottom
      INTEGER,PARAMETER  :: bbot = 1700
      INTEGER            :: bhigh          ! Button highlight colour
      LOGICAL            :: bred           ! Red highlight
      LOGICAL            :: toplin = .TRUE.  !  True - use met data
      LOGICAL            :: botlin = .FALSE. !  False - use index data
      INTEGER :: grl,grr,grt,grb   ! current graf - left,right,top and bottom

! Display style
      TYPE dst                  ! Displayed lines   
        LOGICAL       :: on     ! Line status
        INTEGER       :: th     ! Thickness
        INTEGER       :: st     ! Style
        INTEGER       :: co     ! Colour
        CHARACTER(30) :: lab    ! Description
      END TYPE
      TYPE(dst),DIMENSION(103) :: lins  ! For each line being plotted
      LOGICAL :: oocol  ! Buttons on/off or change line style option
      LOGICAL :: psok   ! Writing a .eps file
! Display values
      INTEGER :: RDfy      ! RCS display first year
      INTEGER :: RDly      ! RCS display last year
      INTEGER :: CDfy      ! Chronology display first year
      INTEGER :: CDly      ! Chronology display last year
      INTEGER :: CDsp      ! Chronology smoothing
      INTEGER :: TDfy      ! Tree display first year
      INTEGER :: TDly      ! Tree display last year
      INTEGER :: RBsp      ! RBar segment length
! From chrondata
      REAL(8),DIMENSION(1:mmx,maxm)    :: met    ! Met file values
      LOGICAL,DIMENSION(1:mmx,maxm)    :: okm    ! Met value exists
      CHARACTER(20),DIMENSION(maxm)    :: mnam   ! Names of met monthly files
      CHARACTER(60),DIMENSION(maxm/12) :: metn   ! Raw met file names
      INTEGER                          :: mfil   ! Number of raw met files
      INTEGER,DIMENSION(maxm)  :: mfy,mly,myr    ! First,last,years rain/temp               
      INTEGER                  :: fm,lm          ! Earliest / latest met years

      CHARACTER(3),DIMENSION(12),PARAMETER :: mth = (/"Jan","Feb", &
        "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

      END MODULE DETDATA
