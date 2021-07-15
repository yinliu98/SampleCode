module SimulationParameters
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif

  !$ use omp_lib
  implicit none
  !****************************************************************
  !
  !  SimulationParameters
  !   Simulation parameters, renormalize methods, initialize methods
  !    Methods:
  !      parameterCheck
  !      renormalize
  !      initialize
  !
  !          Version FP1.1   November, 2016
  !%****************************************************************

  ! With ElectroStatic mode / Without Electrostatic mode
  logical         ,parameter :: HAS_ELECTROSTATIC_COMPONENT = .false.

  integer,parameter          :: NTIME   = 262144      ! Number of time step
  ! Grid parameter
  double precision,parameter :: DX      = 1.0d0       ! mesh size
  double precision,parameter :: DT      = 0.0078125d0 ! Timestep 1/128
  integer         ,parameter :: NX      = 32768       ! Number of step
  double precision,save      :: cv      = 100.0d0     ! Speed of light

  ! Particle parameters
  integer         ,parameter :: NS            = 2          ! Number of spicies
  integer         ,parameter :: NP(NS)        = (/ NX*256, NX*256  /) ! Number of particles !!! should be multiples of 64
  double precision,save      :: omega_p(NS)   = (/ 4.00d0, 4.0d-1  /) ! Plasma Frequency (/ dsqrt(2.0d0),dsqrt(2.0d0) /)
  double precision,parameter :: QM(NS)        = (/ -1.0d0, -1.0d0  /) ! Charge to Mass retio

  double precision,save      :: v_para(NS)    = (/ 1d0, 25d0/) ! Velocity to parallel
  double precision,save      :: v_perp(NS)    = (/ 1d0, 30d0/) ! Velocity to perpendicular

  double precision,save      :: v_d(NS)       = (/ 0d0, 0d0/) ! Drift velocity
  double precision,parameter :: PITCH_DEG(NS) = (/ 0d0, 0d0/) ! Pitch angle
  logical         ,parameter :: SUBTRACT_BI_MAXWELLIAN(NS) =  (/ .false., .true./)
  double precision,parameter :: SUBTRACT_BETA(NS)          =  (/ 0.0d0  , 0.3d0 /) ! A parameter is neglected for the species for SUBTRACT_BI_MAXWELLIAN(NS) = .false..
  double precision,parameter :: SUBTRACT_RHO(NS)           =  (/ 0.0d0  , 1.0d0 /)
  double precision,parameter :: IS_PARABOLIC_DISTRIBUTION(NS) = (/.false., .true. /)

  ! Injection of particles
  logical         ,parameter :: ENABLE_INJECTION           =  .false.
  logical         ,parameter :: IS_INJECTION(NS)           =  (/.false., .false. /)
  integer         ,parameter :: iINJECTION_TIME(NS)        =  (/0      , 0       /)

  double precision,save      :: omega_c       =  -1.0d0        ! Cyclotron Frequency should be minu
  ! magnetic field diraction in (r, THETA_DEG)
  double precision,parameter :: THETA_DEG = 0d0

  ! Particle---------------------------------------------------------------------------------------
  ! Magnetic Mirror motion
  ! CAUTION :: it works ONLY THETA_DEG = 0
  logical         ,parameter :: EXT_MIRROR       = .true.
  integer         ,parameter :: MIRROR_NTIME     = NTIME/8  ! Number of time step for stable state
  double precision,save      :: ext_mirror_a     = 5.0d-10  ! B = Beq(1+EXT_MIRROR_A*x**2)
  !Boundary Elimination
  logical         ,parameter :: PARTICLE_LOSSCORN_ELIMINATE  = .true.
  double precision,parameter :: LOSSCORN_ANGLE_DEG           = 5d0     ! losscorn angle of the system

  ! Field -----------------------------------------------------------------------------------------
  ! External Current
  ! CAUTION: These parameters are neglected in MIRROR INIT program
  logical         ,parameter :: EXT_CURRENT       = .true.
  integer         ,parameter :: iANTENNA_OFFSET   = 0        ! Less than NX/2 - NO NEED TO OFFSET
  integer         ,save      :: ctime             = 2097152  ! Duration of wave, shown in number of cyclotron gyro motions.
  integer         ,save      :: cwidth_begin      = 128  ! time lag to start generating waves
                                                         ! Note: this parameter is also used as CWITDH_END when EXT_SOLITARY_WAVE is true.
  integer         ,save      :: cwidth_end        = 128  ! time lag to quit generateing  wave packet

  double precision,save      :: ext_j_amp         = 12d0
  double precision,save      :: ext_omega_a       = 0.0d0   ! f = a * t + b0
  double precision,save      :: ext_omega_b0      = 0.30d0

  !Wave Number Analysis
  integer         ,parameter :: iX_WAVE_NUMBER_DIV = 2

  ! Charge neutrization
  logical         ,parameter :: LOCALIZED_CHARGE_AVG = .true.

  ! Boundary
  logical         ,parameter :: FIELD_ABSORB_BOUNDARY  = .true.
  integer         ,parameter :: DAMPING_REGION_ND      = 2048
  double precision,parameter :: DAMPING_WAVE_FREQUENCY = 0.30d0 ! Low frequency gives more strict condition (high DANPING COEFFICIENT)
  double precision,parameter :: DAMPING_INTEGRAL_K     = 4d0    ! this parameter should be between 3.0 to 4.0 from numerical experiment.

  ! Record parameters -----------------------------------------------------------------------------
  integer         ,parameter :: ISKIP  = 4096
  integer         ,parameter :: BSKIP  = 262144*8 !SKIPS for backup
  integer         ,parameter :: IXSKIP = 256      !SKIPS to X direction in dynamic spectrum

  ! Histogram Number of division
  ! Histogram Number of division
  integer         ,parameter :: iV_HISTO = 11
  double precision,parameter, dimension(iV_HISTO) :: HISTO_X_POS = (/-12288d0, -8192d0, -4096d0, -2048d0, -1024d0, 0d0, 1024d0, 2048d0, 4096d0, 8192d0, 12288d0/)
  double precision,parameter, dimension(iv_HISTO) :: HISTO_XW    = (/64d0,     64d0,    64d0,    64d0,    64d0,    64d0,  64d0,   64d0,   64d0,   64d0, 64d0/)
  integer         ,parameter :: HISTOGRAM_NUM_V = 32

  !WPIA
  integer         ,parameter :: WPIA_NUM_V_PARA = 200
  integer         ,parameter :: WPIA_NUM_ZETA   = 90
  ! X-zeta parallel velocity division
  integer         ,parameter :: iX_ZETA = 10
  ! MIN or MAX velocity is written as ratio to CV
  double precision,parameter, dimension(iX_ZETA) :: X_ZETA_MIN = (/-0.54822086d0, -0.32837618d0, -0.19498561d0, -0.10788471d0, -0.06246954d0, 0.0d0       , 0.06246954d0, 0.10788471d0, 0.19498561d0, 0.32837618d0/)
  double precision,parameter, dimension(iX_ZETA) :: X_ZETA_MAX = (/-0.32837618d0, -0.19498561d0, -0.10788471d0, -0.06246954d0, 0.0d0        , 0.06246954d0, 0.10788471d0, 0.19498561d0, 0.32837618d0, 0.54822086d0/)

  ! PARTICLE TRACKER
  logical          ,parameter :: PARTICLE_TRACKER = .true. !when we set PARTICLE_TRACKER is true, bucket sort is disabled
  integer          ,parameter :: iZETA_V_PARA      = 12
  integer          ,parameter :: TRACKER_MAX_PARTICLE_PER_PROCESS = 256

  integer          ,parameter, dimension(iZETA_V_PARA) :: TRACKER_ITIME   = (/8449, 8449, 11521, 11521, 25601, 25601, 10241, 10241, 26113, 26113, 44673, 44673/)
  integer          ,parameter, dimension(iZETA_V_PARA) :: TRACKER_SPECIES = 2 
  integer          ,parameter, dimension(iZETA_V_PARA) :: TRACKER_POSITION_NUM = (/2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2/) ! Number defined in iV_HISTO
  ! MIN or MAX velocity is written as ratio from 0 to 2*pi (example) pi/2 => 0.5d0, pi => 1d0, 2*pi =>2d0
  double precision ,parameter, dimension(iZETA_V_PARA) :: TRACKER_ZETA_MIN = (/1.5d0, 0d0   , 1.25d0, 0.5d0, 0.5d0, 0.25d0, 0.25d0, 0.25d0, 1.25d0, 0.5d0, 0.5d0, 0.6d0/)
  double precision ,parameter, dimension(iZETA_V_PARA) :: TRACKER_ZETA_MAX = (/2.0d0, 0.75d0, 1.75d0, 1.0d0, 1.0d0, 1.5d0 , 0.5d0 , 0.5d0 , 1.75d0, 1.5d0, 1.0d0, 1.5d0/)
  ! MIN or MAX velocity is written as ratio to CV
  double precision ,parameter, dimension(iZETA_V_PARA) :: TRACKER_V_PARA_MIN = (/-0.28d0, -0.28d0, -0.22d0, -0.28d0, -0.20d0, -0.26d0, 0.20d0, 0.22d0, -0.24d0, -0.30d0, -0.24d0, -0.16d0/)
  double precision ,parameter, dimension(iZETA_V_PARA) :: TRACKER_V_PARA_MAX = (/-0.24d0, -0.24d0, -0.18d0, -0.24d0, -0.16d0, -0.24d0, 0.22d0, 0.26d0, -0.22d0, -0.26d0, -0.20d0, -0.12d0/)
  !Label is determined by following rule. (itime)_(v_hist_number)_(d:depression or e:enhancement) (example) 8449_002_d
  character(len=30),parameter, dimension(iZETA_V_PARA) :: TRACKER_LABEL =(/'8449_002_d', '8449_002_e', '11521_002_d', '11521_002_e', '25601_002_d', '25601_002_e','10241_002_d', '10241_002_e', &
                                                                          '26113_007_d', '26113_007_e', '44673_007_d', '44673_007_e' /)

  character(len=5),parameter    :: VERSION = '2.0.0'

  integer,parameter,dimension(17) :: ISEED = (/383318173,863273061,392709476,220722442,285831396,&
                                               248143550,770881566, 67591433,866577399,366621657,&
                                               487806116,310673413,881339797,775845601,855232128,&
                                               769700443,988880450/) 
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  private
  !Parameters/Variables
  public :: NTIME
  public :: DX, DT, NX
  public :: cv
  public :: NS, NP, omega_p, QM, v_para, v_perp, v_d, PITCH_DEG, omega_c
  public :: THETA_DEG
  public :: ENABLE_INJECTION, IS_INJECTION, iINJECTION_TIME
  public :: SUBTRACT_BI_MAXWELLIAN, SUBTRACT_BETA, SUBTRACT_RHO
  public :: IS_PARABOLIC_DISTRIBUTION
  public :: HAS_ELECTROSTATIC_COMPONENT
  public :: EXT_MIRROR, MIRROR_NTIME, ext_mirror_a
  public :: PARTICLE_LOSSCORN_ELIMINATE, LOSSCORN_ANGLE_DEG
  public :: EXT_CURRENT, ctime, cwidth_begin, cwidth_end, iANTENNA_OFFSET
  public :: ext_j_amp, ext_omega_a, ext_omega_b0
  public :: ISKIP, BSKIP, IXSKIP
  public :: iX_WAVE_NUMBER_DIV
  public :: LOCALIZED_CHARGE_AVG
  public :: FIELD_ABSORB_BOUNDARY, DAMPING_REGION_ND, DAMPING_WAVE_FREQUENCY, DAMPING_INTEGRAL_K
  public :: iV_HISTO, HISTO_X_POS, HISTO_XW, HISTOGRAM_NUM_V
  public :: WPIA_NUM_V_PARA, WPIA_NUM_ZETA
  public :: iX_ZETA, X_ZETA_MIN, X_ZETA_MAX
  public :: PARTICLE_TRACKER, iZETA_V_PARA, TRACKER_MAX_PARTICLE_PER_PROCESS
  public :: TRACKER_ITIME, TRACKER_SPECIES
  public :: TRACKER_POSITION_NUM, TRACKER_ZETA_MIN, TRACKER_ZETA_MAX, TRACKER_V_PARA_MIN, TRACKER_V_PARA_MAX
  public :: TRACKER_LABEL
  public :: VERSION, ISEED
  !Methods
  public parametersCheck, checkBinaryFileParameters
contains
  subroutine parametersCheck(lerr)
    logical,intent(out) :: lerr
    ! Condition check
    double precision :: cfl_vel  , debye_len,  gyro_len
    double precision :: omega_c_mirror = 0.0d0
    logical ::          lcfl_vel  , ldebye_len, lgyro_len

    lcfl_vel   = .false.
    ldebye_len = .false.
    lgyro_len  = .false.

    cfl_vel   = DX/DT
    debye_len = minval(dsqrt(v_perp**2+v_para**2+v_d**2)/omega_p)

    if (EXT_MIRROR == .true.) omega_c_mirror = omega_c * (1.0d0+ext_mirror_a*0.25d0*dble(NX)**2*DX**2)
    gyro_len  = max(maxval(omega_p),dabs(omega_c),dabs(omega_c_mirror)) * DT

    if (cfl_vel   >  cv      )  lcfl_vel   = .true.
    if (debye_len >= DX*0.5d0)  ldebye_len = .true.
    if (gyro_len  <= 6.25d-2 )  lgyro_len  = .true.

    write (*,'(a16,a30,a2,L1)')  'CFL (dx/dt > cv)',                             '',': ', lcfl_vel
    write (*,'(a31,a15,a2,L1)')  'Debye length (lambda_D >= dx/2)',              '',': ', ldebye_len
    write (*,'(a45,a1 ,a2,L1)')  'Boris arc length (max(omega_p,omega_c) <=0.1)','',': ', lgyro_len

    if ((lcfl_vel .and. ldebye_len .and. lgyro_len) .eqv. .false.) then
       print *, "Error: Invalid simulation parameters."
       lerr = .true.
       return
    else
       lerr = .false.
    endif

    ! DAMPING_REGION
    if (DAMPING_REGION_ND > NX/2 ) then
      print *, "Error: Invalid Damping region ND > NX / 2 or NY / 2 "
       lerr = .true.
       return
    else
      lerr = .false.
    endif
  end subroutine parametersCheck

  subroutine checkBinaryFileParameters(filename)

    character(*), intent(in) :: filename

    character(len=5) :: crank
    character (:),allocatable,save :: datfile_name

    integer :: nx_f, ns_f
    double precision :: cv_f
    integer,allocatable,dimension(:) :: np_f
    double precision,allocatable,dimension(:) :: omega_p_f, qm_F
    double precision,allocatable,dimension(:) :: v_para_f, v_perp_f, v_d_f, pitch_deg_f
    double precision :: omega_c_f, theta_deg_f

    integer :: is

#ifndef _DEBUG
          write(crank,'(i5.5)') irank_mpi
#else
          crank = '00000'
#endif
          allocate(character(len_trim(filename)+10) :: datfile_name)
          datfile_name =  trim(filename) //'_'// crank // '.dat'

          open(99,file=datfile_name, form='unformatted', access='stream', status='old')
          read(99) nx_f, cv_f, ns_f

          if (nx_f /= NX) then
            print *, 'NX','file', nx_f,'program', NX ;goto 999
          endif
          if (cv_f /= cv) then
            print *, 'cv','file', cv_f,'program', cv ;goto 999
          endif
          if (ns_f /= NS) then
            print *, 'NS','file', ns_f,'program', NS ;goto 999
          endif

          allocate(np_f(NS))
          allocate(omega_p_f(NS))
          allocate(qm_f(NS))
          allocate(v_para_f(NS))
          allocate(v_perp_f(NS))
          allocate(v_d_f(NS))
          allocate(pitch_deg_f(NS))

          read(99) np_f, omega_p_f, qm_F, v_para_f, v_perp_f, v_d_f, pitch_deg_f

          do is = 1, NS
            if (np_f(is)        /= NP(is)       ) then
              print *, 'NP',       'file', np_f,       'program', NP        ;goto 999
            endif
            if (omega_p_f(is)   /= omega_p(is)  ) then
              print *, 'omega_p',  'file', omega_p_f,  'program', omega_p   ;goto 999
            endif
            if (qm_f(is)        /= QM(is)       ) then
              print *, 'QM',       'file', qm_f,       'program', QM        ;goto 999
            endif
            if (v_para_f(is)    /= v_para(is)   ) then
              print *, 'v_para',   'file', v_para_f,   'program', v_para    ;goto 999
            endif
            if (v_perp_f(is)    /= v_perp(is)   ) then
              print *, 'v_perp',   'file', v_perp_f,   'program', v_perp    ;goto 999
            endif
            if (v_d_f(is)       /= v_d(is)      ) then
              print *, 'v_d',      'file', v_d_f,      'program', v_d       ;goto 999
            endif
            if (pitch_deg_f(is) /= PITCH_DEG(is)) then
              print *, 'PITCH_DEG','file', pitch_deg_f,'program', PITCH_DEG ;goto 999
            endif
          end do

          deallocate(np_f)
          deallocate(omega_p_f)
          deallocate(qm_F)
          deallocate(v_para_f)
          deallocate(v_perp_f)
          deallocate(v_d_f)
          deallocate(pitch_deg_f)

          read(99) omega_c_f, theta_deg_f

          if (omega_c_f   /= omega_c  ) then
            print *, 'omega_c'  ,'file', omega_c_f,  'program', omega_c   ;goto 999
          endif
          if (theta_deg_f /= THETA_DEG) then
            print *, 'THETA_DEG','file', theta_deg_f,'program', THETA_DEG ;goto 999
          endif

          return

999       print *, "Error: Parameters does not correspond."
          close(99)
          stop
  end subroutine
end module
