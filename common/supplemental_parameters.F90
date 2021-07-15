module SupplementalParameters
  use SimulationParameters, only: NTIME, DX, DT, NX, cv,&
                                  NS, NP, omega_p, QM, v_para, v_perp, v_d, PITCH_DEG,&
                                  omega_c, &
                                  ENABLE_INJECTION, IS_INJECTION, &
                                  ext_mirror_a, &
                                  THETA_DEG, LOSSCORN_ANGLE_DEG, &
                                  ctime, cwidth_begin, cwidth_end, &
                                  ext_j_amp, ext_omega_a, ext_omega_b0, &
                                  LOCALIZED_CHARGE_AVG, & 
                                  FIELD_ABSORB_BOUNDARY, DAMPING_WAVE_FREQUENCY, DAMPING_INTEGRAL_K, DAMPING_REGION_ND, ISKIP
  use Slps, only: PI
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif
  !$ use omp_lib
  implicit none
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! sub parameters
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !Renormalize Parameters-------------------------------------------------------
  double precision,parameter :: RN_X            = DX                  ! Distance
  double precision,parameter :: RN_T            = DT      / 2.0d0     ! Time
  double precision,parameter :: RN_V            = RN_X    / RN_T      ! Velocity
  double precision,parameter :: RN_E_FIELD      = RN_X    / (RN_T**2) ! Electric field
  double precision,parameter :: RN_POTENTIAL    = RN_X**2 / (RN_T**2) ! Electric potential
  double precision,parameter :: RN_B_FIELD      = 1.0d0   / RN_T      ! Magnetic field
  double precision,parameter :: RN_CURRENT_DENS = RN_X    / (RN_T**3) ! Current density
  double precision,parameter :: RN_CHARGE_DENS  = 1.0d0   / (RN_T**2) ! Charge density
  double precision,parameter :: RN_CHARGE       = RN_X**2 / (RN_T**2) ! Electric charge
  double precision,parameter :: RN_MASS         = RN_X**2 / (RN_T**2) ! Mass

  ! Parameters------------------------------------------------------------------
  double precision,parameter :: COSTH = dcos(dble(THETA_DEG*PI/180d0))
  double precision,parameter :: SINTH = dsin(dble(THETA_DEG*PI/180d0))
  double precision,parameter :: LOSSCORN_ANGLE_TAN = dtan(dble(LOSSCORN_ANGLE_DEG*PI/180d0))

  integer,parameter          :: NP_TOTAL = sum(NP(1:NS))
  double precision,parameter :: XLEN = dble(NX)
  double precision,save      :: tcs, csq
  double precision,save,dimension(ns)   :: q
  double precision,save,dimension(ns)   :: mass
  double precision,save,dimension(NX+8) :: rho0
  !DEC$ATTRIBUTES ALIGN: 64:: rho0

  double precision,save                 :: damping_rd,    damping_rr
  double precision,save,dimension(NX+8) :: damping_fm_rd, damping_fm_rr
  !DEC$ATTRIBUTES ALIGN: 64:: damping_fm_rd, damping_fm_rr
  integer,save :: ix_sim_first, ix_sim_last
  double precision,save :: x_sim_first, x_sim_last
  double precision,save :: x_sim_size

  !Integral steps for WPIA and resonant current
  integer,save :: IT_INTEGRAL_STEPS

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  private
  !Parameters/Variables
  public :: RN_X, RN_T, RN_V, RN_E_FIELD, RN_POTENTIAL
  public :: RN_B_FIELD, RN_CURRENT_DENS, RN_CHARGE_DENS, RN_CHARGE, RN_MASS
  public :: COSTH, SINTH, LOSSCORN_ANGLE_TAN
  public :: NP_TOTAL
  public :: XLEN
  public :: tcs, csq
  public :: q, mass, rho0
  public :: damping_rd, damping_rr, damping_fm_rd, damping_fm_rr
  public :: ix_sim_first, ix_sim_last, x_sim_first, x_sim_last, x_sim_size
  public :: IT_INTEGRAL_STEPS
  !Methods
  public initSupplementalParameters, setCounterCharge
contains
  subroutine initSupplementalParameters
    call dampingFactor
    call renormalize
    ! Parameters from speed of light
    tcs = 2.0d0 * cv**2
    csq = cv**2
    IT_INTEGRAL_STEPS = min(ISKIP, floor(PI/abs(omega_c)))
    call setCharge
  end subroutine

  subroutine renormalize
    !! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !! This program OVERRIDES simulation parameters as below.
    !! omega_p,omega_c,omega,cv,v_para,v_perp,v_d
    !! DO NOT excute TWICE in single code.
    !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    omega_p = omega_p * RN_T
    omega_c = omega_c * RN_T
    ext_omega_a  = ext_omega_a  * RN_T
    ext_omega_b0 = ext_omega_b0 * RN_T
    cv      = cv      / RN_V
    v_para  = v_para  / RN_V
    v_perp  = v_perp  / RN_V
    v_d     = v_d     / RN_V

    ctime        = ctime / abs(omega_c) * 0.5d0
    cwidth_begin = cwidth_begin / abs(omega_c) * 0.5d0
    cwidth_end   = cwidth_end / abs(omega_c) * 0.5d0
    ext_mirror_a = ext_mirror_a * RN_X ** 2
    ext_j_amp    = ext_j_amp / RN_CURRENT_DENS
  end subroutine renormalize

  subroutine setCounterCharge(rho)
    double precision,intent(inout),dimension(:) :: rho
    rho0 = - rho
  end subroutine setCounterCharge

  !===============================================================================
  ! Private Subroutines for variableListSetVariables
  !===============================================================================

  subroutine setCharge
    double precision,dimension(ns) :: tmp
    double precision :: num_process

    integer :: is

#ifdef _DEBUG
      num_process = 1d0
#else
      num_process = dble(iprocess)
#endif
    tmp  = 1.0d0/(dble(NP(1:NS))*dble(QM(1:NS))*dble(num_process))
    q    = dble(XLEN)*omega_p**2*tmp
    mass = q/dble(QM(1:NS))

    if (LOCALIZED_CHARGE_AVG == .true.) then
      rho0 = 0d0
    else
      if (ENABLE_INJECTION) then
         rho0 = 0d0
         do is = 1, NS
            if (IS_INJECTION(is) == .true.) cycle
            rho0 = rho0 - q(is)*dble(NP(is))
         end do
         rho0 = rho0 * num_process / XLEN
      else
         rho0   = - sum(q(1:NS)*dble(NP(1:NS))*num_process)/dble(NX)
      endif
      rho0(1) = rho0(NX+1)
      rho0(NX+2) = rho0(2)
    endif
  end subroutine setCharge

  subroutine dampingFactor
    double precision :: ld
    integer :: ix

    double precision :: rr_condition
    double precision :: first_coeff, integral_max, integral_fraction, integral_term

    damping_fm_rd = 1d0
    damping_fm_rr = 1d0

    if (FIELD_ABSORB_BOUNDARY == .false.) then
       ix_sim_first = 2
       ix_sim_last  = NX + 1
       x_sim_first  = 0d0
       x_sim_last   = dble(NX)
       x_sim_size   = dble(NX)
       return
    endif
    
    ! Calculate rr
    rr_condition = cv * PI / DAMPING_WAVE_FREQUENCY * ((2*dble(DAMPING_REGION_ND) -1)/dble(DAMPING_REGION_ND)**2)
    if (DX >= rr_condition) then
       damping_rr = dsqrt(1d0 - DAMPING_WAVE_FREQUENCY*DX/(cv*PI)) * dble(DAMPING_REGION_ND - 1) / dble(DAMPING_REGION_ND)
    else
       damping_rr = 1d0
    endif
    
    ! Calculate rd
    first_coeff = dble(DAMPING_REGION_ND) / ( cv * DT / DX)
    integral_max = dble(DAMPING_REGION_ND - 1) / dble(DAMPING_REGION_ND)
    integral_fraction = (1d0 + damping_rr * integral_max)/(1d0 - damping_rr * integral_max)
    integral_term = dlog(integral_fraction)/(2d0*damping_rr) - integral_max
    
    damping_rd = dsqrt(DAMPING_INTEGRAL_K /(first_coeff * integral_term)) * damping_rr

    ! Set Analysis boundaries
    ix_sim_first = DAMPING_REGION_ND + 1
    ix_sim_last  = NX + 2 - DAMPING_REGION_ND
    x_sim_first  = DAMPING_REGION_ND
    x_sim_last   = NX - DAMPING_REGION_ND
    x_sim_size   = NX - 2*DAMPING_REGION_ND

    ld   = dble(DAMPING_REGION_ND)*dx

    ! x only damping region
    do ix = 1, ix_sim_first - 1
      damping_fm_rd(ix) = 1d0 - (damping_rd*(dble(ix)*dx-ld)/ld)**2
      damping_fm_rr(ix) = 1d0 - (damping_rr*(dble(ix)*dx-ld)/ld)**2
    end do

    do ix = ix_sim_last + 1, NX+2
      damping_fm_rd(ix) = 1d0 - (damping_rd*(dble(ix-NX-3)*dx+ld)/ld)**2
      damping_fm_rr(ix) = 1d0 - (damping_rr*(dble(ix-NX-3)*dx+ld)/ld)**2
    enddo
  end subroutine
end module SupplementalParameters
