module CalcForwardBackwardWaves
  use SimulationParameters,   only: NX
  use SupplementalParameters, only: XLEN
  use calcBfield,             only: bfieldBoundary
  use calcEfield,             only: efieldBoundary
  implicit none
  include 'fftw3.f'

  double precision,dimension(NX+8),save :: ey_fwd, ez_fwd, ey_bwd, ez_bwd
  !DEC$ATTRIBUTES ALIGN: 64:: ey_fwd, ez_fwd, ey_bwd, ez_bwd
  double precision,dimension(NX+8),save :: by_fwd, bz_fwd, by_bwd, bz_bwd
  !DEC$ATTRIBUTES ALIGN: 64:: by_fwd, bz_fwd, by_bwd, bz_bwd

  integer(8) :: ifwd_ey, ifwd_ez, ifwd_by, ifwd_bz
  integer(8) :: ibwd_ey_fwd, ibwd_ez_fwd, ibwd_ey_bwd, ibwd_ez_bwd
  integer(8) :: ibwd_by_fwd, ibwd_bz_fwd, ibwd_by_bwd, ibwd_bz_bwd

  ! Temporal storage of raw electromagnetic field
  double precision,dimension(NX+8) :: ey_1d, ez_1d, by_1d, bz_1d
  !DEC$ATTRIBUTES ALIGN: 64:: ey_1d, ez_1d, by_1d, bz_1d

  ! Wave number array
  complex(kind(0d0)),dimension(NX/2+1) :: ey_k, ez_k, by_k, bz_k
  ! Wave number of Electric field for forward and backward waves
  complex(kind(0d0)),dimension(NX/2+1) :: ey_k_fwd, ez_k_fwd, ey_k_bwd, ez_k_bwd
  ! Wave number of Magnetic field for forward and backward waves
  complex(kind(0d0)),dimension(NX/2+1) :: by_k_fwd, bz_k_fwd, by_k_bwd, bz_k_bwd
  ! Temporal storage of Electric field for forward and backward waves
  double precision,dimension(NX) :: ey_fwd_1d, ez_fwd_1d, ey_bwd_1d, ez_bwd_1d
  ! Temporal storage of Magnetic field for forward and backward waves
  double precision,dimension(NX) :: by_fwd_1d, bz_fwd_1d, by_bwd_1d, bz_bwd_1d

  !DEC$ATTRIBUTES ALIGN: 64:: ey_fwd_1d, ez_fwd_1d, ey_bwd_1d, ez_bwd_1d
  !DEC$ATTRIBUTES ALIGN: 64:: by_fwd_1d, bz_fwd_1d, by_bwd_1d, bz_bwd_1d

  public :: ey_fwd, ez_fwd, ey_bwd, ez_bwd, by_fwd, bz_fwd, by_bwd, bz_bwd
  public initForwardBackwardWaves, destroyForwardBackwardWaves, forwardBackwardWaves
contains
  subroutine initForwardBackwardWaves
    call dfftw_plan_dft_r2c_1d(ifwd_ey, NX, ey_1d, ey_k, FFTW_PATIENT)
    call dfftw_plan_dft_r2c_1d(ifwd_ez, NX, ez_1d, ez_k, FFTW_PATIENT)
    call dfftw_plan_dft_r2c_1d(ifwd_by, NX, by_1d, by_k, FFTW_PATIENT)
    call dfftw_plan_dft_r2c_1d(ifwd_bz, NX, bz_1d, bz_k, FFTW_PATIENT)

    call dfftw_plan_dft_c2r_1d(ibwd_ey_fwd, NX, ey_k_fwd, ey_fwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_ez_fwd, NX, ez_k_fwd, ez_fwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_ey_bwd, NX, ey_k_bwd, ey_bwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_ez_bwd, NX, ez_k_bwd, ez_bwd_1d, FFTW_PATIENT)

    call dfftw_plan_dft_c2r_1d(ibwd_by_fwd, NX, by_k_fwd, by_fwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_bz_fwd, NX, bz_k_fwd, bz_fwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_by_bwd, NX, by_k_bwd, by_bwd_1d, FFTW_PATIENT)
    call dfftw_plan_dft_c2r_1d(ibwd_bz_bwd, NX, bz_k_bwd, bz_bwd_1d, FFTW_PATIENT)
  end subroutine initForwardBackwardWaves

  subroutine destroyForwardBackwardWaves
    call dfftw_destroy_plan(ifwd_ey)
    call dfftw_destroy_plan(ifwd_ez)
    call dfftw_destroy_plan(ifwd_by)
    call dfftw_destroy_plan(ifwd_bz)

    call dfftw_destroy_plan(ibwd_ey_fwd)
    call dfftw_destroy_plan(ibwd_ez_fwd)
    call dfftw_destroy_plan(ibwd_ey_bwd)
    call dfftw_destroy_plan(ibwd_ey_bwd)
    call dfftw_destroy_plan(ibwd_by_fwd)
    call dfftw_destroy_plan(ibwd_bz_fwd)
    call dfftw_destroy_plan(ibwd_by_bwd)
    call dfftw_destroy_plan(ibwd_by_bwd)
  end subroutine

  subroutine forwardBackwardWaves(ey,ez,by,bz)
    double precision,dimension(:),intent(in) :: ey,ez
    double precision,dimension(:),intent(in) :: by,bz

    double precision,parameter :: X_INV = 1d0/XLEN 

    ey_1d(1:NX) = ey(2:NX+1)
    ez_1d(1:NX) = ez(2:NX+1)
    by_1d(1:NX) = by(2:NX+1)
    bz_1d(1:NX) = bz(2:NX+1)

    call dfftw_execute_dft_r2c(ifwd_ey, ey_1d, ey_k)
    call dfftw_execute_dft_r2c(ifwd_ez, ez_1d, ez_k)
    call dfftw_execute_dft_r2c(ifwd_by, by_1d, by_k)
    call dfftw_execute_dft_r2c(ifwd_bz, bz_1d, bz_k)
    
    ez_k_fwd = 0.5d0 * (ey_k * (0,  1d0) + ez_k) ! (0, 1d0) stands for imaginary unit 1i
    ez_k_bwd = 0.5d0 * (ey_k * (0, -1d0) + ez_k)
    ey_k_fwd = ez_k_fwd * (0, -1d0)
    ey_k_bwd = ez_k_bwd * (0, 1d0 )
    
    bz_k_fwd = 0.5d0 * (by_k * (0,  1d0) + bz_k)
    bz_k_bwd = 0.5d0 * (by_k * (0, -1d0) + bz_k)
    by_k_fwd = bz_k_fwd * (0, -1d0)
    by_k_bwd = bz_k_bwd * (0, 1d0 )
    
    call dfftw_execute_dft_c2r(ibwd_ey_fwd, ey_k_fwd, ey_fwd_1d)
    call dfftw_execute_dft_c2r(ibwd_ez_fwd, ez_k_fwd, ez_fwd_1d)
    call dfftw_execute_dft_c2r(ibwd_ey_bwd, ey_k_bwd, ey_bwd_1d)
    call dfftw_execute_dft_c2r(ibwd_ez_bwd, ez_k_bwd, ez_bwd_1d)
    
    call dfftw_execute_dft_c2r(ibwd_by_fwd, by_k_fwd, by_fwd_1d)
    call dfftw_execute_dft_c2r(ibwd_bz_fwd, bz_k_fwd, bz_fwd_1d)
    call dfftw_execute_dft_c2r(ibwd_by_bwd, by_k_bwd, by_bwd_1d)
    call dfftw_execute_dft_c2r(ibwd_bz_bwd, bz_k_bwd, bz_bwd_1d)

    ey_fwd(2:NX+1) = ey_fwd_1d * X_INV
    ez_fwd(2:NX+1) = ez_fwd_1d * X_INV
    ey_bwd(2:NX+1) = ey_bwd_1d * X_INV
    ez_bwd(2:NX+1) = ez_bwd_1d * X_INV

    by_fwd(2:NX+1) = by_fwd_1d * X_INV
    bz_fwd(2:NX+1) = bz_fwd_1d * X_INV
    by_bwd(2:NX+1) = by_bwd_1d * X_INV
    bz_bwd(2:NX+1) = bz_bwd_1d * X_INV

    ! Boundary conditions
    call efieldBoundary(ey_fwd, ez_fwd)
    call efieldBoundary(ey_bwd, ez_bwd)
    call bfieldBoundary(by_fwd, bz_fwd)
    call bfieldBoundary(by_bwd, bz_bwd)
  end subroutine forwardBackwardWaves
end module CalcForwardBackwardWaves
