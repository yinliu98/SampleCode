module CalcWaveNumber
  use SimulationParameters, only: NX, iX_WAVE_NUMBER_DIV
  use SupplementalParameters, only: XLEN
  use CalcForwardBackwardWaves, only:  ey_fwd, ez_fwd, ey_bwd, ez_bwd, &
                                       by_fwd, bz_fwd, by_bwd, bz_bwd
  use OutputHDF5, only: outputKx

! DEBUG mode unables MPI (Compiler option)------------------------------------
! -D_DEBUG    => debug mode ON
! Not defined => debug mode OFF
!-----------------------------------------------------------------------------
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif

  implicit none
  include 'fftw3.f'
  public waveNumber
contains
  subroutine waveNumber(it,ex)
    integer,intent(in) :: it
    double precision,dimension(:),intent(in) :: ex

    double precision,dimension(NX+2) :: var
    !DEC$ATTRIBUTES ALIGN: 64:: var
    integer,parameter :: KX_MAX = NX/iX_WAVE_NUMBER_DIV + 2

    double precision,dimension(KX_MAX) :: k
    !DEC$ATTRIBUTES ALIGN: 64:: k
    double precision,parameter :: XLEN_INV = 1d0/XLEN

    integer(8) :: iforward
    integer :: ix ! Iterator to X,Y direction

    call dfftw_plan_dft_r2c_1d(iforward , NX, ex(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, ex(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Ex', k, .true.)

    call dfftw_plan_dft_r2c_1d(iforward , NX, ey_fwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, ey_fwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Ey_fwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, ey_bwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, ey_bwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Ey_bwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, ez_fwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, ez_fwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Ez_fwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, ez_bwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, ez_bwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do

    call outputKx(it, 'Ez_bwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, by_fwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, by_fwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'By_fwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, by_bwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, by_bwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'By_bwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, bz_fwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, bz_fwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Bz_fwd', k, .false.)

    call dfftw_plan_dft_r2c_1d(iforward, NX, bz_bwd(2:NX+1), var, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iforward, bz_bwd(2:NX+1), var)
    call dfftw_destroy_plan(iforward)

    do ix = 1, KX_MAX
       k(ix) = var(ix)*XLEN_INV
    end do
    call outputKx(it, 'Bz_bwd', k, .false.)

#ifndef _DEBUG
    call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
  end subroutine waveNumber
end module
