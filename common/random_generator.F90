module RandomGenerator
  use SimulationParameters, only: ISEED
  use Slps,                 only: PI
  use mt_stream, mtNew  => new     ,&
                 mtInit => init    ,&
                 mtDelete => delete
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif
  !$ use omp_lib
  implicit none
  type(mt_state), save :: mts
  !$omp threadprivate(mts)
  private
  public mts
  public :: initMTStream, finMTstream, randN, subRandN, randUi
contains
  subroutine initMTStream
    !MT generation
    type(mt_state) :: mts_init

    integer      :: iomp_num_threads
    integer,save :: irank_omp
    !$omp threadprivate(irank_omp)

#ifdef _DEBUG
    integer,parameter :: irank_mpi = 0
#endif

    !Initialize for without openMP
    iomp_num_threads = 1
    irank_omp = 0

    !Mersenne Twister init------------------------------------------------------
    call set_mt19937
    call mtNew(mts_init)
    call mtInit(mts_init,ISEED)

    !$omp parallel shared(iomp_num_threads)
    !$ iomp_num_threads = omp_get_num_threads()
    !$ irank_omp = omp_get_thread_num()
    call create_stream(mts_init,mts,iomp_num_threads*irank_mpi+irank_omp+1)
    !$omp end parallel
  end subroutine initMTstream

  subroutine finMTstream
    !$omp parallel
    call mtDelete(mts)
    !$omp end parallel
  end subroutine finMTstream

  double precision function randN(mts,mu,sigma)
    type(mt_state)  ,intent(inout) :: mts
    double precision,intent(in),optional :: mu    ! Average
    double precision,intent(in),optional :: sigma ! dispersion

    double precision :: tmp1, tmp2, tmp_mu, tmp_sigma

    double precision,save :: rand_cache
    logical,save :: is_cos

    if (present(mu)) then
       tmp_mu = mu
    else
       tmp_mu = 0.0d0
    endif

    if (present(sigma)) then
       tmp_sigma = sigma
    else
       tmp_sigma = 1.0d0
    endif
   
    if (is_cos) then
      tmp1 = genrand_double3(mts)
      tmp2 = genrand_double3(mts)
      randN      = tmp_sigma * dsqrt(-2.0d0*dlog(tmp1)) * dcos(2.0d0*PI*tmp2) + tmp_mu
      rand_cache = dsqrt(-2.0d0*dlog(tmp1)) * dsin(2.0d0*PI*tmp2)
      is_cos = .false.
    else 
      randN = tmp_sigma * rand_cache + tmp_mu
      is_cos = .true.
    endif
  end function

  double precision function subRandN(mts,beta,rho)
    type(mt_state)  ,intent(inout) :: mts
    double precision,intent(in)    :: beta  ! diviation ratio [0,1]
    double precision,intent(in)    :: rho

    double precision,parameter :: range = 10d0

    double precision :: x, f, y

    do
       x = genrand_double3(mts) * range
       f = (1d0/(1d0-rho*beta))*x*(dexp(-x**2/2d0)-rho*dexp(-x**2/2d0/beta))
       y = genrand_double3(mts)

       if ( y < f ) then
          subRandN = x
          return
       else
          continue
       endif
    end do
  end function subRandN

  double precision function randUi(mts)
    type(mt_state)  ,intent(inout) :: mts
    randUi= dble(floor(genrand_double2(mts) + 0.5d0))
  end function
end module RandomGenerator
