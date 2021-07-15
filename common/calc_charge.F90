module CalcCharge
  use SimulationParameters, only: NX, NS, NP
  use SupplementalParameters, only: rho0, q
  ! DEBUG mode unables MPI (Compiler option)------------------------------------
  ! -D_DEBUG    => debug mode ON
  ! Not defined => debug mode OFF
  !-----------------------------------------------------------------------------
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif
  !$ use omp_lib
  implicit none
contains
  subroutine charge(x,iq,rho)
    double precision,dimension(:)  ,intent(in)  :: x
    integer(1),      dimension(:)  ,intent(in)  :: iq
    double precision,dimension(:),intent(inout) :: rho
#ifndef _DEBUG
    ! Particle decomposition =====================================================
    integer         ,parameter            :: COM_LEN_MPI = NX+8
    !=============================================================================
#endif

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2
    integer :: ifirst    ! first index of particle
    integer :: ilast     ! last  index of particle

    integer :: is      ! iteretor for Species
    integer :: ip      ! iteretor for Particle

    rho = 0.0d0

    !$omp parallel private(is, ifirst, ilast),reduction(+:rho)
    ifirst = 0
    ilast  = 0
    do is = 1,NS
      ifirst  = ilast + 1
      ilast   = ilast + NP(is)
      !$omp do private(x_p,ixp,sf1,sf2)
      do ip = ifirst,ilast
        x_p = x(ip) + 2.0d0
        ixp = floor(x_p)
 
        sf2 = x_p - dble(ixp)
        sf1 = 1 - sf2

        rho(ixp  ) = rho(ixp  ) + q(is)*sf1*dble(iq(ip))
        rho(ixp+1) = rho(ixp+1) + q(is)*sf2*dble(iq(ip))
      end do
      !$omp end do nowait
    end do
    !$omp end parallel

#ifndef _DEBUG
    ! Particle decomposition ========================================================================
    call MPI_ALLREDUCE(MPI_IN_PLACE,rho,com_len_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    !================================================================================================
#endif

    rho = rho + rho0
    rho(2) = rho(2) + rho(NX+2) - rho0(NX+2)
    rho(1) = rho(NX+1)
    rho(NX+2) = rho(2)
  end subroutine charge
end module
