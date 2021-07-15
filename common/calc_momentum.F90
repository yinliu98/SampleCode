module CalcMagneticMomentum
  use SimulationParameters,   only: NTIME, NX, NS, NP, cv, &
                                    EXT_MIRROR, ext_mirror_a
  use SupplementalParameters, only: mass, x_sim_first, x_sim_last, x_sim_size
  use InitFieldParticle,      only: bx0
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
  double precision, dimension(:,:),allocatable,save :: magnetic_momentum_total

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  private
  !Parameters/Variables
  public :: magnetic_momentum_total
  !Methods
  public initMagneticMomentum, quitMagneticMomentum, magneticMomentum
contains
  subroutine initMagneticMomentum(isize)
    integer,intent(in) :: isize
    allocate(magnetic_momentum_total(isize,NS))

    magnetic_momentum_total = 0.0d0
  end subroutine initMagneticMomentum

  subroutine quitMagneticMomentum
    deallocate(magnetic_momentum_total)
  end subroutine quitMagneticMomentum

  subroutine magneticMomentum(it,by,bz,x,vx,vy,vz,iq)
    integer,intent(in) :: it
    double precision, dimension(:),intent(in) :: by, bz
    double precision, dimension(:)  ,intent(in) :: x
    double precision, dimension(:)  ,intent(in) :: vx, vy, vz
    integer(1),       dimension(:)  ,intent(in) :: iq

    double precision, dimension(NS) :: mu_thread

    double precision,dimension(NX+8) :: by_avg
    !DEC$ATTRIBUTES ALIGN: 64:: by_avg

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2
    double precision :: p_bx, p_by, p_bz
    double precision :: p_vx, p_vy, p_vz

    double precision :: b0_abs2, v_perp2_b0_abs, gamma
    integer :: is, ip, ix
    integer :: ifirst,ilast

    double precision, parameter :: NXC = dble(NX/2)

    by_avg = 0.0d0

    !==============================================================================================
    ! Field cancelation to prevent from self-force oscillation
    !==============================================================================================
    do ix = 2, NX+1
       by_avg(ix) = (by(ix+1) + by(ix))*0.5d0
    end do
    by_avg(1) = by_avg(NX+1)
    mu_thread = 0d0

    !$omp parallel private(is, ifirst, ilast)
    ifirst = 0
    ilast  = 0
    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)
      !$omp do private(ixp,sf1,sf2),&
      !$omp private(p_bx,p_by,p_bz, p_vx, p_vy, p_vz),&
      !$omp private(b0_abs2, v_perp2_b0_abs,gamma),&
      !$omp reduction(+:mu_thread)
      do ip = ifirst, ilast
         if( x(ip) < x_sim_first .or. x(ip) > x_sim_last ) cycle

          ! Ambient magnetic field
          if (EXT_MIRROR == .true. ) then
            p_bx=bx0*(1d0 + ext_mirror_a * (x(ip)-NXC)**2)
          else
            p_bx=bx0
          endif

          x_p = x(ip) + 1.5d0
          ixp = floor(x_p)
          sf2 = x_p - dble(ixp)
          sf1 = 1d0 - sf2

          p_by=sf1*by_avg(ixp) + sf2*by_avg(ixp+1)
          p_bz=sf1*bz(ixp)     + sf2*bz(ixp+1)

          p_vx = vx(ip)
          p_vy = vy(ip)
          p_vz = vz(ip)

          b0_abs2 = 0.5d0/(p_bx ** 2 + p_by ** 2 + p_bz ** 2)
          v_perp2_b0_abs = (p_vy*p_bz - p_vz*p_by)**2 + (p_vz*p_bx-p_vx*p_bz)**2 + (p_vx*p_by - p_vy*p_bx)**2
          gamma = cv/dsqrt(cv**2 - p_vx**2 - p_vy**2 - p_vz**2)
          mu_thread(is) = mu_thread(is) + mass(is) * gamma * v_perp2_b0_abs * b0_abs2 * dble(iq(ip))
      end do
      !$omp end do nowait
    end do
    !$omp end parallel

#ifndef _DEBUG
      ! Particle decomposition =======================================================================
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, mu_thread,NS,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(mu_thread,    mu_thread,NS,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
      endif
      !==============================================================================================
#endif

    mu_thread = mu_thread / dble(x_sim_size)
    magnetic_momentum_total(it,:) = mu_thread
  end subroutine magneticMomentum
end module CalcMagneticMomentum
