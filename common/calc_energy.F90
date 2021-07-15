module CalcEnergy
  use SimulationParameters,   only: NTIME, NX, NS, NP, cv, ext_mirror_a
  use SupplementalParameters, only: mass, &
                                    ix_sim_first, ix_sim_last, x_sim_first, x_sim_last, &
                                    x_sim_size      
  use InitFieldParticle, only: bx0, by0
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
  double precision, dimension(:,:),allocatable,save :: energy_kinetic
  double precision, dimension(:,:),allocatable,save :: energy_kinetic_para
  double precision, dimension(:,:),allocatable,save :: energy_kinetic_perp
  double precision, dimension(:,:),allocatable,save :: energy_e_field
  double precision, dimension(:,:),allocatable,save :: energy_b_field
  double precision, dimension(:)  ,allocatable,save :: energy_total
contains
  subroutine initEnergy(isize)
    integer,intent(in) :: isize
    allocate(energy_kinetic(isize,NS))
    allocate(energy_kinetic_para(isize,NS))
    allocate(energy_kinetic_perp(isize,NS))
    allocate(energy_e_field(isize,3))
    allocate(energy_b_field(isize,3))
    allocate(energy_total(isize))

    energy_kinetic = 0.0d0
    energy_kinetic_para = 0.0d0
    energy_kinetic_perp = 0.0d0
    energy_e_field = 0.0d0
    energy_b_field = 0.0d0
    energy_total = 0.0d0
  end subroutine initEnergy

  subroutine quitEnergy
    deallocate(energy_kinetic)
    deallocate(energy_kinetic_para)
    deallocate(energy_kinetic_perp)
    deallocate(energy_e_field)
    deallocate(energy_b_field)
    deallocate(energy_total)
  end subroutine quitEnergy

  subroutine energy(it,ex,ey,ez,by,bz,x,vx,vy,vz,iq)
    integer,intent(in) :: it
    double precision, dimension(:),intent(in) :: ex, ey, ez
    double precision, dimension(:),intent(in) :: by, bz
    double precision, dimension(:),intent(in) :: x
    double precision, dimension(:),intent(in) :: vx, vy, vz
    integer(1),       dimension(:),intent(in) :: iq

    double precision, dimension(NS) :: e_k_thread, e_k_tmp
    double precision, dimension(NS) :: e_k_thread_para, e_k_tmp_para
    double precision, dimension(NS) :: e_k_tmp_perp
    double precision, dimension(3)  :: e_e_tmp, e_b_tmp

    double precision, parameter :: NXC = dble(NX/2)  
    double precision :: v_perp_eq_2, v_para_eq_2
    double precision :: tmp

    double precision :: ek_stational

    integer :: is, ip
    integer :: ix

    integer :: ifirst, ilast

    ! Kinetic energy
    !$omp parallel private(is, ifirst, ilast, ek_stational, v_perp_eq_2, tmp, v_para_eq_2),reduction(+:e_k_thread, e_k_thread_para)
    ifirst = 0
    ilast  = 0
    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)
      ek_stational = mass(is)*cv**2
      !$omp do simd linear(ip:1)
      do ip = ifirst, ilast
         if( x(ip) < x_sim_first .or. x(ip) > x_sim_last ) cycle
          e_k_thread(is)      = e_k_thread(is)      + (cv/dsqrt(cv**2 - vx(ip)**2 - vy(ip)**2 - vz(ip)**2)-1.0d0)* ek_stational * dble(iq(ip))
          ! get equatorial parallel velocity
          v_perp_eq_2 = vy(ip)**2 + vz(ip)**2
          tmp =   ext_mirror_a * (x(ip)-NXC)**2 / (1d0 + ext_mirror_a * (x(ip)-NXC)**2)
          v_para_eq_2 = vx(ip)**2 + v_perp_eq_2 * tmp
          e_k_thread_para(is) = e_k_thread_para(is) + (cv/dsqrt(cv**2 - v_para_eq_2)-1.0d0)                      * ek_stational * dble(iq(ip))
      end do
      !$omp end do simd
    end do
    !$omp end parallel

#ifdef _DEBUG
    e_k_tmp      = e_k_thread
    e_k_tmp_para = e_k_thread_para
#else
    call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
    call MPI_REDUCE(e_k_thread     ,e_k_tmp     ,NS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_REDUCE(e_k_thread_para,e_k_tmp_para,NS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
#endif

  e_k_tmp_perp = e_k_tmp - e_k_tmp_para
  
    ! electric field energy & magnetic field energy
    !$omp parallel do private(ix),reduction(+:e_e_tmp),reduction(+:e_b_tmp)
    do ix = ix_sim_last, ix_sim_last
        e_e_tmp(1) = e_e_tmp(1) + ex(ix)**2
        e_e_tmp(2) = e_e_tmp(2) + ey(ix)**2
        e_e_tmp(3) = e_e_tmp(3) + ez(ix)**2

        e_b_tmp(2) = e_b_tmp(2) + (by(ix)-by0)**2
        e_b_tmp(3) = e_b_tmp(3) + bz(ix)**2
    end do
    !$omp end parallel do

    e_k_tmp = e_k_tmp                  / dble(x_sim_size)
    e_k_tmp_para = e_k_tmp_para        / dble(x_sim_size)
    e_k_tmp_perp = e_k_tmp_perp        / dble(x_sim_size)
    e_e_tmp = 0.5d0*           e_e_tmp / dble(x_sim_size)
    e_b_tmp = 0.5d0* (cv**2) * e_b_tmp / dble(x_sim_size)

    energy_kinetic(it,:) = e_k_tmp
    energy_kinetic_para(it,:) = e_k_tmp_para
    energy_kinetic_perp(it,:) = e_k_tmp_perp
    energy_e_field(it,:) = e_e_tmp
    energy_b_field(it,:) = e_b_tmp
    energy_total(it) = sum(e_k_tmp) + sum(e_e_tmp) + sum(e_b_tmp)
  end subroutine
end module
