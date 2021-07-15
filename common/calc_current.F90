module CalcCurrent
  use Slps, only: PI
  use SimulationParameters, only: NX, NS, NP, &
                                  HAS_ELECTROSTATIC_COMPONENT, &
                                  EXT_CURRENT, &
                                  ctime, cwidth_begin, cwidth_end, iANTENNA_OFFSET, &
                                  ext_omega_a, ext_omega_b0, ext_j_amp
  use SupplementalParameters, only: q, rn_t, rn_current_dens, XLEN
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
  public current
contains
  subroutine current(it,x,vx,vy,vz,iq,jx,jy,jz)
    integer,intent(in) :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1),      dimension(:),intent(in)    :: iq
    double precision,dimension(:),intent(inout) :: jx, jy, jz

#ifndef _DEBUG
    ! Particle decomposition =====================================================
    integer,parameter :: com_len_mpi = NX+8
    !=============================================================================
#endif

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2
    !Current density from every single particle
    double precision :: p_jy, p_jz
    !Average current density
    double precision :: jx_avg,  jy_avg,  jz_avg

    ! Charge cancelletion
    double precision :: qhs, abs_vx, x1, x2
    integer :: i1, i2

    double precision,parameter :: XLEN_INV = 1d0/XLEN

    integer :: ifirst    ! first index of particle
    integer :: ilast     ! last  index of particle

    integer :: is      ! iteretor for Species
    integer :: ip      ! iteretor for Particle
    integer :: ix

    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0

    !$omp parallel private(is, ifirst, ilast) reduction(+:jx), reduction(+:jy), reduction(+:jz)
    ifirst = 0
    ilast  = 0
    do is = 1, NS
      ifirst  = ilast + 1
      ilast   = ilast + NP(is)
      !$omp do private(x_p,ixp,sf1,sf2,p_jy,p_jz,qhs,abs_vx,x1,x2,i1,i2)
      do ip = ifirst, ilast
         x_p = x(ip) + 1.5d0
         ixp = floor(x_p)
         
         sf2 = x_p - dble(ixp)
         sf1 = 1d0 - sf2
         
         p_jy = q(is)*vy(ip)*dble(iq(ip))
         p_jz = q(is)*vz(ip)*dble(iq(ip))

         jy(ixp  ) = jy(ixp  ) + p_jy*sf1
         jy(ixp+1) = jy(ixp+1) + p_jy*sf2

         jz(ixp  ) = jz(ixp  ) + p_jz*sf1
         jz(ixp+1) = jz(ixp+1) + p_jz*sf2

         if (HAS_ELECTROSTATIC_COMPONENT == .true.) then
            qhs = q(is) * 0.5d0 * dsign(1d0, vx(ip))
            abs_vx = dabs(vx(ip))

            x1 = x_p + 0.5d0 - abs_vx 
            x2 = x_p + 0.5d0 + abs_vx 
            i1 = floor(x1)
            i2 = floor(x2)

            jx(i1) = jx(i1) + (i2 - x1) * qhs * dble(iq(ip))
            jx(i2) = jx(i2) + (x2 - i2) * qhs * dble(iq(ip))
         endif
      end do
      !$omp end do nowait
    end do
    !$omp end parallel

#ifndef _DEBUG
    ! Particle decomposition =======================================================================
    call MPI_ALLREDUCE(MPI_IN_PLACE,jx,com_len_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    call MPI_ALLREDUCE(MPI_IN_PLACE,jy,com_len_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    call MPI_ALLREDUCE(MPI_IN_PLACE,jz,com_len_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    !==============================================================================================
#endif

    ! For X direction 
    ! periodic boundary conditions-----------------------------------------------------------------
    jx(NX+1) = jx(1) + jx(NX+1) 
    jx(2)    = jx(2) + jx(NX+2)

    jy(NX+1) = jy(1) + jy(NX+1) 
    jy(2)    = jy(2) + jy(NX+2)
    jy(1)    = jy(NX+1)

    jz(NX+1) = jz(1) + jz(NX+1) 
    jz(2)    = jz(2) + jz(NX+2)

    ! Relocation of current -----------------------------------------------------------------------

    !DIR$ SIMD
    do ix = NX+1, 2, -1 !Decrement
       jy(ix) = (jy(ix) + jy(ix-1))*0.5d0
    end do

    ! Cancelation of current-----------------------------------------------------------------------
    jx_avg = sum(jx(2:NX+1))*XLEN_INV
    jy_avg = sum(jy(2:NX+1))*XLEN_INV
    jz_avg = sum(jz(2:NX+1))*XLEN_INV

    do ix = 2,NX+1
       jx(ix) = jx(ix) - jx_avg
       jy(ix) = jy(ix) - jy_avg
       jz(ix) = jz(ix) - jz_avg
    end do

    if (EXT_CURRENT == .false.) return
    call extCurrent_phased_array(it,jy,jz)
  end subroutine current

  subroutine extCurrent_phased_array(it,jy,jz)
    integer,intent(in)          :: it
    double precision,intent(inout),dimension(:) :: jy,jz

    double precision :: j_amp, omega
    double precision,dimension(nx+8) :: jy_ext, jz_ext
    !DEC$ATTRIBUTES ALIGN: 64:: jy_ext, jz_ext

    integer,parameter :: NXC = NX/2+1

    call amplitude(it,j_amp)
    call frequency(it,omega)

    jy_ext = 0.0d0
    jz_ext = 0.0d0

    jy_ext(NXC+iANTENNA_OFFSET)   = j_amp * dcos(2.0d0*omega*dble(it))
    jz_ext(NXC+iANTENNA_OFFSET)   = j_amp * dsin(2.0d0*omega*dble(it))
    
    if (it <= CTIME+CWIDTH_END/2) then
      jy = jy + jy_ext 
      jz = jz + jz_ext
    end if
  end subroutine extCurrent_phased_array

  subroutine frequency(it,omega)
    integer,intent(in) :: it
    double precision,intent(out) :: omega

    double precision :: t

    t = dble(it)
    omega = ext_omega_a * t + ext_omega_b0
  end subroutine frequency

  subroutine amplitude(it,j_amp)
    integer,intent(in) :: it
    double precision,intent(out) :: j_amp

    double precision,parameter :: f = 0.01d0
    double precision :: a_begin_half, a_end_half
    double precision :: tmp_b, tmp_e

    ! Coefficient for begin 
    a_begin_half = 1d0/cwidth_begin*log((1d0-f)/f)
    a_end_half   = 1d0/cwidth_end  *log((1d0-f)/f)

    tmp_b = 0.5d0 * dtanh(a_begin_half*(it - cwidth_begin/2.0) + 1d0)
    tmp_e = 0.5d0 * dtanh(a_end_half*(it - cwidth_begin - cwidth_end/2.0 - ctime) + 1d0)
    j_amp = ext_j_amp * (tmp_b - tmp_e)
  end subroutine amplitude
end module CalcCurrent
