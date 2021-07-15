module CalcPosition
  use SimulationParameters,   only: EXT_MIRROR, PARTICLE_LOSSCORN_ELIMINATE, &
                                    NS, NP, IS_PARABOLIC_DISTRIBUTION
  use SupplementalParameters, only: NP_TOTAL, XLEN, LOSSCORN_ANGLE_TAN
  use Slps,                   only: PI
  use mt_stream,              only: genrand_double2
  use RandomGenerator,        only: mts
  implicit none
  private
  public position
contains
  subroutine position(vx,vy,vz,x,iq)
    double precision,intent(inout),dimension(:) :: vx, vy, vz
    double precision,intent(inout),dimension(:) :: x
    integer(1)      ,intent(inout),dimension(:) :: iq

    integer :: ip

    !Random phase boundary
    double precision :: v_perp, phase
    double precision :: particle_angle_tan
    double precision, parameter :: epsilon = 1d-10

    integer :: is
    integer :: ifirst, ilast

    if (EXT_MIRROR == .true.) then
      !$omp parallel private(is, ifirst, ilast)
      ifirst = 0
      ilast = 0
      do is = 1, NS
        ifirst  = ilast+1
        ilast   = ilast + NP(is)
        !$omp do simd safelen(16), linear(ip:1), &
        !$omp private(v_perp, phase, particle_angle_tan)
        do ip = ifirst, ilast
          x(ip) = x(ip) + vx(ip)
 
          if (x(ip) >= 0d0 .and. x(ip) < XLEN ) cycle
          x(ip) = x(ip) - vx(ip)
          v_perp = dsqrt(vy(ip)**2+vz(ip)**2)
          v_perp = max(v_perp, 0d0)
          phase  = 2.0d0 * PI * genrand_double2(mts)

          vx(ip) = -vx(ip)
          vy(ip) = v_perp * dcos(phase)
          vz(ip) = v_perp * dsin(phase)

          ! LOSS CORN CHECK
          if (PARTICLE_LOSSCORN_ELIMINATE == .false. .and. IS_PARABOLIC_DISTRIBUTION(is) == .true. ) cycle

          particle_angle_tan = v_perp/(dabs(vx(ip))+epsilon)
          if(particle_angle_tan < LOSSCORN_ANGLE_TAN) then
            iq(ip) = 0
          endif
        end do
       !$omp end do simd
      end do
      !$omp end parallel
    else
       !DIR$ SIMD
       do ip = 1, NP_TOTAL
          x(ip) = x(ip) + vx(ip)
          x(ip) = dmod(x(ip)+XLEN,XLEN)
       end do
    endif
  end subroutine position
end module
