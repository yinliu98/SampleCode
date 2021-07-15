module InitFieldParticle
  use SimulationParameters,   only: NTIME, DX, DT, NX, cv,&
                                    NS, NP, omega_p, QM, v_para, v_perp, &
                                    v_d, PITCH_DEG, omega_c, &
                                    ENABLE_INJECTION, IS_INJECTION, &
                                    THETA_DEG, &
                                    SUBTRACT_BI_MAXWELLIAN, SUBTRACT_BETA, &
                                    SUBTRACT_RHO, &
                                    IS_PARABOLIC_DISTRIBUTION, &
                                    PARTICLE_LOSSCORN_ELIMINATE, &
                                    EXT_MIRROR, ext_mirror_a, ISEED
  use Slps, only: PI, elle
  use supplementalParameters, only: COSTH, SINTH, LOSSCORN_ANGLE_TAN, &
                                    NP_TOTAL, XLEN, csq
  use mt_stream,              only: genrand_double2
  use RandomGenerator,        only: mts, randN, subRandN, randUi
  use InjectParticles,        only: initInjection
#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif
  !$ use omp_lib
  implicit none
  double precision,save :: b0
  double precision,save :: bx0
  double precision,save :: by0
  private
  !Parameters/Variables
  public :: b0, bx0, by0

  public initField, initParticle
contains
  subroutine initField(ex,ey,ez,by,bz,jx,jy,jz,rho)
    double precision,intent(out),dimension(:) :: ex, ey, ez
    double precision,intent(out),dimension(:) :: by, bz
    double precision,intent(out),dimension(:) :: jx, jy, jz
    double precision,intent(out),dimension(:) :: rho

    call efield (ex ,ey ,ez)
    call bfield (by ,bz)
    call current(jx ,jy ,jz)
    call charge (rho)
  end subroutine initField

  subroutine initParticle(x,vx,vy,vz,iq)
    double precision,intent(out),dimension(:) :: x
    double precision,intent(out),dimension(:) :: vx, vy, vz
    integer(1)      ,intent(inout),dimension(:) :: iq

    integer :: is, ip
    integer :: ifirst, ilast
 
    ! First touch to particles
    
    !$omp parallel private(ifirst, ilast)
    ifirst = 0
    ilast  = 0
    do is = 1, NS
     ifirst = ilast + 1
     ilast  = ilast + NP(is)
     !$omp do simd linear(ip:1)
     do ip = ifirst, ilast
       x(ip)  = 0d0
       vx(ip) = 0d0
       vy(ip) = 0d0
       vz(ip) = 0d0
       iq(ip) = 1
     end do
     !$omp end do simd
    end do
    !$omp end parallel

    ifirst = 0
    ilast = 0
    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      if (EXT_MIRROR == .true. .and. IS_PARABOLIC_DISTRIBUTION(is) == .true.) then
        call mirrorParticles(is, ifirst,ilast,x,vx,vy,vz)
      else
        call uniformParticles(is, ifirst,ilast,x,vx,vy,vz)
      end if
    end do

    ! Ensure particles are in the region
    do ip=1, NP_TOTAL
      x(ip) = dmod(x(ip)+XLEN,XLEN)
    end do

    call initInjection(iq)
  end subroutine initParticle
  !===============================================================================
  ! Private Subroutines for variableListSetVariables
  !===============================================================================
  subroutine bfield(by,bz)
    double precision,intent(out),dimension(:) :: by, bz
    integer :: ix

    b0  = omega_c/QM(1)
    by = 0d0
    bz = 0d0
    bx0 = 0d0
    by0 = 0d0

    bx0 = b0 * COSTH
    by0 = b0 * SINTH
    by = by0
  end subroutine bfield

  subroutine efield(ex,ey,ez)
    double precision,intent(out),dimension(:) :: ex, ey, ez
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
  end subroutine efield

  subroutine current(jx,jy,jz)
    double precision,intent(out),dimension(:) :: jx, jy, jz
    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0
  end subroutine current

  subroutine charge(rho)
    double precision,intent(out),dimension(:) :: rho
    rho = 0.0d0
  end subroutine charge

  subroutine uniformParticles(is, ifirst, ilast, x, vx, vy, vz)
    integer,         intent(in)               :: is, ifirst, ilast
    double precision,intent(out),dimension(:) :: x
    double precision,intent(out),dimension(:) :: vx, vy, vz

    double precision :: u_perp, tmp

    double precision :: pchr, vdpa, vdpe
    double precision :: phase

    double precision :: uxi, uyi, uzi, ux, uy, uz

    double precision :: g

    integer :: ip

    !$omp parallel do
    do ip = ifirst, ilast
      x(ip) = genrand_double2(mts)
    end do
    !$omp end parallel do

    x = XLEN * x

    pchr  = PI/180.0d0*PITCH_DEG(is)
    vdpa  = v_d(is)*dcos(pchr)
    vdpe  = v_d(is)*dsin(pchr)
    !$omp parallel do default(shared), private(u_perp, tmp,phase,uxi,uyi,uzi, ux,uy,uz,g)
    do ip = ifirst, ilast
      if (SUBTRACT_BI_MAXWELLIAN(is)) then
        u_perp = v_perp(is) * subRandN(mts,SUBTRACT_BETA(is),SUBTRACT_RHO(is))
        tmp = genrand_double2(mts)
        phase = 2.0d0*PI*tmp
        uyi   = u_perp * dcos(phase)
        uzi   = u_perp * dsin(phase)
        tmp = genrand_double2(mts)
        phase = 2.0d0*PI*tmp
        uxi   = v_para(is) * randN(mts) + vdpa
        uyi   = uyi + vdpe*dcos(phase)
        uzi   = uzi + vdpe*dsin(phase)
      else
        tmp = genrand_double2(mts)
        phase = 2.0d0*PI*tmp
        uxi   = v_para(is) * randN(mts) + vdpa
        uyi   = v_perp(is) * randN(mts) + vdpe*dcos(phase)
        uzi   = v_perp(is) * randN(mts) + vdpe*dsin(phase)
      endif

      ! rotation to the direction of the magnetic field
      ux =  uxi*COSTH - uyi*SINTH
      uy =  uxi*SINTH + uyi*COSTH
      uz =  uzi

      g = cv / dsqrt(csq +ux*ux +uy*uy +uz*uz)

      vx(ip) = ux*g
      vy(ip) = uy*g
      vz(ip) = uz*g
    end do
    !$omp end parallel do
  end subroutine uniformParticles

  subroutine mirrorParticles(is, ifirst, ilast, x, vx, vy, vz)
    integer,         intent(in)               :: is, ifirst, ilast
    double precision,intent(out),dimension(:) :: x
    double precision,intent(out),dimension(:) :: vx, vy, vz

    !MT generation
    double precision :: pchr, vdpa, vdpe

    double precision :: u_perp, u_para
    double precision :: alpha_tan

    double precision :: mirror_point, xmax, u_para_xmax
    double precision :: pitch_rad_at_xmax

    double precision :: randpi
    double precision :: pos_phase_rad
    double precision :: pos

    double precision :: u_perp_t, u_para_t
    double precision :: phase
    double precision :: uxi, uyi, uzi, ux, uy, uz
    double precision :: g

    double precision :: u_para_max
    double precision :: max_ellips_len

    double precision :: ellips_a, ellips_b
    double precision :: ellips_kc, ellips_len

    !Vector length for gyrotropic parameters
    integer :: ivlen
    integer :: ip
    integer :: ivfirst, ivlast

    integer, parameter :: NPV = 1024

    double precision, parameter :: X_CENTER = XLEN/2.0d0
    double precision, parameter :: epsilon = 1d-16

    double precision :: beta
    double precision :: rho

    ivfirst = ifirst
    ivlast = ifirst

    pchr  = PI/180.0d0*PITCH_DEG(is)
    vdpa  = v_d(is)*dcos(pchr)
    vdpe  = v_d(is)*dsin(pchr)

    beta = SUBTRACT_BETA(is)
    rho = SUBTRACT_RHO(is)
    
    !STEP1: Set maximum arc length of bounce motion
    u_para_max = v_para(is)*5d0
    ellips_a = max(X_CENTER, u_para_max)
    ellips_b = min(X_CENTER, u_para_max)
    ellips_kc = ellips_b/ellips_a
    max_ellips_len = 4d0 * ellips_a * elle(0.5d0*PI, ellips_kc)

    do
      do
        ! STEP2: generate distribution function at equator
        if (SUBTRACT_BI_MAXWELLIAN(is)) then
          u_perp = v_perp(is) * subRandN(mts,beta,rho)
        else
          uyi = v_perp(is) * randN(mts)
          uzi = v_perp(is) * randN(mts)
          u_perp = dsqrt(uyi ** 2 + uzi ** 2)
        endif

        u_para = v_para(is)* randN(mts) + vdpa
        u_perp = u_perp + vdpe

        ! STEP3: loss cone
        if (PARTICLE_LOSSCORN_ELIMINATE) then
          alpha_tan = u_perp/(dabs(u_para) + epsilon)
          if (alpha_tan > LOSSCORN_ANGLE_TAN) exit
        else
          exit
        endif
      end do

      ! STEP4: Determine arc length in phase space
      mirror_point = dabs(u_para)/(u_perp + epsilon)/dsqrt(ext_mirror_a)
      xmax = min(X_CENTER, mirror_point)
      u_para_xmax = dsqrt(max(0d0, u_para **2 - u_perp **2 * ext_mirror_a * xmax **2))
      pitch_rad_at_xmax = datan(u_para_xmax*mirror_point/(dabs(u_para)+epsilon)/xmax)

      if ( mirror_point > u_para_max) then
        ellips_a = mirror_point
        ellips_b = dabs(u_para)
        ellips_kc = ellips_b/ellips_a
      else
        ellips_a = dabs(u_para)
        ellips_b = mirror_point
        ellips_kc = ellips_a/ellips_b
        pitch_rad_at_xmax = 0.5d0 * PI - pitch_rad_at_xmax
      endif

      ellips_len = 4d0 * ellips_a * (elle(0.5d0*PI, ellips_kc) - elle(pitch_rad_at_xmax, ellips_kc))

      ivlen = floor(NPV*ellips_len/max_ellips_len)
      ivlen = min(ivlen, ilast-ivfirst+1)

      ivlast  = ivlast + ivlen
      !$omp parallel do if(ivlen > 1024), &
      !$omp default(shared), &
      !$omp private(randpi, pos_phase_rad, pos), &
      !$omp private(u_perp_t, u_para_t), &
      !$omp private(phase, uxi, uyi, uzi, ux, uy, uz, g)
      do ip = ivfirst, ivlast
        ! STEP4 : uniform random numver in x-upara phase elipse
        randpi = randUi(mts) * PI
        pos_phase_rad = (PI - 2d0* pitch_rad_at_xmax) * genrand_double2(mts) + pitch_rad_at_xmax
        pos_phase_rad = pos_phase_rad + randpi
        pos = mirror_point * dcos(pos_phase_rad)

        ! STEP5 : adiabatic motion from equator to pos
        u_perp_t = u_perp * dsqrt(1d0+ext_mirror_a*pos**2)
        u_para_t = dsqrt(max(0d0, u_para**2 + u_perp**2 - u_perp_t**2))

        x(ip) = pos + X_CENTER
        u_para_t = u_para_t*(randUi(mts) * 2d0 - 1d0)
        phase = 2d0*PI*genrand_double2(mts)
        uxi = u_para_t
        uyi = u_perp_t * dcos(phase)
        uzi = u_perp_t * dsin(phase)

        ux =  uxi*COSTH - uyi*SINTH
        uy =  uxi*SINTH + uyi*COSTH
        uz =  uzi

        g = cv / dsqrt(csq +ux*ux +uy*uy +uz*uz)
        vx(ip) = ux*g
        vy(ip) = uy*g
        vz(ip) = uz*g
      end do
      !$omp end parallel do

      if (ivfirst >= ilast) then
        exit
      endif

      ivfirst = ivlast + 1
    end do
  end subroutine mirrorParticles
end module InitFieldParticle
