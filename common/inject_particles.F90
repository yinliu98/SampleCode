module InjectParticles
  use SimulationParameters,   only: NX, NS, NP, &
                                    ENABLE_INJECTION, IS_INJECTION, &
                                    iINJECTION_TIME, &
                                    LOCALIZED_CHARGE_AVG, &
                                    HAS_ELECTROSTATIC_COMPONENT
  use SupplementalParameters, only: setCounterCharge, rho0, q, XLEN
  use CalcCharge,             only: charge
  use CalcPotential,          only: potential                           
  !$ use omp_lib

#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif

  implicit none
  private
  !Parameters/Variables
  public initInjection, injection
contains
  subroutine initInjection(iq)
    integer(1),dimension(:),intent(inout) :: iq

    integer :: is, ip 
    integer :: ifirst, ilast

    if (ENABLE_INJECTION == .false.) return 

    ifirst = 0
    ilast = 0

    do is =1, NS
       ifirst = ilast + 1
       ilast  = ilast + NP(is)
       if (IS_INJECTION(is)) then
          do ip = ifirst, ilast
             iq(ip) = 0 
          end do
       end if
    end do
  end subroutine initInjection
  
  subroutine injection(it,x,ex,rho,iq)
    integer   ,intent(in) :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(inout) :: ex
    double precision,dimension(:),intent(inout) :: rho
    integer(1),dimension(:),intent(inout) :: iq
    integer :: num_process

    integer :: is, ip,isi
    integer :: ifirst, ilast

    if (ENABLE_INJECTION == .false. ) return

    ifirst = 0
    ilast = 0
    
    do is = 1,NS
       ifirst = ilast + 1
       ilast  = ilast + NP(is)

       if (IS_INJECTION(is) == .false.) cycle
       if (it  /= iINJECTION_TIME(is) ) cycle

#ifdef _DEBUG
      num_process = 1
#else
      num_process = dble(iprocess)
#endif

       do ip = ifirst, ilast
          iq(ip) = 1 
       end do
       
       if (HAS_ELECTROSTATIC_COMPONENT == .true. ) then 
          if (LOCALIZED_CHARGE_AVG == .true.) then
             rho0 = 0d0
          else
             rho0 = 0d0
             do isi = 1, NS
                if (IS_INJECTION(isi) == .true.) cycle
                rho0 = rho0 - q(isi)*dble(NP(isi))
             end do
             rho0 = rho0 * num_process / XLEN
             rho0(1) = rho0(NX+1)
             rho0(NX+2) = rho0(2)
          endif
          
          call charge   (x,iq,                rho          ) ! MPI_ALLREDUCE
          if (LOCALIZED_CHARGE_AVG == .true.) then
             call setCounterCharge(rho)
             call charge   (x,iq,                rho          ) ! MPI_ALLREDUCE
          endif
          call potential(rho,                 ex           )
       endif
    end do
  end subroutine injection
end module InjectParticles
