module CalcPotential
  use SimulationParameters, only: NX
  use SupplementalParameters, only: XLEN
  implicit none
  public potential
contains

  subroutine potential(rho,ex)
    double precision,dimension(:),intent(in)    :: rho
    double precision,dimension(:),intent(inout) :: ex

    double precision,parameter :: XLEN_INV = 1d0/XLEN
    integer :: ix ! Iterator to X,Y direction
    
    double precision ex_avg

    do ix = 2, NX+1
       ex(ix) = ex(ix-1) + rho(ix)
    end do

    ex_avg = 0d0   
    do ix = 2, NX+1
       ex_avg = ex_avg + ex(ix)*XLEN_INV
    end do

    do ix = 2, NX+1
       ex(ix) = ex(ix) - ex_avg
    end do

    ex(1) = ex(NX+1)
    ex(NX+2) = ex(2)
  end subroutine potential
end module
