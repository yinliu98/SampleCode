module BucketSort
  use SimulationParameters,   only: NX, NS, NP
  use SupplementalParameters, only: NP_TOTAL
  !$ use omp_lib
  implicit none
  private
  public bsort
contains
  subroutine bsort(x,vx,vy,vz,iq)
    double precision, dimension(:), intent(inout) :: x
    double precision, dimension(:), intent(inout) :: vx, vy, vz
    integer(1)      , dimension(:), intent(inout) :: iq

    integer         , dimension(NP_TOTAL) :: key
    call sort_key(x, key)
    call key_to_variables(key, x, vx, vy, vz, iq)
  end subroutine

  subroutine sort_key(x, key)
    double precision,dimension(:), intent(in)  :: x
    integer         ,dimension(:), intent(out) :: key

    integer         ,dimension(NP_TOTAL) :: idx

    integer,parameter :: BOX_NX = min(2, NX)
    integer,parameter :: iBOX = NX/box_NX
    integer,dimension(ibox) :: hist, offset

    integer :: ifirst, ilast
    integer :: is, ip, ib

    integer :: tmp

    ifirst = 0
    ilast = 0

    do is = 1, NS
      ifirst  = ilast + 1
      ilast   = ilast + NP(is)

      !$omp parallel do
      do ip = ifirst, ilast
        idx(ip) = floor(x(ip)/BOX_NX)
        key(ip) = ip
      end do
      !$omp end parallel do
      ! count times
      hist = 0

      !$omp parallel do reduction(+:hist)
      do ip = ifirst, ilast
         hist(idx(ip)+1) = hist(idx(ip)+1) + 1
      end do
      !$omp end parallel do

      ! Start position
      offset(1) = ifirst
      do ib = 2, iBOX
        offset(ib) = offset(ib-1) + hist(ib-1)
      end do

      ! Sort
      do ip = ifirst, ilast
        tmp = idx(ip) + 1
        key(offset(tmp)) = ip
        offset(tmp) = offset(tmp) + 1
      end do
    end do
  end subroutine

  subroutine key_to_variables(key,x,vx,vy,vz,iq)
    integer         ,dimension(:), intent(in)    :: key
    double precision,dimension(:), intent(inout) :: x
    double precision,dimension(:), intent(inout) :: vx, vy, vz
    integer(1)      ,dimension(:), intent(inout) :: iq

    double precision,dimension(NP_TOTAL) :: swap
    integer(1),dimension(NP_TOTAL) :: iq_swap

    integer :: ip

    swap = x
    !$omp parallel do
    do ip = 1, NP_TOTAL
      x(ip) = swap(key(ip))
    end do
    !$omp end parallel do

    swap = vx
    !$omp parallel do
    do ip = 1, NP_TOTAL
      vx(ip) = swap(key(ip))
    end do
    !$omp end parallel do

    swap = vy
    !$omp parallel do
    do ip = 1, NP_TOTAL
      vy(ip) = swap(key(ip))
    end do
    !$omp end parallel do

    swap = vz
    !$omp parallel do
    do ip = 1, NP_TOTAL
      vz(ip) = swap(key(ip))
    end do
    !$omp end parallel do

    iq_swap = iq
    !$omp parallel do
    do ip = 1, NP_TOTAL
      iq(ip) = iq_swap(key(ip))
    end do
    !$omp end parallel do
  end subroutine
end module
