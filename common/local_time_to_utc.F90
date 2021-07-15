module localTimeToUTC
  implicit none
  private
  public get_utc_date_and_time
contains
    subroutine get_utc_date_and_time(date, time)
    character(8),intent(out) :: date
    character(10),intent(out) :: time

    integer,dimension(8) :: value_day
    integer    :: l_year, l_month, l_day, l_hour, l_min
    integer    :: u_year, u_month, u_day, u_hour, u_min
    integer    :: zone_diff_min, zone_hour, zone_min

    call date_and_time(values=value_day)

    l_year = value_day(1)
    l_month = value_day(2)
    l_day = value_day(3)
    zone_diff_min = value_day(4)
    l_hour = value_day(5)
    l_min = value_day(6)

    zone_hour = floor(zone_diff_min / 60.0)
    zone_min = mod(zone_diff_min, 60)

    call lst2utc(l_year,l_month,l_day,l_hour,l_min, u_year, u_month, u_day, u_hour, u_min, zone_hour, zone_min)

    write(date(1:4),'(i4.4)') u_year
    write(date(5:6),'(i2.2)') u_month
    write(date(7:8),'(i2.2)') u_day
    write(time(1:2),'(i2.2)') u_hour
    write(time(3:4),'(i2.2)') u_min
    write(time(5:6),'(i2.2)') value_day(7)
    time(7:7) = '.'
    write(time(8:10),'(i3.3)') value_day(8)
  end subroutine

  subroutine lst2utc(lyr,lmo,ldy, lhr, lmin, uyr, umo, udy, uhr, umin, tdh, tdm)
    ! Description:
    !
    ! Author: am
    !
    ! Host: aofd165.fish.nagasaki-u.ac.jp
    ! Directory: /work2/kunoki/to_manda_sensei/Tools/MGDSST.AlongTrack
    !
    ! Revision history:
    !  This file is created by /usr/local/mybin/nff.sh at 09:19 on 09-12-2014.
    ! Reference
    !   1行で書けるうるう年判別法
    !   http://d.hatena.ne.jp/Kappuccino/20081025/1224906495
    !
    implicit none
    integer,intent(in)    :: lyr,lmo,ldy, lhr, lmin
    integer,intent(inout) :: uyr,umo,udy, uhr, umin
    integer,intent(in):: tdh, tdm
    integer y

    integer,dimension(12)::last_day
    data last_day/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

    uyr=lyr
    umo=lmo
    udy=ldy
    uhr=lhr

    uhr=lhr-tdh
    umin=lmin-tdm

    if(uhr<=0)then
       uhr=uhr+24
       udy=udy-1
    endif

    if(uhr==24)then
       uhr=0
       udy=udy+1
    endif

    if(udy<=0)then
       umo=umo-1
       udy=last_day(umo)
       if(umo == 2)then
          y=uyr
          udy = 28 + (1 / (y - y / 4 * 4 + 1)) * &
               &          (1 - 1 / (y - y / 100 * 100 + 1)) &
               &          +  (1 / (y - y / 400 * 400 + 1));
       endif
    endif

    if(umo<=0)then
       umo=1
       uyr=uyr-1
    endif
  end subroutine lst2utc
end module
