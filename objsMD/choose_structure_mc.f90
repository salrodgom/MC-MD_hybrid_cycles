PROGRAM choose
 implicit none
 real(8)              :: beta,en0,sumw
 real(8), allocatable :: w(:),ener(:)
 integer              :: j,k,n,seed
 call get_seed(seed)
 read(5,*)k
 allocate(ener(k))
 allocate(w(k))
 read(5,*)beta
 beta = 1.0 / beta
 en0 = 0.0
 read_econfig: do j=1,k
   read(5,*)n,ener(n)
 end do read_econfig
 do0: do j=1,k
   !ener(j)=r8_uniform(-814100.1d+0,-814000.0d+0,seed)
   if (ener(j)<=en0) en0 = ener(j)
 end do do0
 do1: do j=1,k
   ener(j) = ener(j) - en0
   w(j) = exp(-beta*ener(j))
   sumw = sumw + w(j)
 end do do1
 write_: do j=1,k
   write(6,*)j,ener(j),w(j),sumw
 end do write_
 call select_config(k,w,sumw,n)
 write(6,*)n,ener(n)+en0
 CONTAINS
 SUBROUTINE select_config(k,w,sumw,n)
  implicit none
  integer              :: seed,k
  real(8),intent(in)   :: w(k),sumw
  integer,intent(out)  :: n
  real(8)              :: cumw,ws
  call get_seed(seed) 
  ws = r8_uniform(0.0d+0,sumw,seed)
  write(6,*) ws
  cumw = w(1)
  n = 1
  wh1: do while (cumw.lt.ws)
     n = n +1
     cumw = cumw + w(n)
  end do wh1
  RETURN
 END SUBROUTINE select_config
 SUBROUTINE get_seed(seed)
  implicit none
  integer day,hour,i4_huge,seed,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
!  Force 0 < TEMP <= 1.
  DO
    IF(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       CYCLE
    ELSE
       EXIT
    ENDIF
  ENDDO
  DO
    IF(1.0D+00<temp)then
       temp=temp-1.0D+00
       CYCLE
    else
       EXIT
    ENDIF
  ENDDO
  seed=int(dble(i4_huge)*temp)
  if(seed==0)       seed = 1
  if(seed==i4_huge) seed = seed-1
  RETURN
 END subroutine get_seed
 REAL(8) function r8_uniform(b1,b2,seed)
  implicit none
  real(8) b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed==0)then
       write(*,'(b1)')' '
       write(*,'(b1)')'R4_UNIFORM - Fatal error!'
       write(*,'(b1)')'Input value of SEED = 0.'
      stop '[ERROR] r4_uniform'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r8_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  RETURN
 end function r8_uniform
END PROGRAM choose
