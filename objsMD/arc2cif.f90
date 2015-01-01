PROGRAM pdb2gin
 IMPLICIT NONE
 REAL,    ALLOCATABLE           :: xcryst(:,:),xinbox(:,:)
 INTEGER, ALLOCATABLE           :: id(:)
 CHARACTER (LEN=4), ALLOCATABLE :: label(:)
 INTEGER                        :: n_atoms = 0
 INTEGER                        :: n_T = 0
 INTEGER                        :: IERR,i,j,k
 REAL, PARAMETER                :: pi=ACOS(-1.0)
 REAL                           :: atom(3),ouratom(3),cell_0(1:6)
 REAL                           :: rv(1:3,1:3),vr(1:3,1:3),average(1:3),x1,x2,x3
 CHARACTER (LEN=80)             :: line,string
 CHARACTER (LEN=4)              :: mol
 CHARACTER (LEN=6)              :: atomc
 CHARACTER (LEN=5)              :: model
 CHARACTER (LEN=6)              :: charct(6)
 CHARACTER (LEN=4)              :: typ(3)
 CHARACTER (LEN=10)             :: spacegroup = "P1"
! open PDB
! input=100
! output=200
!---------
 OPEN(100,file="input.arc",IOSTAT=IERR)
 OPEN(200,file="out.cif")
 IF( IERR /= 0 ) STOP "[ERROR] No se puede abrir input.pdb"
 fileopen_ARC: IF( IERR == 0) THEN
   do1: DO
    READ (100,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT do1
    IF (line(1:4 )=='PBC ') THEN
      n_T = n_T + 1
      n_atoms = 0
    ENDIF
    IF (line(1:2)=='Ge'.OR.line(1:2)=='Si'.or.&
        line(1:2)=='Li'.or.line(1:2)=='Mg'.or.&
        line(1:2)=='O '.or.line(1:2)=='Al'.or.&
        line(1:2)=='Ca'.or.line(1:2)=='Na'.or.&
        line(1:2)=='Cs'.or.line(1:2)=='Sr'.or.&
        line(1:2)=='K '.or.line(1:2)=='Fe') n_atoms = n_atoms + 1
   END DO do1
 ENDIF fileopen_ARC
! alicatamiento
 ALLOCATE(xcryst(3,n_atoms)   ,STAT=IERR)
 ALLOCATE(xinbox(3,n_atoms)   ,STAT=IERR)
 ALLOCATE(id(n_atoms)         ,STAT=IERR)
 ALLOCATE(label(n_atoms)      ,STAT=IERR)
 IF( IERR /= 0 ) STOP "Problema al alicatar variables"
 REWIND(100)
! Read ARC
   READ(100,'(A)', iostat = IERR ) line
   WRITE(6,*)line
   do2: DO
    IF(line(1:4)=='PBC ')THEN
      EXIT do2
    ELSE
      READ (100,'(A)', iostat = IERR ) line
      WRITE(6,*)line
      IF (IERR/=0) EXIT do2
    END IF
   END DO do2
   READ(line(6:13),*) cell_0(1)
   READ(line(16:23),*)cell_0(2)
   READ(line(26:33),*)cell_0(3)
   READ(line(36:43),*)cell_0(4)
   READ(line(46:53),*)cell_0(5)
   READ(line(56:63),*)cell_0(6)
   CALL cell(rv,vr,cell_0)
 read_coor_ARC: DO i=1,n_atoms    
    READ(100, '(A)', iostat = IERR ) line
    WRITE(6,*)line
    IF (IERR/=0) EXIT read_coor_ARC
    id(i)=i
    READ(line,*) typ(1),x1,x2,x3,string
!    READ(line,*) typ(1),x1,x2,x3,string
!    PRINT*,x1,x2,x3,typ(1)
    xinbox(1,i)=x1
    xinbox(2,i)=x2
    xinbox(3,i)=x3
    label(i)=typ(1)
    FORALL ( j=1:3 )
      xcryst(j,id(i)) = mod(vr(j,1)*xinbox(1,id(i))+vr(j,2)*xinbox(2,id(i))+vr(j,3)*xinbox(3,id(i))+100.0,1.0)
    END FORALL
    
 ENDDO read_coor_ARC
 READ (100,'(A)') line
 call outcif(200,xcryst,n_atoms,label,cell_0)
 CLOSE(100)
 CLOSE(200)
 STOP
END PROGRAM pdb2gin
!
subroutine GULP(u,xfrac,n,label,cell_0)
 implicit none
 integer, intent(in)              :: n,u
 real, intent(in)                 :: xfrac(3,n),cell_0(6)
 CHARACTER (LEN=4), intent(in)    :: label(n,2)
 CHARACTER (LEN=80)               :: line
 integer :: k,i
 WRITE(u,'(A)')'single conv'
 WRITE(u,'(A4,1x,a8,1x,i10)')'name','PDB2GIN:',n
 WRITE(u,'(A)')'cell'
 WRITE(u,'(6f11.6)')(cell_0(k), k=1,6 )
 !WRITE(u,'(6(F14.7))')(cell_0(k), k=1,6)
 WRITE(u,'(A)')'frac'
 write_: do k=1,n
!   WRITE(u,*)LABEL(k,1),'core',xfrac(1,k),xfrac(2,k),xfrac(3,k),0.0
  do i = 1,80
    line(i:i) = ' '
  enddo
  IF(label(k,1)/='sh   ')THEN
   line(1:5) = LABEL(k,1)
   line(7:10) = 'core'
  ELSE
   line(1:5) = ' O   '
   line(7:10) = 'shel' 
  END IF
  WRITE(line(12:20),'(f8.7)') xfrac(1,k)
  WRITE(line(22:30),'(f8.7)') xfrac(2,k)
  WRITE(line(32:40),'(f8.7)') xfrac(3,k)
  WRITE(line(42:51),'(f10.8)') 0.0
  WRITE(line(53:67),'(a)')'1.00000 0.00000'
  WRITE(u,'(a)')line
  !WRITE(u,'(a4,1X,a4,1X,4(F14.7,1x))')LABEL(k,1),'core',xfrac(1,k),xfrac(2,k),xfrac(3,k),0.0
 enddo write_
 close(1)
 RETURN
end subroutine
!
subroutine outcif(iout,xfrac,n,label,cell_0)
 implicit none
 integer, intent(in)              :: n,iout
 real, intent(in)                 :: xfrac(3,n),cell_0(6)
 CHARACTER (LEN=4), intent(in)    :: label(n)
 integer :: k,i
!
 write(iout,'(''data_cif'')')
 write(iout,'(''_symmetry_space_group_name_H-M'',5x,a)') "'P 1'"
 write(iout,'(''_symmetry_Int_Tables_number   '',5x,''1'')')
 write(iout,'(''_symmetry_cell_setting        '',5x,''triclinic'',/)')
 write(iout,'(''_cell_length_a'',21x,f9.4)')cell_0(1)
 write(iout,'(''_cell_length_b'',21x,f9.4)')cell_0(2)
 write(iout,'(''_cell_length_c'',21x,f9.4)')cell_0(3)
 write(iout,'(''_cell_angle_alpha'',18x,f8.4)')cell_0(4)
 write(iout,'(''_cell_angle_beta'',19x,f8.4)')cell_0(5)
 write(iout,'(''_cell_angle_gamma'',18x,f8.4)')cell_0(6)
 write(iout,'(/,''loop_'')')
 write(iout,'(''_atom_site_label'')')
 write(iout,'(''_atom_site_fract_x'')')
 write(iout,'(''_atom_site_fract_y'')')
 write(iout,'(''_atom_site_fract_z'')')
 write(iout,'(''_atom_site_occupancy'')')
 do i = 1,n
  write(iout,'(2x,a5,2x,3(f12.5,1x),f7.4)')label(i),(xfrac(k,i),k=1,3),1.0
 end do
end subroutine outcif
!
SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
 pi = ACOS(-1.0)
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3)) 
 call inverse(rv,vr,3)
 RETURN
END SUBROUTINE cell
!
SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  REAL :: radtodeg,PI
  PI=ACOS(-1.0)
  radtodeg=180.0/PI
!
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
!  Avoid round off errors for 90.0 and 120.0 degrees
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
!
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
END SUBROUTINE inverse
!
real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  RETURN
end function
