PROGRAM gin2dlpoly
! {{ ... }}
 IMPLICIT NONE
 INTEGER       :: h,i,j,k,l,m,n
 INTEGER       :: IERR = 0
 INTEGER       :: n_atoms = 0
 REAL          :: dist1 = 0.0
 REAL          :: dist2 = 0.0
 REAL          :: r1(3),r2(3),r3(3),cell_0(6),rv(3,3),vr(3,3),r
 REAL, ALLOCATABLE              :: xcryst(:,:),dist(:,:)
 REAL, ALLOCATABLE              :: xcarte(:,:)
 REAL, ALLOCATABLE              :: charge(:)
 CHARACTER (LEN=4), ALLOCATABLE :: label(:,:,:)
 CHARACTER (LEN=4)              :: label0(2)
 CHARACTER (LEN=1)              :: lf = ' '
 CHARACTER (LEN=8)              :: atmnam
 CHARACTER (LEN=80)             :: line
 INTEGER, ALLOCATABLE           :: adj(:,:)
 INTEGER                        :: tbp(20)
 CHARACTER (LEN=1), ALLOCATABLE :: adj_char(:,:)
!
 OPEN(10,file="CONFIG")
 OPEN(11,file="FIELD.shell")
 OPEN(12,file="FIELD.bends")
 OPEN(13,file="ConnectivityList.txt")
!
! open GULP file:
!
 fileopen_gin: IF( IERR == 0) THEN
  do1: DO
    READ(5,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT do1
    IF (line(1:4)=='name') THEN
       WRITE(10, Fmt='(a72)') line(1:72)
       WRITE(10,'(2i10)')0,3
    END IF
    IF (line(1:4)=='cell') EXIT do1
  ENDDO do1
  READ(5,'(6f11.6)')(cell_0(k), k=1,6 )
  CALL cell(rv,vr,cell_0)
  DO k=1,3
   WRITE(10, Fmt='(3f20.10,a12)') ( rv(k,l), l=1,3 ), Repeat( ' ', 12 )
  END DO
  READ(5,'(A)',IOSTAT=IERR) line ! frac
  DO
    READ(5,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    IF ((line(7:10)=='core'.or.line(7:10)=='shel').and. &
         line(53:67)=='1.00000 0.00000') THEN
         n_atoms = n_atoms + 1
    ENDIF
  END DO
 END IF fileopen_gin
 REWIND(5)
! allocate
 ALLOCATE(xcryst(n_atoms,0:3) ,STAT=IERR)
 ALLOCATE(dist(n_atoms,n_atoms),STAT=IERR)
 ALLOCATE(xcarte(n_atoms,1:3) ,STAT=IERR)
 ALLOCATE(label(n_atoms,1:2,1:2) ,STAT=IERR)
 ALLOCATE(charge(n_atoms) ,STAT=IERR)
 IF(IERR/=0) STOP '[ERROR] allocate variables'
! CONFIG FILE (10)
 do2: DO
    READ(5,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT do2
    IF (line(1:4)=='frac') EXIT do2
 END DO do2
! Read coordinates
 DO k=1,n_atoms
    READ(5,'(A)',IOSTAT=IERR ) line
!  Atom Type: 1, element 2, core/shel
    READ(line(1:4),*)  label(k,1,1)
    READ(line(7:10),*) label(k,2,1)
!  Coordinates (gulp code: wcoord.f90 )
    READ(line(12:20),*) xcryst(k,1)
    READ(line(22:30),*) xcryst(k,2)
    READ(line(32:40),*) xcryst(k,3)
    DO j=1,3
     xcryst(k,j)=MOD(xcryst(k,j)+100.0,1.0)
    END DO
    READ(line(42:51),*) charge(k)
    do l = 1,8
     atmnam(l:l)=' '
    enddo
    atmnam(1:4)=label(k,1,1)(1:4)
    shelllabel: IF(label(k,2,1)=='shel')THEN
      do3: do l=1,8
         if(atmnam(l:l)==' ') then
           atmnam(l:l+1)='sh'
           EXIT do3
         end if
      end do do3
    END IF shelllabel
    label(k,1,2)=atmnam
    WRITE(10,Fmt='(a8,i10,a54,a1)')atmnam,k,Repeat(' ',54),lf
    xcarte(k,1) = rv(1,1)*xcryst(k,1)+rv(1,2)*xcryst(k,1)+rv(1,3)*xcryst(k,1)
    xcarte(k,2) = rv(2,1)*xcryst(k,2)+rv(2,2)*xcryst(k,2)+rv(2,3)*xcryst(k,2)
    xcarte(k,3) = rv(3,1)*xcryst(k,3)+rv(3,2)*xcryst(k,3)+rv(3,3)*xcryst(k,3)
    WRITE(10,Fmt='(3g20.10,a12,a1)')(xcarte(k,l),l=1,3),Repeat(' ',12),lf
    !WRITE(10,Fmt='(3g20.10,a12,a1)')(0.0,l=1,3),Repeat(' ',12),lf
    !WRITE(10,Fmt='(3g20.10,a12,a1)')(0.0,l=1,3),Repeat(' ',12),lf
 END DO
! {{
! distance matrix:
! ===============
  DO i=1,n_atoms
     DO j=i,n_atoms
        IF(i==j)dist(i,j)=0.0
         FORALL ( k=1:3 )
          r1(k)= xcryst(j,k)
          r2(k)= xcryst(i,k)
         END FORALL
         CALL make_distances(cell_0,r1,r2,rv,r) 
         dist(i,j)=r
         dist(j,i)=dist(i,j)
     END DO
  END DO
!                  }}
! {{
! acotar distancias
! =================
  ALLOCATE(adj(n_atoms,n_atoms),STAT=IERR)
  ALLOCATE(adj_char(n_atoms,n_atoms) ,STAT=IERR)
  IF(IERR/=0) STOP '[ERROR] allocate'
  pairs: DO i=1,n_atoms
   DO j=i,n_atoms
      adj(i,j)=0.0
      adj_CHAR(i,j)=' '
      IF(j/=i)THEN
       IF((label(i,1,2)=='O   '.and.label(j,1,2)=='Osh ').or.&
          (label(j,1,2)=='O   '.and.label(i,1,2)=='Osh '))THEN
           dist1 = 0.0 !core ... shell
           dist2 = 1.0
       ELSE
        IF((label(i,1,2)=='Osh '.and.label(j,1,2)=='Si  ').or.&
          (label(j,1,2)=='Osh '.and.label(i,1,2)=='Si  '))THEN
           dist1 = 1.0 !Si ... O
           dist2 = 2.0
        ELSE
         IF((label(i,1,2)=='Osh '.and.label(j,1,2)=='Al  ').or.&
          (label(j,1,2)=='Osh '.and.label(i,1,2)=='Al  '))THEN
           dist1 = 1.0 !Al ... O
           dist2 = 2.0
         ELSE
           dist1 = 999999999999.0 ! cation ... O
           dist2 = 9999999999999.0
         END IF
        END IF
       END IF
! {{ condicion de enlace en un par ( XX ... YY )
       IF(dist(i,j)>=dist1.AND.dist(i,j)<dist2)THEN
         adj(i,j)=1.0
         adj(j,i)=adj(i,j)
         adj_CHAR(j,i)='@'
         adj_CHAR(i,j)=adj_CHAR(j,i)
       ENDIF
      END IF
   END DO
  END DO pairs
  conectivity: DO i=1,n_atoms
     k=0
     DO j=1,n_atoms
        k=k+adj(i,j)
     ENDDO
     xcryst(i,0)=k
  END DO conectivity
  WRITE(13,'(1000(I1))')(int(xcryst(k,0)),k=1,n_atoms) ! }}
!{{ 
!  FIELD FILE (11)
!  core...shell bonds
  DO k=1,n_atoms
      IF(label(k,1,2)=='O   ')THEN
        DO l=1,n_atoms
           IF(k/=l)THEN
            IF(label(l,1,2)=='Osh ')THEN
              IF(adj(k,l)==1)THEN
                WRITE(6,*)k,label(k,1,2),(xcryst(k,i),i=0,3)
                WRITE(6,*)l,label(l,1,2),(xcryst(l,i),i=0,3)
                WRITE(6,*)dist(k,l)
                IF(dist(k,l)>=1.0)STOP 'Non-realistic O...Osh distance'
                WRITE(6,*)'---------------------------------'
                WRITE(11,'(1x,i4,1x,i4,1x,f14.7)')k,l,74.920
              END IF
            END IF
           END IF
        END DO
      END IF
  END DO
! angles
  l=0
  doi: DO i=1,n_atoms
     k=0
     IF(label(i,1,2)=='Si  '.or.label(i,1,2)=='Al  ')THEN
       doj: DO j=1,n_atoms
          IF(i/=j)THEN
            IF(label(j,1,2)=='Osh ')THEN
              IF(adj(i,j)==1)THEN ! 1 2, 1 3, 1 4, 2 3, 2 4, 3 4 
                k=k+1             ! combinaciones de angulos
                tbp(k)=j          ! para el tetraedro
              END IF
            END IF
          END IF
       END DO doj
       DO j=1,k-1
          DO m=j+1,k
           l=l+1
           WRITE(12,'(a4,1x,i4,1x,i4,1x,i4,1x,2(f14.7,1x))')'harm',tbp(j),i,tbp(m),2.0972,109.47000
           WRITE(6,*)label(tbp(j),1,2),label(i,1,2),label(tbp(m),1,2),'(',l,')'
           WRITE(6,*)tbp(j),label(tbp(j),1,2),(xcryst(tbp(j),n),n=1,3)
           WRITE(6,*)i,label(i,1,2),(xcryst(i,n),n=1,3)
           WRITE(6,*)tbp(m),label(tbp(m),1,2),(xcryst(tbp(m),n),n=1,3)
           WRITE(6,*)label(tbp(j),1,2),'...',label(i,1,2),dist(i,tbp(j))
           IF(dist(i,tbp(j))==0.0) STOP 'ERROR: r=0.0 Angles'
           WRITE(6,*)label(i,1,2),'...',label(tbp(m),1,2),dist(i,tbp(m))
           IF(dist(i,tbp(m))==0.0) STOP 'ERROR: r=0.0 Angles'
           WRITE(6,*)label(tbp(m),1,2),'...',label(tbp(j),1,2),dist(tbp(m),tbp(j))
           IF(dist(tbp(m),tbp(j))==0.0) STOP 'ERROR: r=0.0 Angles'
           WRITE(6,*)'---------------------------------' 
          END DO
       END DO
     END IF
  END DO doi
  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
! WRITE(11,'(a6)')'FINISH'
 STOP
 CONTAINS
! {{ main contains :
  SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
   IMPLICIT NONE
   REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
   REAL,    intent(out) :: dist
  !REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
   REAL                 :: d_image(1:125),image(3,125)
   INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
   REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
   k=0
   do l=-2,2
    do m=-2,2
      do n=-2,2
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
    enddo
   enddo
   dist=MINVAL(d_image)
   RETURN
  END SUBROUTINE
!
  REAL FUNCTION DISTANCE(atom,ouratom,rv)
   IMPLICIT NONE
   INTEGER :: j
   REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
   REAL :: rv(3,3)
   FORALL ( j=1:3 )
    o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
    o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
    dist(j) = o_ouratom(j) - o_atom(j)
   END FORALL
   DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
  END FUNCTION
!
 SUBROUTINE cell(rv,vr,cell_0)
! ======================
! GULP
! ======================
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
! ======================
! GULP
! ======================
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
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
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
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
!==========================================================
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
! }}
END PROGRAM gin2dlpoly
