PROGRAM Modelo
IMPLICIT NONE
!REAL, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
!REAL, PARAMETER:: DBL = SELECTED_REAL_KIND(p=14, r=200)
INTEGER::i, j, Nx, Nz, zr, zc
REAL, ALLOCATABLE, DIMENSION(:,:)::vel

Nx=602
Nz=302

ALLOCATE(vel(Nz,Nx))

vel=1500.0
!Interface 1
DO j=1,Nx
  zr=NINT(-0.1495*j+140)
  DO i=1,zr
    vel(i,j)=2000.0
  ENDDO
ENDDO


!Interface 2
DO j=1,Nx
  IF (j.ge.245 .and. j.le.345) THEN
    zc=230+NINT(SQRT(2500.0-(j-295.0)**2))
  ELSE
    zc=230
  ENDIF
  DO i=1,zc
    vel(i,j)=2100.0
  ENDDO
ENDDO




!Modelo de duas camadas
! Este laço gera a matriz de velocidades camada1
!DO j=1,601
!  DO i=1,201
!    vel(i,j)=1500.0
!  ENDDO
!ENDDO


! Este laço gera a matriz de velocidades camada2
!DO j=1,601
!  DO i=202,401
!    vel(i,j)=2000.0
!  ENDDO
!ENDDO



!PRINT*,vel

! Salva a matriz de velicidades em um arquivo binário
OPEN(10,FILE='VELOCIDADES.bin', STATUS='UNKNOWN',ACCESS='DIRECT',&
FORM='UNFORMATTED',RECL=4*Nx*Nz)
WRITE(10,REC=1)((vel(i,j),i=1,Nz),j=1,Nx)! Isto é um laço interno
CLOSE(10)


END PROGRAM Modelo
