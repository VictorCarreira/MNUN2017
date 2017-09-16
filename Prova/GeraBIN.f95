PROGRAM Modelo
IMPLICIT NONE
REAL, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
REAL, PARAMETER:: DBL = SELECTED_REAL_KIND(p=14, r=200)
INTEGER::i, j, k, Nx, Nz, zr
REAL, ALLOCATABLE, DIMENSION(:,:)::vel

Nx=602
Nz=302

ALLOCATE(vel(Nz,Nx))

vel=1500




DO i=1,Nx
  zr=NINT(0.07*i+619,57)
  DO j=1,zr
    vel(i,j)=2000
  ENDDO
ENDDO
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
