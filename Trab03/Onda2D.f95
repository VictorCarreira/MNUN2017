PROGRAM Onda2D
! Declaração de Variáveis
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER(KIND=SGL)::i, j, k, nx, nz, nt, snap_passos ! k = loop temporal
INTEGER(KIND=SGL)::nfonte
REAL(KIND=SGL)::DeltaT, DeltaX, DeltaZ, rho, c, alfa2, alfa4, beta, inicial, final, custocomputacional
REAL(KIND=SGL)::  xfonte, zfonte, t, t0,td ,fc ,fcorte, ampl_fonte
REAL(KINDSGL), PARAMETER::pi=3.1416
!REAL(KIND=SGL), ALLOCATABLE, DIMENSION(:,:,:):: P
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:)::P1, P2, P3
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:,:)::S

!INPUT
nx=10
nz=10
nt=10
x=10.0
z=10.0
t=100.0
c=1500
alfa2=5!segunda ordem no espaço
alfa4=10!quarta ordem no espaço
beta=4!segunda ordem no tempo

c=c**2 ! Escrever a eq da Onda2D

CALL Fonte()

DeltaX=x/nx
DeltaZ=z/nz
DeltaT=t/nt

!Posicão da fonte no grid
xi_fonte=NINT((sfonte/h)+1)

!ALLOCATE(P1(nx,nz), P2(nx,nz), Paux(nx,nz), S(nx,nt)) !Alocando as matrizes
ALLOCATE(P(nx,nz,nt), S(nx,nz,nt))

!Condição inicial
 P1=0.0
 P2=0.0
 P3=0.0




!Condição de contorno (Oneway):

!Borda Direita
DO k=1,nt
    DO i=1,nx
        P(i,j,k+1)= - (c(i,j)*DeltaT/DeltaX) * (P(i,j,k)-P(i-1,j,k)) +P(i,j,k)
        !P2(i,j) = -(c(i,j)*DeltaT/DeltaX)*(P1(i,j)-P1(i-1,j)) + P1(i,j)
    ENDDO
ENDDO


!Borda Esquerda
DO k=1,nt
    DO i=1,nx
        P(i,j,k+1)= - (c(i,j)*DeltaT/DeltaX) * (P(i,j,k)-P(i-1,j,k)) +P(i,j,k)
        !P2(i,j)= + (c(i,j)*DeltaT/DeltaX) * (P1(i,j)-P1(i-1,j)) +P1(i,j)
    ENDDO
ENDDO


!Fundo

DO k=1,nt
    DO j=1,nz
        P(i,j,k+1)= - (c(i,j)*DeltaT/DeltaZ) * (P(i,j,k)-P(i,j-1,k)) +P(i,j,k)
        !P2(i,j)= - (c(i,j)*DeltaT/DeltaZ) * (P1(i,j)-P1(i,j-1)) +P1(i,j)
    ENDDO
ENDDO


!Contição de Contorno (Cerjan):
DO k=1,nt
    DO j=1,nz
        DO i=1, nx
            P(i,j,k+1)=P(i,j,k)+DeltaT*P(i,j,k+1/2)!OLHAR NO EXCEL
            !P2(i,j)=P1(i,j)+DeltaT*Paux(i,j)
        ENDDO
    ENDDO
ENDDO


!Critério da Não-dispersão de segunda ordem:
IF h .LE. (MIN(C)/alfa2*fcorte) THEN
    PRINT*, 'Não-dispersa'
ELSE
    PRINT*,'Dispersa'
END IF

!Critério da Não-dispersão de quarta ordem:
IF h .LE. (MIN(C)/alfa4*fcorte) THEN
    PRINT*, 'Não-dispersa'
ELSE
    PRINT*,'Dispersa'
ENDIF


!Critério de Estabilidade:
IF DeltaT .LE. (h/beta*MAX(c)) THEN
    PRINT*,'ESTÁVEL'
ELSE
    PRINT*,'INSTÁVEL'
ENDIF


!Fonte





!Cálculo das DF (Stencil)
DO n=1,nt
  P2(zfonte,xfonte)=P1(zfonte,xfonte)-fonte(n)
        DO i-1,Nx
            DO j=1,Nz
                !P2(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P1(i-2,j,k)+P1(i,j-2,k)+P1(i+2,j,k)+P1(i,j+2,k)) + 16*(P1(i-1,j,k)+P1(i,j-1,k)+P1(i+1,j,k)+P1(i,j+1,k))-60*P1(i,j+1,k)] + 2* P1(i,j,k) - P1(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k)
                !P3(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P2(i-2,j,k)+P2(i,j-2,k)+P2(i+2,j,k)+P2(i,j+2,k)) + 16*(P2(i-1,j,k)+P2(i,j-1,k)+P2(i+1,j,k)+P2(i,j+1,k))-60*P2(i,j+1,k)] + 2* P2(i,j,k) - P2(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k)
                P3(i,j)=2*P2(i,j)-P1(i,j)+c(i,j)*(P2())
            ENDDO
            CALL Oneway()
            CALL Cerjan()
          ENDDO
    P1(i,j,k)= P3(i,j,k) - source(sa,(n-1)*dt-st0, sf2) ! que source é esse?
ENDDO


CONTAINS

  SUBROUTINE Fonte()
    IMPLICIT NONE
    t0 = 2*SQRT(pi)/fcorte
    td = t-t0
    fc = fcorte/3*SQRT(pi)
    ALLOCATE(fonte(nt))
    fonte=0.0
    nfonte=NINT(2*t0/dt)


    DO i=1,nfonte
      fonte=ampl_fonte*(2*pi*(pi)*fc*(i*dt-t0))**2-1.0)EXP(-pi*(pi)*fc*(i*dt-t0))**2)
    ENDDO

  ENDSUBROUTINE Fonte

  SUBROUTINE Oneway()
    IMPLICIT NONE

  ENDSUBROUTINE Oneway

SUBROUTINE Cerjan()
  IMPLICIT NONE

ENDSUBROUTINE Cerjan


END PROGRAM Onda2D
