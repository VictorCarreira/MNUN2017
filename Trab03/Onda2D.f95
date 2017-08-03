PROGRAM Onda2D
! Declaração de Variáveis
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER(KIND=SGL)::i, j, k, nx, nz, nt ! k = loop temporal
REAL(KIND=SGL)::DeltaT, DeltaX, DeltaZ, rho, c, pi, a, fc, fcorte, alfa2, alfa4, beta, inicial, final, custocomputacional
REAL(KIND=SGL)::  x, z, t, t0
REAL(KIND=SGL), ALLOCATABLE, DIMENSION(:,:,:):: P
!REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:)::P1, P2, Paux
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:,:)::S

!INPUT
nx=10
nz=10
nt=10
pi=3.1416
x=10.0
z=10.0
t=100.0
c=1500
alfa2=5!segunda ordem no espaço
alfa4=10!quarta ordem no espaço
beta=4!segunda ordem no tempo

DeltaX=x/nx
DeltaZ=z/nz
DeltaT=t/nt

!ALLOCATE(P1(nx,nz), P2(nx,nz), Paux(nx,nz), S(nx,nt)) !Alocando as matrizes
ALLOCATE(P(nx,nz,nt), S(nx,nz,nt))

!Condição inicial
 !P1(nx,nz)=0
 !P2(nx,nz)=0
 P(nx,nz,nt)=0

 !P4(Nx,Nz,Nt=0)=0

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
            P(i,j,k+1)=P(i,j,k)+DeltaT*P(i,j,k+1/2)
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
fc = fcorte/3*SQRT(pi)
t0 = 2*SQRT(pi)/fcorte
td = t-t0
S=a*[(2*pi*(pi*fc*td)**2)-1] * EXP(-pi*(pi*fc*td)**2) ! Pergutar que tipo de fonte é essa?


!Cálculo das DF (Stencil)
DO k=2,Nt
        DO i-1,Nx
            DO j=1,Nz
                P2(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P1(i-2,j,k)+P1(i,j-2,k)+P1(i+2,j,k)+P1(i,j+2,k)) + 16*(P1(i-1,j,k)+P1(i,j-1,k)+P1(i+1,j,k)+P1(i,j+1,k))-60*P1(i,j+1,k)] + 2* P1(i,j,k) - P1(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k)
                P3(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P2(i-2,j,k)+P2(i,j-2,k)+P2(i+2,j,k)+P2(i,j+2,k)) + 16*(P2(i-1,j,k)+P2(i,j-1,k)+P2(i+1,j,k)+P2(i,j+1,k))-60*P2(i,j+1,k)] + 2* P2(i,j,k) - P2(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k)
            ENDDO
        ENDDO
    P4(i,j,k)= P3(i,j,k) - source(sa,(n-1)*dt-st0, sf2) ! que source é esse?
ENDDO




!OUTPUT




END PROGRAM Onda2D
