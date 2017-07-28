PROGRAM Onda2D
! Declaração de Variáveis
IMPLICIT NONE
REAL, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
REAL, PARAMETER:: DBL = SELECTED_REAL_KIND(p=14, r=200)
INTEGER::i, j, k, h, Nx, Nz, Nt, x, z, n, t, t0, ! k = loop temporal
REAL(KIND=DBL)::DeltaT, DeltaX, DeltaZ, rho, c, pi, a, fc, fcorte, alfa2, alfa4, beta, inicial, final, custocomputacional
!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:):: P
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:,:,:)::P1, P2, P3
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:,:)::S

!INPUT
Nx=10
Nz=10
Nt=10
pi=3.1416
n=10
x=10
z=10
t=100
alfa2=5!segunda ordem no espaço
alfa4=10!quarta ordem no espaço
beta=4!segunda ordem no tempo

DeltaX=x/Nx
DeltaZ=z/Nz
DeltaT=t/Nt

ALLOCATE(P1(Nx,Nz,Nt), P2(Nx,Nz,Nt), P3(Nx,Nz,Nt), S(Nx,Nt)) !Alocando as matrizes

!Condição inicial 
 P1(Nx,Nz,Nt=0)=0
 P2(Nx,Nz,Nt=0)=0
 P3(Nx,Nz,Nt=0)=0
 !P4(Nx,Nz,Nt=0)=0

!Condição de contorno (Oneway):

!Borda Direita 
DO k= 1,n !n=Nt??
    DO i=1,Nx   
    P(i,j,k+1)= - (c(i,j)*DeltaT/DeltaX) * (P(i,j,k)-P(i-1,j,k)+P(i,j,k)  
    ENDDO        
ENDDO

!Borda Esquerda
DO k= 1,n !n=Nt??
    DO i=1,Nx   
    P(i,j,k+1)= + (c(i,j)*DeltaT/DeltaX) * (P(i,j,k)-P(i-1,j,k)+P(i,j,k)   
    ENDDO        
ENDDO

!Fundo
DO k= 1,n !n=Nt??
    DO j=1,Nz   
    P(i,j,k+1)= - (c(i,j)*DeltaT/DeltaZ) * (P(i,j,k)-P(i,j-1,k)+P(i,j,k) 
    ENDDO        
ENDDO

!Contição de Contorno (Cerjan):





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


!Critério da Não-dispersão de segunda ordem:
IF h .LE. (MIN(C)/alfa2*fcorte) 
    PRINT*, 'Não-dispersa'
    THEN 
    PRINT*,'Dispersa' 
ENDIF

!Critério da Não-dispersão de quarta ordem:
IF h .LE. (MIN(C)/alfa4*fcorte) 
    PRINT*, 'Não-dispersa'
    THEN 
    PRINT*,'Dispersa' 
ENDIF


!Critério de Estabilidade:
IF DeltaT .LE. (h/beta*MAX(c))
    PRINT*,'ESTÁVEL'
    THEN
    PRINT*,'INSTÁVEL'
ENDIF


!OUTPUT




END PROGRAM Onda2D