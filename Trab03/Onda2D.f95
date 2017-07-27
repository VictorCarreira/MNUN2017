PROGRAM Onda2D
! Declaração de Variáveis
IMPLICIT NONE
REAL, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
REAL, PARAMETER:: DBL = SELECTED_REAL_KIND(p=14, r=200)
INTEGER::i, j, k, h, Nx, Nz, Nt,x, z ! k = loop temporal
REAL(KIND=DBL)::DeltaT, rho, c, pi, a, fc, fcorte, t, t0, inicial, final, custocomputacional
!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:):: P
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:,:)::P1, P2, P3, S



!Condição inicial 
 P1(i,j,k=0)=0
 P2(i,j,k=0)=0
 P3(i,j,k=0)=0


!Condição de contorno



!Fonte

pi = 3.1416
fc = fcorte/3*SQRT(pi)
t0 = 2*SQRT(pi)/fcorte
td = t-t0

S=a*[(2*pi*(pi*fc*td)**2)-1] EXP(-pi*(pi*fc*td)**2) ! Pergutar que tipo de fonte é essa?

!Cálculo das DF (Stencil)

DO k=2,Nt
    P(z,x,k)= P(z,x,k) - source(sa,(n-1)*dt-st0, sf2) ! que source é esse?
        DO i-1,Nx
            DO j=1,Nz
                P2(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P1(i-2,j,k)+P1(i,j-2,k)+P1(i+2,j,k)+P1(i,j+2,k)) + 16*(P1(i-1,j,k)+P1(i,j-1,k)+P1(i+1,j,k)+P1(i,j+1,k))-60*P1(i,j+1,k)] + 2* P1(i,j,k) - P1(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k) 
                P3(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P2(i-2,j,k)+P2(i,j-2,k)+P2(i+2,j,k)+P2(i,j+2,k)) + 16*(P2(i-1,j,k)+P2(i,j-1,k)+P2(i+1,j,k)+P2(i,j+1,k))-60*P2(i,j+1,k)] + 2* P2(i,j,k) - P2(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k) 
            ENDDO
        ENDDO
ENDDO







END PROGRAM Onda2D