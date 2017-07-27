PROGRAM Onda2D
! Declaração de Variáveis
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER, PARAMETER:: DBL = SELECTED_REAL_KIND(p=1, r=200)
INTEGER(KIND=SGL)::i, j, k,np, kt ! kt= loop temporal
REAL(KIND=DBL)::x, dx, dt, lambda, l, inicial, final, custocomputacional
!REAL(KIND=4),DIMENSION(0:5000,0:5):: T !A expressão "0:" inicia o contador a partir do zero!!! Caso
!contrário ele se inicia a partir do 1!!!
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:):: aux
!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:):: T
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:)::T1,T2














END PROGRAM Onda2D