MODULE Unidimensional
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER, PARAMETER:: DBL = SELECTED_REAL_KIND(p=1, r=200)
INTEGER(KIND=SGL)::i, l, np, Lt ! Lt= loop temporal
REAL(KIND=DBL)::x, dx, dt, lambda, k, inicial, final, custocomputacional
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:):: aux
REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:)::T1,T2


CONTAINS

SUBROUTINE DF (lp,nt,lambda,T1,T2)
INTEGER, INTENT(IN)::lp,nt 
REAL, INTENT(IN)::lambda
REAL, INTENT(INOUT)::T1
REAL, INTENT(OUT)::T2

!Cáculo Numérico por Diferenças Finitas progressivas de maneira explícita
DO l=0,Lt       ! loop temporal
  DO i=1,np-2   ! loop só do miolo da malha (0,np-1) são as bordas
    T2(i)=T1(i)+lambda*(T1(i+1)-2.0*T1(i)+T1(i-1))
  ENDDO
  T1=T2! Atualização da variável temperatura. Mapeia a variação temporal
  WRITE(1,*) T2 !Registra as temperaturas ao longo de lt instantes de tempo
ENDDO


END SUBROUTINE DF



END MODULE Unidimensional