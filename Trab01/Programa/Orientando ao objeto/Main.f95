PROGRAM Main
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Aluno: Victor Ribeiro Carreira
!Professor: Leandro Di Bartolo
!Este programa visa cumprir os requisitos da disciplina MNUM,
!tópico 2.1- formulação explícita.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Declaração de módulos
USE Unidimensional

! Declaração de Variáveis
IMPLICIT NONE
!INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
!INTEGER, PARAMETER:: DBL = SELECTED_REAL_KIND(p=14, r=200)
!INTEGER::i!, l, np, Lt ! Lt= loop temporal
!REAL(KIND=DBL)::x, dx, dt, lambda, k, inicial, final, custocomputacional
!REAL(KIND=4),DIMENSION(0:5000,0:5):: T !A expressão "0:" inicia o contador a partir do zero!!! Caso
!contrário ele se inicia a partir do 1!!!
!REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:):: aux
!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:):: T
!REAL(KIND=DBL),ALLOCATABLE, DIMENSION(:)::T1,T2

!                           EXPERIMENTO
!
! Definição da malha unidimensional de 6 pontos geometricamente variando no tempo:
!
!l=n|======0======0======0======0=======|
!
!l=3|======0======0======0======0=======|
!
!l=2|======0======0======0======0=======|
!         /|\    /|\    /|\    /|\
!l=1|======0======0======0======0=======|
!  i=0    i=1    i=2    i=3    i=4     i=5
!
!(T0=100°C)                         (T5=50°C) => Condição de contorno
!
!---|-------|------|------|------|------|------|-------> x (Comprimento, cm)
!   0      2.0    4.0    6.0    8.0    10.0   12.0
!
! Variáveis do problema:
! x, tamanho da barra(cm)
! k, condutividade térmica(W/mK)
! alfa², difusividade térmica
! lambda, número de Fourier
! t, tempo(s)
! T, Temperatura (°C)

CALL cpu_time(inicial)

! Leitura do input file:
OPEN(UNIT=2,FILE='input.txt')
READ(2,*) x     ! comprimento da barra (cm)
READ(2,*) np    ! discretizacao da barra
READ(2,*) Lt    ! Iteracoes temporais
READ(2,*) dt    ! discretizacao temporal (tempo fisico do experimento)
READ(2,*) k     ! condutividade termica do material

!ALLOCATE( T(0:Lt,0:np-1), aux(0:np-1) )
ALLOCATE(T1(0:np-1), T2(0:np-1), aux(0:np-1)) !np-1 indica o número de pontos da barra menos o 0 inicial que eu comecei a contagem

! Condições de valor inicial da barra de alumínio:
 dx =  x / ( REAL(np,KIND=DBL) - 1.0 ) ! conversao de inteiro em real 4 somente para o cálculo de dx
                                     ! diferenca entre numero de passos e numero de pontos!!!!!



!k=0.835
!lambda=0.020875
! Calculo do numero de fourier:
lambda= k*dt/dX**2

!Definindo o valor inicial:
!T = 0.0
T1=0.0 ! Zerando as variáveis
T2=0.0


!Condições de Contorno de Dirichlet!
!T(0:Lt,0)    = 100.0
!T(0:Lt,np-1) = 50.0

T1(0)=100.0
T1(np-1)=50.0

T2(0)=100.0
T2(np-1)=50.0


! Criando um vetor auxiliar para dimensionar o arquivo de saida:
aux(0:np-1) = 0.0
DO i=0,np-2
  aux(i+1) = aux(i) + dX
ENDDO
!Registra os resultados de Temperatura em um arquivo *txt
OPEN(UNIT=1,FILE='output.txt')! arquivo de saida
WRITE(1,*) aux

CALL DF (lp,nt,lambda,T1,T2)

!Cáculo Numérico por Diferenças Finitas progressivas de maneira explícita
!DO l=0,Lt       ! loop temporal
!  DO i=1,np-2   ! loop só do miolo da malha (0,np-1) são as bordas
    !T(l+1,i)=T(l,i)+lambda*(T(l,i+1)-2.0*T(l,i)+T(l,i-1))
!    T2(i)=T1(i)+lambda*(T1(i+1)-2.0*T1(i)+T1(i-1))
!  ENDDO
!  T1=T2! Atualização da variável temperatura. Mapeia a variação temporal
!  WRITE(1,*) T2 !Registra as temperaturas ao longo de lt instantes de tempo
!ENDDO





 !OPEN(UNIT=1,FILE='Trab01b2.txt')
!  WRITE(1,*) aux
  !DO l=1,Lt
  !  WRITE(1,*)T(l,:)
  !ENDDO


!Cálculo do Custo Computacional
CALL cpu_time(final)
custocomputacional=final-inicial
PRINT*, 'Custo Computacional=',custocomputacional, 'segundos'
PRINT*,'*********************** FIM ***************************'

END PROGRAM Main
