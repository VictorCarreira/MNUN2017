PROGRAM Trab01
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Aluno: Victor Ribeiro Carreira
!Professor: Leandro Di Bartolo
!Este programa visa cumprir os requisitos da disciplina MNUM,
!tópico 2.1- formulação explícita.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Declaração de Variáveis
IMPLICIT NONE
INTEGER(KIND=4)::i, l
REAL(KIND=4)::x,dX, dt, lambda, alfa, k, inicial, final, custocomputacional
REAL(KIND=4),DIMENSION(0:5000,0:5):: T !A expressão "0:" inicia o contador a partir do zero!!! Caso contrário ele se inicia a partir do 1!!!

!                           CASO I
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

! Condições de valor inicial da barra de alumínio:
dX=2.0
dt=0.1
k=0.835
!lambda=0.020875

lambda= k*dt/dx**2

!Definindo o valor inicial:
T=0.0

!Condições de Contorno de Dirichlet!
T(:,0)=100.0
T(:,5)=50.0

!Cáculo Numérico por Diferenças Finitas progressivas de maneira explícita

DO l=1,5000
  DO i=1,4
    T(l+1,i)=T(l,i)+lambda*(T(l,i+1)-2*T(l,i)+T(l,i-1))
  ENDDO
ENDDO


!Registra os resultados de Temperatura em um arquivo *txt

 OPEN(UNIT=1,FILE='Trab01b.txt')
 WRITE(1,*)'         0                 1                 2               3               4              5'      
 
  DO l=1,5000
    WRITE(1,*)T(l,:)
  ENDDO
!  WRITE(1,*)'++++++++++++++ INFORMAÇÕES ADICIONAIS +++++++++++++'
!  WRITE(1,*)"Tempo de Máquina =",custocomputacional, 'Segundos'
!  WRITE(1,*),T

 !Cálculo do Custo Computacional
CALL cpu_time(final)
custocomputacional=final-inicial
PRINT*, 'Custo Computacional=',custocomputacional, 'segundos'
PRINT*,'*********************** FIM ***************************'


END PROGRAM Trab01
