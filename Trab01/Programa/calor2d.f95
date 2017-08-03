PROGRAM calor2d
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Aluno: Victor Ribeiro Carreira
!Professor: Leandro Di Bartolo
!Este programa visa cumprir os requisitos da disciplina MNUM,
!tópico 2.1- formulação explícita.
!O problema do calor 2D
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Declaração de Variáveis
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER(KIND=SGL)::i, j, l, nx, ny, nz ! Contadores(i,j,l) e Parametros numéricos[nx(dimensão), ny(dimensão), nz(tempo real devido a instabilidade numérica) e nl(é um índice auxiliar de tempo inteiro)]
INTEGER(KIND=SGL)::nsnap,csnap=0!Número de snapshots para impressão e o contador de snapshots.
INTEGER(KIND=SGL)::dx,dy !Parâmetros numéricos distância nos três eixos.
REAL(KIND=SGL)::x, y, z, dz, Tesq, Tdir, Tinf, Tsup, inicial, final, custocomputacional
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:)::T1,T2, alfa !Vetores com as temperaturas
CHARACTER(LEN=3)::num_snap!Contador de snaps para impressão

!                           EXPERIMENTO
!
! Definição do Stencil bidimensional de 6x6 pontos variando somente no tempo futuro:
!
!
!  y(cm)
!   ^    (T0=100°C)                         (T5=50°C) => Condição de contorno
!   |       |-----------------------------------|
!   |       |       l=5                         |
!   |       |       /                           |
!  12.0  j=5|------0------0------0------0-------|
!   |       |      l=4                          |
!   |       |       /                           |
!  10.0  j=4|------0------0------0------0-------|
!   |       |       l=3                         |
!   |       |       /                           |
!  8.0   j=3|------0------0------0------0-------|
!   |       |       l=2                         |
!   |       |       /                           |
!  6.0   j=2|------0------0------0------0-------|
!   |       |       l=1                         |
!   |       |       /                           |
!  4.0   j=1|------0------0------0------0-------|
!   |       |      l=0                          |
!   |       |       /                           |
!  2.0   j=0|------0------0------0------0-------|
!   |      i=0    i=1    i=2    i=3    i=4     i=5
!   |   /   |-----------------------------------|
!   |  /z=(tempo)
!   | /
!---|/-------|------|------|------|------|------|-------> x (cm)
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
READ(2,*) x     ! comprimento da placa em x (cm)
READ(2,*) y     ! dimensao da placa em cm (y);
READ(2,*) dx    ! subdivisoes horizontais na placa (dx);
READ(2,*) dy    ! subdivisoes verticais na placa (dy);
READ(2,*) z     ! Iteracoes temporais
READ(2,*) dz    ! discretizacao temporal (tempo fisico do experimento)
READ(2,*) alfa  ! difusividade termica do material
READ(2,*) Tesq  ! temperatura na borda esquerda (°C)
READ(2,*) Tdir  ! temperatura na borda direita (°C)
READ(2,*) Tsup  ! temperatura na borda superior (°C)
READ(2,*) Tinf  ! temperatura na borda inferior (°C)
READ(2,*) nsnap ! número de snapshots

! Condições de valor inicial da barra de alumínio:
nx =  INT(x/dx)+1 ! conversao de inteiro em real 4 somente para o cálculo de dx ! diferenca entre numero de passos e numero de pontos!!!!!
ny =  INT(y/dy)+1
nz =  INT(z/dz)+1

PRINT*,'ny=', ny
PRINT*,'nx=', nx

ALLOCATE(T1(ny,nx), T2(ny,nx), alfa(ny,nx)) !n...-1 indica o número de pontos da barra menos o 0 inicial que eu comecei a contagem na direção x ou y. E !A expressão "0:" inicia o contador a partir do zero!!! Caso contrário ele se inicia a partir do 1!!!




!Definindo o valor inicial:
T1=0.0 ! Zerando as variáveis
T2=0.0


!Condições de Contorno de Dirichlet!

T1(:,1)=Tesq
T1(:,nx)=Tdir
T1(1,:)=Tsup
T1(ny,:)=Tinf

T2=T1 !As condições de contorno não variam tanto no tempo presente quanto no tempo futuro.


!Difusividade Térmica ao longo da placa
alfa=alfa*dz/(dx*dy)



!Cáculo Numérico por Diferenças Finitas progressivas de maneira explícita
DO l=1,nz ! laço temporal (A marcha no tempo)
  IF(MOD(l,50)==0) WRITE(*,*)'passo',l
  DO j=2,nx-1   ! laço só do miolo da malha (2,nx-1) são as bordas. O cálculo das DFS (espaço).
    DO i=2,ny-1   ! laço só do miolo da malha (2,nz-1) são as bordas. O cálculo das DFS (espaço).
    T2(i,j)=T1(i,j)+alfa(i,j)*(T1(i+1,j)-2.0*T1(i,j)+T1(i-1,j)+ T1(i,j+1)-2.0*T1(i,j)+T1(i,j-1))
    ENDDO
  ENDDO
  CALL Snap()
  T1(2:ny-1,2:nx-1)=T2(2:ny-1,2:nx-1)! Atualização da variável temperatura. Mapeia a variação temporal
ENDDO
WRITE(*,*)'Final da marcha temporal'

DEALLOCATE(T1,T2,alfa)

PRINT*,'Tamanho de T2',SIZE(T2)

!Cálculo do Custo Computacional
CALL cpu_time(final)
custocomputacional=final-inicial
PRINT*, 'Custo Computacional=',custocomputacional, 'segundos'
PRINT*,'*********************** FIM ***************************'



CONTAINS
!Criando os snapshots e gerando o arquivo de saída
SUBROUTINE Snap()
IMPLICIT NONE
IF(MOD(l,nsnap)==0)THEN !imprime os snapsshots de nsnap em nsnaps passos! MOD é uma intrínseca que devolve o resto de l na divisão por nsnap, ambos os argumentos sendo do mesmo tipo (l-INT(l/nsnap) * nsnap)
  csnap=csnap+1
  WRITE(num_snap,'(I3.3)')csnap
  WRITE(*,*) 'Imprimindo snap', num_snap, '...'
  OPEN(3,FILE='snap'//num_snap//'.bin',STATUS='replace', ACCESS='direct', FORM='unformatted', RECL=SGL*nx*ny)
  WRITE(3,REC=1) ((T2(i,j),i=1,ny), j=1,nx)
  CLOSE(3)
ENDIF
END SUBROUTINE Snap

END PROGRAM calor2d
