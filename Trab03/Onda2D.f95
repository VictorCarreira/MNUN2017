PROGRAM Onda2D
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !Aluno: Victor Ribeiro Carreira
  !Professor: Leandro Di Bartolo
  !Este programa visa cumprir os requisitos da disciplina MNUM,
  !simulando uma Onda 2D acústica variando no tempo.
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Declaração de Variáveis Globais
IMPLICIT NONE
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=6, r=10)
INTEGER(KIND=SGL)::i, j, n, nx, nz, nt, x, z, snap_passos, xifonte, zifonte, nsnap, csnap=0
INTEGER(KIND=SGL):: ncer
REAL(KIND=SGL)::DeltaT, DeltaX, DeltaZ, inicial, final, custocomputacional, dt
REAL(KIND=SGL)::  xfonte, zfonte, t,fcorte, ampl_fonte, fcer, h,c_const
REAL(KIND=SGL)::alfa2, alfa4, beta
REAL(KIND=SGL), PARAMETER::pi=3.1416
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:):: fonte
REAL(KIND=SGL),ALLOCATABLE, DIMENSION(:,:):: P1, P2, P3, c
REAL(KIND=SGL), DIMENSION(:), ALLOCATABLE:: g
CHARACTER(LEN=30)::file, vp_nome
CHARACTER(LEN=3):: num_snap
LOGICAL:: existe_arq

!########################################################################################
CALL cpu_time(inicial)

CALL Entrada()

 c=c*c*(dt*dt)/(h*h)

CALL Dispersao(c,alfa2,alfa4,fcorte,h)

CALL Estabilidade(h,beta,c,DeltaT)

CALL Explosiva(fcorte, ampl_fonte, t, nt, fonte)

!Cálculo das DF (Stencil)
DO n=1,nt
  P1(zifonte,xifonte)=P1(zifonte,xifonte)-fonte(n)
  IF(mod(n,50)==0) WRITE(*,*) "passo",n !imprime na tela o avan�o no tempo de 50 em 50 passos
    DO j=2,nx-1!dimensão lenta
      DO i=2,nz-1!dimensão rápida
       !P3(i,j,k)=[(DeltaT**2 * c(i,j)**2)/12*h]*[-(P2(i-2,j,k)+P2(i,j-2,k)+P2(i+2,j,k)+P2(i,j+2,k)) + 16*(P2(i-1,j,k)+P2(i,j-1,k)+P2(i+1,j,k)+P2(i,j+1,k))-60*P2(i,j+1,k)] + 2* P2(i,j,k) - P2(i,j,k-1) + (DeltaT**2) * [c(i,j)**2] * rho(i,j) * s(i,j,k) !Quarta ordem
        P3(i,j)=2*P2(i,j)-P1(i,j)+c(i,j)*(P2(i-1,j)-2*P2(i,j)+P2(i+1,j)+P2(i,j-1)-2*P2(i,j)+P2(i,j+1))
      ENDDO
    CALL Oneway(nx, nz, nt, x, z, t, c, P2, P3)
    CALL Cerjan(g,ncer,P2,P3)
   ENDDO
  CALL Snap()
  P1(2:nz-1,2:nx-1)=P2(2:nz-1,2:nx-1) !atualiza as temperaturas para continuar a marcha
  P2(2:nz-1,2:nx-1)=P3(2:nz-1,2:nx-1)
ENDDO

OPEN(5,file="P3.txt")
WRITE(5,*)P3

CALL cpu_time(final)
custocomputacional=final-inicial
PRINT*, 'Custo Computacional=',custocomputacional, 'segundos'
PRINT*,'*********************** FIM ***************************'

!################################################################################################################################################


CONTAINS
SUBROUTINE Explosiva(fcorte, ampl_fonte, t, nt, fonte)
  IMPLICIT NONE
  REAL(KIND=SGL),INTENT(IN)::fcorte, ampl_fonte, t
  INTEGER(KIND=SGL), INTENT(IN)::  nt
  INTEGER(KIND=SGL)::nfonte
  REAL(KIND=SGL), DIMENSION(:), ALLOCATABLE, INTENT(OUT):: fonte
  REAL(KIND=SGL):: t0, dt, fc
  REAL(KIND=SGL), PARAMETER::pi=3.141593


  t0 = 2*SQRT(pi)/fcorte
  fc = fcorte/3*SQRT(pi)

  ALLOCATE(fonte(nt))

  fonte=0.0
  nfonte=NINT(2*t0/dt)

  DO i=1, nfonte
    fonte(i)=ampl_fonte*(2*pi*(pi*fc*(i*dt-t0))**2-1.0)*EXP(-pi*(pi*fc*(i*dt-t0))**2)
  ENDDO

  OPEN(20,file="fonte.txt")
  DO i=1,nt
    WRITE(20,*)i*dt,fonte(i)
  ENDDO

ENDSUBROUTINE Explosiva

!------------------------------------------------------------------------------------
SUBROUTINE Entrada() ! Leitura dos dados de entrada


WRITE(*,*) "Entre com os arquivos de entrada. 1,  para entrada 'input.txt' "
READ(*,*) file
IF(file=="1") file="input.txt"
INQUIRE(file=file, exist=existe_arq)
IF(existe_arq) THEN



  OPEN(UNIT=2,FILE='input.txt')
  READ(2,*) x     ! comprimento Stencil na direção x (cm)
  READ(2,*) z     ! comprimento Stencil na direção z (cm);
  READ(2,*) t    ! Iteracoes temporais (s);
  READ(2,*) c_const,vp_nome ! Velocidade vp
  READ(2,*) h     ! Espaçamento entre os pontos
  READ(2,*) dt    ! discretizacao temporal (tempo fisico do experimento)
  !DADOS DE ENTRADA DA FONTE
  READ(2,*) xfonte ! Posição da fonte na direção horizontal
  READ(2,*) zfonte ! Posição da Fonte em profundidade
  READ(2,*) ampl_fonte ! Amplitude da fonte
  READ(2,*) fcorte ! Frequencia de corte
  READ(2,*) fcer,ncer  ! Fator de atenuação
 ! Critérios de Estabilidade e Dispersao
  READ(2,*) alfa2  ! Alfa de segunda ordem
  READ(2,*) alfa4  ! Alfa de quarta ordem
  READ(2,*) beta  ! Número de iterações para levar de 1 para 2
! Número de arquivos de saída
  READ(2,*) nsnap ! número de snapshots

  !**********************************************
  WRITE(*,*)'Parametos lidos:'
  WRITE(*,*) x     ! comprimento Stencil na direção x (cm)
  WRITE(*,*) z     ! comprimento Stencil na direção z (cm);
  WRITE(*,*) t    ! Iteracoes temporais (s);
  WRITE(*,*) c_const,vp_nome  ! Velocidade de onda p
  WRITE(*,*) h     ! Espaçamento entre os pontos
  WRITE(*,*) dt    ! discretizacao temporal (tempo fisico do experimento)
  !DADOS DE ENTRADA DA FONTE
  WRITE(*,*) xfonte ! Posição da fonte na direção horizontal
  WRITE(*,*) zfonte ! Posição da Fonte em profundidade
  WRITE(*,*) ampl_fonte ! Amplitude da fonte
  WRITE(*,*) fcorte ! Frequencia de corte
  WRITE(*,*) fcer,ncer ! Fator de atenuação
  ! Critérios de Estabilidade e Dispersao
  WRITE(*,*) alfa2  ! Alfa de segunda ordem
  WRITE(*,*) alfa4  ! Alfa de quarta ordem
  WRITE(*,*) beta  ! temperatura na borda direita (°C)

  ! Número de arquivos de saída
  WRITE(2,*) nsnap ! número de snapshots
        !c=1500
        !c=c**2 ! Escrever a eq da Onda2D

      !discretizacao da malha
      nx=nint(x/h)+1
      nz=nint(z/h)+1
      nt=int(t/dt)+1

      !posição da fonte na malha
      xifonte=nint(xfonte/h)+1
      zifonte=nint(zfonte/h)+1

      !calculo do passo de tempo entre cada snap
      snap_passos=int(nt/nsnap)

      WRITE(*,*) ""
      WRITE(*,*) "Parametros numericos:"
      WRITE(*,*) "nx,nz,nt = ",nx,nz,nt
      WRITE(*,*) ""

      ALLOCATE(P1(nz,nx),P2(nz,nx),P3(nz,nx),c(nz,nx))


  ELSE
      WRITE(*,*)'ERRO: arquivo de entrada ', trim(file),' não encontrado!'
      STOP
  ENDIF


 !Condição inicial. Zerando variáveis
        P1=0.0
        P2=0.0
        P3=0.0

     !função de Cerjan
      g=Fcerj(g)
ENDSUBROUTINE Entrada

!------------------------------------------------------------------------------------


SUBROUTINE Oneway(nx, nz, nt, x, z, t, c, P2, P3)
  IMPLICIT NONE
  INTEGER(KIND=SGL),INTENT(IN)::nx,nz,nt,x,Z
  INTEGER(KIND=SGL)::i,j
  REAL(KIND=SGL), INTENT(IN)::t
  REAL(KIND=SGL), DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P2, P3, c

!ALLOCATE(P2(nz,nx),P3(nz,nx), c(nz,nx))


  DeltaX=x/nx
  DeltaZ=z/nz
  DeltaT=t/nt

!Condição de contorno (Oneway):

!Borda Direita
j=Nx !indices i e j trocados em todas. conferir sinais e usar sqrt(c)
DO i=1,nz
  P3(i,j)= + (c(i,j)*DeltaT/DeltaX) * (P2(i,j)-P2(i-1,j)) +P2(i,j)
ENDDO

!Borda Esquerda
j=1
DO i=1,nz
  P3(i,j)= - (c(i,j)*DeltaT/DeltaX) * (P2(i,j)-P2(i-1,j)) +P2(i,j)
ENDDO

!Fundo
i=nz
DO j=1,nx
  P3(i,j)= + (c(i,j)*DeltaT/DeltaZ) * (P2(i,j)-P2(i,j-1)) +P2(i,j)
ENDDO

ENDSUBROUTINE Oneway

!-----------------------------------------------------------------------------------

SUBROUTINE Cerjan(g,ncer,P2,P3)
   IMPLICIT NONE
    INTEGER(KIND=SGL)::j
    INTEGER(KIND=SGL),INTENT(IN):: ncer
    REAL(KIND=SGL),DIMENSION(:),INTENT(IN)::g
    REAL(KIND=SGL),DIMENSION(:,:),INTENT(INOUT)::P2, P3

    !Borda esquerda
    DO j= 1, ncer-1
      DO i= 1, nz
        P2(i,j)=g(i)*P2(i,j)
        P3(i,j)=g(i)*P3(i,j)
      ENDDO
    ENDDO


ENDSUBROUTINE Cerjan

!-----------------------------------------------------------------------------------


SUBROUTINE Dispersao(c, alfa2, alfa4, fcorte, h)
  IMPLICIT NONE
  REAL(KIND=SGL),INTENT(IN), ALLOCATABLE, DIMENSION(:,:)::c
  REAL(KIND=SGL),INTENT(IN)::alfa2, alfa4, fcorte
  REAL(KIND=SGL),INTENT(INOUT)::h


  !Critério da Não-dispersão de segunda ordem:
  IF (h .LE. (MIN(nz,nx)/alfa2*fcorte)) THEN
      PRINT*, 'Não-dispersa'
  ELSE
      PRINT*,'Dispersa'
      STOP "Checar o critério de dispersão de segunda ordem"
  ENDIF

  !Critério da Não-dispersão de quarta ordem:
  IF (h .LE. (MIN(nz,nx)/alfa4*fcorte)) THEN
      PRINT*, 'Não-dispersa'
  ELSE
      PRINT*,'Dispersa'
      STOP"Checar o critério de dispersão de quarta ordem"
  ENDIF

ENDSUBROUTINE Dispersao


!----------------------------------------------------------------------------------

SUBROUTINE Estabilidade(h, beta, c, DeltaT)
  IMPLICIT NONE
  REAL(KIND=SGL),INTENT(IN), ALLOCATABLE, DIMENSION(:,:)::c
  REAL(KIND=SGL),INTENT(IN)::beta, h
  REAL(KIND=SGL),INTENT(OUT)::DeltaT

  !Critério de Estabilidade:
  IF (DeltaT .LE. (h/beta*MAX(nz,nx))) THEN
      PRINT*,'ESTÁVEL'
  ELSE
      PRINT*,'INSTÁVEL'
      STOP"Checar estabilidade numérica"
  ENDIF

ENDSUBROUTINE Estabilidade

!-----------------------------------------------------------------------------------

SUBROUTINE Snap()

  !calculo do passo de tempo entre cada snap
  snap_passos=int(nt/nsnap)

IF(MOD(n,snap_passos)==0)THEN !imprime os snapsshots de nsnap em nsnaps passos! MOD é uma intrínseca que devolve o resto de l na divisão por nsnap, ambos os argumentos sendo do mesmo tipo (l-INT(l/nsnap) * nsnap)
  csnap=csnap+1
  WRITE(num_snap,'(I3.3)')csnap
  WRITE(*,*) 'Imprimindo snap', num_snap, '...'
  OPEN(3,FILE='snap'//num_snap//'.bin',STATUS='replace', ACCESS='direct', FORM='unformatted', RECL=SGL*nx*nz)
  WRITE(3,REC=1) ((P2(i,j),i=1,nz), j=1,nx)
  CLOSE(3)
ENDIF

ENDSUBROUTINE Snap

!---------------------------------------------------------------------------------------

SUBROUTINE Modelo()

  INQUIRE(file=vp_nome,exist=existe_arq)

  IF(existe_arq) then
  OPEN(10,FILE=vp_nome, STATUS='UNKNOWN',ACCESS='DIRECT', FORM='UNFORMATTED',RECL=4*nz*nx)
  READ(10,REC=1)((c(i,j),i=1,nz),j=1,nx)
  ELSE
    WRITE(*,*)'Rodando na agua'
    c=1500.0
  ENDIF


ENDSUBROUTINE Modelo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNÇÕES UTILIZADAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL FUNCTION Fcerj(g)
     IMPLICIT NONE
     REAL(KIND=SGL), DIMENSION(:), ALLOCATABLE:: g
      REAL(KIND=SGL)::fat
      INTEGER(KIND=SGL), PARAMETER::ncer=100
      INTEGER(KIND=SGL)::i
       ALLOCATE(g(ncer-1))
        !Cálculo do função de atenuação de Cerjan:
      DO i=1,ncer-1
        !Fcerj(i)=EXP(-(fat*(ncerj-i))**2)
        g(i)=EXP(-(fat*(ncer-i))**2)!OLHAR NO EXCEL
      ENDDO

  RETURN
END FUNCTION Fcerj


END PROGRAM Onda2D
