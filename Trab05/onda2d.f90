program onda2D

implicit none

integer :: nx,nz,nt,snap_passos !parametros numericos
real    :: h,dt,lbd !parametros numericos
real    :: x,z,t !parametros fisicos
real    :: t0,td,ampl_fonte,fcorte,fc,zfonte,xfonte
integer :: nfonte,zi_fonte,xi_fonte
integer :: nsnap,count_snap=0 !numero de snapshots e contador
character(len=3) :: num_snap  !contador de snap para impress�o
character(len=30) :: vp_nome, in_nome
real, dimension(:,:), allocatable :: P1,P2,P3,vp !vetores com as temperaturas
real, dimension(:), allocatable :: fonte
integer :: i,j,n !contadores
logical :: existe_arq

call inicio() !Leitura de variaveis e calculos iniciais

call cria_fonte() !cria vetor de fonte

!inquire(file=vp_nome,exist=existe_arq)

!IF(existe_arq) then
!OPEN(10,FILE=vp_nome, STATUS='UNKNOWN',ACCESS='DIRECT', FORM='UNFORMATTED',RECL=4*nz*nx)
!READ(10,REC=1)((vp(i,j),i=1,nz),j=1,nx)
!ELSE
!  write(*,*)'Rodando na agua'
  vp=1500.0
!ENDIF

vp=vp*vp*(dt*dt)/(h*h)

write(*,*) "Inicio da marcha no tempo..."
do n=1,nt!marcha no tempo
   P1(zi_fonte,xi_fonte)=P1(zi_fonte,xi_fonte)-fonte(n)
   if(mod(n,50)==0) write(*,*) "passo",n !imprime na tela o avan�o no tempo de 50 em 50 passos
   do j=2,nx-1 !c�lculo do operador de DF
      do i=2,nz-1
         P3(i,j)=2*P2(i,j)-P1(i,j)+vp(i,j)*(P2(i-1,j)-2*P2(i,j)+P2(i+1,j)+P2(i,j-1)-2*P2(i,j)+P2(i,j+1))
      enddo
   enddo
   call imprime_snap() !imprime snapshot de temperatura
   P1(2:nz-1,2:nx-1)=P2(2:nz-1,2:nx-1) !atualiza as temperaturas para continuar a marcha
   P2(2:nz-1,2:nx-1)=P3(2:nz-1,2:nx-1)
enddo
write(*,*) "Fim da marcha no tempo."

contains !subrotinas internas

subroutine inicio()

implicit none

write(*,*) "Entre com o nome do arquivo de entrada (0 para 'in_onda.txt')"
read(*,*) in_nome
if(in_nome=="0") in_nome="in_onda.txt"

inquire(file=in_nome,exist=existe_arq)

if(existe_arq) then
    !lendo par�metros de arquivo
    open(20,file=in_nome)
    read(20,*) x,z,t,nsnap
    read(20,*) h,dt
    read(20,*) xfonte,zfonte
    read(20,*) ampl_fonte,fcorte
    read(20,*) vp_nome

    write(*,*) "Parametros lidos:"
    write(*,*) x,z,t,nsnap
    write(*,*) h,dt
    write(*,*) xfonte,zfonte
    write(*,*) ampl_fonte,fcorte
    write(*,*) vp_nome

!    lbd=alfa*dt/(h*h)

    nx=nint(x/h)+1
    nz=nint(z/h)+1
    nt=int(t/dt)+1

    !posição da fonte no grid
    xi_fonte=nint(xfonte/h)+1
    zi_fonte=nint(zfonte/h)+1

    !calculo do passo de tempo entre cada snap
    snap_passos=int(nt/nsnap)

    write(*,*) ""
    write(*,*) "Parametros numericos:"
    write(*,*) "nx,nz,nt = ",nx,nz,nt
!    write(*,*) "lambda",lbd
    write(*,*) ""

 !   if(lbd>0.5) then
 !       write(*,*)"Lambda maior que 0.5: ESQUEMA INSTAVEL!"
 !       write(*,*)"Escolha parametros adequado e rode novamente."
 !       stop
 !   endif

    write(*,*) ""

    allocate(P1(nz,nx),P2(nz,nx),P3(nz,nx),vp(nz,nx))

    P1=0.0
    P2=0.0
    P3=0.0

else
    write(*,*)'ERRO: arquivo de entrada ', trim(in_nome),' nao encontrado!'
    stop
endif

end subroutine

subroutine cria_fonte()

implicit none

real,parameter :: pi=3.141593

t0=2*SQRT(pi)/fcorte
fc=fcorte/(3*SQRT(pi))

allocate(fonte(nt))

fonte=0.0

nfonte=nint(2*t0/dt)
write(*,*) nfonte

do i=1,nfonte
    fonte(i)=ampl_fonte*(2*pi*(pi*fc*(i*dt-t0))**2-1.)*exp(-pi*(pi*fc*(i*dt-t0))**2)
enddo

open(20,file="fonte.txt")
do i=1,nt
    write(20,*)i*dt,fonte(i)
enddo

end subroutine cria_fonte

subroutine imprime_snap()

implicit none

if (mod(n,snap_passos)==0) then !imprime snapshots de nsnap em nsnap passos
   count_snap=count_snap+1
   write(num_snap,"(I3)")count_snap
   write(*,*) "Imprimindo snap ", num_snap,"..."
   open(10,file="snap"//trim(adjustl(num_snap))//".bin",status="replace",access="direct",form="unformatted",recl=4*nx*nz)
   write(10,rec=1) ((P2(i,j),i=1,nz),j=1,nx)
   close(10)
endif
end subroutine

 end program onda2D
