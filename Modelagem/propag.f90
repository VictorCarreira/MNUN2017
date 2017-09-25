Program MDF_acustico
! **************************************************
!  PROGRAMA EXEMPLO para gerar seção empilhada simples
!	pelo modelo do refletor explosivo
! **************************************************

use modulo_propag

implicit none

integer :: k_snap,k_sismo,k1,zfonte
real    :: pi,nf,tf,fc,cput2,cput1,fon_d2g,aux

parameter(pi=4*atan(1.0)) !Pi

integer :: i,j,k,l      !variaveis auxiliares
! *****************************************************

call leitura()

call alocacao()

write(*,*) 'Ns = ',Ns

!Tempo de processamento
 call CPU_TIME(cput1)

!CÁLCULOS INICIAIS
!Cálculos para o termo fonte
Nf = 4*sqrt(pi)/(fcorte*DT)
Tf = 2*sqrt(pi)/fcorte
fc = fcorte/(3.*sqrt(pi))

!Cálculos dos fatores de amortecimento
do i=1,cerjan_n
  g(i)=exp( - ( cerjan_fat*(cerjan_n-i) )*( cerjan_fat*(cerjan_n-i) )  )
enddo

write(*,*) "Duracao da aplicacao da fonte: Nf=",Nf
write(*,*) "Numero de amostras temporais do sismograma Ns = ",Ns

!Leitura do grid de velocidades
write(*,*) ' '
write(*,*) 'Salvando modelo de velocidade...'
write(*,*) ' '

open(10,file=trim(ADJUSTL(dir_entrada))//vel_arq,status="unknown",access="direct",&
& form="unformatted",recl=4*Nx*Nz) !Abre o arquivo
read(10,rec=1)((vel(i,j),i=1,Nz),j=1,Nx) !Armazena o modelo de velocidade em vel
close(10)


vel=vel/2.0


!Inicializacao dos campos (condição inicial)
do j=1,Nx
  do i=1,Nz
    p2(i,j)=0.0
    p1(i,j)=0.0
  enddo
enddo

k_snap = 0 !Contador para snapshots
k_sismo = 0 !Contador de amostras do sismograma

write(*,*) "Inicio"

vel=vel*Dt/h
vel=vel*vel
vel=vel/12.0

do k=1,Nt !Inicio do loop da marcha no tempo
    
  !imprime o andamento da marcha no tempo na tela
  if(mod(k,50).eq.0) write(*,*)"tempo = ",k

!Aplicação das fontes do modelo do refletor explosivo
  do j=1,Nx
    if (j>= 270 .and. j<= 390) then
      zfonte =nint( 100 + sqrt(3600. - (j - 330.)**2))
    else
     zfonte = 100
    endif
    if(k.le.Nf) p1(zfonte,j) = p1(zfonte,j) - fon_d2g((k-1)*Dt-Tf,fc,famp)
  enddo
  
    do j=1,Nx
    if (j>= 151 .and. j<= 601) then
      zfonte =nint( (131./450.)*j + 56714./450.)
    else
     zfonte = 170
    endif
    if(k.le.Nf) p1(zfonte,j) = p1(zfonte,j) - fon_d2g((k-1)*Dt-Tf,fc,famp)
  enddo
  
    do j=255,Nx
     zfonte = 200
    if(k.le.Nf) p1(zfonte,j) = p1(zfonte,j) - fon_d2g((k-1)*Dt-Tf,fc,famp)
  enddo
  
  
  !Loop do operador do MDF (envolvendo toda a malha, exceto bordas)
  !Operador de quarta ordem espacial
  do j=3,Nx-2	
    do i=3,Nz-2
      p1(i,j) = vel(i,j)*											&

      (	-( p2(i+2,j)+p2(i-2,j)+p2(i,j+2)+p2(i,j-2) )+16*( p2(i+1,j)+p2(i-1,j)+p2(i,j+1)+p2(i,j-1) ) - 60*p2(i,j)	&

      )															&

        +   2*p2(i,j) - p1(i,j)

    enddo
  enddo

  !IMPLEMENTAR A SEGUNDA ORDEM NOS PONTOS PRÓXIMOS À BORDA
  call borda_2O()
  
  !Aplicacao das condicoes de contorno nao-reflexivas nas bordas
  call CCNR_Reinolds()
    
  !Aplicacao de zona de amortecimento próximo às bordas
  call AB_Cerjan()
    
  !Salva de amostra do sismograma  cada "dt_sismo" passos de tempo           
  if(mod(k,dt_sismo).eq.0) then
    k_sismo=k_sismo + 1
    do i=1,Nx
      sismo(k,i)=p1(iprof,i)
    enddo
  endif

  !Impressao dos snapshots a cada "dt_snap" intervalos de tempo           
  if(mod(k,dt_snap).eq.0) then
    write(varsnap,'(i3)')k_snap !'(i3.3)'
    write(*,*)"snapshot = ",trim(varsnap)
    !Impressao do campo de pressoes em binario
    open(10,file=trim(ADJUSTL(dir_saida))//"snap"//trim(ADJUSTL(varsnap))//".bin",status="replace",access="direct",&
    & form="unformatted",recl=4*Nx*Nz) !Abre o arquivo
    write(10,rec=1)((p1(i,j),i=1,Nz),j=1,Nx)
    close(10)
    k_snap=k_snap + 1
  endif	
  
 
  !Atualizacao dos campos para o próximo passo da marcha no tempo
  do j=1,Nx
    do i=1,Nz
      aux=p2(i,j)
      p2(i,j)=p1(i,j)
      p1(i,j)=aux
    enddo
  enddo
  
enddo !Fim do loop da marcha no tempo

!Impressao do sismograma em arquivo binário
write(*,*) ' ' ; write(*,*) 'Salvando sismograma...' ; write(*,*) ' '
open(10,file=trim(ADJUSTL(dir_saida))//"sismo.bin",status="replace",access="direct", &
& form="unformatted",recl=4*Ns*Nx)
write(10,rec=1)((sismo(j,i),j=1,Ns),i=1,Nx)
close(10)

!Calcula o tempo de processamento e imprime na tela
call CPU_TIME(cput2)
WRITE(*,*) "Tempo de execucao ", cput2-cput1, "segundos"

! **************************************************
! SUBROTINAS INTERNAS
! **************************************************
  CONTAINS

! **************************************************
    
  subroutine CCNR_Reinolds()
      
  implicit none
  
  !Superior
  do j=2,Nx-1
    p1(1,j)=0.0 !Condição de superfície livre
  enddo
  
  !Esquerda
  do i=1,Nz
    p1(i,1)=p2(i,1)+sqrt(vel(i,1)*12.0)*(p2(i,2)-p2(i,1))
  enddo

  !Direita
  do i=1,Nz
    p1(i,Nx)=p2(i,Nx)-sqrt(vel(i,Nx)*12.0)*(p2(i,Nx)-p2(i,Nx-1))
  enddo

  !Inferior
  do j=2,Nx-1
    p1(Nz,j)=p2(Nz,j)-sqrt(vel(Nz,j)*12.0)*(p2(Nz,j)-p2(Nz-1,j))
  enddo
  
  end subroutine CCNR_Reinolds
    
! **************************************************
    
  subroutine AB_Cerjan()
    
  implicit none
    
  ! ATENUAÇÃO DO CAMPO DE PESSÕES ACÚSTICAS PASSADOS
  !esquerda
  do j=1,cerjan_n-1!
    do i=1,Nz
      p1(i,j)=g(j)*p1(i,j)
      p2(i,j)=g(j)*p2(i,j)
    enddo
  enddo

  !direita
  k1=cerjan_n
  do j=Nx-cerjan_n+2,Nx
    k1=k1-1
    do i=1,Nz
      p1(i,j)=g(k1)*p1(i,j)
      p2(i,j)=g(k1)*p2(i,j)
    enddo
  enddo
  
  !inferior
  k1=cerjan_n
  do i=Nz-cerjan_n+2,Nz
    k1=k1-1
    do j=1,Nx
      p1(i,j)=g(k1)*p1(i,j)
      p2(i,j)=g(k1)*p2(i,j)
    enddo
  enddo
  
  end subroutine AB_Cerjan

  subroutine borda_2O()
    
  implicit none
  
 !Superior
  do j=2,Nx-1
    p1(2,j) = vel(2,j)*12.0*			&

      (	p2(2+1,j)-2*p2(2,j)+p2(2-1,j) +		&

      p2(2,j+1)-2*p2(2,j)+p2(2,j-1) 		&

      )						&

        +   2*p2(2,j) - p1(2,j)
  enddo
  
  !Esquerda
  do i=1,Nz
    p1(i,2) = vel(i,2)*12.0*			&

      (	p2(i+1,2)-2*p2(i,2)+p2(i-1,2) +		&

      p2(i,2+1)-2*p2(i,2)+p2(i,2-1) 		&

      )						&

        +   2*p2(i,2) - p1(i,2)
  enddo

  !Direita
  do i=1,Nz
    p1(i,Nx-1) = vel(i,j)*12.0*			&

      (	p2(i+1,Nx-1)-2*p2(i,Nx-1)+p2(i-1,Nx-1) +&

      p2(i,Nx-1+1)-2*p2(i,Nx-1)+p2(i,Nx-1-1) 	&

      )						&

        +   2*p2(i,Nx-1) - p1(i,Nx-1)
  enddo

  !Inferior
  do j=2,Nx-1
    p1(Nz-1,j) = vel(Nz-1,j)*12.0*		&

      (	p2(Nz-1+1,j)-2*p2(Nz-1,j)+p2(Nz-1-1,j) +&

      p2(Nz-1,j+1)-2*p2(Nz-1,j)+p2(Nz-1,j-1)	&

      )						&

        +   2*p2(Nz-1,j) - p1(Nz-1,j)
  enddo
   
  end subroutine borda_2O
  
end program MDF_acustico
