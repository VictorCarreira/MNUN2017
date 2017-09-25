module modulo_propag

integer :: Nx,Nz,Nt,Ns,dt_sismo,dt_snap
integer :: ixf,izf,iprof,cerjan_n

real    :: Dt,h,fcorte,cerjan_fat,famp

character(len=100) :: dir_entrada,dir_saida
character(len=50)  :: vel_arq,secao_arq, mtt_arq !nome dos arquivos de velocidade
character(len=3)   :: varsnap !nome dos arquivos de snapshot

real,dimension(:,:),allocatable  :: vel,sismo,p1,p2
real, dimension(:), allocatable :: g

! *****************************************************
CONTAINS
! *****************************************************

subroutine leitura()

write(*,*) "Entre com os parâmetros de entrada pedidos ou rode novamente com './proag <  script' (i.e., utilise um script &
   & direcionado a antrada padrão para ele):"

write(*,*) "Número de pontos Nx,Nz,Nt"
read(*,*) Nx,Nz,Nt
write(*,*) "Número o espaçamento e o passo de tempo"
read(*,*) h,Dt
write(*,*) "Amplitude e frequencia de corte da fonte"
read(*,*) famp,fcorte
write(*,*) "Posição da fonte (x,z)"
read(*,*) ixf,izf
write(*,*) "Passos de tempo entre snapshots"
read(*,*) dt_snap
write(*,*) "Número de pontos para pular e salvar o sismograma reduzido (dt_sismo) e profundidade"
read(*,*) dt_sismo, iprof
write(*,*) "Número de pontos da camada de Cerjan"
read(*,*) cerjan_n,cerjan_fat
write(*,*) "Nome do modelo de velocidades"
read(*,*) vel_arq
write(*,*) "Nome da seção sísmica (empilhada)"
read(*,*) secao_arq
write(*,*) "Nome da matriz de tempos de transito"
read(*,*) mtt_arq
write(*,*) "Diretório de entrada"
read(*,*) dir_entrada
write(*,*) "Diretório de saida"
read(*,*) dir_saida

Ns=Nt/dt_sismo
!close(10)

end subroutine leitura

! *****************************************************

subroutine alocacao()

allocate(vel(Nz,Nx),sismo(Ns,Nx),p1(Nz,Nx),p2(Nz,Nx))

allocate(g(cerjan_n))

end subroutine alocacao


end module modulo_propag

