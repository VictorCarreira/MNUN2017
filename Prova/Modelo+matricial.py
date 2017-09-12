
# coding: utf-8

# # Transformando imagem em um arquivo de dados

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# In[8]:

############Pacotes################
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import misc
from PIL import Image, ImageFilter
###################################


# A função misc serve para transformar qualquer imagem em um arquivo de geométrica de dados que carece apenas de propriedade física.

# In[28]:

ma = misc.imread('Modelo.png',mode='P') #carrega o modelo da hipótese, em imagem, e transforma em uma matriz de dados. Função misc.imread()
print(ma.shape) #mostra a matriz

# no 'mode' usar:
#* 'L' (8-bit pixels, black and white)
#* 'P' (8-bit pixels, mapped to any other mode using a color palette)
#* 'RGB' (3x8-bit pixels, true color)
#* 'RGBA' (4x8-bit pixels, true color with transparency mask)
#* 'CMYK' (4x8-bit pixels, color separation)
#* 'YCbCr' (3x8-bit pixels, color video format)
#* 'I' (32-bit signed integer pixels)
#* 'F' (32-bit floating point pixels)

plt.matshow(ma,vmin=0.0,vmax=2.0)#Cria o gráfico com a matriz da imagem (cria um valor mínimo e máximo)
plt.colorbar(mappable=None, cax=None, ax=None)#Define o tipo de cor da barra de dados
plt.savefig('modelo.pdf')
#plt.grid(color='r', linestyle='-', linewidth=2)
plt.show()#mostra o gráfico



# In[10]:

#Analisando as dimensões da matriz de dados:
tamanho=len(ma)
ma-np.array(ma)
print(ma)


# In[11]:

#Atribuindo propriedade física as camadas do modelo

#Vel = [[0.0]*np.shape(ma)[1]]*np.shape(ma)[0]
#Vel = np.array(vel,float)

#for i in range(np.shape(ma)[0]-1):
#    for j in range(np.shape(ma)[1]-1):
#        if ma[i,j] == 0:
#            Vel[i,j] = 2000
#        if ma[i,j] == 1:
#            Vel[i,j] = 1500
#        if ma[i,j] == 2:
#            Vel[i,j] = 2100


#print(np.shape(Vel))

Vel = [[0]*np.shape(ma)[1]]*np.shape(ma)[0]
Vel = np.array(Vel,float)

for i in range(np.shape(ma)[0]-1):
    for j in range(np.shape(ma)[1]-1):
        if ma[i,j] == 0:
            Vel[i,j] = 2000
        if ma[i,j] == 1:
            Vel[i,j] = 1500
        if ma[i,j] == 2:
            Vel[i,j] = 2100



print(np.shape(Vel))
print(Vel[0,0])
plt.imshow(Vel)
plt.colorbar()
plt.show()



# In[12]:

#Salvando a matriz de dados:
np.savetxt('Modelo.txt', Vel)
plt.close()


# ___________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

# In[ ]:
