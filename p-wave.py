import matplotlib.pyplot as plt
import numpy as np
from math import pow, cos, sin, sqrt, pi
import scipy as sp
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
#se establece el tamaño de las arrays de los momentos y los valores de las constantes
Nx=25
Ny=25
a=1
V0=1
Vn=1
t=1
#se definen las arrays de momentos kx y ky
kx=np.zeros(Nx)
ky=np.zeros(Ny)
for rx in range(0,Nx):
    kx[rx] = -pi/a + 2*pi*(rx)/(a*(Nx-1))
for ry in range(0,Ny):
    ky[ry] = -pi/a + 2*pi*(ry)/(a*(Ny-1))
    
#se construye la matriz que se usará para relacionar los indices de V(i,j) con sus valores dados
#se obtiene una matriz q de dimensión Nx X Ny enumerada desde el 1 al Nx*Ny, a esta matriz nos referiremos cuando queramos usar las indices del operador V(i,j)
p=np.zeros(Nx*Ny)
b=0
for rx in range(0,Nx):
    for ry in range(0,Ny):
        p[b]=b+1
        b=b+1
q = np.zeros((Nx,Ny))
c=0
for nx in range(0,Nx):
    for ny in range(0,Ny):
        q[nx][ny] = p[c]
        c=c+1

#se define las funciones i y j de dimension NxN=N**2 que dan las componentes de la matriz q al darle los indices que tiene esa componente
def i(rx,ry):
    return q[rx,ry]
def j(rx,ry):
    return q[rx,ry]




#relacionamos cada elemento de la matriz con sus indices
#esta parte pide una matriz y un elemento, si encuentra tal elemento en esa matriz da como resultado los indices de su ubicacion (fila, columna)
def find(element, matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == element:
                return (i, j)


#mostramos la función V que solo depende de i y j, buscamos tales indices en la matriz q y relacionamos tales indices con los de los momentos kx y ky
def V(i,j):
    k1=find(i,q)
    k2=find(j,q)
    
    kx1=kx[k1[0]]
    ky1=ky[k1[1]]
    
    kx2=kx[k2[0]]
    ky2=ky[k2[1]]
    return V0+Vn*(cos((kx1-kx2)*a)+cos((ky1-ky2)*a))




#se muestra la función V en forma matricial para poder graficarla
Vm = np.zeros((Nx**2,Ny**2))
for m in range(0,Nx*Nx):
    for n in range(0,Ny*Ny):
        Vm[m][n] = V(m+1,n+1)
        
        
#mostramos la función de energía cinética E
def E(i):
    k1=find(i,q)
    
    kx1=kx[k1[0]]
    ky1=ky[k1[1]]
    return -2*t*(cos(kx1*a)+cos(ky1*a))



#mostramos la función K, que depende de la energía w
def K(i,j,w):
    return (1/(w-2*E(j)))*V(i,j)






U = np.zeros((Nx**2, Ny**2))
#se da a función K en forma matricial que depende de w
def Km(w):
    for nx in range(0,Nx**2):
        for ny in range(0,Ny**2):
            U[nx][ny] = K(nx+1,ny+1,w)
    return U


#se define la matriz identidad de dimensión Nx**2 X Ny**2
I = np.zeros((Nx**2,Ny**2))
for nx in range(0,Nx**2):
    I[nx][nx] = 1.0
    
    
#se define la matriz T como una función, usando la definición T = (I-K)^{-1}*V
def T(w):
    K2 = (I-Km(w))
    K3 = np.linalg.inv(K2) #matriz inversa
    K4 = np.dot(K3, Vm) #producto matricial
    return K4

#sea un valor dado para w = g, para este caso indiqué w=2.0
g=2.0

# hagamos un enmallado que mida Nx**2 y Ny**2 para poder graficar los componentes de la matriz T(w)
x = range(Nx*Nx)
y = range(Ny*Ny)

data = T(g)

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x, y)  # `plot_surface` espera `x` e `y` 


mycmap = plt.get_cmap('gist_earth')
ha.set_title('T( 2.0) matrix') #titulo de la gráfica
surf1 = ha.plot_surface(X, Y, data, cmap=mycmap) #tipo de grafica que se ofrecerá
hf.colorbar(surf1, ax=ha, shrink=0.5, aspect=5) # colores que tendrá la gráfica segun los valores que ofrezca

#nombres para los ejes
ha.set_xlabel('i')

ha.set_ylabel('j')

ha.set_zlabel('T(w= 2.0)')
#se ejecuta la gráfica
plt.show()


#tiempo de espera: en mi caso fue aprox. 2min para un arreglo inicial de 25x25 y final de 625x625


