import matplotlib.pyplot as plt
import numpy as np
from math import pow, cos, sin, sqrt, pi
import scipy as sp
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
Nx=25
Ny=25
a=1
V0=1
Vn=1
t=1

kx=np.zeros(Nx)
ky=np.zeros(Ny)
for rx in range(0,Nx):
    kx[rx] = -pi/a + 2*pi*(rx)/(a*(Nx-1))
for ry in range(0,Ny):
    ky[ry] = -pi/a + 2*pi*(ry)/(a*(Ny-1))
    
#se construye la matriz que se usará para relacionar los indices de V(i,j) con sus valores dados
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

#se define la matriz i y j de dimension NxN=N**2
def i(rx,ry):
    return q[rx,ry]
def j(rx,ry):
    return q[rx,ry]




#relacionamos cada elemento de la matriz con sus indices
def find(element, matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == element:
                return (i, j)


#mostramos la matriz N**2 X N**2
def V(i,j):
    k1=find(i,q)
    k2=find(j,q)
    
    kx1=kx[k1[0]]
    ky1=ky[k1[1]]
    
    kx2=kx[k2[0]]
    ky2=ky[k2[1]]
    return V0+Vn*(cos((kx1-kx2)*a)+cos((ky1-ky2)*a))




#se muestra la función V en forma matricial
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



#mostramos la función K, que depende de w
def K(i,j,w):
    return (1/(w-2*E(j)))*V(i,j)





#se da a función K en forma matricial que depende de w
U = np.zeros((Nx**2, Ny**2))
def Km(w):
    for nx in range(0,Nx**2):
        for ny in range(0,Ny**2):
            U[nx][ny] = K(nx+1,ny+1,w)
    return U


#se define la matriz identidad
I = np.zeros((Nx**2,Ny**2))
for nx in range(0,Nx**2):
    I[nx][nx] = 1.0
    
    
#se define la matriz T
def T(w):
    K2 = (I-Km(w))
    K3 = np.linalg.inv(K2) #matriz inversa
    K4 = np.dot(K3, Vm) #producto matricial
    return K4

#sea
g=2.0

# Set up grid and test data
x = range(Nx*Nx)
y = range(Ny*Ny)

data = T(g)

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
#ha.plot_surface(X, Y, data)

mycmap = plt.get_cmap('gist_earth')
ha.set_title('T( 2.0) matrix')
surf1 = ha.plot_surface(X, Y, data, cmap=mycmap)
hf.colorbar(surf1, ax=ha, shrink=0.5, aspect=5)


ha.set_xlabel('i')

ha.set_ylabel('j')

ha.set_zlabel('T(w= 2.0)')

plt.show()





