import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

##### Graficas punto 2 #####
#trayectoria 45 grados
proyectil = np.loadtxt('proyectil.txt')
x_45 = proyectil[:,1]
y_45 = proyectil[:,2]
plt.plot(x_45, y_45)
plt.ylim(0)
plt.title('Trayectoria proyectil a 45$^\circ$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.savefig('trayectoria_45.png')
plt.close()

#otras trayectorias
proyectil2 = np.loadtxt('proyectil2.txt')
angulos = [10, 20, 30, 40, 50, 60, 70]
indexes = np.where(proyectil2[:,0]==0)[0]
plt.figure()
for i in range(len(angulos)):
	i1 = indexes[i]
	if i != len(angulos)-1:
		i2 = indexes[i+1]
		plt.plot(proyectil2[i1:i2,1], proyectil2[i1:i2,2], label='${}^\circ$'.format(angulos[i]))
	else:
		plt.plot(proyectil2[i1:,1], proyectil2[i1:,2], label='${}^\circ$'.format(angulos[i]))
plt.ylim(0)
plt.title('Trayectoria proyectil a diferentes angulos')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend()
plt.savefig('trayectoria_varios_angulos.png')
plt.close()

##### Graficas punto 3 #####
#caso fixed boundaries
temp_fixed = np.loadtxt('temp_fixed_boundaries.txt')
mean_fixed = np.loadtxt('mean_fixed_boundaries.txt')
#Se convierte el arreglo 2D en uno 3D
a, b = temp_fixed.shape
temp_fixed = temp_fixed.reshape((a/(b-1), b-1, b))
print temp_fixed.shape
x_fixed = temp_fixed[0,:,0]
y_fixed = temp_fixed[0,:,0]
x, y = np.meshgrid(x_fixed, x_fixed)
#Graficas
#Temperatura promedio
plt.plot(mean_fixed[:,0], mean_fixed[:,1])
plt.show()
plt.close()
#condiciones iniciales
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_fixed[0,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
