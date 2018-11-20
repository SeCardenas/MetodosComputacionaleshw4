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
plt.xlabel('$x [m]$')
plt.ylabel('$y [m]$')
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
plt.xlabel('$x [m]$')
plt.ylabel('$y [m]$')
plt.legend()
plt.savefig('trayectoria_varios_angulos.png')
plt.close()

##### Graficas punto 3 #####
#caso fixed boundaries
temp_fixed = np.loadtxt('temp_fixed_boundaries.txt')
mean_fixed = np.loadtxt('mean_fixed_boundaries.txt')
#Se convierte el arreglo 2D en uno 3D
a, b = temp_fixed.shape
print temp_fixed.shape
temp_fixed = temp_fixed.reshape((a/(b-1), b-1, b))
print temp_fixed.shape
x_fixed = temp_fixed[0,:,0]
y_fixed = temp_fixed[0,:,0]
x, y = np.meshgrid(x_fixed, x_fixed)
#Graficas
#Temperatura promedio
plt.plot(mean_fixed[:,0], mean_fixed[:,1])
plt.title('Temperatura promedio fronteras fijas')
plt.xlabel('$t [s]$')
plt.ylabel('$T\ promedio\ [^\circ C]$')
plt.savefig('mean_fixed_boundaries.png')
plt.close()
#condiciones iniciales
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_fixed[0,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras fijas t=0')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_fixed_1.png')
plt.close()
#estados intermedios
a = temp_fixed.shape[0]
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_fixed[1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras fijas t={}s'.format(mean_fixed[125,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_fixed_2.png')
plt.close()

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_fixed[2,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras fijas t={}s'.format(mean_fixed[250,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_fixed_3.png')
plt.close()

#estado final
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_fixed[-1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras fijas t={}s'.format(mean_fixed[-1,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_fixed_4.png')
plt.close()


#caso periodic boundaries
temp_periodic = np.loadtxt('temp_periodic_boundaries.txt')
mean_periodic = np.loadtxt('mean_periodic_boundaries.txt')
#Se convierte el arreglo 2D en uno 3D
a, b = temp_periodic.shape
print temp_periodic.shape
temp_periodic = temp_periodic.reshape((a/(b-1), b-1, b))
print temp_periodic.shape
x_periodic = temp_periodic[0,:,0]
y_periodic = temp_periodic[0,:,0]
x, y = np.meshgrid(x_periodic, x_periodic)
#Graficas
#Temperatura promedio
plt.plot(mean_periodic[:,0], mean_periodic[:,1])
plt.title('Temperatura promedio fronteras periodicas')
plt.xlabel('$t [s]$')
plt.ylabel('$T\ promedio\ [^\circ C]$')
plt.savefig('mean_periodic_boundaries.png')
plt.close()
#condiciones iniciales
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_periodic[0,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras periodicas t=0')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_periodic_1.png')
plt.close()
#estados intermedios
a = temp_periodic.shape[0]
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_periodic[1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras periodicas t={}s'.format(mean_periodic[1000,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_periodic_2.png')
plt.close()

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_periodic[2,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras periodicas t={}s'.format(mean_periodic[2000,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_periodic_3.png')
plt.close()

#estado final
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_periodic[-1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras periodicas t={}s'.format(mean_periodic[-1,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_periodic_4.png')
plt.close()


#caso periodic boundaries
temp_open = np.loadtxt('temp_open_boundaries.txt')
mean_open = np.loadtxt('mean_open_boundaries.txt')
#Se convierte el arreglo 2D en uno 3D
a, b = temp_open.shape
print temp_open.shape
temp_open = temp_open.reshape((a/(b-1), b-1, b))
print temp_open.shape
x_open = temp_open[0,:,0]
y_open = temp_open[0,:,0]
x, y = np.meshgrid(x_open, x_open)
#Graficas
#Temperatura promedio
plt.plot(mean_open[:,0], mean_open[:,1])
plt.title('Temperatura promedio fronteras abiertas')
plt.xlabel('$t [s]$')
plt.ylabel('$T\ promedio\ [^\circ C]$')
plt.savefig('mean_open_boundaries.png')
plt.close()
#condiciones iniciales
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_open[0,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras abiertas t=0')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_open_1.png')
plt.close()
#estados intermedios
a = temp_open.shape[0]
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_open[1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras abiertas t={}s'.format(mean_open[1000,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_open_2.png')
plt.close()

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_open[2,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras abiertas t={}s'.format(mean_open[2000,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_open_3.png')
plt.close()

#estado final
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, temp_open[-1,:,1:], rstride=1, cstride=1, cmap=cm.autumn_r, linewidth=0, vmin=10, vmax=100)
ax.set_zlim(10,100)
ax.set_xlabel('$x [m]$')
ax.set_ylabel('$y [m]$')
ax.set_zlabel('$T [^\circ C]$')
ax.set_title('Temperatura fronteras abiertas t={}s'.format(mean_open[-1,0]))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('temp_open_4.png')
plt.close()
