import numpy as np
import matplotlib.pyplot as plt

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
plt.show()
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
plt.show()
plt.close()
