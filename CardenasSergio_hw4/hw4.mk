#
trayectoria_45.png trayectoria_varios_angulos.png mean_fixed_boundaries.png mean_periodic_boundaries.png mean_open_boundaries.png temp_fixed_1.png temp_fixed_2.png temp_fixed_3.png temp_fixed_4.png temp_periodic_1.png temp_periodic_2.png temp_periodic_3.png temp_periodic_4.png temp_open_1.png temp_open_2.png temp_open_3.png temp_open_4.png : proyectil.txt proyectil2.txt mean_fixed_boundaries.txt mean_open_boundaries.txt mean_periodic_boundaries.txt temp_fixed_boundaries.txt temp_open_boundaries.txt temp_periodic_boundaries.txt
	python Plots_hw4.py

mean_open_boundaries.txt mean_periodic_boundaries.txt temp_fixed_boundaries.txt temp_open_boundaries.txt temp_periodic_boundaries.txt: PDE.cpp
	g++ PDE.cpp
	/a.out
	rm ./a.out

proyectil.txt proyectil2.txt : ODE.cpp
	g++ ODE.cpp
	/a.out
	rm ./a.out

clean : 
	rm *.txt *.png