
FC = gfortran
FFLAGS = -fbounds-check -fopenmp
LFLAGS = $(FFLAGS)

.PHONY: plots, clean, clobber

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

jacobi2d.exe: jacobi2d_main_omp.o jacobi2d_sub_omp.o
	$(FC) $(LFLAGS) jacobi2d_main_omp.o jacobi2d_sub_omp.o -o jacobi2d.exe

solution.txt: jacobi2d.exe
	@echo 
	@echo Running code...
	./jacobi2d.exe

plots: solution.txt
	@echo 
	@echo Plotting results...
	python plot_solution.py

clean:
	rm -f *.o *.exe

clobber: clean
	rm -f solution.txt *.png

