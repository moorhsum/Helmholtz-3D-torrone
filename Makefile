CC = g++
CFLAGS = -O3 -I.
PFLAG = -fopenmp
DEPS = distance_application3D.h finite_difference3D.h interpolate3D.h tools.h BICGSTAB.h functions.h HelmholtzKernel.h printslice.h weno.h frame.h PDE.h constants.h
OBJ = distance_application3D.o finite_difference3D.o interpolate3D.o tools.o BICGSTAB.o functions.o HelmholtzKernel.o printslice.o weno.o large3Dexample2.o

large3Dexample2: $(OBJ)
	$(CC) $(CFLAGS) $(PFLAG) -o large3Dexample2 $(OBJ)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(PFLAG) -c $< -o $@

.PHONY: clean

clean:
	rm *.o *~
