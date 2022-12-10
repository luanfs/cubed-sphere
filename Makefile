#--------------------------------------
# Makefile
#--------------------------------------

#Use gfortran
F90 = gfortran
FFLAG := -O3 -fopenmp

#Include path - for modules created
IMOD=-Ibin 

#Objects
OBJ=bin/constants.obj \
bin/datastruct.obj \
bin/miscellaneous.obj \
bin/linear_algebra.obj \
bin/sphgeo.obj \
bin/ppm_reconstruction.obj \
bin/ppm_flux.obj \
bin/advection_vars.obj \
bin/discrete_operators.obj \
bin/allocation.obj \
bin/deallocation.obj \
bin/output.obj \
bin/input.obj \
bin/cubed_sphere.obj \
bin/diagnostics.obj \
bin/advection_ic.obj \
bin/advection_timestep.obj \
bin/simulpack.obj \

#Compile and build all
all: header config bin/main ending

#Make all and run executable
run: all
	./main

#Print heading
header:
	@echo --------------------------------------------
	@echo Compiling and building the software   
	@echo --------------------------------------------
	@echo 
	@echo 

#Configure Enviroment (directories)
config:
	chmod +x sh/*.sh
	. sh/dirs.sh

#Constants
bin/constants.obj: src/constants.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv constants.mod bin/.

#Data struct
bin/datastruct.obj: src/datastruct.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv datastruct.mod bin/.

#Miscellaneous
bin/miscellaneous.obj: src/miscellaneous.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv miscellaneous.mod bin/.

#Numerical linear algebra
bin/linear_algebra.obj: src/linear_algebra.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv linear_algebra.mod bin/.

#Spherical geometry
bin/sphgeo.obj: src/sphgeo.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv sphgeo.mod bin/.

#Data allocation
bin/allocation.obj: src/allocation.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv allocation.mod bin/.

#Data deallocation
bin/deallocation.obj: src/deallocation.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv deallocation.mod bin/.

#Output
bin/output.obj: src/output.f90
	$(F90) $(FFLAG) $(NFFLAG) -c  $^ -o $@ $(IMOD)
	mv output.mod bin/.

#Input
bin/input.obj: src/input.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv input.mod bin/.

#Cubed-sphere generation
bin/cubed_sphere.obj: src/cubed_sphere.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv cubed_sphere.mod bin/.

#Diagnostics
bin/diagnostics.obj: src/diagnostics.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv diagnostics.mod bin/.

#PPM flux
bin/ppm_flux.obj: src/ppm_flux.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv ppm_flux.mod bin/.

#PPM reconstruction
bin/ppm_reconstruction.obj: src/ppm_reconstruction.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv ppm_reconstruction.mod bin/.

#Discrete operators
bin/discrete_operators.obj: src/discrete_operators.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv discrete_operators.mod bin/.

#Advection vars
bin/advection_vars.obj: src/advection_vars.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv advection_vars.mod bin/.

#Advection initial condition
bin/advection_ic.obj: src/advection_ic.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv advection_ic.mod bin/.

#Advection timestep
bin/advection_timestep.obj: src/advection_timestep.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv advection_timestep.mod bin/.

# Simulation package
bin/simulpack.obj: src/simulpack.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv simulpack.mod bin/.

#Main executable
bin/main: src/main.f90 $(OBJ)
	$(F90) $(FFLAG)  src/main.f90 $(OBJ) -o $@ $(IMOD)

#Creates a link for executable and prints ending
ending: 
	chmod +x sh/link.sh
	sh/link.sh
	@echo End of compilation
	@echo
	@echo "Set parameter files (pars / *.par )" 
	@echo "   and then run 'main'"
	@echo "------------------------------------------------------------------"


# Create tarball and backup
archive: 
        #Backup all important files
	chmod +x sh/backup.sh
	./sh/backup.sh

#Clean targets
clean: 
	rm -rf bin/*.obj bin/*.o bin/*.mod	
	rm -rf bin/main*
	rm -rf *~

cleandata: clean
	rm -rf data/
	rm -rf graphs/
	rm -rf bin/
	rm main

cleangrids: clean
	rm -rf graphs/
	rm -rf grids/

cleanall: clean cleandata cleangrids
