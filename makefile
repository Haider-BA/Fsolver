objects = settings.o partitioner.o io.o communicator.o schemes.o boundaryConditions.o solver.o cylinderIBMForce.o fcm.o channel.o main.o


FC = $(HOME)/apps/bin/mpif90
FFLAGS = -ffree-line-length-none -O3 -Wno-all
# for profiling: -pg 
#
# intel: (for profiling: -pg)
# flags = -fc=ifort
#FFLAGS = -no-wrap-margin -O1 -debug -heap-arrays
###FFLAGS = -fast -heap-arrays -wrap-margin-


channel: $(objects)
	$(FC) $(FFLAGS) -o channel $(objects)

#platelet_count.exe: platelet_count.f90
#	$(FC) $(FFLAGS) -o platelet_count.exe platelet_count.f90

settings.o : settings.f90
	$(FC) $(FFLAGS) -c settings.f90

io.o : io.f90
	$(FC) $(FFLAGS) -c io.f90

partitioner.o : partitioner.f90
	$(FC) $(FFLAGS) -c partitioner.f90

communicator.o : communicator.f90
	$(FC) $(FFLAGS) -c communicator.f90

schemes.o : schemes.f90
	$(FC) $(FFLAGS) -c schemes.f90

boundaryConditions.o : boundaryConditions.f90
	$(FC) $(FFLAGS) -c boundaryConditions.f90 

solver.o : solver.f90
	$(FC) $(FFLAGS) -c solver.f90

cylinderIBMForce.o : cylinderIBMForce.f90
	$(FC) $(FFLAGS) -c cylinderIBMForce.f90

fcm.o : fcm.f90
	$(FC) $(FFLAGS) -c fcm.f90

channel.o : channel.f90
	$(FC) $(FFLAGS) -c channel.f90

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90

# platelet_count.exe: platelet_count.f90
# 	ifort -o platelet_count.exe platelet_count.f90

clean:
	rm *.o *.obj *.mod *.exe *.objut *.dat *.plt *~ *pbs.e* *pbs_run.obj* *pbs_run.e* core.* *.pdb *.ilk channel
