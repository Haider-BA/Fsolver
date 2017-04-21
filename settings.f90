module settings

use mpi

    implicit none
    
    double precision, parameter :: pi = 4*atan(1.0)
    
    integer :: nx,ny,nz
    double precision :: dx,dy,dz,dt
    double precision :: rho,mu,nu
    double precision :: Dt_fbg, Dt_fn, Dt_th
    integer :: nt,nit
    double precision :: omega_P, omega_U, omega_V, omega_W, tol
    double precision :: Lx,Ly,Lz,RADIUS
    double precision :: U0
    double precision :: a		! particle radius
    integer:: nLog, nWrite		! nLog: write log to terminal every nLog timesteps
    							! write results to file every nwrite timesteps
    
    integer :: Np				! number of particles
    double precision :: x1, x2  ! Initial placement of particles between x1 and x2
	double precision :: x0		! location of injury (x0,y0,z0)
    
    
    contains
        subroutine setup(filename)

	        character(len=*) :: filename
	        character(len=30) :: temp

            integer :: id, Nproc, ierr

            call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)

            if (id==0) then
	        
                open(unit=10,file=filename,status='old',action='read')
                
                read (10,*) temp, Lx
                read (10,*) temp, Ly
                read (10,*) temp, Lz
                read (10,*) temp, RADIUS
                read (10,*) temp, dx
                read (10,*) temp, dy
                read (10,*) temp, dz
                read (10,*) temp, dt
                read (10,*) temp, nt
                read (10,*) temp, nit
                read (10,*) temp, rho
                read (10,*) temp, nu
                read (10,*) temp, Dt_fbg
                read (10,*) temp, Dt_fn
                read (10,*) temp, Dt_th
                read (10,*) temp, omega_P
                read (10,*) temp, omega_U
                read (10,*) temp, omega_V
                read (10,*) temp, omega_W
                read (10,*) temp, tol
                read (10,*) temp, U0
                read (10,*) temp, Np
                read (10,*) temp, x1
                read (10,*) temp, x2
				read (10,*) temp, x0
                read (10,*) temp, a
                read (10,*) temp, nLog
                read (10,*) temp, nWrite
                
			    close(10)

            end if

            ! broadcast to all
            call mpi_bcast(Lx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Ly,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Lz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(RADIUS,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(rho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Dt_fbg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Dt_fn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Dt_th,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_P,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_U,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_V,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(omega_W,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(U0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(Np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(x1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(x2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			call mpi_bcast(x0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nLog,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call mpi_bcast(nWrite,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

            !call mpi_barrier(MPI_COMM_WORLD,ierr)

			nx = floor(Lx/dx)
			ny = floor(Ly/dy)
			nz = floor(Lz/dz)

            !print*, Lx, Ly, Lz
			
			mu = nu*rho
			
            
        end subroutine setup
        
        
        
        subroutine setTimestep(a)
        	double precision :: a
        	dt = a
        end subroutine
        
        
        
        subroutine setParticleVariables(Yp,Vp,Vp0)
            ! use ifport
        	
        	double precision, parameter :: pi = 4*atan(1.0)
        	double precision ,dimension(Np,3) :: Yp, Vp, Vp0
        	
        	double precision :: r_max, r_rand, theta_rand, temp_rand
        	double precision :: dist
        	logical :: place_found
        	integer :: i,j

            integer :: error, id

            call mpi_comm_rank(MPI_COMM_WORLD,id,error)
			!call random_seed()								! <-- Will generate different random numbers for every run
        	
        	!Yp(1,1)=Lx/4; Yp(2,1)=Lx/4; Yp(3,1)=Lx/2; Yp(4,1)=3*Lx/4; Yp(5,1)=3*Lx/4
        	!Yp(1,2)=Ly/3; Yp(2,2)=2*Ly/3; Yp(3,2)=Ly/2; Yp(4,2)=Ly/3; Yp(5,2)=2*Ly/3
        	!Yp(1,3)=Lz/2; Yp(2,3)=Lz/2; Yp(3,3)=Lz/2; Yp(4,3)=Lz/2; Yp(5,3)=Lz/2;
        	
		! particles will be placed in a cylindrical volume of length
		! x2-x1 and radius r_max

            if (id==0) then

                r_max = 0.9*RADIUS
                
                ! -> Read from setup.txt file
                !x1 = 5e-6
                !x2 = 50e-6
                
                print*, ''
                print '(A,I0,A)', 'Finding locations for ', Np, ' particles...'
                
                do i=1,Np
                
                    ! print '(A20,I3)', 'Placing particle ', i
                    
                    place_found = .false.
                    
                    do while (place_found .eqv. .false.)
                        
						!call random_number(temp_rand)
                        !r_rand = temp_rand*r_max			! pseudo-random: same random numbers every time
                        r_rand = rand(0)*r_max
						
						!call random_number(temp_rand)
                        !theta_rand = temp_rand*2*pi
                        theta_rand = rand(0)*2*pi
                        
						!call random_number(temp_rand)
                        !Yp(i,1) = x1 + (x2-x1)*temp_rand
                        Yp(i,1) = x1 + (x2-x1)*rand(0)
                        Yp(i,2) = Ly/2 + r_rand*cos(theta_rand)	
                        Yp(i,3) = Lz/2 + r_rand*sin(theta_rand)
                        
                        !print*, Yp(i,:)
                        
                        place_found = .true.
                        
                        if (i>1) then
                            do j=1,i-1
                                dist = sqrt(sum((Yp(i,:)-Yp(j,:))**2))
                                if (dist<a) then
                                    place_found = .false.
                                    exit
                                end if
                            end do
                        end if
                        
                    end do
                
                end do

                ! Place three activated particles between 15 and 20 um
                !Yp(1,1) = 15e-6
                !Yp(2,1) = 17.5e-6
                !Yp(3,1) = 20e-6
                
                print '(A,I0,A)', 'Placed ', Np, ' particles.'
                print*, ''

            end if
            
            ! Broadcast to all processes
        	call mpi_bcast(Yp,Np*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)

        	Vp = 0
        	Vp0 = 0
        	
        end subroutine setParticleVariables
        
        
        
        

end module settings
