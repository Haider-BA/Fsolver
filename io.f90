module io

use mpi
use settings, only: nx, ny, nz, RADIUS, Lx, dt, Np, a
use partitioner, only: mesh, m
implicit none

contains
	subroutine printHeader
	
		character(len=100), parameter :: hline = '-----------------------------------------------------------------'
		
		print*, hline
		print*, 'PLATELET SIMULATION PROGRAM'
		print*, hline

		100 format(A20,'(',I3, ' x ',I3, ' x ',I3, ')')
		print 100, 'Mesh size: ', nx, ny, nz
		101 format(A20,ES10.2)
		print 101, 'Timestep: ',dt
		102 format(A20, ES10.2)
		print 102, 'Channel radius: ', RADIUS
		103 format(A20, ES10.2)
		print 103, 'Channel length: ', Lx
		104 format(A20, I5)
		print 104, 'Np : ',Np
	
		print*, hline
		
		print*, ''
		print*, 'Starting Run'
		print*, ''
	
	end subroutine printHeader
		
	
	
	subroutine writeFlowDataToFile(u,v,w,P,timestamp)
	
	double precision, dimension(:,:,:) :: u,v,w,P
	integer :: timestamp
	
	character(len=100) :: u_file, v_file, w_file, p_file
	integer :: n
	
	n = size(u)
	
	100 format(A2,I0.10,A4)
	write(u_file,100) 'u_',timestamp,'.dat'
	write(v_file,100) 'v_',timestamp,'.dat'
	write(w_file,100) 'w_',timestamp,'.dat'
	write(p_file,100) 'p_',timestamp,'.dat'
	
	
	open(unit=10,file=u_file,status='replace',action='write')
	open(unit=20,file=v_file,status='replace',action='write')
	open(unit=30,file=w_file,status='replace',action='write')
	open(unit=40,file=p_file,status='replace',action='write')
	
	write(10,*) reshape(u,(/n,1/))
	write(20,*) reshape(v,(/n,1/))
	write(30,*) reshape(w,(/n,1/))
	write(40,*) reshape(p,(/n,1/))
	
	close(10)
	close(20)
	close(30)
	close(40)
	
	write(*,*) '--> Data written to file.'
	write(*,*)
	
	end subroutine writeFlowDataToFile
	
	

	subroutine writeToTecplot(x,y,z,u,v,w,P,Yp,timestamp)
		
		double precision, dimension(nz,ny,nx) :: x,y,z,u,v,w,P
		double precision, dimension(Np,3) :: Yp
		integer :: timestamp

		integer, parameter :: nvar = 7					! number of variables to be written
	
		integer :: imax, jmax, kmax
		double precision, dimension(nvar,nz,ny,nx) :: var
		integer i,j,k,iv
		character(30) :: filename, particle_file
	
		imax = nz; jmax=ny; kmax=nx
	
		! set up var matrix
		var(1,:,:,:) = x
		var(2,:,:,:) = y
		var(3,:,:,:) = z
		var(4,:,:,:) = u
		var(5,:,:,:) = v
		var(6,:,:,:) = w
		var(7,:,:,:) = P



		! WRITE FLOW DATA	
		100 format(A5,I0.10,A4)
		write(filename,100) 'data_',timestamp,'.plt'
	
		open(unit=10,file=filename,status='replace',action='write')
		write(10,*) 'VARIABLES = "X","Y","Z","U","V","W","P"'
		write(10,*) 'ZONE DATAPACKING=BLOCK, I=', imax, ', J=', jmax, ', K=', kmax
		write(10,*) 'STRANDID=1  SOLUTIONTIME=',timestamp*dt
		do iv=1,nvar
			write(10,*) var(iv,:,:,:)
		end do
	
		close(10)
	
		
	
		! ALSO WRITE PARTICLE POSITIONS TO SEPARATE FILE
		200 format(A9,I0.10,A4)
		write(particle_file,200) 'particle_',timestamp,'.plt'
		
		open(unit=20,file=particle_file,status='replace',action='write')
		write(20,*) 'VARIABLES = "X","Y","Z"'
		do i=1,3
			write(20,*) Yp(:,i)
		end do
		close(20)
		
		
		
		write(*,*)		
		write(*,'(A35,A20,A2,A25)') '--> Data written to files: ',filename, ', ', particle_file
		write(*,*)
		
	
	end subroutine writeToTecplot



	subroutine writeToTecplot_MPI(x,y,z,u,v,w,P,c_fbg,c_fn,c_th,Yp,Vp,timestamp)
		
		double precision, dimension(nz,ny,1:m) :: x,y,z
		double precision, dimension(nz,ny,0:m+1) :: u,v,w,P
		double precision, dimension(nz,ny,0:m+1) :: c_fbg,c_fn,c_th
		double precision, dimension(Np,3) :: Yp
		double precision, dimension(Np,3) :: Vp
		integer :: timestamp

		integer :: imax, jmax, kmax
		double precision, dimension(:,:,:), allocatable ::x_all, y_all, z_all, u_all, v_all, w_all, P_all, c_fbg_all, c_fn_all, c_th_all ! all gathered data (only at root)
        integer, dimension(:), allocatable :: recvcounts, displs ! only at root
		integer :: i,j,k
		character(30) :: filename, particle_file

        integer :: id, Nproc, ierr

        call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)
        call mpi_comm_size(MPI_COMM_WORLD,Nproc,ierr)
	
		imax = nz; jmax=ny; kmax=nx

        if (id==0) then
            allocate(x_all(nz,ny,nx))
            allocate(y_all(nz,ny,nx))
            allocate(z_all(nz,ny,nx))
            allocate(u_all(nz,ny,nx))
            allocate(v_all(nz,ny,nx))
            allocate(w_all(nz,ny,nx))
            allocate(P_all(nz,ny,nx))
			allocate(c_fbg_all(nz,ny,nx))
			allocate(c_fn_all(nz,ny,nx))
			allocate(c_th_all(nz,ny,nx))

            allocate(recvcounts(Nproc))
            allocate(displs(Nproc))
        end if

        call mpi_gather(m,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (id==0) then
            recvcounts = nz*ny*recvcounts
            displs(1) = 0
            do i=2,Nproc
                displs(i) = sum(recvcounts(1:i-1))
            end do

        !print*, recvcounts
        !print*, displs

        end if

        call mpi_gatherv(x,nz*ny*m,MPI_DOUBLE_PRECISION,x_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(y,nz*ny*m,MPI_DOUBLE_PRECISION,y_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(z,nz*ny*m,MPI_DOUBLE_PRECISION,z_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(u,nz*ny*m,MPI_DOUBLE_PRECISION,u_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(v,nz*ny*m,MPI_DOUBLE_PRECISION,v_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(w,nz*ny*m,MPI_DOUBLE_PRECISION,w_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call mpi_gatherv(P,nz*ny*m,MPI_DOUBLE_PRECISION,P_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call mpi_gatherv(c_fbg,nz*ny*m,MPI_DOUBLE_PRECISION,c_fbg_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call mpi_gatherv(c_fn,nz*ny*m,MPI_DOUBLE_PRECISION,c_fn_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call mpi_gatherv(c_th,nz*ny*m,MPI_DOUBLE_PRECISION,c_th_all,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


        if (id==0) then
		
            ! WRITE FLOW DATA	
            100 format(A5,I0.10,A4)
            write(filename,100) 'data_',timestamp,'.tec'
        
            open(unit=10,file=filename,status='replace',action='write')
            write(10,*) 'VARIABLES = "X","Y","Z","U","V","W","P","C_fbg","C_fn","C_th"'
            write(10,*) 'ZONE DATAPACKING=BLOCK, I=', imax, ', J=', jmax, ', K=', kmax
            write(10,*) 'STRANDID=1  SOLUTIONTIME=',timestamp*dt
            
            write(10,*) x_all
            write(10,*) y_all
            write(10,*) z_all
            write(10,*) u_all
            write(10,*) v_all
            write(10,*) w_all
            write(10,*) P_all
			write(10,*) C_fbg_all
			write(10,*) C_fn_all
			write(10,*) C_th_all

        
            close(10)
        
            
        
            ! ALSO WRITE PARTICLE POSITIONS TO SEPARATE FILE
            200 format(A9,I0.10,A4)
            write(particle_file,200) 'particle_',timestamp,'.csv'
            
            open(unit=20,file=particle_file,status='replace',action='write')
            write(20,*) 'X  Y  Z  Ux Uy Uz RADIUS'
            do i=1,Np
                !write(20,'(ES10.2,A1,ES10.2,A1,ES10.2)') Yp(i,1),',',Yp(i,2),',',Yp(i,3)
                write(20,'(7ES12.2)') Yp(i,1),Yp(i,2),Yp(i,3),Vp(i,1),Vp(i,2),Vp(i,3),a
            end do
            close(20)
            
            
            
            write(*,*)		
            write(*,'(A35,A20,A2,A25)') '--> Data written to files: ',filename, ', ', particle_file
            write(*,*)
            
        end if
	
	end subroutine writeToTecplot_MPI






    !--------------------------------------------------------------
    ! writeToBOV_MPI:
    !--------------------------------------------------------------
	subroutine writeToBOV_MPI(u,v,w,P,Yp,timestamp)
		
		double precision, dimension(nz,ny,0:m+1) :: u,v,w,P
		double precision, dimension(Np,3) :: Yp
		integer :: timestamp

		character(30) :: filename, particle_file

        integer :: id,Nproc,ierr

        call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

        ! write BOV header file



		! WRITE FLOW DATA	
		
        
        100 format(A,I0,A)
		write(filename,100) 'u_',id,'.dat'
	
		open(unit=10,file=filename,status='replace',action='write', form='unformatted')
        write(10) u
        close(10)
	
		write(filename,100) 'v_',id,'.dat'
	
		open(unit=10,file=filename,status='replace',action='write', form='unformatted')
        write(10) v
        close(10)

		write(filename,100) 'w_',id,'.dat'
	
		open(unit=10,file=filename,status='replace',action='write', form='unformatted')
        write(10) w
        close(10)

	
		! ALSO WRITE PARTICLE POSITIONS TO SEPARATE FILE
        !200 format(A9,I0.10,A4)
		!write(particle_file,200) 'particle_',timestamp,'.plt'
		!
		!open(unit=20,file=particle_file,status='replace',action='write')
		!write(20,*) 'VARIABLES = "X","Y","Z"'
		!do i=1,3
		!	write(20,*) Yp(:,i)
		!end do
		!close(20)
		!
		
		
		!write(*,*)		
		!write(*,'(A35,A20,A2,A25)') '--> Data written to files: ',filename, ', ', particle_file
		!write(*,*)
		!

        print*, 'data written from ', id
	
	end subroutine writeToBOV_MPI


    !-----------------------------------------------------------------------------------
    ! writeParticleData:
    ! write particle data only to file
    !-----------------------------------------------------------------------------------
    subroutine writeParticleData(Yp,timestamp)

        double precision, dimension(Np,3) :: Yp
        integer :: timestamp

        character(30) :: particle_file
        integer :: i

        integer :: Nproc, id, ierr

        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
	
        if (id==0) then

            ! WRITE PARTICLE POSITIONS ONLY
            200 format(A9,I0.10,A4)
            write(particle_file,200) 'particle_',timestamp,'.plt'
            
            open(unit=20,file=particle_file,status='replace',action='write')
            write(20,*) 'VARIABLES = "X","Y","Z"'
            do i=1,3
                write(20,*) Yp(:,i)
            end do
            close(20)
            
            
            
            write(*,*)		
            write(*,'(A,A)') '--> Data written to file: ', particle_file
            write(*,*)
        end if
            

    end subroutine writeParticleData

	

    !---------------------------------------------------------------------------------------------
    ! writeLogToTerminal:
    ! write requested data to terminal (serial version)
    !---------------------------------------------------------------------------------------------
	subroutine writeLogToTerminal(timestep,p_res,p_err,NITS,Umax,Vmax,Wmax,Pmax,Yp,Vp,V1_name,V1,V2_name,V2,V3_name,V3,V4_name,V4,V5_name,V5)
	
	integer :: timestep
	double precision :: p_res,p_err
	integer :: NITS									! number of iterations for pressure equation to converge
	double precision :: Umax,Vmax,Wmax,Pmax				! flow variables
	double precision, dimension(Np,3) :: Yp,Vp			! particle variables

	! (optional additional variables to write to terminal)
	character(len=*),optional :: V1_name, V2_name, V3_name, V4_name, V5_name
	double precision, optional :: V1,V2,V3,V4,V5
	
	integer :: i
	
	character(len=100), parameter :: hline = '-----------------------------------------------------------------'
	
	100 format (A1,A25,A1,I10,A30)
	101 format (A1,A25,A1,ES10.2,A30)
	102 format (A1,I10,A2,3ES10.2,A2,3ES10.2,A2)
	
	print *, hline
	print '(A1,A15,I8,A25,ES10.2,A2,A6)', '|','Timestep = ',timestep,' Simulation time = ', timestep*dt, ' s', '|'
	print *, hline
	print 100, '|','Pressure Iterations', '|',NITS, '|'
	print 101, '|','Pressure Error', '|',p_err, '|'
	print *, hline

	! optional variables to be written to terminal
	if (present(V1)) then
		print 101, '|',V1_name, '|',V1, '|'
	end if
	if (present(V2)) then
		print 101, '|',V2_name, '|',V2, '|'
	end if
	if (present(V3)) then
		print 101, '|',V3_name, '|',V3, '|'
	end if
	if (present(V4)) then
		print 101, '|',V4_name, '|',V4, '|'
	end if
	if (present(V5)) then
		print 101, '|',V5_name, '|',V5, '|'
	end if

	print *, hline
	print 101, '|','U_max', '|', Umax, '|'
	print 101, '|','V_max', '|', Vmax, '|'
	print 101, '|','W_max', '|', Wmax, '|'
	print 101, '|','P_max', '|', Pmax, '|'
	print *, hline
	print '(A1,A10,A2,A30,A2,A30,A2)', '|','Particle', '|', 'Position', '|', 'Velocity','|'
	do i=1,Np
		print 102, '|', i, '|', Yp(i,:), '|', Vp(i,:), '|'
	end do
	print *, hline
	!print*, ''
	
	
	end subroutine writeLogToTerminal



    !------------------------------------------------------------
    ! writeLogToTerminal_MPI:
    ! write requested data to terminal (parallel)
    !-----------------------------------------------------------
    subroutine writeLogToTerminal_MPI(timestep, &
    V1_name,V1_type,V1, &
    V2_name,V2_type,V2, &
    V3_name,V3_type,V3, &
    V4_name,V4_type,V4, &
	V5_name,V5_type,V5, &
	V6_name,V6_type,V6, &
	V7_name,V7_type,V7, &
	V8_name,V8_type,V8, &
	V9_name,V9_type,V9, &
    V10_name,V10_type,V10)
	
	integer :: timestep

	! (optional variables to write to terminal)
	character(len=*),optional :: V1_name, V2_name, V3_name, V4_name, V5_name, V6_name, V7_name, V8_name, V9_name, V10_name
	character(len=*),optional :: V1_type, V2_type, V3_type, V4_type, V5_type, V6_type, V7_type, V8_type, V9_type, V10_type ! must be 'integer' or 'double'
	double precision, optional :: V1,V2,V3,V4,V5,V6,V7,V8,V9,V10
	
	integer :: i
	
    character(len=100), parameter :: error_msg = 'Error in writeLogToTerminal_MPI(): value type must be "integer" or "double"'
	character(len=100), parameter :: hline = '-----------------------------------------------------------------'
	
	101 format (A1,A25,A1,I10,A30)
	102 format (A1,A25,A1,ES10.2,A30)
	
	print *, hline
	print '(A1,A15,I8,A25,ES10.2,A2,A6)', '|','Timestep = ',timestep,' Simulation time = ', timestep*dt, ' s', '|'
	print *, hline

	! optional variables to be written to terminal
	if (present(V1)) then
        if (V1_type .eq. 'integer') then
		    print 101, '|',V1_name, '|',int(V1), '|'
        else if (V1_type .eq. 'double') then
		    print 102, '|',V1_name, '|',V1, '|'
        else
            print *, error_msg
        end if
	end if
    
    if (present(V2)) then
        if (V2_type .eq. 'integer') then
		    print 101, '|',V2_name, '|',int(V2), '|'
        else if (V2_type .eq. 'double') then
		    print 102, '|',V2_name, '|',V2, '|'
        else
            print *, error_msg
        end if
	end if
	
    if (present(V3)) then
        if (V3_type .eq. 'integer') then
		    print 101, '|',V3_name, '|',int(V3), '|'
        else if (V3_type .eq. 'double') then
		    print 102, '|',V3_name, '|',V3, '|'
        else
            print *, error_msg
        end if
	end if
    
    if (present(V4)) then
        if (V4_type .eq. 'integer') then
		    print 101, '|',V4_name, '|',int(V4), '|'
        else if (V4_type .eq. 'double') then
		    print 102, '|',V4_name, '|',V4, '|'
        else
            print *, error_msg
        end if
	end if
    
    if (present(V5)) then
        if (V5_type .eq. 'integer') then
		    print 101, '|',V5_name, '|',int(V5), '|'
        else if (V5_type .eq. 'double') then
		    print 102, '|',V5_name, '|',V5, '|'
        else
            print *, error_msg
        end if
	end if
	
	if (present(V6)) then
        if (V6_type .eq. 'integer') then
		    print 101, '|',V6_name, '|',int(V6), '|'
        else if (V6_type .eq. 'double') then
		    print 102, '|',V6_name, '|',V6, '|'
        else
            print *, error_msg
        end if
	end if
	
	if (present(V7)) then
        if (V7_type .eq. 'integer') then
		    print 101, '|',V7_name, '|',int(V7), '|'
        else if (V7_type .eq. 'double') then
		    print 102, '|',V7_name, '|',V7, '|'
        else
            print *, error_msg
        end if
	end if
	
	if (present(V8)) then
        if (V8_type .eq. 'integer') then
		    print 101, '|',V8_name, '|',int(V8), '|'
        else if (V8_type .eq. 'double') then
		    print 102, '|',V8_name, '|',V8, '|'
        else
            print *, error_msg
        end if
	end if
	
	if (present(V9)) then
        if (V9_type .eq. 'integer') then
		    print 101, '|',V9_name, '|',int(V9), '|'
        else if (V9_type .eq. 'double') then
		    print 102, '|',V9_name, '|',V9, '|'
        else
            print *, error_msg
        end if
	end if
	
	if (present(V10)) then
        if (V10_type .eq. 'integer') then
		    print 101, '|',V10_name, '|',int(V10), '|'
        else if (V10_type .eq. 'double') then
		    print 102, '|',V10_name, '|',V10, '|'
        else
            print *, error_msg
        end if
	end if
	

    ! draw a line if at least one additional info was printed
    if (present(V1)) print *, hline

    print*, ''
		
	end subroutine writeLogToTerminal_MPI


end module io
