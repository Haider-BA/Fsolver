module boundaryConditions
use settings
use partitioner
use communicator

implicit none

contains

	subroutine setIC(u,v,w,P,c_fbg, c_fn, c_th, x,y,z)
		double precision, dimension(nz,ny,0:m+1) :: u,v,w,P
		double precision, dimension(nz,ny,0:m+1) :: c_fbg, c_fn, c_th
		double precision, dimension(nz,ny,1:m) :: x,y,z
		
		double precision, dimension(nz,ny,1:m) :: r
		
		

		u = U0
		v = 0
		w = 0
		P = 0

		P(3,1,1) = 1

        !print*, 'setIC: all set to zero done'
		
		r = sqrt( (y-Ly/2)**2 + (z-Lz/2)**2 )
		where (r<RADIUS) u(:,:,1:m) = U0 * ( 1 - (r/RADIUS)**2 ) 

        !print*, 'done 2'
            ! only the non-ghost points are set here; call sync function
            ! to set the ghost cells as well
        call communicate(u)
        call communicate(v)
        call communicate(w)
        call communicate(P)
        
        
        ! (1) Setting a constant value of fbg at inlet and a constant value of 
        ! Thrombin at the sticky patch and simulate the evolution of concs with flow
        ! (no plats in this case)
        
        c_fbg = 0; c_fn = 0; c_th = 0;

        call communicate(c_fbg)
        call communicate(c_fn)
        call communicate(c_th)

        !print*, 'done 3'
	
	end subroutine setIC
	
	
	
	subroutine setPressureBC(P_new,P,u,v,w)
	! calculates pressure at boundary points

	! pressure eqn: 
	! DEL2(P) = (rho/dt) * DIV(U)
	
		double precision, dimension(nz,ny,nx) :: u,v,w,P,P_new
		
		double precision :: ap,ae,aw,an,as,at,ab
		double precision, dimension(nz-2,ny-2) :: DuDx,DvDy,DwDz,b

		ap = 2/dx**2+2/dy**2+2/dz**2
		ae = 1/dx**2
		aw = 1/dx**2
		an = 1/dy**2
		as = 1/dy**2
		at = 1/dz**2
		ab = 1/dz**2

		! at x = 0 (excluding edge nodes)
		DuDx = (u(2:nz-1,2:ny-1,2)-u(2:nz-1,2:ny-1,nx))/(2*dx)
		DvDy = (v(2:nz-1,3:ny,1)-v(2:nz-1,1:ny-2,1))/(2*dy)
		DwDz = (w(3:nz,2:ny-1,1)-w(1:nz-2,2:ny-1,1))/(2*dz)

		b = -rho*(DuDx+DvDy+DwDz)/dt

		P_new(2:nz-1,2:ny-1,1) = (an*P(2:nz-1,3:ny,1)+as*P(2:nz-1,1:ny-2,1) &
			+ae*P(2:nz-1,2:ny-1,2)+aw*P(2:nz-1,2:ny-1,nx) &
			+at*P(3:nz,2:ny-1,1)+ab*P(1:nz-2,2:ny-1,1)+b)/ap


		! at x = L (excluding edge nodes)
		DuDx = (u(2:nz-1,2:ny-1,1)-u(2:nz-1,2:ny-1,nx-1))/(2*dx)
		DvDy = (v(2:nz-1,3:ny,nx)-v(2:nz-1,1:ny-2,nx))/(2*dy)
		DwDz = (w(3:nz,2:ny-1,nx)-w(1:nz-2,2:ny-1,nx))/(2*dz)

		b = -rho*(DuDx+DvDy+DwDz)/dt

		P_new(2:nz-1,2:ny-1,nx) = (an*P(2:nz-1,3:ny,nx)+as*P(2:nz-1,1:ny-2,nx) &
			+ae*P(2:nz-1,2:ny-1,1)+aw*P(2:nz-1,2:ny-1,nx-1) &
			+at*P(3:nz,2:ny-1,nx)+ab*P(1:nz-2,2:ny-1,nx)+b)/ap


		! dP/dy = 0 at y=0,H (including edge nodes)
		P_new(:,1,:) = P_new(:,2,:)
		P_new(:,ny,:) = P_new(:,ny-1,:)

		! dP/dz = 0 at z=0,W (including edge nodes)
		P_new(1,:,:) = P_new(2,:,:)
		P_new(nz,:,:) = P_new(nz-1,:,:)
		
		! fix at least one point (scarsborough criterion)
		!P_new(nz/2,ny/2,nx-5) = 0.0

	end subroutine setPressureBC


    !-------------------------------------------------------------
    ! setPressureBC_MPI:
    ! set dP/dn=0 at boundaries on y and z directions
    !-------------------------------------------------------------
    subroutine setPressureBC_MPI(P)
        double precision, dimension(nz,ny,0:m+1) :: P
    
        ! dP/dy = 0 at y=0,H (including edge nodes)
		P(:,1,:) = P(:,2,:)
		P(:,ny,:) = P(:,ny-1,:)

		! dP/dz = 0 at z=0,W (including edge nodes)
		P(1,:,:) = P(2,:,:)
		P(nz,:,:) = P(nz-1,:,:)
		

    end subroutine setPressureBC_MPI



    !---------------------------------------------------------------
    ! setVelocityBC_MPI:
    ! set velocities at no-slip walls in y and z directions to 0
    !---------------------------------------------------------------
    subroutine setVelocityBC_MPI(u,v,w)
        double precision, dimension(nz,ny,0:m+1) :: u,v,w

        call setToZero(u)
        call setToZero(v)
        call setToZero(w)

    end subroutine setVelocityBC_MPI

    subroutine setToZero(u)
        double precision, dimension(nz,ny,0:m+1) :: u

        u(1,:,:) = 0
        u(nz,:,:) = 0
        u(:,1,:) = 0
        u(:,ny,:) = 0
    end subroutine setToZero
	
	
	!~subroutine setNormalFluxZero(c)
	!~	double precision, dimension(nz,ny,0:m+1) :: c
	!~	
	!~end subroutine setNormalFluxZero
    
	
    !--------------------------------------------------------------------------
	! setConcBC_MPI:
	! set boundary conditions for concs of fibrinogen, thrombin and fibrin
	!--------------------------------------------------------------------------
    subroutine setConcBC_MPI(c_fbg, c_fn, c_th,x,y,z)
        double precision, dimension(nz,ny,0:m+1) :: c_fbg, c_fn, c_th
        double precision, dimension(nz,ny,1:m) :: x,y,z
		double precision, dimension(nz,ny,1:m) :: r
  		integer :: id, Nproc, ierr
  		
  		double precision, parameter :: a_inj = 1e-6	 		! radius of injury patch
  		!double precision, parameter :: x0 = 13e-6			! location of injury patch center (x0, 0, Lz/2)
  				
		call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    
    	r = sqrt( (y-Ly/2)**2 + (z-Lz/2)**2 )
		
		! (1) Setting a constant value of fbg at inlet and a constant value of 
        ! Thrombin at the sticky patch and simulate the evolution of concs with flow
        ! (no plats in this case)
        
		! inlet for fbg
        if (id==0) then
        	where (r(:,:,1)<3e-6) c_fbg(:,:,1) = 1
        endif
		
        ! injury patch for thrombin
        where ( r>RADIUS*0.9 .and. sqrt((x-x0)**2+y**2+(z-Lz/2)**2)<a_inj ) c_th(:,:,1:m) = 1
        
		!print*, 'end of concBC'
        call communicate(c_fbg)
        call communicate(c_fn)
        call communicate(c_th)

    end subroutine setConcBC_MPI


end module boundaryConditions
