module fcm

use mpi
use settings
use partitioner, only: mesh,m
implicit none

! Delta function = (2*pi*sigma**2)^(-3/2) * exp( (x-Yn)**2/(2*sigma**2) )
! to be calculated once at each timestep for all particles
double precision, dimension(:,:,:,:), allocatable :: delta_fn

contains


subroutine calc_delta_fn(Yp)
    implicit none
    double precision, dimension(Np,3) :: Yp
    double precision :: temp2
    double precision, dimension(nz,ny,m) :: dist_sqr, temp1
    double precision :: sigma
    integer :: i

    if (.not. allocated(delta_fn)) allocate(delta_fn(nz,ny,1:m,Np))

    sigma = a/sqrt(pi)

    temp2 = (2*pi*sigma**2)**(-1.5)
    
    do i=1,Np
        dist_sqr = (mesh%x-Yp(i,1))**2 + (mesh%y-Yp(i,2))**2 + (mesh%z-Yp(i,3))**2
	    temp1 = -dist_sqr/(2*sigma**2)
        delta_fn(:,:,:,i) = temp2*exp(temp1)
    end do
end subroutine calc_delta_fn




subroutine updateParticles(Yp,Vp,Vp0,u,v,w)

	implicit none
	
	double precision, dimension(nz,ny,0:m+1) :: u,v,w
	double precision, dimension(Np,3) :: Yp, Vp, Vp0

    double precision :: dv
    double precision, dimension(Np,3) :: V_loc
	
	integer :: i

    integer :: id, Nproc, ierr
	
    !call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    !call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    if (.not. allocated(delta_fn)) then
        !print*, 'Error (forceCoupling): call subroutine calc_delta_fn() before calling forceCoupling()'
        call calc_delta_fn(Yp)
    end if


	Vp0 = Vp

    dv = dx*dy*dz

    do i=1,Np    
	!do i=id*(Np/Nproc)+1,min((id+1)*(Np/Nproc),Np)
	
		!dist_sqr = (x-Yp(i,1))**2 + (y-Yp(i,2))**2 + (z-Yp(i,3))**2
	
		!temp1 = -dist_sqr/(2*sigma**2)

		V_loc(i,1) = sum(u(:,:,1:m)*delta_fn(:,:,:,i))*dv
		V_loc(i,2) = sum(v(:,:,1:m)*delta_fn(:,:,:,i))*dv
		V_loc(i,3) = sum(w(:,:,1:m)*delta_fn(:,:,:,i))*dv
	
	end do

    !call mpi_bcast(Vp(id*(Np/Nproc)+1:min((id+1)*(Np/Nproc),Np),:),3*(min((id+1)*(Np/Nproc)-id*(Np/Nproc)+1,Np)),MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)

    call mpi_allreduce(V_loc(:,1),Vp(:,1),Np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(V_loc(:,2),Vp(:,2),Np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(V_loc(:,3),Vp(:,3),Np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    !print*, V_loc(1,1), V_loc(1,2), V_loc(1,3)

    Yp = Yp + dt * Vp

	!Yp(4:Np,:) = Yp(4:Np,:) + dt*Vp(4:Np,:)


    ! Don't move three fixed particles


end subroutine updateParticles





subroutine forceCoupling(Gx,Gy,Gz,Yp)
	
	double precision, dimension(nz,ny,1:m) :: Gx,Gy,Gz
	double precision, dimension(Np,3) :: Yp
	
	double precision, dimension(Np,3) :: F, F_stokes, F_part_coll, F_wall_coll, F_vwf, F_pivkin
	double precision, dimension(Np,3) :: F_fbg, F_fn, F_collagen			! Fbg and fn are interparticle forces and collagen in particle-sticky patch attraction
	
	integer :: i
	
	!integer :: error, Nproc, id

    !print*, 'entered forceCoupling sub'

    call calc_delta_fn(Yp)

	!call mpi_comm_size(MPI_COMM_WORLD, Nproc, error)
    !call mpi_comm_rank(MPI_COMM_WORLD, id, error)
    
    !if (id==0) then
    !    print*, 'nprocs=', Nproc
    !end if
    !print*, id

	! sigma = a/sqrt(pi)


	! assume particles are neutrally buoyant. therefore no inertia term

	!F = 
	!F = 1e-15;
	
	!temp2 = (2*pi*sigma**2)**(-1.5)

	! Stokes force to ensure terminal velocity of U0
	!F_stokes(:,1) = 6*pi*mu*a*U0
	F_stokes(:,1) = 0
	F_stokes(:,2:3) = 0
	
	
	
	! Initialize all forces
	F_part_coll = 0
	F_wall_coll = 0
	F_vwf = 0
	
	F_fbg = 0
	F_fn = 0
	F_collagen = 0
	
	F_pivkin = 0
	
	
	! CALCULATE FORCES
	
	call particleCollisionForce(F_part_coll,Yp)
	call wallCollisionForce(F_wall_coll,Yp)
	!call vwfForce(F_vwf,Yp)
	
	!call pivkinForceField(F_pivkin,Yp) 
    !F_pivkin = -F_pivkin
	
	call collagenForce(F_collagen,Yp)
	
	
	! Add all forces
	F = F_stokes + F_part_coll + F_wall_coll + F_vwf + F_pivkin + F_fbg + F_fn + F_collagen
	!print*, F, -6*pi*mu*a*U0 
	!F(:,1) = -6*pi*mu*a*U0; F(:,2:3) = 0
	!F = 0
		
	Gx = 0; Gy  =0; Gz = 0;

    !print*, id*(Np/Nproc)+1,min((id+1)*(Np/Nproc),Np)

    ! local sum
	!do i=id*(Np/Nproc)+1,min((id+1)*(Np/Nproc),Np)
	do i=1,Np
	
		!dist_sqr = (x-Yp(i,1))**2 + (y-Yp(i,2))**2 + (z-Yp(i,3))**2
	
		!temp1 = -dist_sqr/(2*sigma**2)

		Gx = Gx + F(i,1)*delta_fn(:,:,:,i)/rho
		Gy = Gy + F(i,2)*delta_fn(:,:,:,i)/rho
		Gz = Gz + F(i,3)*delta_fn(:,:,:,i)/rho
	
	end do

    !print *, 'before allreduce'

    !call mpi_allreduce(Gx_loc,Gx,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    !call mpi_allreduce(Gy_loc,Gy,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    !call mpi_allreduce(Gz_loc,Gz,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)

    !print*, 'after allreduce'


end subroutine forceCoupling


subroutine pivkinForceField(F,Yp)

	double precision, dimension(Np,3) :: Yp
	double precision, dimension(Np,3) :: F
	
	double precision, dimension(3) :: r, r_wall

	double precision :: r_mag, r_wall_mag

	integer :: i,j

	! parameters
	double precision, parameter :: R_piv = 1.0
	double precision, parameter :: L_piv = 1.5
	double precision, parameter :: B_piv = 3.5
	double precision, parameter :: M_piv = (L_piv+B_piv)/2
	double precision, parameter :: alpha_piv = 4e-12

    ! activation zone
    double precision, parameter :: x0 = 15e-6
    double precision, parameter :: x1 = 20e-6


	do i=1,Np
		F(i,:) = 0
	
		! Wall
		r_wall_mag = a - sqrt(Yp(i,2)**2+Yp(i,3)**2)
		r_wall(1)=0; r_wall(2)=Yp(i,2);	r_wall(3)=Yp(i,3);
		r_wall = -r_wall_mag*r_wall/sqrt(sum(r_wall**2))	! directed away from wall
		if (r_wall_mag/a <= R_piv) then
			F(i,:) = F(i,:) + 2*alpha_piv*(1-(r_wall_mag/a)/R_piv)*r_wall/r_wall_mag
		end if

		
		! other particles
		do j=1,Np
			if (i/=j) then
				r = Yp(i,:)-Yp(j,:)	! directed j -> i
				r_mag = sqrt(sum(r**2))
				if (r_mag/a <= R_piv) then
					F(i,:) = F(i,:) + 2*alpha_piv*(1-(r_mag/a)/R_piv)*r/r_mag ! repulsive (accounting for collision)
				else if ( (r_mag/a >= R_piv) .and. (r_mag/a <= L_piv) ) then
					F(i,:) = F(i,:) + 0
				else if ( (r_mag/a >= L_piv) .and. (r_mag/a <= M_piv) ) then
                    if (Yp(i,1)>x0 .and. Yp(i,1)<x1) then
    					F(i,:) = F(i,:) - alpha_piv*((r_mag/a)/L_piv-1)*r/r_mag ! attractive (incresing with distance)
                    end if
				else if ( (r_mag/a >= M_piv) .and. (r_mag/a <= B_piv) ) then
                    if (Yp(i,1)>x0 .and. Yp(i,1)<x1) then
					    F(i,:) = F(i,:) - alpha_piv*(M_piv/L_piv-1)*((B_piv-(r_mag/a))/(B_piv-M_piv))*r/r_mag ! attractive (decreasing with distance)
                    end if
				end if
			end if
		end do
		
		
	end do

end subroutine pivkinForceField




subroutine particleCollisionForce(F,Yp)

	double precision, dimension(Np,3) :: Yp
	
	double precision :: R_ref
	double precision :: F_ref
	double precision, dimension(3) :: rij, r			! rij: vector from j to i, r: vector from center of i to surface of j
	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision :: rij_mag, r_mag
	
	integer :: i, j
	
	! Parameters
	R_ref = 1.0
	F_ref = 1e-9
	
	do i=1,Np
		F(i,:) = 0
		
		! other particles
		do j=1,Np
			if (i/=j) then
				rij = Yp(i,:)-Yp(j,:)
				rij_mag = sqrt(sum(rij**2))
				r = (1-a/rij_mag)*rij
				r_mag = sqrt(sum(r**2))
				if (r_mag/a < R_ref) then
					!F(i,:) = F(i,:) + F_ref*( (R_ref**2 - r_mag**2)/(R_ref**2 - 4*a**2) )**2 * r / r_mag
					F(i,:) = F(i,:) + F_ref*( 1 - (r_mag/a)/R_ref )* r / r_mag
				end if
			end if
		end do
		
		
	end do
				
end subroutine particleCollisionForce



subroutine wallCollisionForce(F,Yp)

	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision, dimension(Np,3) :: Yp
	
	double precision, dimension(3) :: ri, r					! ri: axis -> Yp, r: Yp -> wall
	double precision :: ri_mag, r_mag
	
	double precision :: R_ref		! close to channel RADIUS
	double precision :: F_ref		! force
	
	integer :: i
	
	! Parameters
	R_ref = 1.0
	F_ref = 1e-7
	
	F = 0
	
	do i=1,Np
		ri(1) = 0; ri(2) = Yp(i,2)-Ly/2; ri(3) = Yp(i,3)-Lz/2	
		ri_mag = sqrt(sum(ri**2))
		r = (1-RADIUS/ri_mag)*ri
		r_mag = sqrt(sum(r**2))
		if (r_mag/a < R_ref) then
			!F(i,:) = F_ref*( (r_mag**2 - R_ref**2)/(RADIUS**2 - R_ref**2) )**2 * r / r_mag
			F(i,:) = F_ref*( 1 - (r_mag/a)/R_ref )* r / r_mag
		end if
	end do

end subroutine wallCollisionForce


! Is vWF inter-particle or particle-wound ?
subroutine vwfForce(F,Yp)
	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision, dimension(Np,3) :: Yp

	double precision :: x0, y0, z0 	! location of injury
	
	double precision, dimension(3) :: r			! r vector from j to i
	double precision :: r_mag
	double precision :: vwfRange
	
	double precision :: F_ref
	
	integer :: i
	
	x0 = 30e-6; y0=0.0; z0 = Lz/2;
	vwfRange = 5e-6
	F_ref = 1e-10
	
	F = 0
	
	do i=1,Np
	
		r(1) = Yp(i,1) - x0; r(2) = Yp(i,2) - y0; r(3) = Yp(i,3) - z0;
		r_mag = sqrt(sum(r**2))
		
		if (r_mag<vwfRange) then
			F(i,:) = -F_ref* r
		end if
	end do
	
end subroutine vwfForce


subroutine fbgForce(F,Yp)

	! depends on local fbg concentration
	! To be added when biochem is implemented
	
	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision, dimension(Np,3) :: Yp
	
	double precision, dimension(3) :: r			! r vector from j to i
	double precision :: r_mag
	
	integer :: i, j
	
	do i=1,Np
		F(i,:) = 0
		
	end do
	
end subroutine fbgForce



subroutine fnForce(F,Yp)
	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision, dimension(Np,3) :: Yp
end subroutine fnForce



subroutine collagenForce(F,Yp)
	double precision, dimension(Np,3) :: F		! force vector for each particle
	double precision, dimension(Np,3) :: Yp

	double precision :: y0, z0 	! location of injury (x0,y0,z0), x0 read from setup.txt
	
	double precision, dimension(3) :: ri, r					! ri: sticky-point -> Yp, r: sticky-point -> Yp surface
	double precision :: ri_mag, r_mag
	double precision :: collagenRange
	
	double precision :: F_ref
	
	integer :: i
	
	integer :: nPart
	
	! PARAMETERS FOR THIS COLLAGEN MODEL
	y0=Ly/2-RADIUS; z0 = Lz/2;
	collagenRange = 5e-6
	F_ref = 1.0e-9
	
	F = 0
	
	nPart=0
	do i=1,Np
	
		ri(1) = Yp(i,1) - x0; ri(2) = Yp(i,2) - y0; ri(3) = Yp(i,3) - z0;
		ri_mag = sqrt(sum(ri**2))
		r = (1-a/ri_mag)*ri
		r_mag = sqrt(sum(r**2))
		
		if (r_mag<collagenRange) then
			F(i,:) = -F_ref* r/collagenRange
			nPart = nPart+1
		end if
	end do
	
	print*, 'Particles in range = ', nPart
end subroutine collagenForce



end module fcm
