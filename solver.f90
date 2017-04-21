module solver

use mpi
use settings, only: nx, ny, nz, dx, dy, dz, dt, rho, nit, omega_P, omega_U, omega_V, omega_W, tol
use schemes
use boundaryConditions
use communicator

implicit none

contains
	subroutine PressureSolver(P,R_rms,E_rms,NITS_final,u,v,w)

		! Returns pressure and rms residual

		double precision, dimension(nz,ny,nx) :: u,v,w
		double precision, dimension(nz,ny,nx) :: P
		double precision :: R_rms, E_rms
		integer :: NITS_final
		
		
		double precision, dimension(nz-2,ny-2,nx-2) :: DIVU, b, R
		double precision, dimension(nz,ny,nx) :: E
		double precision, dimension(nz,ny,nx) :: P_new
		double precision :: ap,ae,aw,an,as,at,ab,P_max
		integer :: i,n
		
		
		call DIV(DIVU,u,v,w)

		ap = 2/dx**2 + 2/dy**2 + 2/dz**2;
		ae = 1/dx**2;
		aw = 1/dx**2;
		an = 1/dy**2;
		as = 1/dy**2;
		at = 1/dz**2;
		ab = 1/dz**2;
		b = -rho*(DIVU)/dt; ! dt affects possibility of convergence

		n = nx*ny*nz

		!P_new = P

		do i=1,nit
			!print*, i
			P_new(2:nz-1,2:ny-1,2:nx-1) = (1-omega_P)*P(2:nz-1,2:ny-1,2:nx-1) &
				+omega_P*(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
				+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
				+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b)/ap;

			call setPressureBC(P_new,P,u,v,w)

			if (mod(i,50)==0) then

				E = P_new - P
				E_rms = sqrt(sum(E**2)/n)

				if (E_rms<tol) then
					exit
				end if

			end if

			P = P_new

		end do
		
		NITS_final = i
		
		R = ap*P(2:nz-1,2:ny-1,2:nx-1)-(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
			+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
			+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b);
		R_rms = sqrt(sum(R**2)/n)

	end subroutine PressureSolver


	
    subroutine PressureSolver_MPI(P,E_rms,NITS_final,u,v,w)
        double precision, dimension(nz,ny,0:m+1) :: u,v,w,P
        double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
        
        double precision, dimension(nz-2,ny-2,1:m) :: DIVU, b
        double precision, dimension(nz,ny,1:m) :: E
        double precision, dimension(nz,ny,0:m+1) :: P_new
        double precision :: ap,ae,aw,an,as,at,ab,P_max
        integer :: i,n

        integer :: ierr

        call DIV(DIVU,u,v,w)

		ap = 2/dx**2 + 2/dy**2 + 2/dz**2;
		ae = 1/dx**2;
		aw = 1/dx**2;
		an = 1/dy**2;
		as = 1/dy**2;
		at = 1/dz**2;
		ab = 1/dz**2;
		b = -rho*(DIVU)/dt; ! dt affects possibility of convergence

		n = nx*ny*nz

		!P_new = P

		do i=1,nit
			!print*, i
			P_new(2:nz-1,2:ny-1,1:m) = (an*P(2:nz-1,3:ny,1:m)+as*P(2:nz-1,1:ny-2,1:m) &
				+ae*P(2:nz-1,2:ny-1,2:m+1)+aw*P(2:nz-1,2:ny-1,0:m-1) &
				+at*P(3:nz,2:ny-1,1:m)+ab*P(1:nz-2,2:ny-1,1:m)+b)/ap;

			call setPressureBC_MPI(P_new)

            ! COMMUNICATE(P)
            call communicate(P_new)

			if (mod(i,50)==0) then

				E = P_new(:,:,1:m) - P(:,:,1:m)
				E_sqr_loc = sum(E**2)

                call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                E_rms = sqrt(E_sqr/n)

				if (E_rms<tol) then
					exit
				end if

			end if

            ! MPI_Barrier ??

			P = P + omega_P * (P_new - P)

		end do
		
		NITS_final = i
		
		! R = ap*P(2:nz-1,2:ny-1,2:nx-1)-(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
		!	+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
		!	+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b);
		! R_rms = sqrt(sum(R**2)/n)



    end subroutine PressureSolver_MPI

	
	
		!----------------------------------------------
    	! PressureSolver_MPI_BiCGStab
    	!----------------------------------------------    
    subroutine PressureSolver_MPI_BiCGStab(X,E_rms,NITS_final,u,v,w)
        use settings, only: nx, ny, nz, dx, dy, dz, dt, rho, nit, tol
        double precision, dimension(nz,ny,0:m+1) :: u,v,w,X
        double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
        
        double precision, dimension(nz,ny,0:m+1) :: DIVU, b ! includes boundary nodes and flaps
        double precision, dimension(nz,ny,1:m) :: E
        double precision, dimension(nz,ny,0:m+1) :: X_new
        double precision :: ap,ae,aw,an,as,at,ab,X_max
        integer :: i,n

        ! BiCGStab Variables
        double precision, dimension(nz,ny,0:m+1) :: r, r0_star, p, s, r_new, p_new,  AX_temp
        double precision :: temp_loc, temp1, temp2
        double precision :: alpha, beta, omega


        integer :: id, Nproc, ierr

        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

		!call setPressureBC_MPI(X)
        !call communicate(X)
        call DIV_All(DIVU,u,v,w)

		!ap = 2/dx**2 + 2/dy**2
		!ae = 1/dx**2;
		!aw = 1/dx**2;
		!an = 1/dy**2;
		!as = 1/dy**2;
		b = -rho*(DIVU)/dt; ! dt affects possibility of convergence

        ! for BC's
	b(:,1,:) = 0; b(:,ny,:) = 0 ! lateral sides, where dPdn=0 -> b=0
        b(1,:,:) = 0; b(nz,:,:) = 0 ! 
        ! for cavity, left and right boundaries
        !if (id==0) b(:,1) = 0
        !if (id==Nproc-1) b(:,m) = 0


		n = nz*nx*ny

		!P_new = P

        r = b - AX(X)
        r0_star = r
        p = r
		do i=1,nit
			!print*, i
            
            AX_temp = AX(p)
            temp_loc = dot_product(pack(r(:,:,1:m),.true.),pack(r0_star(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(AX_temp(:,:,1:m),.true.),pack(r0_star(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            alpha = temp1/temp2

            !print*, temp_loc, temp1, temp2, alpha
            !if(i==50) stop

            s = r - alpha * AX(p)
            call communicate(s)

            AX_temp = AX(s)
            temp_loc = dot_product(pack(AX_temp(:,:,1:m),.true.),pack(s(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(AX_temp(:,:,1:m),.true.),pack(AX_temp(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            omega = temp1/temp2

            ! update X
            X_new = X + alpha*p + omega*s
            call communicate(X_new)

            ! update r
            r_new = s - omega*AX(s)
            call communicate(r_new)

            temp_loc = dot_product(pack(r_new(:,:,1:m),.true.),pack(r0_star(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            temp_loc = dot_product(pack(r(:,:,1:m),.true.),pack(r0_star(:,:,1:m),.true.))
            call mpi_allreduce(temp_loc, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            beta = (temp1/temp2)*(alpha/omega)

            ! update p
            p_new = r_new + beta * (p-omega*AX(p))
            call communicate(p_new)




            ! COMMUNICATE(P)
            !call communicate(X_new)

			!call setPressureBC_MPI(X_new)


			if (mod(i,1)==0) then

				E = X_new(:,:,1:m) - X(:,:,1:m)
				E_sqr_loc = sum(E**2)

                call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                E_rms = sqrt(E_sqr/n)

                !print*, 'Pressure Error:',E_rms

				if (E_rms<tol) then
					exit
				end if

			end if

            ! MPI_Barrier ??

			!X = X_new 
            !r = r_new
            !p = p_new

			X = X + omega_P*(X_new-X)
			r = r + omega_P*(r_new-r)
			p = p + omega_P*(p_new-p)

		end do
		
		NITS_final = i
		
		! R = ap*P(2:nz-1,2:ny-1,2:nx-1)-(an*P(2:nz-1,3:ny,2:nx-1)+as*P(2:nz-1,1:ny-2,2:nx-1) &
		!	+ae*P(2:nz-1,2:ny-1,3:nx)+aw*P(2:nz-1,2:ny-1,1:nx-2) &
		!	+at*P(3:nz,2:ny-1,2:nx-1)+ab*P(1:nz-2,2:ny-1,2:nx-1)+b);
		! R_rms = sqrt(sum(R**2)/n)



    end subroutine PressureSolver_MPI_BiCGStab



    function AX(X)
        use settings, only: nx, ny, nz, dx, dy, dz, dt

        double precision, dimension(nz,ny,0:m+1) :: X
        double precision, dimension(nz,ny,0:m+1) :: AX
        double precision :: ap, ae, aw, an, as, at, ab

        integer :: id, Nproc, ierr

        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    
        ap = 2/dx**2 + 2/dy**2 + 2/dz**2
	ae = 1/dx**2;
	aw = 1/dx**2;
	an = 1/dy**2;
	as = 1/dy**2;
	at = 1/dz**2
	ab = 1/dz**2

        AX(2:nz-1,2:ny-1,1:m) = ap*X(2:nz-1,2:ny-1,1:m) - ( an*X(2:nz-1,3:ny,1:m) + as*X(2:nz-1,1:ny-2,1:m) + aw*X(2:nz-1,2:ny-1,0:m-1) + ae*X(2:nz-1,2:ny-1,2:m+1) + at*X(3:nz,2:ny-1,1:m) + ab*X(1:nz-2,2:ny-1,1:m) )
        AX(2:nz-1,1,1:m) = X(2:nz-1,1,1:m) - X(2:nz-1,2,1:m)
        AX(2:nz-1,ny,1:m) = X(2:nz-1,ny,1:m) - X(2:nz-1,ny-1,1:m)
	AX(1,2:ny-1,1:m) = X(1,2:ny-1,1:m) - X(2,2:ny-1,1:m)
	AX(nz,2:ny-1,1:m) = X(nz,2:ny-1,1:m) - X(nz-1,2:ny-1,1:m)

        ! for cavity case, bc at left and right edges
        !if (id==0) AX(:,1) = X(:,1) - X(:,2)
        !if (id==Nproc-1) AX(:,m) = X(:,m) - X(:,m-1)

        call communicate(AX)

    end function AX



    subroutine calcUstar_AB(u_star,v_star,w_star,u,v,w,u_last,v_last,w_last)
        use settings, only: nx, ny, nz, dx, dy, dz, dt, rho, nit, omega_P, tol

        double precision, dimension(nz,ny,0:m+1) :: u,v,w,u_star,v_star,w_star,u_last,v_last, w_last
        double precision, dimension(2:nz-1,2:ny-1,1:m) :: SU,SV,SW,SU_last,SV_last,SW_last

        double precision, dimension(2:nz-1,2:ny-1,1:m) :: D_CONV_U,D_CONV_V,D_CONV_W,D_CONV_U_last, D_CONV_V_last,D_CONV_W_last
        double precision, dimension(2:nz-1,2:ny-1,1:m) :: DEL2U,DEL2V,DEL2W,DEL2U_last, DEL2V_last, DEL2W_last

        call CONV_U_UPWIND(D_CONV_U,u,v,w)
        call CONV_V_UPWIND(D_CONV_V,u,v,w)
        call CONV_W_UPWIND(D_CONV_W,u,v,w)
        call CONV_U_UPWIND(D_CONV_U_last,u_last,v_last,w_last)
        call CONV_V_UPWIND(D_CONV_V_last,u_last,v_last,w_last)
        call CONV_W_UPWIND(D_CONV_W_last,u_last,v_last,w_last)

        call DEL2(DEL2U,u)
        call DEL2(DEL2V,v)
        call DEL2(DEL2W,w)
        call DEL2(DEL2U_last,u_last)
        call DEL2(DEL2V_last,v_last)
        call DEL2(DEL2W_last,w_last)

        SU = -D_CONV_U + nu*DEL2U
        SV = -D_CONV_V + nu*DEL2V
        SW = -D_CONV_W + nu*DEL2W
        SU_last = -D_CONV_U_last + nu*DEL2U_last
        SV_last = -D_CONV_V_last + nu*DEL2V_last
        SW_last = -D_CONV_W_last + nu*DEL2W_last

        u_star(2:nz-1,2:ny-1,1:m) = u(2:nz-1,2:ny-1,1:m) + dt * (1.5*SU-0.5*SU_last)
        v_star(2:nz-1,2:ny-1,1:m) = v(2:nz-1,2:ny-1,1:m) + dt * (1.5*SV-0.5*SV_last)
        w_star(2:nz-1,2:ny-1,1:m) = w(2:nz-1,2:ny-1,1:m) + dt * (1.5*SW-0.5*SW_last)

        call setVelocityBC_MPI(u,v,w)

    end subroutine
	

	
	subroutine calcUstar_Implicit(u_star,v_star,w_star,E_rms,NITS_final,u,v,w,Bx,By,Bz)
	
		double precision, dimension(nz,ny,0:m+1) :: u, v, w, u_star,v_star,w_star
		double precision, dimension(nz,ny,0:m+1) :: u_star_new,v_star_new,w_star_new
		
		double precision ,dimension(nz-2,ny-2,1:m) :: D_CONV_U, D_CONV_V, D_CONV_W, DEL2U, DEL2V, DEL2W
		
		double precision ,dimension(nz,ny,1:m) :: Bx, By, Bz 
		
		double precision, dimension(nz,ny,1:m) :: E
		double precision :: E_rms, E_sqr_loc, E_sqr
        integer :: NITS_final
		integer :: i, n, Nproc, ierr, id

		call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)
		
		
		n = nx*ny*nz
		
		u_star = u
		v_star = v
		w_star = w
		
		u_star_new = 0
		v_star_new = 0
		w_star_new = 0
		
		do i=1,nit
			call CONV_U_UPWIND(D_CONV_U,u_star,v_star,w_star)
			call CONV_V_UPWIND(D_CONV_V,u_star,v_star,w_star)
			call CONV_W_UPWIND(D_CONV_W,u_star,v_star,w_star)
	
			call DEL2(DEL2U,u_star)
			call DEL2(DEL2V,v_star)
			call DEL2(DEL2W,w_star)

			u_star_new(2:nz-1,2:ny-1,1:m) = u(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_U + nu*DEL2U + Bx(2:nz-1,2:ny-1,:))
			v_star_new(2:nz-1,2:ny-1,1:m) = v(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_V + nu*DEL2V + By(2:nz-1,2:ny-1,:))
			w_star_new(2:nz-1,2:ny-1,1:m) = w(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_W + nu*DEL2W + Bz(2:nz-1,2:ny-1,:))

			call setVelocityBC_MPI(u_star_new,v_star_new,w_star_new)
			
			! COMMUNICATE (u_star)
			call communicate(u_star_new)
			call communicate(v_star_new)
			call communicate(w_star_new)
			
			E = u_star_new(:,:,1:m) - u_star(:,:,1:m)
			E_sqr_loc = sum(E**2)

			call mpi_allreduce(E_sqr_loc,E_sqr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
			E_rms = sqrt(E_sqr/n)

			!if (id==0) print*, 'U_star error:', E_rms
			
			if (E_rms<tol) then
				exit
			end if
			
			u_star = u_star + omega_U * (u_star_new - u_star)
			v_star = v_star + omega_V * (v_star_new - v_star)
			w_star = w_star + omega_W * (w_star_new - w_star)

			
		end do
		
		NITS_final = i
	
	end subroutine calcUstar_Implicit


end module solver
