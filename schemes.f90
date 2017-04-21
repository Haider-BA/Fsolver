module schemes
use settings, only: nx, ny, nz, dx, dy, dz
use partitioner, only : m
use communicator

implicit none
    
contains
    subroutine CONV_U_UPWIND(D,u,v,w)
    
		double precision, dimension(nz,ny,0:m+1) :: u,v,w
		double precision, dimension(nz-2,ny-2,1:m) :: D 
		
		double precision, dimension(nz-2,ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus,w_plus,w_minus
		double precision, dimension(nz-2,ny-2,1:m) :: Du2Dx,DuvDy,DuwDz
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;w_plus=0;w_minus=0

		where (u(2:nz-1,2:ny-1,1:m)>0) u_plus = 1; where (v(2:nz-1,2:ny-1,1:m)>=0) v_plus = 1; where (w(2:nz-1,2:ny-1,1:m)>=0) w_plus = 1
		where (u(2:nz-1,2:ny-1,1:m)<0) u_minus = 1; where (v(2:nz-1,2:ny-1,1:m)<0) v_minus = 1; where (w(2:nz-1,2:ny-1,1:m)<0) w_minus = 1

		Du2Dx = u_plus*(u(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m)-u(2:nz-1,2:ny-1,0:m-1)*u(2:nz-1,2:ny-1,0:m-1))/dx &
		    + u_minus*(u(2:nz-1,2:ny-1,2:m+1)*u(2:nz-1,2:ny-1,2:m+1)-u(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m))/dx
		DuvDy = v_plus*(u(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m)-u(2:nz-1,1:ny-2,1:m)*v(2:nz-1,1:ny-2,1:m))/dy &
		    + v_minus*(u(2:nz-1,3:ny,1:m)*v(2:nz-1,3:ny,1:m)-u(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m))/dy
		DuwDz = w_plus*(u(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m)-u(1:nz-2,2:ny-1,1:m)*w(1:nz-2,2:ny-1,1:m))/dz &
		    + w_minus*(u(3:nz,2:ny-1,1:m)*w(3:nz,2:ny-1,1:m)-u(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m))/dz

		D = Du2Dx + DuvDy + DuwDz

    end subroutine CONV_U_UPWIND


    subroutine CONV_V_UPWIND(D,u,v,w)
    
		double precision, dimension(nz,ny,0:m+1) :: u,v,w
		double precision, dimension(nz-2,ny-2,1:m) :: D 
		
		double precision, dimension(nz-2,ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus,w_plus,w_minus
		double precision, dimension(nz-2,ny-2,1:m) :: DvuDx,Dv2Dy,DvwDz
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;w_plus=0;w_minus=0

		where (u(2:nz-1,2:ny-1,1:m)>0) u_plus = 1; where (v(2:nz-1,2:ny-1,1:m)>=0) v_plus = 1; where (w(2:nz-1,2:ny-1,1:m)>=0) w_plus = 1
		where (u(2:nz-1,2:ny-1,1:m)<0) u_minus = 1; where (v(2:nz-1,2:ny-1,1:m)<0) v_minus = 1; where (w(2:nz-1,2:ny-1,1:m)<0) w_minus = 1

		DvuDx = u_plus*(v(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m)-v(2:nz-1,2:ny-1,0:m-1)*u(2:nz-1,2:ny-1,0:m-1))/dx &
		    + u_minus*(v(2:nz-1,2:ny-1,2:m+1)*u(2:nz-1,2:ny-1,2:m+1)-v(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m))/dx
		Dv2Dy = v_plus*(v(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m)-v(2:nz-1,1:ny-2,1:m)*v(2:nz-1,1:ny-2,1:m))/dy &
		    + v_minus*(v(2:nz-1,3:ny,1:m)*v(2:nz-1,3:ny,1:m)-v(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m))/dy
		DvwDz = w_plus*(v(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m)-v(1:nz-2,2:ny-1,1:m)*w(1:nz-2,2:ny-1,1:m))/dz &
		    + w_minus*(v(3:nz,2:ny-1,1:m)*w(3:nz,2:ny-1,1:m)-v(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m))/dz

		D = DvuDx + Dv2Dy + DvwDz

    end subroutine CONV_V_UPWIND
    
    
    subroutine CONV_W_UPWIND(D,u,v,w)
    
		double precision, dimension(nz,ny,0:m+1) :: u,v,w
		double precision, dimension(nz-2,ny-2,1:m) :: D 
		
		double precision, dimension(nz-2,ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus,w_plus,w_minus
		double precision, dimension(nz-2,ny-2,1:m) :: DwuDx,DwvDy,Dw2Dz
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;w_plus=0;w_minus=0

		where (u(2:nz-1,2:ny-1,1:m)>0) u_plus = 1; where (v(2:nz-1,2:ny-1,1:m)>=0) v_plus = 1; where (w(2:nz-1,2:ny-1,1:m)>=0) w_plus = 1
		where (u(2:nz-1,2:ny-1,1:m)<0) u_minus = 1; where (v(2:nz-1,2:ny-1,1:m)<0) v_minus = 1; where (w(2:nz-1,2:ny-1,1:m)<0) w_minus = 1

		DwuDx = u_plus*(w(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m)-w(2:nz-1,2:ny-1,0:m-1)*u(2:nz-1,2:ny-1,0:m-1))/dx &
		    + u_minus*(w(2:nz-1,2:ny-1,2:m+1)*u(2:nz-1,2:ny-1,2:m+1)-w(2:nz-1,2:ny-1,1:m)*u(2:nz-1,2:ny-1,1:m))/dx
		DwvDy = v_plus*(w(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m)-w(2:nz-1,1:ny-2,1:m)*v(2:nz-1,1:ny-2,1:m))/dy &
		    + v_minus*(w(2:nz-1,3:ny,1:m)*v(2:nz-1,3:ny,1:m)-w(2:nz-1,2:ny-1,1:m)*v(2:nz-1,2:ny-1,1:m))/dy
		Dw2Dz = w_plus*(w(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m)-w(1:nz-2,2:ny-1,1:m)*w(1:nz-2,2:ny-1,1:m))/dz &
		    + w_minus*(w(3:nz,2:ny-1,1:m)*w(3:nz,2:ny-1,1:m)-w(2:nz-1,2:ny-1,1:m)*w(2:nz-1,2:ny-1,1:m))/dz

		D = DwuDx + DwvDy + Dw2Dz

    end subroutine CONV_W_UPWIND
    
    
    subroutine CONV_C_UPWIND(D,u,v,w,c)
    
		double precision, dimension(nz,ny,0:m+1) :: u,v,w,C
		double precision, dimension(nz-2,ny-2,1:m) :: D 
		
		double precision, dimension(nz-2,ny-2,1:m) :: u_plus,u_minus,v_plus,v_minus,w_plus,w_minus
		double precision, dimension(nz-2,ny-2,1:m) :: DucDx,DvcDy,DwcDz
		
		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;w_plus=0;w_minus=0

		where (u(2:nz-1,2:ny-1,1:m)>0) u_plus = 1; where (v(2:nz-1,2:ny-1,1:m)>=0) v_plus = 1; where (w(2:nz-1,2:ny-1,1:m)>=0) w_plus = 1
		where (u(2:nz-1,2:ny-1,1:m)<0) u_minus = 1; where (v(2:nz-1,2:ny-1,1:m)<0) v_minus = 1; where (w(2:nz-1,2:ny-1,1:m)<0) w_minus = 1

		DucDx = u_plus*(u(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m)-u(2:nz-1,2:ny-1,0:m-1)*c(2:nz-1,2:ny-1,0:m-1))/dx &
		    + u_minus*(u(2:nz-1,2:ny-1,2:m+1)*c(2:nz-1,2:ny-1,2:m+1)-u(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m))/dx
		DvcDy = v_plus*(v(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m)-v(2:nz-1,1:ny-2,1:m)*c(2:nz-1,1:ny-2,1:m))/dy &
		    + v_minus*(v(2:nz-1,3:ny,1:m)*c(2:nz-1,3:ny,1:m)-v(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m))/dy
		DwcDz = w_plus*(w(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m)-w(1:nz-2,2:ny-1,1:m)*c(1:nz-2,2:ny-1,1:m))/dz &
		    + w_minus*(w(3:nz,2:ny-1,1:m)*c(3:nz,2:ny-1,1:m)-w(2:nz-1,2:ny-1,1:m)*c(2:nz-1,2:ny-1,1:m))/dz

		D = DucDx + DvcDy + DwcDz

    end subroutine CONV_C_UPWIND



	subroutine DEL2(D,T)

		double precision, dimension(nz,ny,0:m+1) :: T
		double precision, dimension(nz-2,ny-2,1:m) :: D
		
		double precision, dimension(nz-2,ny-2,1:m) :: D2TDx2, D2TDy2, D2TDz2
		
		D2TDx2 = (T(2:nz-1,2:ny-1,2:m+1)-2*T(2:nz-1,2:ny-1,1:m)+T(2:nz-1,2:ny-1,0:m-1))/dx**2
		D2TDy2 = (T(2:nz-1,3:ny,1:m)-2*T(2:nz-1,2:ny-1,1:m)+T(2:nz-1,1:ny-2,1:m))/dy**2
    D2TDz2 = (T(3:nz,2:ny-1,1:m)-2*T(2:nz-1,2:ny-1,1:m)+T(1:nz-2,2:ny-1,1:m))/dz**2

		D = D2TDx2 + D2TDy2 + D2TDz2

	end subroutine DEL2
	
	
	subroutine DIV(D,u,v,w)
	
		double precision, dimension(nz,ny,0:m+1) :: u,v,w
		double precision, dimension(nz-2,ny-2,1:m) :: D 
		
		double precision, dimension(nz-2,ny-2,1:m) :: DuDx,DvDy,DwDz

        ! Central difference

		! DuDx
		DuDx = (u(2:nz-1,2:ny-1,2:m+1)-u(2:nz-1,2:ny-1,0:m-1))/(2*dx)

		! DvDy
		DvDy = (v(2:nz-1,3:ny,1:m)-v(2:nz-1,1:ny-2,1:m))/(2*dy)

		! DwDz
		DwDz = (w(3:nz,2:ny-1,1:m)-w(1:nz-2,2:ny-1,1:m))/(2*dz)

		D = DuDx + DvDy + DwDz

	end subroutine DIV
	
	
    !------------------------------------------------------
    ! Divergence at all points including boundary nodes
    !------------------------------------------------------
    subroutine DIV_All(D,u,v,w)
	
		double precision, dimension(nz,ny,0:m+1) :: u,v,w
		double precision, dimension(nz,ny,0:m+1) :: D 
		
		double precision, dimension(nz,ny,0:m+1) :: DuDx,DvDy,DwDz

        ! Central difference at internal nodes and forward/backward difference
        ! at boundary nodes

		! DuDx
		DuDx(:,:,1:m) = (u(:,:,2:m+1)-u(:,:,0:m-1))/(2*dx)

		! DvDy
		DvDy(:,2:ny-1,1:m) = (v(:,3:ny,1:m)-v(:,1:ny-2,1:m))/(2*dy)
        	DvDy(:,1,1:m) = (v(:,2,1:m)-v(:,1,1:m))/dy
        	DvDy(:,ny,1:m) = (v(:,ny,1:m)-v(:,ny-1,1:m))/dy

		! DwDz
		DwDz(2:nz-1,:,1:m) = (w(3:nz,:,1:m)-w(1:nz-2,:,1:m))/(2*dz)
        	DwDz(1,:,1:m) = (w(2,:,1:m)-w(2,:,1:m))/dz
        	DwDz(nz,:,1:m) = (w(nz,:,1:m)-w(nz-1,:,1:m))/dz


		D = DuDx + DvDy + DwDz

        call communicate(D)

	end subroutine DIV_All


	subroutine GRAD(DpDx,DpDy,DpDz,P)
	
		double precision, dimension(nz,ny,0:m+1) :: P
		double precision, dimension(nz-2,ny-2,1:m) :: DpDx,DpDy,DpDz

        ! Central difference

		DpDx = (P(2:nz-1,2:ny-1,2:m+1)-P(2:nz-1,2:ny-1,0:m-1))/(2*dx)
		DpDy = (P(2:nz-1,3:ny,1:m)-P(2:nz-1,1:ny-2,1:m))/(2*dy)
		DpDz = (P(3:nz,2:ny-1,1:m)-P(1:nz-2,2:ny-1,1:m))/(2*dz)

	end subroutine GRAD
	
	
	
	subroutine GRAD_UPWIND(DpDx,DpDy,DpDz,P,u,v,w)
	
		double precision, dimension(nz,ny,nx) :: P,u,v,w
		double precision, dimension(nz-2,ny-2,nx-2) :: DpDx,DpDy,DpDz
		
		double precision, dimension(nz-2,ny-2,nx-2) :: u_plus,u_minus,v_plus,v_minus,w_plus,w_minus

		u_plus=0;u_minus=0;v_plus=0
		v_minus=0;w_plus=0;w_minus=0

		where (u(2:nz-1,2:ny-1,2:nx-1)>0) u_plus = 1; where (v(2:nz-1,2:ny-1,2:nx-1)>=0) v_plus = 1; where (w(2:nz-1,2:ny-1,2:nx-1)>=0) w_plus = 1
		where (u(2:nz-1,2:ny-1,2:nx-1)<0) u_minus = 1; where (v(2:nz-1,2:ny-1,2:nx-1)<0) v_minus = 1; where (w(2:nz-1,2:ny-1,2:nx-1)<0) w_minus = 1

		DpDx = u_plus*(P(2:nz-1,2:ny-1,2:nx-1)-P(2:nz-1,2:ny-1,1:nx-2))/dx &
			+ u_minus*(P(2:nz-1,2:ny-1,3:nx)-P(2:nz-1,2:ny-1,2:nx-1))/dx
		DpDy = v_plus*(P(2:nz-1,2:ny-1,2:nx-1)-P(2:nz-1,1:ny-2,2:nx-1))/dy &
			+ v_minus*(P(2:nz-1,3:ny,2:nx-1)-P(2:nz-1,2:ny-1,2:nx-1))/dy
		DpDz = w_plus*(P(2:nz-1,2:ny-1,2:nx-1)-P(1:nz-2,2:ny-1,2:nx-1))/dz &
			+ w_minus*(P(3:nz,2:ny-1,2:nx-1)-P(2:nz-1,2:ny-1,2:nx-1))/dz

	end subroutine GRAD_UPWIND
	
	
	
	


end module schemes
