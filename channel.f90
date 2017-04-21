subroutine channel()
    use settings
    use partitioner
	use io
	use schemes
	use solver
	use boundaryConditions
	use fcm

	implicit none

	
    double precision ,dimension(Np,3) :: Yp, Vp, Vp0		! Particle variables

	double precision ,dimension(nz-2,ny-2,1:m) :: D_CONV_U, D_CONV_V, D_CONV_W, DEL2U, DEL2V, DEL2W, DpDx, DpDy, DpDz
	double precision ,dimension(nz-2,ny-2,1:m) :: D_CONV_C_FBG, D_CONV_C_FN, D_CONV_C_TH, DEL2C_FBG, DEL2C_FN, DEL2C_TH
	
    double precision :: R, E_P, E_U
	integer :: NITS_final_P, NITS_final_U

	integer :: i, j

    integer :: id, Nproc, ierr

    ! for io
    double precision :: Umax, Pmax, Umax_loc, Pmax_loc

    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)

    
	! Initialize Flow Variables
	field%u=0; field%v=0; field%w=0; field%P=0
	field%c_fbg = 0; field%c_fn = 0; field%c_th = 0
	field%u_star=0; field%v_star=0; field%w_star=0
	field%u_star2=0; field%v_star2=0; field%w_star2=0

    !print*, 'variables intialized'

    !if(id==0) print*, 'flow initialized.'

    call setParticleVariables(Yp,Vp,Vp0)

    !print*, 'particle vars initialized'

	! Pressure driven force
	forces%F = 4*nu*U0/RADIUS**2

	! construct x,y,z matrices (mesh)
	do i=1,m; mesh%x(:,:,i) = (iStartx+i-1)*dx;	end do
	do i=1,ny; mesh%y(:,i,:) = (i-1)*dy;	end do
	do i=1,nz; mesh%z(i,:,:) = (i-1)*dz;	end do

    !print*, 'mesh set'

	! Velocity IC
	call setIC(field%u,field%v,field%w,field%P,field%c_fbg,field%c_fn,field%c_th,mesh%x,mesh%y,mesh%z)

    !print*, 'IC set'

	![u,v] = setVelocityBC1(u,v,dx,dy,dt,nu);

	! write data before starting simulation
	call writeToTecplot_MPI(mesh%x,mesh%y,mesh%z,field%u,field%v,field%w,field%P,field%c_fbg,field%c_fn,field%c_th,Yp,Vp,0)
	
	do i=1,nt
		
    !print*, 'Timeloop', i, 'started'

		! Calculate force for cylinder IBM
		call cylinderIBMForce(forces%Hx,forces%Hy,forces%Hz,field%u,field%v,field%w,mesh%x,mesh%y,mesh%z)
		
		! Calculate force for FCM
		call forceCoupling(forces%Gx,forces%Gy,forces%Gz,Yp)
		!forces%Gx=0; forces%Gy=0; forces%Gz=0;

		! Total body force
		forces%Bx = forces%F + forces%Hx + forces%Gx
		forces%By = forces%Hy + forces%Gy
		forces%Bz = forces%Hz + forces%Gz
		
		call calcUstar_Implicit(field%u_star,field%v_star,field%w_star,E_U,NITS_final_U,field%u,field%v,field%w,forces%Bx,forces%By,forces%Bz)

		!call CONV_U_UPWIND(D_CONV_U,field%u,field%v,field%w)
		!call CONV_V_UPWIND(D_CONV_V,field%u,field%v,field%w)
		!call CONV_W_UPWIND(D_CONV_W,field%u,field%v,field%w)
		
		
	
		!call DEL2(DEL2U,field%u)
		!call DEL2(DEL2V,field%v)
		!call DEL2(DEL2W,field%w)
		
		

		!field%u_star(2:nz-1,2:ny-1,1:m) = field%u(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_U + nu*DEL2U)
		!field%v_star(2:nz-1,2:ny-1,1:m) = field%v(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_V + nu*DEL2V)
		!field%w_star(2:nz-1,2:ny-1,1:m) = field%w(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_W + nu*DEL2W)
		
		
        
        !print*, 'u_star done'

		! call setVelocityBC_MPI(u,v,w)
		
		!call calcUstar(u_star,v_star,w_star,u,v,w)

        ! COMMUNICATE (u_star)
        !call communicate(field%u_star)
        !call communicate(field%v_star)
        !call communicate(field%w_star)

        !if(id==0) print*, 'communication done'

		call PressureSolver_MPI(field%P,E_P,NITS_final_P,field%u_star,field%v_star,field%w_star)
        !print*, 'preconditioner done'
		call PressureSolver_MPI_BiCGStab(field%P,E_P,NITS_final_P,field%u_star,field%v_star,field%w_star)
        !print*, 'BiCGStab solver done'

        !if(id==0) print*, 'Pressure solver done: ', NITS_final
		call GRAD(DpDx,DpDy,DpDz,field%P)

        !if(id==0) print*, 'GRAD calculated'

		field%u_star2(2:nz-1,2:ny-1,1:m) = field%u_star(2:nz-1,2:ny-1,1:m) - dt/rho*DpDx
		field%v_star2(2:nz-1,2:ny-1,1:m) = field%v_star(2:nz-1,2:ny-1,1:m) - dt/rho*DpDy
		field%w_star2(2:nz-1,2:ny-1,1:m) = field%w_star(2:nz-1,2:ny-1,1:m) - dt/rho*DpDz

        !if(id==0) print*, 'u_star2 calculated'

		! call setVelocityBC2(u_star2,v_star2,w_star2,u_star,v_star,w_star,P)

		

        !if(id==0) print*, 'IBM done'

    	

        !if(id==0) print*, 'fcm done'

		! march in time using force
		! (Entire x dimension is used, since periodic BC in that direction)
		!field%u(2:nz-1,2:ny-1,1:m) = field%u_star2(2:nz-1,2:ny-1,1:m) + dt*( forces%F(2:nz-1,2:ny-1,:) + forces%Hx(2:nz-1,2:ny-1,:) + forces%Gx(2:nz-1,2:ny-1,:) )
		!field%v(2:nz-1,2:ny-1,1:m) = field%v_star2(2:nz-1,2:ny-1,1:m) + dt*( forces%Hy(2:nz-1,2:ny-1,:) + forces%Gy(2:nz-1,2:ny-1,:) )
		!field%w(2:nz-1,2:ny-1,1:m) = field%w_star2(2:nz-1,2:ny-1,1:m) + dt*( forces%Hz(2:nz-1,2:ny-1,:) + forces%Gz(2:nz-1,2:ny-1,:) )
		
		field%u(2:nz-1,2:ny-1,1:m) = field%u_star2(2:nz-1,2:ny-1,1:m)
		field%v(2:nz-1,2:ny-1,1:m) = field%v_star2(2:nz-1,2:ny-1,1:m)
		field%w(2:nz-1,2:ny-1,1:m) = field%w_star2(2:nz-1,2:ny-1,1:m)
        
        ! COMMUNICATE (u)
        call communicate(field%u)
        call communicate(field%v)
        call communicate(field%w)


        ! Calculate Concentration fields
        !call CONV_C_UPWIND(D_CONV_C_FBG,field%u,field%v,field%w,field%c_fbg)
		!call CONV_C_UPWIND(D_CONV_C_FN,field%u,field%v,field%w,field%c_fn)
		!call CONV_C_UPWIND(D_CONV_C_TH,field%u,field%v,field%w,field%c_th)
		
		!call DEL2(DEL2C_FBG,field%c_fbg)
		!call DEL2(DEL2C_FN,field%c_fn)
		!call DEL2(DEL2C_TH,field%c_th)
		
		!field%c_fbg(2:nz-1,2:ny-1,1:m) = field%c_fbg(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_C_FBG + Dt_fbg*DEL2C_FBG)
		!field%c_fn(2:nz-1,2:ny-1,1:m) = field%c_fn(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_C_FN + Dt_fn*DEL2C_FN)
		!field%c_th(2:nz-1,2:ny-1,1:m) = field%c_th(2:nz-1,2:ny-1,1:m) + dt*(-D_CONV_C_TH + Dt_th*DEL2C_TH)

		! COMMUNICATE (c)
        !call communicate(field%c_fbg)
        !call communicate(field%c_fn)
        !call communicate(field%c_th)
        
		!if(id==0) print*, 'Before BC reached'
		
        ! Set Boundary conditions
        !call setConcBC_MPI(field%c_fbg, field%c_fn, field%c_th, mesh%x, mesh%y, mesh%z)

        !if(id==0) print*, 'field updated'

		call updateParticles(Yp,Vp,Vp0,field%u,field%v,field%w)
	
        !if(id==0) print*, 'update particles done'

        Umax_loc = maxval(field%u(:,:,1:m))
        Pmax_loc = maxval(field%P(:,:,1:m))
        call mpi_reduce(Umax_loc,Umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call mpi_reduce(Pmax_loc,Pmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        !if(id==0) print*, 'reduction for io done'

        if (id==0) then       
            if (mod(i,nLog)==0) then
                call writeLogToTerminal_MPI(i, &
                'Umax','double',Umax, &
                'Pmax','double',Pmax, &
                'Presure Error','double',E_P, &
                'Pressure Iterations','integer',dble(NITS_final_P),&
				'Velocity Error','double',E_U, &
				'Velocity Iterations','integer',dble(NITS_final_U),&
                'particle 1 position','double',Yp(1,1))
            end if
            
            
        end if

        if (mod(i,nWrite)==0) then
            !call writeToBOV_MPI(field%u,field%v,field%w,field%P,Yp,i)
            !call writeParticleData(Yp,i)
            call writeToTecplot_MPI(mesh%x,mesh%y,mesh%z,field%u,field%v,field%w,field%P,field%c_fbg,field%c_fn,field%c_th,Yp,Vp,i)
        end if

        call mpi_barrier(MPI_COMM_WORLD,ierr)
		
	end do

	!call writeFlowDataToFile(u,v,w,P,i)
	!call writeToTecplot(x,y,z,u,v,w,P,i)

end subroutine channel
