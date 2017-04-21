subroutine cylinderIBMForce(Hx,Hy,Hz,u,v,w,x,y,z)

	use settings, only : nz,ny,Lz,Ly,RADIUS,dt
    use partitioner, only : m

	implicit none

	double precision, dimension(nz,ny,m) :: Hx,Hy,Hz, x,y,z
	double precision, dimension(nz,ny,0:m+1) :: u,v,w
	
	integer, dimension(nz,ny,m) :: eta
	double precision :: yc,zc


	!yc=(maxval(y(1,:,1))-minval(y(1,:,1)))/2
	!zc=(maxval(z(:,1,1))-minval(z(:,1,1)))/2

    yc = Ly/2
    zc = Lz/2

	eta = 0
	where ((y-yc)**2 + (z-zc)**2 > RADIUS**2 ) eta=1

	Hx = -eta*u(:,:,1:m)/(dt)
	Hy = -eta*v(:,:,1:m)/(dt)
	Hz = -eta*w(:,:,1:m)/(dt)

end subroutine cylinderIBMForce
