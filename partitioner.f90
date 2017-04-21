module partitioner 

use mpi
use settings
implicit none


! Data structures
type :: mesh_t
    double precision, dimension(:,:,:), allocatable :: x, y, z
end type mesh_t

type :: field_t
    double precision, dimension(:,:,:), allocatable :: u, v, w, P,   c_fbg, c_fn, c_th
    double precision, dimension(:,:,:), allocatable :: u_star, v_star, w_star
    double precision, dimension(:,:,:), allocatable :: u_star2, v_star2, w_star2
end type field_t

type :: forces_t
    double precision, dimension(:,:,:), allocatable :: F            ! Pressure driven force
    double precision, dimension(:,:,:), allocatable :: Hx, Hy, Hz   ! Cylinder IBM force
    double precision, dimension(:,:,:), allocatable :: Gx, Gy, Gz   ! Force coupling for platelets
	double precision, dimension(:,:,:), allocatable :: Bx, By, Bz	! Total body force (sum of all other forces)
end type forces_t


! global data variables
type(mesh_t) :: mesh
type(field_t) :: field
type(forces_t) :: forces

! Number of nodes in x-direction on this processor
integer :: m

! global node number of first vertex on this processor in x direction
integer :: iStartx

contains

subroutine setPartitions()

    integer :: Nproc, id, ierr

    call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    

    m = Nx/Nproc
    if (id<(Nx-m*Nproc)) then
        iStartx = id*(m+1)+1
        m = m+1
    else
        iStartx = id*m + (Nx-m*Nproc) +1
    end if

    allocate(mesh%x(Nz,Ny,m))
    allocate(mesh%y(Nz,Ny,m))
    allocate(mesh%z(Nz,Ny,m))

    allocate(field%u(Nz,Ny,0:m+1))
    allocate(field%v(Nz,Ny,0:m+1))
    allocate(field%w(Nz,Ny,0:m+1))
    allocate(field%P(Nz,Ny,0:m+1))
    
    allocate(field%c_fbg(Nz,Ny,0:m+1))
    allocate(field%c_fn(Nz,Ny,0:m+1))
    allocate(field%c_th(Nz,Ny,0:m+1))
    
    allocate(field%u_star(Nz,Ny,0:m+1))
    allocate(field%v_star(Nz,Ny,0:m+1))
    allocate(field%w_star(Nz,Ny,0:m+1))
    
    allocate(field%u_star2(Nz,Ny,0:m+1))
    allocate(field%v_star2(Nz,Ny,0:m+1))
    allocate(field%w_star2(Nz,Ny,0:m+1))

    allocate(forces%F(Nz,Ny,m))
    allocate(forces%Gx(Nz,Ny,m))
    allocate(forces%Gy(Nz,Ny,m))
    allocate(forces%Gz(Nz,Ny,m))
    allocate(forces%Hx(Nz,Ny,m))
    allocate(forces%Hy(Nz,Ny,m))
    allocate(forces%Hz(Nz,Ny,m))
	allocate(forces%Bx(Nz,Ny,m))
    allocate(forces%By(Nz,Ny,m))
    allocate(forces%Bz(Nz,Ny,m))

    if (id==0) then
        write(*,'(A,I0,A)') 'Memory allocated on ', Nproc, ' processors'
    end if

end subroutine setPartitions


end module partitioner
