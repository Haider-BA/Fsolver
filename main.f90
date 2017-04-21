program main

    use mpi
    use settings
    use partitioner
    use io, only : printHeader

	implicit none

    integer :: ierr, Nproc, id

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)

    if (Nproc==1) then
        print*, ''
        print*, 'Error: This is a parallel program. Select at least two processors!'
        print*, ''
        call exit()
    end if

	call setup('setup.txt')

    !call mpi_barrier(MPI_COMM_WORLD,ierr)

    if (id==0) then
        call printHeader()
    end if


    call setPartitions()
	
	call channel()

    call mpi_finalize(ierr)

end program main
