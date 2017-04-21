program main
    implicit none

    integer, parameter :: PI = 4*atan(1.0)
    integer, parameter :: conc = 250000     ! /uL
    integer :: N                            ! number of platelets in domain

    double precision :: radius, length

    write(*,'(A)',advance='no') 'Domain Radius: '
    read(*,*) radius
    write(*,'(A)',advance='no') 'Domain Length: '
    read(*,*) length

    N = 250000 * 1e9 * (PI*radius**2*length)

    write(*,'(A,I0)') 'Platelet count in domain = ', N


end program main
