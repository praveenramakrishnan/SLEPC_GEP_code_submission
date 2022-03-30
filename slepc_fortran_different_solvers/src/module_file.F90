!This is module file contating the subroutine "CompMatrixM" to read in the
!complex matrix from data files contating the real and imaginary parts
module module_file
implicit none

contains

subroutine CompMatrixM(n, inputfilename_real, inputfilename_imag, M)
    double complex, allocatable, dimension(:, :) :: M
    double precision, allocatable, dimension(:, :) :: real_M, imag_M

    integer :: n, iloop, jloop
    character*32 :: inputfilename_real, inputfilename_imag

    call number_of_lines_of_file(n, inputfilename_real)
    allocate(M(n, n))
    allocate(real_M(n, n))
    allocate(imag_M(n, n))

    open(unit=1, file=inputfilename_real)
    open(unit=2, file=inputfilename_imag)
    do iloop=1, n
        read(unit=1, fmt=*) real_M(:, iloop) 
        read(unit=2, fmt=*) imag_M(:, iloop) 
    end do

    do iloop=1, n
        do jloop=1, n
            M(iloop, jloop) = complex(real_M(jloop, iloop), imag_M(jloop, iloop))
        end do
    end do
    close(unit=1)
    close(unit=2)

    deallocate(real_M)
    deallocate(imag_M)

end subroutine    

subroutine number_of_lines_of_file(N_lines, inputfilename)
    integer :: N_lines
    character*32 :: inputfilename

    N_lines = 0
    open(unit=1, file=inputfilename)
    do
        read(unit=1, fmt=*, end=10)
        N_lines = N_lines + 1
    end do

    10 close(unit=1)

end subroutine

end module
