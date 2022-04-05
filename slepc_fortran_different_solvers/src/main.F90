!Author: Praveen Kalarickel Ramakrishnan
!This program solves the generalized eigenvalue problem (GEP) for the matrix pencil (A, B)
!The programs calls the subroutine "CompMatrixM" to read the input data files:
!"real_parts.txt", "imag_parts.txt", "weight_matrix.txt", "zero_matrix.txt"
!The output contains the eigen vectors and eigenvalues 
!along with some additional information about the solver
!Many parameters can be set from the command line as mention in the SLEPc manual
!(https://slepc.upv.es/)

!Meanings of relevant variable names
!A, B are the input matrices of the GEP
!n is the size of the problem
!nev is the number of requested eigenvalues
!ncv is the size of the maximum size of the solution subspace
!maxit is maximum number of iterations for convergence
!tol is the value of tolerance level for convergence
!nconv is the number of converged eigenvalues
!kr and ki are used to retrieve the real and imaginary parts of an eigenvalue respectively
!xr and xi are used to retrieve the real and imaginary parts of an eigenvector respectively

#include <slepc/finclude/slepceps.h>
program main
    use slepceps
    use module_file
    implicit none

    Mat            A, B
    EPS            eps
    EPSType        eps_type
    PetscReal      error,tol,re,im;
    PetscScalar    kr,ki
    Vec            xr,xi
    PetscBool      flg
    PetscInt       n,i,j,Istart,Iend,nev,ncv,mpd,maxit,its,nconv, rank_comm
    PetscErrorCode ierr

    !Read the input matrices
    character*32 :: inputfileA_real,inputfileA_imag, inputfileB_real, inputfileB_imag
    double complex, allocatable, dimension(:, :) :: valuesA, valuesB
    double complex, parameter :: mycomplexvalue = complex(21, 42)

    inputfileA_real = "real_parts.txt" 
    inputfileA_imag = "imag_parts.txt" 
    call CompMatrixM(n, inputfileA_real, inputfileA_imag, valuesA)
    inputfileB_real = "weight_matrix.txt" 
    inputfileB_imag = "zero_matrix.txt" 
    call CompMatrixM(n, inputfileB_real, inputfileB_imag, valuesB)
    print*, n

    read*, nev
    ncv = 2*nev+1
    mpd=ncv

    call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank_comm, ierr)

    call MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, PETSC_NULL_SCALAR, A, ierr)
    call MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, PETSC_NULL_SCALAR, B, ierr)
    call MatGetOwnershipRange(A, Istart, Iend, ierr)

    do i=Istart, Iend-1
        do j=Istart, Iend-1
            call MatSetValue(A, i, j, (valuesA(i+1, j+1)-valuesB(i+1, j+1))/complex(0.0,1.0), INSERT_VALUES, ierr)
            call MatSetValue(B, i, j, valuesB(i+1, j+1), INSERT_VALUES, ierr)
        end do
    end do
    if (rank_comm .eq. 0) then
        deallocate(valuesA)
    end if

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY, ierr)

    call MatCreateVecs(A,PETSC_NULL_VEC,xr, ierr)
    call MatCreateVecs(A,PETSC_NULL_VEC,xi, ierr)

!  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                Create the eigensolver and set various options
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
!  /*
!     Create eigensolver context
!  */
    call EPSCreate(PETSC_COMM_WORLD,eps, ierr)
!
!  /*
!     Set operators. In this case, it is a generalized eigenvalue problem
!  */
  call EPSSetOperators(eps,A,B, ierr)
  call EPSSetProblemType(eps,EPS_GNHEP, ierr) !GNHEP is for generalized non hermician eigenvalue problem
!
!  /*
!     Set solver parameters at runtime
!  */
    call EPSSetFromOptions(eps, ierr)
    call EPSSetDimensions(eps, nev,ncv, mpd, ierr)
!
!  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                      Solve the eigensystem
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
!
    call EPSSolve(eps, ierr)

!  /*
!     Optional: Get some information from the solver and display it
!  */
    call EPSGetIterationNumber(eps,its, ierr)
    if (rank_comm .eq. 0) then
        print*, 'Number of iterations of the method: ', its
    end if

    call EPSGetType(eps,eps_type, ierr)
    if (rank_comm .eq. 0) then
        print*, 'Solution method: ', eps_type
    end if

    call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, ierr)
    if (rank_comm .eq. 0) then
        print*, 'Number of requested eigenvalues: ', nev
    end if

    call EPSGetTolerances(eps,tol,maxit, ierr)
    if (rank_comm .eq. 0) then
        print*, 'Stop condition: tol = ', tol, ', maxit = ', maxit
    end if
!
!  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                    Display solution and clean up
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
!  /*
!     Get number of converged approximate eigenpairs
!  */
    call EPSGetConverged(eps,nconv, ierr)
    if (rank_comm .eq. 0) then
        print*, 'Number of converged eigenpairs: ', nconv
    end if
!
    if (nconv>0) then 
!    /*
!       Display eigenvalues and relative errors
!    */
        if (rank_comm .eq. 0) then
            print*, "   Real(lambda)            Imag(lambda)                error=||Ax-kBx||/||kBx||"
            print*, "   ----------------- -----------------------------------------------------------"
        end if 
!
        do i=0, nconv-1 
!      /*
!        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
!        ki (imaginary part)
!      */
            call EPSGetEigenpair(eps,i,kr,ki,xr,xi, ierr)
!      /*
!         Compute the relative error associated to each eigenpair
!      */
            call EPSComputeError(eps,i,EPS_ERROR_RELATIVE,error, ierr)
!
            re = PetscRealPart(kr);
            im = PetscImaginaryPart(kr);
            if(rank_comm .eq. 0) then
                print*, re, im, error
            end if 
        end do
    end if
!
!/*
!   Free work space
!*/
    call EPSDestroy(eps, ierr)
    call MatDestroy(A, ierr)
    call VecDestroy(xr, ierr)
    call VecDestroy(xi, ierr)

    call SlepcFinalize(ierr)

end program main
