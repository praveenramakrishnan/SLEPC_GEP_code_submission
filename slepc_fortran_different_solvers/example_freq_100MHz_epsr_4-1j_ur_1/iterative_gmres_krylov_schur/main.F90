!Author: Praveen Kalarickel Ramakrishnan
!Date: 2022-03-21

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
