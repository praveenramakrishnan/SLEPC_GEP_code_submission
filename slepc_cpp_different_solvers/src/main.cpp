/*!Author: Praveen Kalarickel Ramakrishnan
 This program solves the generalized eigenvalue problem (GEP) for the matrix pencil (A, B)
 The programs calls the subroutine "CompMatrixM" to read the input data files:
 "real_parts.txt", "imag_parts.txt", "weight_matrix.txt", "zero_matrix.txt"
 The output contains the eigen vectors and eigenvalues 
 along with some additional information about the colver
 Many parameters can be set from the commandline as mention in the SLEPc manual
 (https://slepc.upv.es/)
 
 Meanings of relevant variable names
 A, B are the input matrices of the GEP
 n is the size of the problem
 nev is the number of requested eigenvalues
 ncv is the size of the maximum size of the solution subspace
 maxit is maximum number of iterations for convergence
 tol is the value of tolerance level for convergence
 nconv is the number of converged eigenvalues
 kr and ki are used to retrieve the real and imaginary parts of an eigenvalue respectively
 xr and xi are used to retrieve the real and imaginary parts of an eigenvector respectively
 */

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>
#include <iostream>
#include "read_dense_matrix_complex.h"
#include <cstring>
#include <complex>
#include <vector>

int main(int argc,char **argv)
{
  Mat            A, B;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       n,i,j,Istart,Iend,nev,ncv,mpd,maxit,its,nconv;
  PetscErrorCode ierr;

  //Read the input matrix from the file
  int N;
  std::vector<double> valuesA, valuesB;
  std::string inputfileA_real = "real_parts.txt";
  std::string inputfileA_imag = "imag_parts.txt";
  CompMatrixM(N, inputfileA_real, inputfileA_imag, valuesA);
  std::string inputfileB_real = "weight_matrix.txt";
  std::string inputfileB_imag = "zero_matrix.txt";
  CompMatrixM(N, inputfileB_real, inputfileB_imag, valuesB);

  n = N;
  //nev = 30;
  std::cin >> nev;
  ncv = 2*nev+1;
  mpd = ncv;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &A);CHKERRQ(ierr);
  ierr = MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &B);CHKERRQ(ierr);

  std::complex<double> jcomplex(0.0, 1.0);

int count = 0;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
   for (i=Istart;i<Iend;i++) {
    for (j=Istart;j<Iend;j++) {
        double realA = valuesA[count];
        double imagA = valuesA[count+1];
        std::complex<double> complexA = realA + jcomplex*imagA;
        double realB = valuesB[count];
        double imagB = valuesB[count+1];
        std::complex<double> complexB = realB + jcomplex*imagB;
        complexA = (complexA - complexB)/jcomplex;
        ierr = MatSetValue(A,i, j, complexA, INSERT_VALUES);CHKERRQ(ierr);
        ierr = MatSetValue(B,i, j, complexB, INSERT_VALUES);CHKERRQ(ierr);
        count +=2;
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  //ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, nev,ncv, mpd);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %.15e%+.15ei %.15e\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %.15e       %.15e\n",(double)re,(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

/*
   Free work space
*/
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}
