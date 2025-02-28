//
// Created by maksbh on 5/19/20.
//

#include <PETSc/Solver/NonLinearSolver.h>

namespace PETSc{
NonlinearSolver::NonlinearSolver(DA* da,feMat<DIM> * mat, feVec<DIM> * vec,const DistTREE & dtree, const Point<DIM> & domainMin,const Point<DIM> & domainMax, const DENDRITE_UINT dof, const bool mfree, const bool skipBCApplication,const bool oneShotMatrixAssembly)
:Solver(da,mat,vec,dtree,domainMin,domainMax,dof,mfree,skipBCApplication,oneShotMatrixAssembly){

}
PetscErrorCode NonlinearSolver::init() {
#ifdef PROFILING
  Timer timer("NonLinearInit");
#endif
  assert(m_vecSolution == NULL);
  assert(m_vecRHS == NULL);
  assert(m_uiDof != 0);

  int ierr;

  // Allocate memory for working vectors
  ierr = m_octDA->petscCreateVector(m_vecSolution, false, false, m_uiDof); CHKERRQ(ierr);
  ierr = m_octDA->petscCreateVector(m_vecRHS, false, false, m_uiDof); CHKERRQ(ierr);




  if(m_matrixFree) {
    //Matrix free
    PetscInt matsize;
    ierr = VecGetLocalSize(m_vecRHS, &matsize);
    CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD, matsize, matsize, PETSC_DETERMINE, PETSC_DETERMINE, this,
                          &m_matJacobian);
    CHKERRQ(ierr);
    ierr = MatShellSetOperation(m_matJacobian, MATOP_MULT, (void (*)()) (ShellMatMult));
    CHKERRQ(ierr);
  }
  else{
    ierr =  m_octDA->createMatrixWithRows(m_matJacobian, MATBAIJ,m_distTree, m_uiDof); CHKERRQ(ierr);
    CHKERRQ(ierr);
    ierr = MatSetBlockSize(m_matJacobian,(PetscInt) m_uiDof);
    CHKERRQ(ierr);
  }

  // avoids new nonzero errors when boundary conditions are used
  ierr = MatSetOption(m_matJacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

  // Create a KSP context to solve  @ every timestep
  ierr = SNESCreate(m_octDA->getCommActive(), &m_snes); CHKERRQ(ierr);
  ierr = SNESSetFunction(m_snes, m_vecRHS, FormFunction, this); CHKERRQ(ierr);
  ierr = SNESSetJacobian(m_snes, m_matJacobian, m_matJacobian, FormJacobian, this); CHKERRQ(ierr);

  KSP ksp;
  SNESGetKSP(m_snes,&ksp);


  ierr = SNESSetFromOptions(m_snes); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode NonlinearSolver::solve() {
  if (m_octDA->isActive()) {
    int ierr;
    if(not m_skipBCApplication) {
      updateBoundaries(m_octDA);

      ierr = boundary_conditions_.applyVecBC(m_octDA, m_vecSolution);

      CHKERRQ(ierr);


      ierr = VecAssemblyBegin(m_vecSolution);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(m_vecSolution);
      CHKERRQ(ierr);
    }
    ierr = SNESSolve(m_snes, PETSC_NULL, m_vecSolution);

    CHKERRQ(ierr);
    SNESConvergedReason converged_reason;
    SNESGetConvergedReason(m_snes, &converged_reason);

    if (converged_reason < 0) {
      // diverged
      TALYFEMLIB::PrintWarning("Non-linear solve diverged.");

    }
  }
  return 0;

}

PetscErrorCode NonlinearSolver::jacobianMatMult(Vec in, Vec out) {
  assert(m_matrixFree);
  assert(m_Mat != NULL);
  int ierr;
  ierr = VecZeroEntries(out); CHKERRQ(ierr);

  m_Mat->setPlaceholder(m_guess);
  m_Mat->matVec(in, out);
  ierr = getBC().applyMatrixFreeBC(m_octDA, in, out); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode NonlinearSolver::FormJacobian(SNES snes, Vec sol, Mat jac, Mat precond_matrix, void *ctx) {
#ifdef PROFILING
  PetscLogEventBegin(Profiling::matAssembly,0,0,0,0);
#endif
  int ierr;
  NonlinearSolver *nl = (NonlinearSolver *) ctx;

  // note: jacobianMatMult depends on this staying set
  nl->m_guess = sol;

  ierr = MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);
  if(!nl->m_matrixFree and not(nl->isAssembled)){
    ierr = MatZeroEntries(jac); CHKERRQ(ierr);
    nl->m_Mat->setPlaceholder(nl->m_guess);
#ifdef PROFILING
    PetscLogEventBegin(Profiling::matElementalAssembly,0,0,0,0);
    {
      Timer timer("GetAssmebledMatrix");
      DOLLAR("GetAssembledMatrix")
#endif

      nl->m_Mat->getAssembledMatrix(&jac, 0);

#ifdef PROFILING
    }

#endif

#ifdef PROFILING
    {
      Timer timer("AssemblyBegin-End");
      DOLLAR("Mat-Assembly-Begin-End")
#endif
#ifdef PROFILING
    {
      PetscInt nstash[4],globalStash[4];
      MatStashGetInfo(jac,&nstash[0],&nstash[1],&nstash[2],&nstash[3]);
#ifdef PETSC_USE_64BIT_INDICES
      MPI_Reduce(nstash,globalStash,4,MPI_INT64_T,MPI_MAX,0,MPI_COMM_WORLD);
#else
      MPI_Reduce(nstash,globalStash,4,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
#endif
      if(!TALYFEMLIB::GetMPIRank()) {
        std::cout << "Stash Max : " <<  globalStash[0] << " " << globalStash[1] << " " << globalStash[2] << " " << globalStash[3] << "\n";
      }
    }
#endif
      ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
      ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
#ifdef PROFILING
    }
#endif

#ifdef PROFILING
    PetscLogEventEnd(Profiling::matElementalAssembly,0,0,0,0);
#endif

#ifdef PROFILING
    PetscLogEventBegin(Profiling::matBC,0,0,0,0);
#endif
    if(not nl->m_skipBCApplication) {
      ierr = nl->getBC().applyMatBC(nl->m_octDA, jac);
      CHKERRQ(ierr);
    }
#ifdef PROFILING
    PetscLogEventEnd(Profiling::matBC,0,0,0,0);
#endif
    if(nl->m_oneShotMatrixAssembly){
      nl->isAssembled = true;
    }

  }

  nl->m_Mat->setPlaceholder(NULL);
#ifdef PROFILING
  PetscLogEventEnd(Profiling::matAssembly,0,0,0,0);
#endif
  return 0;
}
PetscErrorCode NonlinearSolver::FormFunction(SNES snes, Vec in, Vec out, void *ctx) {
#ifdef PROFILING
  PetscLogEventBegin(Profiling::vecAssembly,0,0,0,0);
#endif
  int ierr;
  NonlinearSolver *nl = (NonlinearSolver *) ctx;

  nl->m_guess = in;

  nl->m_Vec->setPlaceholder(nl->m_guess);
  VecZeroEntries(out);
#ifdef PROFILING
  PetscLogEventBegin(Profiling::vecElementalAssembly,0,0,0,0);
  {
    Timer timer("FormFunction");
#endif
    nl->m_Vec->computeVec(in, out);
#ifdef PROFILING
  }

  PetscLogEventEnd(Profiling::vecElementalAssembly,0,0,0,0);
#endif
#ifdef PROFILING
  PetscLogEventBegin(Profiling::vecBC,0,0,0,0);
#endif
  if(not nl->m_skipBCApplication) {
    ierr = nl->getBC().applyResidualBC(nl->m_octDA, out);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(out);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(out);
    CHKERRQ(ierr);
  }
  nl->m_Vec->setPlaceholder( NULL);
#ifdef PROFILING
  PetscLogEventEnd(Profiling::vecBC,0,0,0,0);
#endif

#ifdef PROFILING
  PetscLogEventEnd(Profiling::vecAssembly,0,0,0,0);
#endif
  return 0;
}

NonlinearSolver::~NonlinearSolver() {
  cleanup();
}

PetscErrorCode NonlinearSolver::cleanup() {
#ifdef PROFILING
  Timer timer("cleanup");
#endif
  int ierr;
  ierr = SNESDestroy(&m_snes); CHKERRQ(ierr);
  ierr = MatDestroy(&m_matJacobian); CHKERRQ(ierr);
  ierr = VecDestroy(&m_vecSolution); CHKERRQ(ierr);
  ierr = VecDestroy(&m_vecRHS); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroyAll(); CHKERRQ(ierr);

  m_guess = NULL;  // managed by petsc
  return 0;
}






}
