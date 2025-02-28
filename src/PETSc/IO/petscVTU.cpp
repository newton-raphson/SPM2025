//
// Created by maksbh on 9/19/19.
//

#include <PETSc/IO/petscVTU.h>
//#include <talyfem/utils/.h>
namespace PETSc{
void petscVectopvtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                    const Vec & vec,
                    const char *fprefix,
                    const char **varName,
                    const DomainExtents & domain,
                    const bool isElemental,
                    const bool isGhosted ,
                    const unsigned int ndof){
#if(DIM == 4)
  TALYFEMLIB::PrintWarning("4D file writing not supported\n");
  return;
#endif

  const DENDRITE_REAL * array;
  VecGetArrayRead(vec, & array);
  IO::writeVecTopVtu(da,treePart,array,fprefix,varName,domain,isElemental,isGhosted,ndof);
  VecRestoreArrayRead(vec,&array);

}
void petscVectopvtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                    const Vec & vec,
                    const char *foldername,
                    const char *fprefix,
                    const char **varName,
                    const DomainExtents & domain,
                    const bool isElemental,
                    const bool isGhosted ,
                    const unsigned int ndof){
#if(DIM == 4)
  TALYFEMLIB::PrintWarning("4D file writing not supported\n");
  return;
#endif
  if((da->isActive())) {

    const DENDRITE_REAL *array;
    VecGetArrayRead(vec, &array);
    IO::writeVecTopVtu(da, treePart,array, foldername, fprefix, varName, domain, isElemental, isGhosted, ndof);
    VecRestoreArrayRead(vec, &array);
  }
  MPI_Barrier(da->getCommActive());
}
}