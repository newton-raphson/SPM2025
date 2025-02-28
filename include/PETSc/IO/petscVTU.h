//
// Created by maksbh on 9/19/19.
//

#ifndef DENDRITEKT_PETSCVTU_H
#define DENDRITEKT_PETSCVTU_H
#include <IO/VTU.h>
namespace PETSc{
 void petscVectopvtu(ot::DA<DIM> *da,
                     const std::vector<TREENODE> & treePart,
                     const Vec & vec,
                     const char *fprefix,
                     const char **varName,
                     const DomainExtents & domain,
                     const bool isElemental = false,
                     const bool isGhosted = false,
                     const unsigned int ndof = 1);

void petscVectopvtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                    const Vec & vec,
                    const char *foldername,
                    const char *fprefix,
                    const char **varName,
                    const DomainExtents & domain,
                    const bool isElemental = false,
                    const bool isGhosted = false,
                    const unsigned int ndof = 1);


}
#endif //DENDRITEKT_PETSCVTU_H
