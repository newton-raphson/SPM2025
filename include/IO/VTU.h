////
//// Created by maksbh on 9/17/19.
////

/**
 * @author maksbh: Write octant to pvtu file.
 *
 */
#ifndef DENDRITEKT_VTU_H
#define DENDRITEKT_VTU_H

#include <oda.h>
#include <DataTypes.h>

#include <sys/stat.h>
namespace IO {

/**
 * Temporary function to get number of local cells.
 * @param da present DA
 * @return number of local cells.
 */
DENDRITE_UINT getNumLocalCells(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart);

/**
 * @brief write to vtu. Each processor write its own vtu
 * @param da present da;
 * @param in in vector
 * @param fprefix  prefix for file
 * @param isElemental if the vector is elemental
 * @param isGhosted whether the vector is ghosted or not.
 * @param ndof degree of freedom
 */
void writeVecToVtu(ot::DA<DIM> *da,
                   const std::vector<TREENODE> & treePart,
                   const VECType *in,
                   const char *fprefix,
                   const char **varName,
                   const DomainInfo & domainInfo,
                   const bool isElemental = false,
                   const bool isGhosted = false,
                   const unsigned int ndof = 1);

/**
 * @brief write to pvtu. Written by rank 0 only
 * @param da present da;
 * @param in in vector
 * @param fprefix  prefix for file
 * @param isElemental if the vector is elemental
 * @param isGhosted whether the vector is ghosted or not.
 * @param ndof degree of freedom
 */
void writeVecTopVtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                    const VECType *in,
                    const char *fprefix,
                    const char **varName,
                    const DomainExtents & domain,
                    const bool isElemental = false,
                    const bool isGhosted = false,
                    const unsigned int ndof = 1);

void writeBoundaryElements(ot::DA<DIM> *da,const std::vector<TREENODE> & treePart,const char *foldername, const char *fprefix,const DomainExtents & domain);

void writeVecTopVtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                   const VECType *in,
                   const char *foldername,
                   const char *fprefix,
                   const char **varName,
                   const DomainExtents & domainInfo,
                   const bool isElemental = false,
                   const bool isGhosted = false,
                   const unsigned int ndof = 1);
}


#endif //DENDRITEKT_VTU_H
