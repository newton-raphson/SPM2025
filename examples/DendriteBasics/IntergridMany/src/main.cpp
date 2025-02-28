//
// Created by maksbh on 2/28/22.
//

#include <DendriteUtils.h>
#include <PETSc/IO/petscVTU.h>
#include <SubDA/SubDomain.h>
#include <DARefine.h>
#include <OctToPhysical.h>
using namespace PETSc;
int main(int argc, char * argv[]) {
  dendrite_init(argc, argv);

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " " << "level eleOrder\n";
    exit(EXIT_FAILURE);
  }

  /// Domain Extents
  DomainInfo cubeDomain, physDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(4.0);
  physDomain.min.fill(0.0);
  physDomain.max.fill(4.0);
  physDomain.max[0] = 4.0;
  DomainExtents domainExtents(cubeDomain, physDomain);

  const DENDRITE_UINT level = static_cast<DENDRITE_UINT>(std::atoi(argv[1]));
  const DENDRITE_UINT eleOrder = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));

  DistTREE distTree;
  SubDomain subDomain(domainExtents);
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };

  /// Create DA
  DA * octDA = createSubDA(distTree,functionToRetain,level,eleOrder);
  const int amount = 2;
  const bool refine_all = false;
  OctToPhysical octToPhysical(domainExtents);

  std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
    TALYFEMLIB::ZEROPTV coords;
    std::memcpy(coords.data(),x, sizeof(double )*DIM);
    octToPhysical.convertCoordsToPhys(coords.data());
    var[0] = coords[0] + coords[1] + coords[2];
    var[1] = 5*(coords[0] + coords[1] + coords[2]);
  };
  Vec initialVector;
  octDA->petscCreateVector(initialVector,false,false,2);
  octDA->petscSetVectorByFunction(initialVector,initial_condition,false,false,2);

  DARefine daRefine(octDA,distTree.getTreePartFiltered(),domainExtents);
  DA * newDA = daRefine.getRefineSubDA(distTree);
  daRefine.petscIntergridTransfer(newDA,distTree,initialVector,2);
  std::swap(octDA,newDA);
  delete newDA;
  static const char * varname[]{"a1","a2"};
  petscVectopvtu(octDA,distTree.getTreePartFiltered(),initialVector,"initial","initial",varname,domainExtents,false,false,2);
  IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"subDARefined","boundary",domainExtents);
  dendrite_finalize(octDA);
}