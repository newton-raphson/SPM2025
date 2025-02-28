//
// Created by maksbh on 6/18/21.
//



#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <feMatrix.h>
#include <IO/VTU.h>
#include <chrono>
#include <SDARefine.h>
#include <BenchDataTypes.h>
#include "dollar_stat.h"
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  if (argc < 7) {
    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " level boundaryLevel eleOrder fractionRefine RefinementStrategy CaseType");
    return EXIT_FAILURE;
  }
  MPI_Barrier(MPI_COMM_WORLD);


  int level = std::atoi(argv[1]);
  int blevel = std::atoi(argv[2]);
  int eleOrder = std::atoi(argv[3]);
  DENDRITE_REAL fracRefine = std::atof(argv[4]);
  RefinementStrategy refineStrategy = static_cast<RefinementStrategy>(static_cast<bool>(std::atoi(argv[5])));
  CaseType caseType =  static_cast<CaseType>(std::atoi(argv[6]));



  DomainInfo fullDomain, physDomain;
  fullDomain.max.fill(8.0);
  fullDomain.min.fill(0.0);
  physDomain.max.fill(8.0);
  physDomain.min.fill(0.0);
  if(caseType == CaseType::CHANNEL) {
    physDomain.max[1] = 1.0;
    physDomain.max[2] = 1.0;
  }
  DomainExtents domain(fullDomain, physDomain);

  SubDomain subDomain(domain);
  if(caseType == CaseType::SPHERE){
    DENDRITE_REAL  center[DIM]{3.0,4.0,4.0};
    VOXEL::Sphere sphere(center,0.5);
    subDomain.addObject(sphere);
  }
  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *physCoords,
                                                                               double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };
  DistTREE dTree;
  DA *octDA = createSubDA(dTree, functionToRetain, level, eleOrder);
  subDomain.finalize(octDA,dTree.getTreePartFiltered(),domain);
  double remesh_time = 0;
  auto start  = std::chrono::high_resolution_clock::now();
  if(caseType == CaseType::CHANNEL) {
    for (int i = 0; i < 4; i++) {

      Refine refine(octDA, dTree.getTreePartFiltered(), domain,caseType, &subDomain,blevel, level, fracRefine);
      auto start_remesh = std::chrono::high_resolution_clock::now();
      DA *newDA = refine.getRefineSubDA(dTree, 0.1, refineStrategy);
      auto end_remesh = std::chrono::high_resolution_clock::now();
      remesh_time += (static_cast<std::chrono::duration<DENDRITE_REAL>>(end_remesh - start_remesh)).count();;

      if (newDA == nullptr) {
        break;
      }
      std::swap(newDA, octDA);
      delete newDA;
      TALYFEMLIB::PrintStatus(i, "iteration of refinement done. Global node Size = ", octDA->getGlobalNodeSz());

      subDomain.finalize(octDA, dTree.getTreePartFiltered(), domain);
    }
  }
  else{
    int counter= 0;
    while(true) {
      Refine refine(octDA, dTree.getTreePartFiltered(), domain, caseType, &subDomain, blevel - 1, level, fracRefine);
      auto start_remesh = std::chrono::high_resolution_clock::now();
      DA *newDA = refine.getRefineSubDA(dTree, 0.1, refineStrategy);
      auto end_remesh = std::chrono::high_resolution_clock::now();
      remesh_time += (static_cast<std::chrono::duration<DENDRITE_REAL>>(end_remesh - start_remesh)).count();;

      if (newDA == nullptr) {
        break;
      }
      std::swap(newDA, octDA);
      delete newDA;
      TALYFEMLIB::PrintStatus(++counter, "iteration of refinement done. Global node Size = ", octDA->getGlobalNodeSz());
    }
    subDomain.finalize(octDA, dTree.getTreePartFiltered(), domain);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  auto end  = std::chrono::high_resolution_clock::now();
  auto totalTime =  (static_cast<std::chrono::duration<DENDRITE_REAL>>(end - start)).count(); ;
  TALYFEMLIB::PrintStatus("Total remesh time  = ", remesh_time);
  TALYFEMLIB::PrintStatus("TotalDA Creation Time = ", totalTime);
  remesh_time = 0;
  {DOLLAR("SphereRefine")

  if(caseType == CaseType::SPHERE) {
    int counter = 0;
    while (true) {
      Refine refine(octDA, dTree.getTreePartFiltered(), domain, caseType, &subDomain, blevel , level, fracRefine);
      auto start_remesh = std::chrono::high_resolution_clock::now();
      DA *newDA = refine.getRefineSubDA(dTree, 0.1, refineStrategy);
      auto end_remesh = std::chrono::high_resolution_clock::now();
      remesh_time += (static_cast<std::chrono::duration<DENDRITE_REAL>>(end_remesh - start_remesh)).count();;

      if (newDA == nullptr) {
        break;
      }
      std::swap(newDA, octDA);
      delete newDA;
      TALYFEMLIB::PrintStatus(++counter, "iteration of refinement done. Global node Size = ", octDA->getGlobalNodeSz());
    }
  }
  TALYFEMLIB::PrintStatus("Final remesh time  = ", remesh_time);
  }
  //IO::writeBoundaryElements(octDA,dTree.getTreePartFiltered(),"bnd","bnd",domain);

  dollar::DollarStat dollar_stat(MPI_COMM_WORLD);
  dollar::clear();
  // Collect mean, min, and max timings over all processes.
  dollar::DollarStat reduce_mean = dollar_stat.mpi_reduce_mean();
  dollar::DollarStat reduce_min = dollar_stat.mpi_reduce_min();
  dollar::DollarStat reduce_max = dollar_stat.mpi_reduce_max();
  if (TALYFEMLIB::GetMPIRank() == 0)
  {
    std::ofstream file("mean_chrome.json");
    reduce_mean.chrome(file);

    std::cout << "\n" << "[Mean]\n";
    reduce_mean.tsv(std::cout);
    std::cout << "\n" << "[Min]\n";
    reduce_min.tsv(std::cout);
    std::cout << "\n" << "[Max]\n";
    reduce_max.tsv(std::cout);

    /// reduce_mean.print(std::cout);
    /// reduce_min.print(std::cout);
    /// reduce_max.print(std::cout);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  dendrite_finalize(octDA);

}