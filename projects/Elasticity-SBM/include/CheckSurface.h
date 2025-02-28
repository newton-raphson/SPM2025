//
// Created by maksbh on 12/5/22.
//

#ifndef DENDRITEKT_CHECKSURFACE_H
#define DENDRITEKT_CHECKSURFACE_H
#include "Traversal/Traversal.h"
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"

class CheckSurface:public Traversal{

    static constexpr  int numPoints = 1u << DIM;
    static constexpr  int numFaces = 2*DIM;
    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    std::vector<std::bitset<numFaces>> faceMarker_;
    SubDomainBoundary * subDomainBoundary_;
    
    std::vector<TALYFEMLIB::ZEROPTV> PosLocal;
    std::vector<TALYFEMLIB::ZEROPTV> PosGlobal;

public:
    CheckSurface(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
                 std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
                 const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary);
    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    void correctCycles();

    void WriteFalseNodeToFile();

};

CheckSurface::CheckSurface(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
                           std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary)
        :eleMarkers_(eleMarkers), Traversal(octDA,treePart,v,domain),subdomain_(subDomain),subDomainBoundary_(subDomainBoundary){

    std::bitset<numFaces> mark;
    mark.reset();

    faceMarker_.resize(octDA->getLocalElementSz());
    this->traverse();
}

void CheckSurface::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
#if (DIM == 3)
    static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0,1,2,3,4,5,6,7},{0,2,6,8,18,20,24,26}};
  static constexpr DENDRITE_UINT numNodesPerFace = 4;
  static constexpr DENDRITE_UINT faceID[2*DIM][4]
  {
      {0,2,4,6}, // Left
      {1,3,5,7}, // Right
      {0,1,4,5}, // Bottom
      {2,3,6,7}, // Top
      {0,1,2,3}, // Back
      {4,5,6,7},  // Front
  };
#elif(DIM == 2)
    static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0,1,2,3},{0,2,6,8}};
    static constexpr DENDRITE_UINT numNodesPerFace = 2;
    static constexpr DENDRITE_UINT faceID[2*DIM][2]
            {
                    {2,0}, // Left
                    {1,3}, // Right
                    {0,1}, // Bottom
                    {3,2}, // Top
            };

#else
    throw std::runtime_error("Not implemented\n");
#endif
    const int eleOrder = m_octDA->getElementOrder();
    // We need to eliminate only intercepted elements.
    // Only Neighbors of intercepted cannot have cycles unless you are doing something stupid with very poor resolution.
    if(eleMarkers_[localElemID].test(ElementMarker::INTERCEPTED_ELEMENT) or eleMarkers_[localElemID].test(ElementMarker::OUT_ELEMENT))
    {
        std::vector<TALYFEMLIB::ZEROPTV> coords(m_octDA->getNumNodesPerElement());
        std::vector<bool> isBoundary(m_octDA->getNumNodesPerElement(),false);
        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            std::memcpy(coords[i].data(), &this->m_coords[i*DIM],sizeof(double)*DIM);
        }
        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            unsigned int id;
            subDomainBoundary_->generateBoundaryFlags(coords[i],id);
        }

        for(int i = 0; i < numFaces; i++) {
            bool isBoundaryFaceinterceptedElement = true;
            for(int j = 0; j < numNodesPerFace; j++) {
                unsigned int cornerNodesOnFace = faceID[i][j];
                unsigned int nodeID = cornerMap[eleOrder - 1][cornerNodesOnFace];
                isBoundaryFaceinterceptedElement = isBoundaryFaceinterceptedElement and isBoundary[nodeID];
            }

            bool SBMCheck = true;
            for (int j = 0; j < numNodesPerFace; j++) {
                SBMCheck = SBMCheck and ( (subdomain_->functionToRetainPhysNodes(&this->m_coords[cornerMap[eleOrder - 1][faceID[i][j]]*DIM]) == ibm::IN)// This node belong to inactive region
                                          or ((values[cornerMap[eleOrder - 1][faceID[i][j]]] == NodeMarker::SBM_FALSE_INTERCEPTED_NODES)));
            }

            bool isBoundaryFace = isBoundaryFaceinterceptedElement or SBMCheck;
            if(isBoundaryFace) {
                faceMarker_[localElemID].set(i,true);
            }
        }

    }
    localElemID++;


}

void CheckSurface::correctCycles(){

    int localElemSize = this->m_octDA->getLocalElementSz();
    for(int i = 0; i < localElemSize; i++){

        for(int nFace = 0; nFace < DIM; nFace++){
            if((faceMarker_[i].test(2*nFace)) and (faceMarker_[i].test(2*nFace + 1))){ // parallel GP insertions
                eleMarkers_[i].reset();
                eleMarkers_[i].set(ElementMarker::SBM_FALSE_INTERCEPTED);
                //CountRemove++;
            }

        }

    }
}

void CheckSurface::WriteFalseNodeToFile() {
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = PosLocal.size();
    std::vector<int> eachProcData(nProc);
    MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> disp(nProc, 0);
    for (int i = 1; i < disp.size(); i++)
    {
        disp[i] = disp[i - 1] + eachProcData[i - 1];
    }

    int totalProcData = 0;
    for (int i = 0; i < nProc; i++)
    {
        totalProcData += eachProcData[i];
    }

    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        PosGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(PosLocal.data(), PosLocal.size(), ZEROPTVtype, PosGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    std::string fPrefix = "trueGP.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "x0,y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "x0,y0,z0\n");
#endif

        for (int i = 0; i < PosGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e\n",
                    PosGlobal[i](0), PosGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e\n",
                    PosGlobal[i](0), PosGlobal[i](1), PosGlobal[i](2));
#endif
        }

        fclose(fp);
    }

}




#endif //DENDRITEKT_CHECKSURFACE_H
