//
// Created by maksbh on 3/3/22.
//

#include <oda.h>
#include <checkPoint.h>
#include <sfcTreeLoop_matvec_io.h>
#include <intergridTransfer.h>

static constexpr int DIM = 2;

typedef ot::TreeNode<unsigned int, DIM> TREENODE;
using DistTree = ot::DistTree<uint, DIM>;
using DA = ot::DA<DIM>;
using DENDRITE_UINT = unsigned int;

void readDACheckpoint(const std::string &foldername, DistTree &dtree, const int rank) {
    std::vector<TREENODE> reloadTreeNode;
    std::string fname = foldername + "/oct_rank_" + std::to_string(rank);
    io::checkpoint::readOctFromFile(fname.c_str(), reloadTreeNode);

    DistTree dtree_(reloadTreeNode, MPI_COMM_WORLD);
    std::swap(dtree, dtree_);
}

void getRefineFlags(DA * octDA, DistTree & dtree, std::vector<ot::OCT_FLAGS::Refine> & refineFlags){
    int numLocalElements = octDA->getLocalElementSz();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int totalSize;
    std::ifstream fin("Coarse"+std::to_string(rank),std::ios::in);
    fin >> totalSize;
    assert(totalSize == numLocalElements);
    refineFlags.resize(totalSize);
    int x;
    for(int i = 0; i < totalSize; i++){
        fin >> x;
        refineFlags[i] = static_cast<ot::OCT_FLAGS::Refine>(x);
    }
    fin.close();

}


void intergrid_fine_to_coarse(
        const double * const input, const DA *oldDA,  const DistTree & oldDistTree,  // fine by fine
        const DA *surrDA, const DistTree & surrDistTree, // coarse by fine
        double * const output, const DA *newDA,  const DistTree & newDistTree,  // coarse by coarse
        int ndof = 1)
{
    // fine -> intergrid surrogate -> shift coarse

    static std::vector<VECType> fineGhosted;
    static std::vector<VECType> surrGhosted;
    oldDA->template createVector(fineGhosted, false, true, ndof);
    surrDA->template createVector<VECType>(surrGhosted, false, true, ndof);
    oldDA->nodalVecToGhostedNodal(input, fineGhosted.data(), ndof);
    oldDA->readFromGhostBegin(fineGhosted.data(), ndof);
    oldDA->readFromGhostEnd(fineGhosted.data(), ndof);
    std::fill(surrGhosted.begin(), surrGhosted.end(), 0);

    fem::MeshFreeInputContext<VECType, TREENODE>
            inctx = fem::mesh_free_input_context(fineGhosted.data(), oldDA, oldDistTree);

    fem::MeshFreeOutputContext<VECType, TREENODE>
            outctx = fem::mesh_free_output_context(surrGhosted.data(), surrDA, surrDistTree);

    // Hack: Only to writeToGhosts, when not all owned nodes touch owned cells.
    // (When fixed, i.e. owned nodes touch owned cells, can use readFromGhosts).
    static std::vector<char> outDirty;
    const char zero = 0;
    surrDA->template createVector<char>(outDirty, false, true, 1);
    surrDA->setVectorByScalar(outDirty.data(), &zero, false, true, 1);

    const RefElement *refel = newDA->getReferenceElement();
    fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

    surrDA->template writeToGhostsBegin<VECType>(surrGhosted.data(), ndof, &(*outDirty.cbegin()));
    surrDA->template writeToGhostsEnd<VECType>(surrGhosted.data(), ndof, false, &(*outDirty.cbegin()));

    ot::distShiftNodes(*surrDA, surrGhosted.data() + ndof * surrDA->getLocalNodeBegin(),
                       *newDA, output,
                       ndof);
}



void checkIntergridTransfer(
        const double * const array,
        const ot::DA<DIM> * octDA,
        const ot::DistTree<unsigned int, DIM> &distTree,
        const unsigned int ndof)
{
    double *ghostedArray = nullptr;
    octDA->nodalVecToGhostedNodal(array,ghostedArray, false,ndof);
    octDA->readFromGhostBegin(ghostedArray,ndof);
    octDA->readFromGhostEnd(ghostedArray,ndof);
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    const unsigned int nPe = octDA->getNumNodesPerElement();
    ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, octDA->getElementOrder(), tnCoords, ghostedArray, &(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront, *partBack);
    bool testPassed = true;
    while (!treeloop.isFinished())
    {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf())
        {
            const double * nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            const PetscScalar * nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
            for(int i = 0; i < nPe; i++){
                double correctValue = 0;
                for(int dim = 0; dim < DIM; dim++){
                    correctValue += nodeCoordsFlat[DIM*i+dim];
                }
                for(int dof = 0; dof < ndof; dof++) {
                    double interpolatedValue = nodeValsFlat[i * ndof + dof];
                    if (fabs(interpolatedValue - correctValue) > 1E-6) {
                        std::cout << "Value at (" << nodeCoordsFlat[DIM * i + 0] << " ," << nodeCoordsFlat[DIM * i + 1]
                                  << ") = " << interpolatedValue << "\n";
                        testPassed = false;
                    }
                }

            }
            treeloop.next();
        }
        else
            treeloop.step();
    }
    bool gtestPassed;
    MPI_Reduce(&testPassed,&gtestPassed,1,MPI_CXX_BOOL,MPI_LAND,0,MPI_COMM_WORLD);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(!rank) {
        if (gtestPassed) {
            std::cout << GRN << "TEST passed" << NRM << "\n";
        } else {
            std::cout << RED << "TEST failed" << NRM << "\n";
        }
    }
    if (ghostedArray)
        delete [] ghostedArray;
    MPI_Barrier(MPI_COMM_WORLD);
}


void testIntergrid() {


    std::string folderName = "Coarse_Before";
    _InitializeHcurve(DIM);
    m_uiMaxDepth = 30;
    std::array<double,DIM> physMin, physMax;
    physMin.fill(0.0);
    physMax.fill(1.0);
    DistTree::BoxDecider boxDecider(physMin,physMax);
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank, comm_size;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);


    DistTree coarse_dtree, fine_dtree;
    readDACheckpoint(folderName, fine_dtree, comm_rank);
    fine_dtree.filterTree(boxDecider);

    double sfc_tol = 0.001;
    int ndof = 1;
    int eleOrder =1;
    DA *fineDA = new DA(fine_dtree, MPI_COMM_WORLD, eleOrder, 100, sfc_tol);

    double * fineVec;
    fineDA->template createVector<VECType>(fineVec,false,false,ndof);
    std::function<void(const double *, double *)> functionPointer = [&](const double *x, double *var) {
        double sum = 0.0;
        for (int d = 0; d < DIM; ++d)
            sum += x[d];
        var[0] = sum;
    };
    fineDA->setVectorByFunction(fineVec,functionPointer,false,false,ndof);
    std::vector<ot::OCT_FLAGS::Refine> octFlags;
    getRefineFlags(fineDA,fine_dtree,octFlags);
    ot::DistTree<unsigned int, DIM> coarseDistTree;
    ot::DistTree<unsigned int, DIM> surrDistTree;
    {
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> coarseTree;
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> surrTree;
        ot::SFC_Tree<DENDRITE_UINT , DIM>::distRemeshWholeDomain(fine_dtree.getTreePartFiltered(), octFlags, coarseTree, sfc_tol, comm);
        surrTree = ot::SFC_Tree<DENDRITE_UINT , DIM>::getSurrogateGrid(ot::RemeshPartition::SurrogateOutByIn, fine_dtree.getTreePartFiltered(), coarseTree, comm); // Line 156
        coarseDistTree = ot::DistTree<unsigned int, DIM>(coarseTree, comm);
        surrDistTree = ot::DistTree<unsigned int, DIM>(surrTree, comm,ot::DistTree<unsigned int,DIM>::NoCoalesce);
    }
    ot::DA<DIM> *coarseDA = new ot::DA<DIM>(coarseDistTree, comm, eleOrder);
    ot::DA<DIM> *surrDA = new ot::DA<DIM>(surrDistTree, comm, eleOrder);
//    std::cout << "Number of elements in fineDA " << fineDA->getLocalElementSz() << "\n";
//    std::cout << "Number of elements in coarseDA " << coarseDA->getLocalElementSz() << "\n";
//    std::cout << "Number of elements in surrDA " << surrDA->getLocalElementSz() << "\n";
    if(!comm_rank){
        std::cout << "Fine DA size = " <<  fineDA->getGlobalNodeSz() << "\n";
        std::cout << "Coarse DA size = " <<  coarseDA->getGlobalNodeSz() << "\n";
    }
    double * coarseVec;
    coarseDA->template createVector<VECType>(coarseVec,false,false,ndof);

    intergrid_fine_to_coarse(fineVec, fineDA, fine_dtree, surrDA, surrDistTree, coarseVec, coarseDA, coarseDistTree, 1);

    checkIntergridTransfer(coarseVec, coarseDA, coarseDistTree, 1);

    delete coarseDA;
    delete fineDA;
    delete surrDA;
    delete coarseVec;
    delete fineVec;
    PetscFinalize();





}

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, NULL);
    DendroScopeBegin() ;

    testIntergrid();

    DendroScopeEnd();
    PetscFinalize();
}