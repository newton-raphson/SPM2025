/*
 * testMovingBall.cpp
 *   Output to vtu files the tree generated by the points of a moving sphere.
 *
 * Masado Ishii  --  UofU SoC, 2019-01-24
 */


#include "treeNode.h"
#include "tsort.h"
#include "octUtils.h"

#include "hcurvedata.h"

#include "octUtils.h"
#include "oct2vtk.h"
#include <vector>

#include <assert.h>
#include <mpi.h>
#include <stdio.h>


// ...........................................................................
template <typename T>
std::vector<ot::TreeNode<T,4>> generateMovingSphere(int numPoints);
// ...........................................................................

//-------------------------------
// test_distOutputTreeBalancing()
//-------------------------------
void test_distOutputTreeBalancing(int numPoints, MPI_Comm comm = MPI_COMM_WORLD)
{
  int nProc, rProc;
  MPI_Comm_size(comm, &nProc);
  MPI_Comm_rank(comm, &rProc);

  using T = unsigned int;
  const unsigned int dim = 4;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = generateMovingSphere<T>(numPoints);
  std::vector<TreeNode> treePart;

  const unsigned int maxPtsPerRegion = 16;

  const double loadFlexibility = 0.2;

  const T leafLevel = m_uiMaxDepth;
  const T firstVariableLevel = 1;      // Not sure about this whole root thing...

  printf("[%d] Starting distTreeBalancing()...\n", rProc);
  ot::SFC_Tree<T,dim>::distTreeBalancing(points, treePart, maxPtsPerRegion, loadFlexibility, comm);
  printf("[%d] Finished distTreeBalancing().\n", rProc);
  printf("[%d] treePart.size == %ld\n", rProc, (long) treePart.size());

  // Since dim==4 and we can only output octrees of dim 3 (or less?) to vtu,
  // use the slicing operator.
  /// unsigned int treeDepth = 0;
  /// for (const TreeNode &tn : treePart)
  /// {
  ///   treeDepth = (tn.getLevel() > treeDepth ? tn.getLevel() : treeDepth);
  /// }
  /// unsigned int timeStep = 1u << (m_uiMaxDepth - treeDepth);
  const char dimNames[] = "XYZT";
  std::cout << "Output files will be placed in the directory '_output'.\n";
  for (unsigned int d = 0; d < dim; d++)
  {
    // Many slices, like a time series.
    unsigned int t = 0;
    while (t < (1u << m_uiMaxDepth))
    {
      constexpr bool RM_DUPS_AND_ANC = false;
      constexpr bool RM_DUPS_ONLY = true;

      std::vector<ot::TreeNode<T,3>> slice3D;
      projectSliceKTree(&(*treePart.begin()), slice3D, (unsigned int) treePart.size(), d, (T) t);
      ot::SFC_Tree<T,3>::distRemoveDuplicates(slice3D, loadFlexibility, RM_DUPS_AND_ANC, comm);
      // Note that the partitioning of the slice is not related to partition of original tree.

      printf("[%d] slice %c size == %ld\n", rProc, dimNames[d], (long) slice3D.size());

      // Output to file with oct2vtu().
      char fPrefix[] =  "                                                             ";  // beware buffer overflow.
      sprintf(fPrefix,  "_output/testSlice-%c-t%u", dimNames[d], t);
      io::vtk::oct2vtu(&(*slice3D.begin()), (unsigned int) slice3D.size(), fPrefix, comm);

      // Advance slice dimension by coarseness of the slice.
      unsigned int sliceDepth = 0;
      for (const ot::TreeNode<T,3> &tn : slice3D)
        sliceDepth = (tn.getLevel() > sliceDepth ? tn.getLevel() : sliceDepth);

      unsigned int sliceThickness = 1u << (m_uiMaxDepth - sliceDepth);
      t += sliceThickness;
    }
  }

  _DestroyHcurve();
}


//
// main()
//
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int ptsPerProc = 200;
  if (argc > 1)
    ptsPerProc = strtol(argv[1], NULL, 0);

  test_distOutputTreeBalancing(ptsPerProc, MPI_COMM_WORLD);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}




//
// generateMovingSphere()
//
template <typename T>
std::vector<ot::TreeNode<T,4>> generateMovingSphere(int numPoints)
{
        std::vector<ot::TreeNode<T,4>> points;
        std::array<T,4> uiCoords;
        std::array<double,4> fCoords;

        const double pi = acos(-1.0);

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<double> dist_theta(0, 2.0*pi);
        std::uniform_real_distribution<double> dist_zu(-1.0, 1.0);
        std::uniform_real_distribution<double> dist_time(0.0, 1.0);

        const double clampL = 0;
        const double clampH = 1.0;//(1u << m_uiMaxDepth);

        for (int ii = 0; ii < numPoints; ii++)
        {
          double theta = dist_theta(gen);
          double zu = dist_zu(gen);
          double time = dist_time(gen);

          double u = fabs(zu);
          double r = sqrt(u);
          double h = sqrt(1-u);

          fCoords[0] = r*cos(theta);
          fCoords[1] = r*sin(theta);
          fCoords[2] = copysign(h, zu);
          
          fCoords[3] = time;

          fCoords[0] = 0.25*fCoords[0] + 0.25;// + 0.5*time;
          fCoords[1] = 0.25*fCoords[1] + 0.25 + 0.5*time;
          fCoords[2] = 0.25*fCoords[2] + 0.25;// + 0.5*time;

          #pragma unroll(4)
          for (unsigned int d = 0; d < 4; d++)
          {
            fCoords[d] = (fCoords[d] < clampL ? clampL : fCoords[d] > clampH ? clampH : fCoords[d]);
            uiCoords[d] = (T) ((1u << m_uiMaxDepth) * fCoords[d]);
          }

          ot::TreeNode<T,4> tn(uiCoords, m_uiMaxDepth);
          points.push_back(tn);

          /// printf("(%f, %f, %f, %f)  ", fCoords[0], fCoords[1], fCoords[2], fCoords[3]);
          /// std::cout << tn.getBase32Hex().data() << "\n";
        }

        return points;
}
