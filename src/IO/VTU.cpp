////
//// Created by maksbh on 9/17/19.
////
//
#include <IO/VTU.h>
#include <oct2vtk.h>
#include <point.h>
#include <OctToPhysical.h>
#include <sfcTreeLoop_matvec_io.h>
#include <talyfem/utils/utils.h>

#define VTK_QUAD 9
namespace IO {

void writeVecToVtu(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart,
                   const VECType *in,
                   const char *fprefix,
                   const char **varName,
                   const DomainExtents &domain,
                   const bool isElemental,
                   const bool isGhosted,
                   const unsigned int ndof) {
  if (not(da->isActive())) {
    return;
  }

  VECType *vec_;
  Point<DIM> domainMin, domainMax;
  OctToPhysical octToPhysical(domain);
  octToPhysical.getMinAndMaxPoints(domainMin, domainMax);
  const DENDRITE_UINT nPe = da->getNumNodesPerElement();
  const DENDRITE_UINT sz = da->getTotalNodalSz();
  // Hack:
  DENDRITE_UINT num_cells;
  if (not(isElemental)) {
    num_cells = treePart.size();
    da->nodalVecToGhostedNodal(in, vec_, false, ndof);
    da->readFromGhostBegin(vec_, ndof);
    da->readFromGhostEnd(vec_, ndof);
  } else {
    num_cells = treePart.size();
    vec_ = new double[num_cells * ndof];
    std::memcpy(vec_, in, sizeof(double) * num_cells * ndof);
  }
  const DENDRITE_UINT num_vertices = nPe * num_cells;
  const DENDRITE_UINT eleOrder = da->getElementOrder();

//  num_cells = num_cells*eleOrder*eleOrder*eleOrder;
  char fname[FNAME_LENGTH];
  char str[2048];
  FILE *fp = NULL;

  auto partFront = da->getTreePartFront();
  auto partBack = da->getTreePartBack();
  const auto tnCoords = da->getTNCoords();

  DENDRITE_UINT multiplicativeFactor = 1;
  for (DENDRITE_UINT d = 0; d < DIM; d++) {
    multiplicativeFactor *= eleOrder;
  }
  int retval;

  DENDRITE_UINT rank = da->getRankActive();
  DENDRITE_UINT npes = da->getNpesActive();
  sprintf(fname, "%s_%d_%d.vtu", fprefix, rank, npes);
  fp = fopen(fname, "w+");
  if (fp == NULL) {
    std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
    return;
  }
//
  // VTK header
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
  fprintf(fp, "<UnstructuredGrid >\n");
  fprintf(fp, "<Piece NumberOfPoints=\" %d \" NumberOfCells=\" %d \" >\n", num_vertices,
          num_cells * multiplicativeFactor);


  {                                   /** Points data **/
    float *coords = new float[nPe * num_vertices * 3];
    /** Parview expects point to be written in 3D**/
    fprintf(fp, "<Points>\n");
#ifdef DENDRITE_VTU_ASCII
    fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"ascii\">\n");
#else
    fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"binary\">\n");
#endif
    ot::MatvecBaseCoords<DIM> loop(sz, eleOrder, false, 0, tnCoords, &(*treePart.cbegin()), treePart.size(),*partFront, *partBack);
    int cellID = 0;
    while (!loop.isFinished()) {
      if (loop.isPre() && loop.subtreeInfo().isLeaf()) {

#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < nPe; i++) {
#if (DIM == 2)
          fprintf(fp, "%f %f %f ", loop.subtreeInfo().getNodeCoords()[i * DIM + 0],
                  loop.subtreeInfo().getNodeCoords()[i * DIM + 1], 0.0);
#endif
#if (DIM == 3)
          fprintf(fp,"%f %f %f ",
              loop.subtreeInfo().getNodeCoords()[i * DIM + 0],
              loop.subtreeInfo().getNodeCoords()[i * DIM + 1],
              loop.subtreeInfo().getNodeCoords()[i * DIM + 2]);
#endif
        }
#else
        for (int i = 0; i < nPe; i++) {
          coords[cellID * nPe * 3 + i * 3 + 0] = static_cast<float>(loop.subtreeInfo().getNodeCoords()[i * DIM + 0]);
          coords[cellID * nPe * 3 + i * 3 + 1] = static_cast<float>(loop.subtreeInfo().getNodeCoords()[i * DIM + 1]);
#if (DIM == 2)
          coords[cellID * nPe * 3 + i * 3 + 2] = 0.0;
#endif
#if (DIM == 3)
          coords[cellID*nPe*3 + i*3 + 2] = static_cast<float>(loop.subtreeInfo().getNodeCoords()[i*DIM + 2]);
#endif
          octToPhysical.convertCoordsToPhys(&coords[cellID * nPe * 3 + i * 3]);
        }

#endif
        cellID++;
        loop.next();
      } else {
        loop.step();
      }
    }
#ifndef DENDRITE_VTU_ASCII
    retval = io::vtk::vtk_write_binary(fp, (char *) coords, sizeof(*coords) * 3 * num_vertices);
    if (retval) {
      std::cout << rank << ": [VTU Error]: " << "base64 encode point data failed" << std::endl;
      fclose(fp);
    }
#endif
    fprintf(fp, "\n");
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    delete[] coords;
  }

  {                                  /** Write connectivity **/

    fprintf(fp, "<Cells>\n");
#ifdef DENDRITE_VTU_ASCII
    fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n");

    for (unsigned int ele = 0; ele < num_cells; ele++) {
#if (DIM == 3)
      for (unsigned int ek = 0; ek < (eleOrder); ek++)
        for (unsigned int ej = 0; ej < (eleOrder); ej++)
          for (unsigned int ei = 0; ei < (eleOrder); ei++) {
            fprintf(fp,"%d %d %d %d %d %d %d %d ",
                    (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + ei),
                    (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + (ei + 1)),
                    (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + (ei + 1)),
                    (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + ei),
                    (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + ei),
                    (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + (ei + 1)),
                    (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + (ei + 1)),
                    (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + ei));

          }
#endif
#if (DIM == 2)
      for (unsigned int ej = 0; ej < (eleOrder); ej++)
        for (unsigned int ei = 0; ei < (eleOrder); ei++) {
          fprintf(fp, "%d %d %d %d ",
                  (ele * nPe + ej * (eleOrder + 1) + ei),
                  (ele * nPe + ej * (eleOrder + 1) + (ei + 1)),
                  (ele * nPe + (ej + 1) * (eleOrder + 1) + (ei + 1)),
                  (ele * nPe + (ej + 1) * (eleOrder + 1) + (ei)));
        }
#endif
    }
#else
    fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"binary\">\n");

    uint64_t *connectivityID = nullptr;
#if (DIM == 3)
    connectivityID = new uint64_t[num_cells * multiplicativeFactor*8];
    for (unsigned int ele = 0; ele < num_cells; ele++) {
      for (unsigned int ek = 0; ek < (eleOrder); ek++)
        for (unsigned int ej = 0; ej < (eleOrder); ej++)
          for (unsigned int ei = 0; ei < (eleOrder); ei++) {
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 0] = (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + ei);
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 1] = (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + (ei + 1));
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 2] = (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + (ei + 1));
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 3] = (ele * nPe + ek * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + ei);
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 4] = (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + ei);
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 5] = (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + ej * (eleOrder + 1) + (ei + 1));
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 6] = (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + (ei + 1));
            connectivityID[(ele * multiplicativeFactor + ek * (eleOrder) * (eleOrder) + ej * (eleOrder) + ei) * 8 + 7] = (ele * nPe + (ek + 1) * (eleOrder + 1) * (eleOrder + 1) + (ej + 1) * (eleOrder + 1) + ei);
          }
    }
    retval = io::vtk::vtk_write_binary (fp, (char *) connectivityID, sizeof (*connectivityID) * 8 * num_cells*multiplicativeFactor);
    if (retval) {

      std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
      fclose(fp);

    }

#endif
#if (DIM == 2)
    connectivityID = new uint64_t[num_cells * multiplicativeFactor * 4];
    for (unsigned int ele = 0; ele < num_cells; ele++) {
      for (unsigned int ej = 0; ej < (eleOrder); ej++)
        for (unsigned int ei = 0; ei < (eleOrder); ei++) {
          connectivityID[(ele * multiplicativeFactor + ej * (eleOrder) + ei) * 4 + 0] = (ele * nPe +
                                                                                         ej * (eleOrder + 1) + ei);
          connectivityID[(ele * multiplicativeFactor + ej * (eleOrder) + ei) * 4 + 1] = (ele * nPe +
                                                                                         ej * (eleOrder + 1) +
                                                                                         (ei + 1));
          connectivityID[(ele * multiplicativeFactor + ej * (eleOrder) + ei) * 4 + 2] = (ele * nPe +
                                                                                         (ej + 1) * (eleOrder + 1) +
                                                                                         (ei + 1));
          connectivityID[(ele * multiplicativeFactor + ej * (eleOrder) + ei) * 4 + 3] = (ele * nPe +
                                                                                         (ej + 1) * (eleOrder + 1) +
                                                                                         (ei));
        }
    }
    retval = io::vtk::vtk_write_binary(fp, (char *) connectivityID,
                                       sizeof(*connectivityID) * 4 * num_cells * multiplicativeFactor);
    if (retval) {

      std::cout << rank << ": [VTU Error]: " << "base64 encode connectivity data failed" << std::endl;
      fclose(fp);

    }
#endif
    delete[] connectivityID;
#endif
    fprintf(fp, "</DataArray>\n");

  }

  {                                 /** Write offsets **/

#ifdef DENDRITE_VTU_ASCII
    fprintf(fp, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");

    for (unsigned int il = 1, sk = 1; il <= num_cells * multiplicativeFactor; ++il, ++sk) {
#if(DIM == 3)
      fprintf(fp," %u ",(8*il));
      if (!(sk % 8) && il != num_cells)
        fprintf(fp,"\n");
    }
#endif
#if(DIM == 2)
      fprintf(fp, " %u ", (4 * il));
      if (!(sk % 4) && il != num_cells)
        fprintf(fp, "\n");
    }
#endif
#else
    fprintf(fp, "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"binary\">\n");
    uint64_t *offsetID = new uint64_t[num_cells * multiplicativeFactor];
    for (unsigned int il = 1; il <= num_cells * multiplicativeFactor; il++) {
#if(DIM == 3)
      offsetID[il-1] = 8*il;
#endif
#if(DIM == 2)
      offsetID[il - 1] = 4 * il;
#endif
    }
    retval = io::vtk::vtk_write_binary(fp, (char *) offsetID, sizeof(*offsetID) * (num_cells * multiplicativeFactor));
    if (retval) {

      std::cout << rank << ": [VTU Error]: " << "base64 encode offset data failed" << std::endl;
      fclose(fp);
    }
    delete[] offsetID;
#endif
    fprintf(fp, "</DataArray>\n");

  }

  {                                   /** Write cellTypes **/
#ifdef DENDRITE_VTU_ASCII
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (unsigned int il = 0, sk = 1; il < num_cells * multiplicativeFactor; ++il, ++sk) {
#if (DIM == 3)
      fprintf(fp,"%d ",VTK_HEXAHEDRON);
#endif
#if (DIM == 2)
      fprintf(fp, "%d ", VTK_QUAD);
#endif
      if (!(sk % 20) && il != (num_cells - 1))
        fprintf(fp, "\n");
    }
#else
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n");
    uint8_t *cellType = new uint8_t[num_cells * multiplicativeFactor];
    for (unsigned int il = 0; il < num_cells * multiplicativeFactor; ++il) {
#if (DIM == 3)
      cellType[il] = VTK_HEXAHEDRON;
#endif
#if (DIM == 2)
      cellType[il] = VTK_QUAD;
#endif
    }
    retval = io::vtk::vtk_write_binary(fp, (char *) cellType, sizeof(*cellType) * num_cells * multiplicativeFactor);
    if (retval) {

      std::cout << rank << ": [VTU Error]: " << "base64 encode element type data failed" << std::endl;
      fclose(fp);
    }

    delete[] cellType;
#endif
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
  }
  {
    if (not(isElemental)) {
      /** Write local value node data **/
      ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, eleOrder, tnCoords, vec_,&(*treePart.cbegin()), treePart.size(), *partFront, *partBack);


      fprintf(fp, "<PointData>\n");
#ifdef DENDRITE_VTU_ASCII
      for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
        fprintf(fp, "<DataArray type=\"Float64\" Name=\" %s \" format=\"ascii\">\n", varName[dof]);
        treeloop.reset();
        while (!treeloop.isFinished()) {
          if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
            for (int i = 0; i < nPe; i++) {
              fprintf(fp, "%f ", treeloop.subtreeInfo().readNodeValsIn()[i * ndof + dof]);
            }
            treeloop.next();
          } else {
            treeloop.step();
          }
        }
        fprintf(fp, "</DataArray>\n");
      }
#else
      double *localVal = new double[num_cells * nPe];
      for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
        fprintf(fp, "<DataArray type=\"Float64\" Name=\" %s \" format=\"binary\">\n", varName[dof]);

        int cellCount = 0;
        treeloop.reset();
        while (!treeloop.isFinished()){
          if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
            for (int i = 0; i < nPe; i++) {
              localVal[cellCount * nPe + i] = treeloop.subtreeInfo().readNodeValsIn()[i * ndof + dof];
            }
            cellCount++;
            treeloop.next();
          }
          else {
            treeloop.step();
          }
        }
        retval = io::vtk::vtk_write_binary(fp, (char *) localVal, sizeof(*localVal) * num_vertices);

        if (retval) {

          std::cout << rank << ": [VTU Error]: " << "base64 encode point vars data failed" << std::endl;
          fclose(fp);
        }
        fprintf(fp, "</DataArray>\n");
      }
      delete[] localVal;
#endif

      fprintf(fp, "</PointData>\n");
    }
    fprintf(fp, "<CellData> \n");


    /*** Elemental cell data **/

    if (isElemental) {
#ifdef DENDRITE_VTU_ASCII
      for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
        fprintf(fp, "<DataArray type=\"Float64\" Name=\" %s \" format=\"ascii\">\n", varName[dof]);
        for (int i = 0; i < num_cells; i++) {
          for (int k = 0; k < multiplicativeFactor; k++) {
            fprintf(fp,"%f ",vec_[i]);
          }
        }
      }
#else
      double *localCellData = new double[num_cells * multiplicativeFactor];
      for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
        fprintf(fp, "<DataArray type=\"Float64\" Name=\" %s \" format=\"binary\">\n", varName[dof]);
        for (int i = 0; i < num_cells; i++) {
          for (int k = 0; k < multiplicativeFactor; k++) {
            localCellData[i * multiplicativeFactor + k] = vec_[i*ndof+dof];
          }
        }
        retval = io::vtk::vtk_write_binary(fp, (char *) localCellData,
                                           sizeof(*localCellData) * num_cells * multiplicativeFactor);
        if (retval) {

          std::cout << rank << ": [VTU Error]: " << "base64 encode element rank data failed" << std::endl;
          fclose(fp);
        }
        fprintf(fp, "\n</DataArray>\n");
      }

      delete[] localCellData;
#endif

    }

    {
#ifdef DENDRITE_VTU_ASCII
      fprintf(fp, "<DataArray type=\"UInt32\" Name=\"rank\" format=\"ascii\">\n");
      for (int i = 0; i < num_cells * multiplicativeFactor; i++) {
        fprintf(fp, "%d ", da->getRankActive());
      }
#else
      fprintf(fp, "<DataArray type=\"UInt32\" Name=\"rank\" format=\"binary\">\n");
      uint32_t *rankLocal = new uint32_t[num_cells * multiplicativeFactor];
      for (int i = 0; i < num_cells * multiplicativeFactor; i++) {
        rankLocal[i] = da->getRankActive();
      }
      retval = io::vtk::vtk_write_binary(fp, (char *) rankLocal,
                                         sizeof(*rankLocal) * num_cells * multiplicativeFactor);
      if (retval) {

        std::cout << rank << ": [VTU Error]: " << "base64 encode element rank data failed" << std::endl;
        fclose(fp);
      }
      delete[] rankLocal;
#endif
      fprintf(fp, "</DataArray>\n");
      fprintf(fp, "</CellData>\n");
    }
  }
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  delete[] vec_;


}

void writeVecTopVtu(ot::DA<DIM> *da,  const std::vector<TREENODE> & treePart,
                    const VECType *in,
                    const char *fprefix,
                    const char **varName,
                    const DomainExtents &domain,
                    const bool isElemental,
                    const bool isGhosted,
                    const unsigned int ndof) {
  DENDRITE_UINT rank = da->getRankActive();
  DENDRITE_UINT npes = da->getNpesActive();
  if (!rank) {
    char pfname[FNAME_LENGTH];

    sprintf(pfname, "%s.pvtu", fprefix);
    std::ofstream file(pfname);
    file << "<?xml version=\"1.0\"?> " << std::endl;
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" >" << std::endl;
    file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
    file << "<PPoints>\n";
    file << "<PDataArray type=\"Float32\" NumberOfComponents=\"" << 3 << "\"/>\n";
    file << "</PPoints>\n";
    if (not(isElemental)) {
      file << "<PPointData>\n";
#ifdef DENDRITE_VTU_ASCII
      for (int dof = 0; dof < ndof; dof++) {
        file << "<PDataArray type=\"Float64\" Name=\"" << varName[dof] << "\" format=\"ascii\"/>\n";
      }
#else
      for (int dof = 0; dof < ndof; dof++) {
        file << "<PDataArray type=\"Float64\" Name=\"" << varName[dof] << "\" format=\"binary\"/>\n";
      }
#endif
      file << "</PPointData>\n";
    }
    file << "<PCellData>\n";
#ifdef DENDRITE_VTU_ASCII
    if (isElemental) {
      for (int dof = 0; dof < ndof; dof++) {
        file << "<PDataArray type=\"Float64\" Name=\"" << varName[dof] << "\" format=\"ascii\"/>\n";
      }
    }
    file << "<PDataArray type=\"UInt32\" Name=\"rank\" format=\"ascii\"/>\n";
#else
    if (isElemental) {
      for (int dof = 0; dof < ndof; dof++) {
        file << "<PDataArray type=\"Float64\" Name=\"" << varName[dof] << "\" format=\"binary\"/>\n";
      }
    }
    file << "<PDataArray type=\"UInt32\" Name=\"rank\" format=\"binary\"/>\n";
#endif
    file << "</PCellData>\n";

    for (int i = 0; i < npes; i++) {
      std::string fname(fprefix);
      auto const pos= fname.find_last_of('/');
      std::string leaf=fname.substr(pos+1) + "_"+ std::to_string(i) + "_"+std::to_string(npes) + ".vtu";
      file << "<Piece Source=\"" << leaf << "\" />\n";
    }
    file << "</PUnstructuredGrid>\n";
    file << "</VTKFile>\n";
    file.close();
  }
  writeVecToVtu(da, treePart,in, fprefix, varName, domain, isElemental, isGhosted, ndof);
}


DENDRITE_UINT getNumLocalCells(ot::DA<DIM> *da, const std::vector<TREENODE> & treePart) {
  const DENDRITE_UINT sz = da->getTotalNodalSz();
// Hack:
  auto partFront = da->getTreePartFront();
  auto partBack = da->getTreePartBack();
  const auto tnCoords = da->getTNCoords();
  DENDRITE_UINT num_cells = 0;
  ot::MatvecBaseCoords<DIM> loop(sz, da->getElementOrder(), false, 0, tnCoords, &(*treePart.cbegin()), treePart.size(),*partFront, *partBack);
  while (!loop.isFinished()) {
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      num_cells++;
      loop.next();
    } else {
      loop.step();
    }
  }
  return num_cells;

}

void writeBoundaryElements(DA *octDA,const std::vector<TREENODE> & treePart,const char *foldername,const char *fprefix,const DomainExtents & domain) {
  if(not(octDA->isActive())) {
    return;
  }
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const unsigned int npe = octDA->getNumNodesPerElement();
  int counter = 0;
  std::vector<double> bElement(octDA->getLocalElementSz(),0.0);
  ot::MatvecBaseCoords <DIM> loop(sz,octDA->getElementOrder(), false,0,tnCoords,&(*treePart.cbegin()), treePart.size(),*partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      bElement[counter++] = loop.subtreeInfo().isElementBoundary();
      loop.next();
    }
    else{
      loop.step();
    }
  }
  static const char*varname[]{"bElement"};
  writeVecTopVtu(octDA,treePart,bElement.data(),foldername,fprefix,varname,domain,true);

}

void writeVecTopVtu(ot::DA<DIM> *da,
                    const std::vector<TREENODE> & treePart,
                   const VECType *in,
                   const char *foldername,
                   const char *fprefix,
                   const char **varName,
                   const DomainExtents & domain,
                   const bool isElemental,
                   const bool isGhosted ,
                   const unsigned int ndof){
  if(not(da->isActive())) {
    return;
  }
  if (not(TALYFEMLIB::GetMPIRank())) {
    int ierr = mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr != 0 && errno != EEXIST) {
      TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
      return;
    }
  }
  MPI_Barrier(da->getCommActive());
  char fname[PATH_MAX];
  snprintf(fname, sizeof(fname), "%s/%s", foldername, fprefix);
  writeVecTopVtu(da,treePart,in,fname,varName,domain,isElemental,isGhosted,ndof);
}


}
