//
// Created by maksbh on 12/20/21.
//

#ifndef DENDRITEKT_CARTESIANCLOUD_H
#define DENDRITEKT_CARTESIANCLOUD_H

#include <DataTypes.h>
namespace CartesianCloud {

  template<int ndof>
  struct valIDPair {
    DENDRITE_REAL value[ndof];
    PetscInt id;

    bool operator==(valIDPair const &other) const {
      return (this->id == other.id);
    }

    bool operator<(valIDPair const &other) const {
      return (this->id < other.id);
    }

    bool operator>(valIDPair const &other) const {
      return not(*this < other);
    }
  };

  template<int ndof>
  class CartesianData {

    PetscInt getGlobalID(const DENDRITE_REAL *position) const;

    PetscInt getGlobalID(const int *id) const;

    void getNeighborIDs(const int *pointID, int neighborsID[][DIM]) const;

    void getValueAtID(const int * id,valIDPair<ndof>  & val) const;

    bool isLocalPointsFinalized = false;
    /// Cartesian Spacing
    DENDRITE_REAL DX_[DIM];
    /// Number of Grids in each direction
    const std::array<int, DIM> &numGrids_;
    /// Domain Information
    const DomainExtents &m_domain;
    std::vector<valIDPair<ndof>> localValIDs;

  public:
    CartesianData(const std::array<int, DIM> &nGrids, const DomainExtents &domain);

    ~CartesianData();

    void addPointInProc(const DENDRITE_REAL *position);

    void finalizeLocalPoints();

    void addValues(const int * id,const DENDRITE_REAL *values);

    void getLocalID(const DENDRITE_REAL * position, int * id)const;

    void getValue(const DENDRITE_REAL * position,DENDRITE_REAL *value ) const;

    void finalize();



  };

  template<int ndof>
  CartesianData<ndof>::CartesianData(const std::array<int, DIM> &nGrids, const DomainExtents &domain)
    :numGrids_(nGrids), m_domain(domain) {
    for (int dim = 0; dim < DIM; dim++) {
      DX_[dim] =
        (domain.physicalDADomain.max[dim] - domain.physicalDADomain.min[dim]) / ((DENDRITE_REAL) (numGrids_[dim] - 1));
    }
  }

  template<int ndof>
  PetscInt CartesianData<ndof>::getGlobalID(const DENDRITE_REAL *position) const {
    DENDRITE_UINT id[DIM];
#pragma unroll (DIM)
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      id[dim] = (position[dim] - m_domain.physicalDADomain.min[dim]) / DX_[dim];
    }

#if(DIM == 3)
    PetscInt globalID = id[2] * numGrids_[1] * numGrids_[0] + id[1] * numGrids_[0] + id[0];
#elif(DIM == 2)
    PetscInt globalID =  id[1]*numGrids_[0] + id[0];
#else
  throw TALYFEMLIB::TALYException() << " Not implemented for dim = " << DIM << "\n"
#endif
    return globalID;
  }

  template<int ndof>
  PetscInt CartesianData<ndof>::getGlobalID(const int *id) const {
#if(DIM == 3)
    PetscInt globalID = id[2] * numGrids_[1] * numGrids_[0] + id[1] * numGrids_[0] + id[0];
#elif(DIM == 2)
    PetscInt globalID =  id[1]*numGrids_[0] + id[0];
#else
    throw TALYFEMLIB::TALYException() << " Not implemented for dim = " << DIM << "\n"
#endif
    return globalID;
  }

  template<int ndof>
  void CartesianData<ndof>::addPointInProc(const DENDRITE_REAL *position) {
    int id[DIM];
    this->getLocalID(position,id);
    assert(not(isLocalPointsFinalized));
    static constexpr DENDRITE_UINT numNeighbors = 1u << DIM;
    int idNeighbors[numNeighbors][DIM];
    this->getNeighborIDs(id, idNeighbors);
    for (int i = 0; i < numNeighbors; i++) {
      PetscInt lID = this->getGlobalID(idNeighbors[i]);
      valIDPair<ndof> val;
      val.id = lID;
      localValIDs.push_back(val);
    }

  }

  template<int ndof>
  void CartesianData<ndof>::getNeighborIDs(const int *pointID, int neighborsID[][DIM]) const {
    static constexpr DENDRITE_UINT numNeighbors = 1u << DIM;
#if (DIM == 2)
    static constexpr DENDRITE_UINT  scale[numNeighbors][DIM] {{0,0},{1,0},{0,1},{1,1}};
#endif

#if (DIM == 3)
    static constexpr DENDRITE_UINT scale[numNeighbors][DIM]{{0, 0, 0},
                                                            {1, 0, 0},
                                                            {0, 1, 0},
                                                            {1, 1, 0},
                                                            {0, 0, 1},
                                                            {1, 0, 1},
                                                            {0, 1, 1},
                                                            {1, 1, 1}};
#endif
    for (DENDRITE_UINT i = 0; i < numNeighbors; i++) {
      for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        DENDRITE_UINT _id = pointID[dim] + scale[i][dim];
        neighborsID[i][dim] = _id < numGrids_[dim] ? _id : numGrids_[dim] - 1;
      }
    }
  }

  template<int ndof>
  void CartesianData<ndof>::finalizeLocalPoints() {
    std::sort(localValIDs.begin(), localValIDs.end());
    localValIDs.erase(std::unique(localValIDs.begin(), localValIDs.end()), localValIDs.end());
    isLocalPointsFinalized = true;
//    std::cout << "[INFO] Rank: " << TALYFEMLIB::GetMPIRank() << ", size: " << localValIDs.size() << "\n";
  }
  template<int ndof>
  void CartesianData<ndof>::addValues(const int * id,const DENDRITE_REAL *values) {
    assert(isLocalPointsFinalized);
    valIDPair<ndof> temp;
    std::memcpy(temp.value,values, sizeof(DENDRITE_REAL)*ndof);
    PetscInt globalID = this->getGlobalID(id);
    temp.id = globalID;
    bool found = std::binary_search(localValIDs.begin(),localValIDs.end(),temp);
    if(not found){
      return;
    }
    const auto lowerbound = std::lower_bound(localValIDs.begin(),localValIDs.end(),temp);
    auto dist = std::distance(localValIDs.begin(),lowerbound);
    std::memcpy(&localValIDs[dist],&temp, sizeof(localValIDs));
  }

  template<int ndof>
  void CartesianData<ndof>::getLocalID(const DENDRITE_REAL * position, int * id)const{
#pragma unroll (DIM)
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      id[dim] = (position[dim] - m_domain.physicalDADomain.min[dim]) / DX_[dim];
    }
  }

  template<int ndof>
  void CartesianData<ndof>::getValue(const DENDRITE_REAL *position, DENDRITE_REAL *value) const{
    int localID[DIM];
    this->getLocalID(position, localID);
    DENDRITE_REAL startPosGrid[DIM];

    for(DENDRITE_UINT dim = 0; dim < DIM ; dim++){
      startPosGrid[dim] = m_domain.physicalDADomain.min[dim] + localID[dim]*DX_[dim];
    }

    static constexpr DENDRITE_UINT numNeighbors = 1u <<DIM;
    int idNeighbors[numNeighbors][DIM];

    this->getNeighborIDs(localID, idNeighbors);
    DENDRITE_REAL xd = (position[0] - startPosGrid[0])/DX_[0];
    DENDRITE_REAL yd = (position[1] - startPosGrid[1])/DX_[1];
#if DIM==3
    DENDRITE_REAL zd = (position[2] - startPosGrid[2])/DX_[2];
#endif
    std::array<valIDPair<ndof>,numNeighbors> neighborValues;
    for(int i = 0; i < numNeighbors; i++){
      this->getValueAtID(idNeighbors[i],neighborValues[i]);
    }
    for(DENDRITE_UINT dof = 0; dof < ndof; dof++){
      DENDRITE_REAL _data[numNeighbors];
      for(int i = 0; i < numNeighbors; i++){
        _data[i] = neighborValues[i].value[dof];
      }
#if (DIM == 3)
      DENDRITE_REAL c00 = (1 - xd) * (_data[0]) + xd*_data[1];
      DENDRITE_REAL c01 = (1 - xd) * (_data[2]) + xd*_data[3];
      DENDRITE_REAL c10 = (1 - xd) * (_data[4]) + xd*_data[5];
      DENDRITE_REAL c11 = (1 - xd) * (_data[6]) + xd*_data[7];

      DENDRITE_REAL c0  = (1 - yd) * (c00) + yd*c01;
      DENDRITE_REAL c1  = (1 - yd) * (c10) + yd*c11;

      value[dof] = (1 - zd)*c0 + zd*c1;
#elif (DIM == 2)
      DENDRITE_REAL valx1 = (1 - xd) * (_data[0]) + xd*_data[1];
      DENDRITE_REAL valx2 = (1 - xd) * (_data[2]) + xd*_data[3];

      value[dof] = (1 - yd)*valx1 + valx2*yd;
#endif
    }

  }


  template<int ndof>
  void CartesianData<ndof>::getValueAtID(const int * id, valIDPair<ndof>  & val) const{
    PetscInt globalID = this->getGlobalID(id);
    val.id = globalID;
    bool found = std::binary_search(localValIDs.begin(),localValIDs.end(),val);
    assert(found);
    const auto lowerbound = std::lower_bound(localValIDs.begin(),localValIDs.end(),val);
    auto dist = std::distance(localValIDs.begin(),lowerbound);
    std::memcpy(&val,&localValIDs[dist], sizeof(localValIDs));
  }

  template<int ndof>
  void CartesianData<ndof>::finalize() {
    localValIDs.clear();
    std::vector<valIDPair<ndof>> temp;
    std::swap(localValIDs,temp);
  }

  template<int ndof>
  CartesianData<ndof>::~CartesianData(){
    this->finalize();
  }
}

#endif //DENDRITEKT_CARTESIANCLOUD_H
