//
// Created by chenghau on 10/10/22.
//

#ifndef LE_KT_SBMCALC_H
#define LE_KT_SBMCALC_H

#include "util.h"
#include <IMGA/IMGA.h>

class SBMCalc
{
private:

    double d[DIM];
  double DirichletBCValue;
  TALYFEMLIB::FEMElm fe_;
  const IMGA *imga_;
  LEInputData *idata_;
  ZEROPTV shift_;
  ZEROPTV normal;
  const double min_domain = 0.0;
  const double max_domain = 1.0 - min_domain;
  const my_kd_tree_t *kd_tree_;

  /**
   * @brief function to check whether the shited point on true boundary is within the corresponding triangle
   * @param pt the carved-out based gauss point position
   * @param d the distance function between true boundary and surrogate boundary
   * @param m_triangle the corresponding triangle to check whether the point is inside
   * @param shift the initial displacement of the stl geometry
   * @return [bool] whether the point is on the plane of the triangle
   */
  bool CheckInside3DTriangle(const ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift);

  /**
   * @brief calculate the shortest distance from octree-based Gauss Points to triangle edges
   * @param pt the carved-out based gauss point position
   * @param m_triangle the corresponding triangle to find the distance
   * @param shift the initial displacement of the stl geometry
   * @param d [out] the distance function from octree-based Gauss Points to triangle edges
   */
  void ShortestDist2TriEdge(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM]);

  void ShortestDist2TriVertex(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM]);

public:

    enum class BCTypes
    {
        DIRICHLET=0,
        NEUMANN=1
    };
  /**
   * @brief constructor
   * @param fe the element we use to access gauss points
   * @param idata input Data
   * @param imga imga context, and we use this to access geometry
   */
  SBMCalc(const TALYFEMLIB::FEMElm &fe, LEInputData *idata, const IMGA *imga);

  SBMCalc(const FEMElm &fe, LEInputData *idata, const IMGA *imga, const my_kd_tree_t *kd_tree);


    /**
   * @brief calculate the distance function for different kinds of geometries
   * @param d distance function
   * @param geom_ID here we also return geometry ID because we do not have 4side_carved
   */
  void Dist2Geo(double (&d)[DIM], int &geom_ID);

  /**
   * @brief calculate the true normal of the SBM geometry
   * @param normal true normal of SBM geometry
   * @param d distance function
   */
  void NormalofGeo(ZEROPTV &normal, const double (&d)[DIM]);

  /**
   * returns the boundary condition for different cases
   * @param d distance function
   * @param DirichletBCValue [out] the calculated BC value
   */
  void GetDirichletBC(const double (&d)[DIM], double *DirichletBCValue, bool &DirichletHaveSet);

  /**
   * returns either dirichlet or neumann boundary condition for different cases and the corresponding value
   * @param d distance function
   * @param BCValue [out] the calculated BC value
   * @param BCType [out] whether dirichlet or neumann BC (the enum class BCTypes)
   */
    void GetBC(const double (&d)[DIM], double *BCValue, BCTypes &BCType);

};

SBMCalc::SBMCalc(const TALYFEMLIB::FEMElm &fe, LEInputData *idata, const IMGA *imga)
    : fe_(fe), imga_(imga), shift_(idata->ibm_geom_def[0].InitialDisplacement)
{
  idata_ = idata;
}
SBMCalc::SBMCalc(const TALYFEMLIB::FEMElm &fe, LEInputData *idata, const IMGA *imga, const my_kd_tree_t *kd_tree)
        : fe_(fe), imga_(imga), shift_(idata->ibm_geom_def[0].InitialDisplacement), kd_tree_(kd_tree)
{
    idata_ = idata;
}



bool SBMCalc::CheckInside3DTriangle(const ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift)
{

  const ZEROPTV ptMove{pt(0) + d[0] - shift(0), pt(1) + d[1] - shift(1), pt(2) + d[2] - shift(2)};

  const ZEROPTV ptA{m_triangle.triangleCoord[0][0], m_triangle.triangleCoord[0][1], m_triangle.triangleCoord[0][2]};
  const ZEROPTV ptB{m_triangle.triangleCoord[1][0], m_triangle.triangleCoord[1][1], m_triangle.triangleCoord[1][2]};
  const ZEROPTV ptC{m_triangle.triangleCoord[2][0], m_triangle.triangleCoord[2][1], m_triangle.triangleCoord[2][2]};

  ZEROPTV u;
  ZEROPTV v;
  ZEROPTV w;

  u.crossProduct(ptB - ptA, ptMove - ptA);
  v.crossProduct(ptC - ptB, ptMove - ptB);
  w.crossProduct(ptA - ptC, ptMove - ptC);

  if (u.innerProduct(v) < 0.0f)
  {
    return false;
  }
  else if (u.innerProduct(w) < 0.0f)
  {
    return false;
  }
  else
  {
    return true;
  }
}

void SBMCalc::ShortestDist2TriEdge(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM])
{
  std::vector<ZEROPTV> ShortestVector2Line(3);
  std::vector<ZEROPTV> ptri(3);
  for (int trinum = 0; trinum < DIM; trinum++)
  {
    ptri[trinum] = {m_triangle.triangleCoord[trinum][0] + shift[0], m_triangle.triangleCoord[trinum][1] + shift[1], m_triangle.triangleCoord[trinum][2] + shift[2]};
  }

  ShortestVector2Line[0] = (ptri[1] - ptri[0]) * (((pt - ptri[0]).innerProduct(ptri[1] - ptri[0])) / (ptri[1] - ptri[0]).norm()) - (pt - ptri[0]);
  ShortestVector2Line[1] = (ptri[2] - ptri[1]) * (((pt - ptri[1]).innerProduct(ptri[2] - ptri[1])) / (ptri[2] - ptri[1]).norm()) - (pt - ptri[1]);
  ShortestVector2Line[2] = (ptri[0] - ptri[2]) * (((pt - ptri[2]).innerProduct(ptri[0] - ptri[2])) / (ptri[0] - ptri[2]).norm()) - (pt - ptri[2]);

  int pickNumber = 0;
  double mindist = 100;
  for (int trinum = 0; trinum < DIM; trinum++)
  {
    if (ShortestVector2Line[trinum].norm() < mindist)
    {
      pickNumber = trinum;
      mindist = ShortestVector2Line[trinum].norm();
    }
  }

  for (int trinum = 0; trinum < DIM; trinum++)
  {
    d[trinum] = ShortestVector2Line[pickNumber](trinum);
  }
}

void SBMCalc::ShortestDist2TriVertex(const ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const ZEROPTV &shift, double (&d)[DIM])
{
    std::vector<ZEROPTV> ShortestVector2Vertex(3);
    std::vector<ZEROPTV> ptri(3);
    for (int trinum = 0; trinum < DIM; trinum++)
    {
        ptri[trinum] = {m_triangle.triangleCoord[trinum][0] + shift[0], m_triangle.triangleCoord[trinum][1] + shift[1], m_triangle.triangleCoord[trinum][2] + shift[2]};
    }

    ShortestVector2Vertex[0] = ptri[0] - pt;
    ShortestVector2Vertex[1] = ptri[1] - pt;
    ShortestVector2Vertex[2] = ptri[2] - pt;

    int pickNumber = 0;
    double mindist = 100;
    for (int trinum = 0; trinum < DIM; trinum++)
    {
        if (ShortestVector2Vertex[trinum].norm() < mindist)
        {
            pickNumber = trinum;
            mindist = ShortestVector2Vertex[trinum].norm();
        }
    }

    for (int trinum = 0; trinum < DIM; trinum++)
    {
        d[trinum] = ShortestVector2Vertex[pickNumber](trinum);
    }
}

//void SBMCalc::Dist2Geo(double (&d)[DIM], int &geom_ID)
//{
//  const ZEROPTV pt = fe_.position();
//  double x = pt.x();
//  double y = pt.y();
//
//  /// initialize it as zero. It is because in most SBM test case, we only have one geometry.
//  geom_ID = 0;
//
//#if (DIM == 2)
//
//  switch (idata_->SbmGeo)
//  {
//      case LEInputData::SBMGeo::CIRCLE:
//      {
//
//          double rad = 0.5;
//          double radius_gp = sqrt((x - shift_.x()) * (x - shift_.x()) + (y - shift_.y()) * (y - shift_.y()));
//          d[0] = (rad * (x - shift_.x()) / radius_gp + shift_.x()) - x;
//          d[1] = (rad * (y - shift_.y()) / radius_gp + shift_.y()) - y;
//          break;
//      }
//
//  case LEInputData::SBMGeo::RING:
//  {
//
//    auto &geo_tmp = idata_->ibm_geom_def.at(0);
//    double radius_tmp = sqrt(pow(x - geo_tmp.InitialDisplacement.x(), 2) + pow(y - geo_tmp.InitialDisplacement.y(), 2));
//
//    DENDRITE_REAL R2;
//    double rad;
//
//    if (fabs(radius_tmp-0.25)> fabs(radius_tmp-1.0))
//    {
//      geom_ID = 0;
//      R2 = 1 * 1;
//      rad = 1;
//    }
//    else
//    {
//      geom_ID = 1;
//      R2 = 0.25 * 0.25;
//      rad = 0.25;
//    }
//
//    auto &geo = idata_->ibm_geom_def.at(geom_ID);
//
//    /// calculate d vector ======================================================================
//    double x_mid = geo.InitialDisplacement.x();
//    double y_mid = geo.InitialDisplacement.y();
//
//    double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
//
//    double radius_gp = sqrt((x - x_mid) * (x - x_mid) + (y - y_mid) * (y - y_mid));
//    d[0] = (rad * (x - x_mid) / radius_gp + x_mid - x);
//    d[1] = (rad * (y - x_mid) / radius_gp + y_mid - y);
//
//    break;
//  }
//
//  case LEInputData::SBMGeo::ROTATE:
//  {
//    double x_true = cos(M_PI / 4) * (x - 0.5 - shift_[0]) - sin(M_PI / 4) * (y - 0.5 - shift_[1]) + 0.5;
//    double y_true = sin(M_PI / 4) * (x - 0.5 - shift_[0]) + cos(M_PI / 4) * (y - 0.5 - shift_[1]) + 0.5;
//
//    double d0, d1;
//
//    if ((x_true - min_domain) * (x_true - max_domain) < 0)
//    {
//      d0 = 0;
//      d1 = (fabs(y_true - min_domain) > fabs(y_true - max_domain)) ? max_domain - y_true : min_domain - y_true;
//    }
//    else if ((y_true - min_domain) * (y_true - max_domain) < 0)
//    {
//      d1 = 0;
//      d0 = (fabs(x_true - min_domain) > fabs(x_true - max_domain)) ? max_domain - x_true : min_domain - x_true;
//    }
//    else
//    {
//      d0 = (fabs(x_true - min_domain) > fabs(x_true - max_domain)) ? max_domain - x_true : min_domain - x_true;
//      d1 = (fabs(y_true - min_domain) > fabs(y_true - max_domain)) ? max_domain - y_true : min_domain - y_true;
//    }
//
//    d[0] = cos(M_PI / 4) * d0 + sin(M_PI / 4) * d1;
//    d[1] = -sin(M_PI / 4) * d0 + cos(M_PI / 4) * d1;
//    break;
//  }
//
//  case LEInputData::SBMGeo::BUNNY:
//  case LEInputData::SBMGeo::STAR:
//  {
//    auto m_lines = imga_->getGeometries()[0]->getMSH()->getLines();
//    int msh_size = m_lines.size();
//
//    double MinDist = 1000.0;
//    for (DENDRITE_UINT i = 0; i < msh_size; i++)
//    {
//      if (sqrt((x - m_lines[i].lineCoord[0][0] - shift_[0]) * (x - m_lines[i].lineCoord[0][0] - shift_[0]) + (y - m_lines[i].lineCoord[0][1] - shift_[1]) * (y - m_lines[i].lineCoord[0][1] - shift_[1])) < MinDist)
//      {
//        MinDist = sqrt((x - m_lines[i].lineCoord[0][0] - shift_[0]) * (x - m_lines[i].lineCoord[0][0] - shift_[0]) + (y - m_lines[i].lineCoord[0][1] - shift_[1]) * (y - m_lines[i].lineCoord[0][1] - shift_[1]));
//        d[0] = m_lines[i].lineCoord[0][0] + shift_[0] - x;
//        d[1] = m_lines[i].lineCoord[0][1] + shift_[1] - y;
//      }
//    }
//    break;
//  }
//
//  /// returning zero distance when it is not specific SBM geometries
//  default:
//  {
//    d[0] = 0.0;
//    d[1] = 0.0;
//    break;
//  }
//  }
//
//#endif
//
//#if (DIM == 3)
//
//  double z = pt.z();
//
//  switch (idata_->SbmGeo)
//  {
//
//  case LEInputData::SBMGeo::cantilever:
//  {
//      d[0] = 0.0;
//      d[1]  = 0.0;
//      d[2] = 0.0;
//
//      for (int dim= 0; dim <DIM;dim++) {
//          if (fabs(fe_.surface()->normal()(dim)) > 0.0){
//              double diff1 = idata_->minmax_cantilever[dim].first - fe_.position()(dim);
//              double diff2 = idata_->minmax_cantilever[dim].second - fe_.position()(dim);
//
//              d[dim] = std::min({diff1, diff2}, [](double a, double b) {
//                  return std::fabs(a) < std::fabs(b);
//              });
//          }
//      }
//
//#ifndef NDEBUG
//      if (fabs(fe_.surface()->normal()(0)) > 0.0) {
//          assert(d[0] < 0);
//      } else if (fabs(fe_.surface()->normal()(1)) > 0.0) {
//          assert(d[1] > 0);
//      }
//#endif
//
//    break;
//  }
//
//  case LEInputData::SBMGeo::SPHERE:
//  {
//    double rad = shift_[0];
//    double radius_gp = sqrt((x - rad) * (x - rad) + (y - rad) * (y - rad) + (z - rad) * (z - rad));
//    d[0] = (rad * (x - rad) / radius_gp + rad) - x;
//    d[1] = (rad * (y - rad) / radius_gp + rad) - y;
//    d[2] = (rad * (z - rad) / radius_gp + rad) - z;
//    break;
//  }
//
//
//  case LEInputData::SBMGeo::ARBITRARY:
//  case LEInputData::SBMGeo::BUNNY3D:
//  case LEInputData::SBMGeo::MOAI:
//  case LEInputData::SBMGeo::ARM:
//  {
//
//    double MinDist = 100.0;
//    ZEROPTV OnePointVector;
//    ZEROPTV PickNormalVector;
//    int PickGeomID = 0;
//    int PickTrianleID = 0;
//
//    switch (idata_->DistCalcType)
//    {
//
//    case LEInputData::typeDistCalc::GP_BASED:
//    {
//        PrintError("Not support GP_BASED distance function calculation!");
//        exit(EXIT_FAILURE);
//        break;
//    }
//
//    case LEInputData::typeDistCalc::NORMAL_BASED:
//    {
//      for (int geoID = 0; geoID < imga_->getGeometries().size(); geoID++)
//      {
//        std::vector<GEOMETRY::Triangles> m_triangles = imga_->getGeometries()[geoID]->getSTL()[0].getTriangles();
//
//        for (int i = 0; i < m_triangles.size(); i++)
//        {
//          if (sqrt(pow(x - (m_triangles[i].triangleCoord[0][0] + m_triangles[i].triangleCoord[1][0] + m_triangles[i].triangleCoord[2][0]) / 3 - shift_[0], 2) + pow(y - (m_triangles[i].triangleCoord[0][1] + m_triangles[i].triangleCoord[1][1] + m_triangles[i].triangleCoord[2][1]) / 3 - shift_[1], 2) + pow(z - (m_triangles[i].triangleCoord[0][2] + m_triangles[i].triangleCoord[1][2] + m_triangles[i].triangleCoord[2][2]) / 3 - shift_[2], 2)) < MinDist)
//          {
//            MinDist = sqrt(pow(x - (m_triangles[i].triangleCoord[0][0] + m_triangles[i].triangleCoord[1][0] + m_triangles[i].triangleCoord[2][0]) / 3 - shift_[0], 2) + pow(y - (m_triangles[i].triangleCoord[0][1] + m_triangles[i].triangleCoord[1][1] + m_triangles[i].triangleCoord[2][1]) / 3 - shift_[1], 2) + pow(z - (m_triangles[i].triangleCoord[0][2] + m_triangles[i].triangleCoord[1][2] + m_triangles[i].triangleCoord[2][2]) / 3 - shift_[2], 2));
//
//            for (int dim = 0; dim < DIM; dim++)
//            {
//              OnePointVector(dim) = m_triangles[i].triangleCoord[0][dim] + shift_[dim] - pt(dim);
//              PickNormalVector(dim) = m_triangles[i].normal[dim];
//              PickTrianleID = i;
//              PickGeomID = geoID;
//            }
//          }
//        }
//      }
//
//      // scaling of vector
//      double scale = 0.0;
//      for (int dim = 0; dim < DIM; dim++)
//      {
//        scale += OnePointVector(dim) * PickNormalVector(dim);
//      }
//
//      for (int dim = 0; dim < DIM; dim++)
//      {
//        d[dim] = scale * PickNormalVector(dim);
//      }
//
//      GEOMETRY::Triangles m_triangle = imga_->getGeometries()[PickGeomID]->getSTL()[0].getTriangles()[PickTrianleID];
//      if (!CheckInside3DTriangle(pt, d, m_triangle, shift_))
//      {
//        ShortestDist2TriEdge(pt, m_triangle, shift_, d);
//        if (!CheckInside3DTriangle(pt, d, m_triangle, shift_)){
//            ShortestDist2TriVertex(pt, m_triangle, shift_, d);
//        }
//      }
//      break;
//    }
//
//    case LEInputData::typeDistCalc::KD_TREE:
//        {
//            size_t num_results = 1;
//            std::vector<uint32_t> ret_index(num_results);
//            std::vector<double> out_dist_sqr(num_results);
//
//
//
//            const double query_pt[3] = {x - shift_[0], y - shift_[1], z - shift_[2]}; // shift_
//
//            num_results = kd_tree_->knnSearch(
//                    &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
//
//            const std::vector<GEOMETRY::Triangles> *m_triangles = &imga_->getGeometries()[0/*TODO: fix hard code*/]->getSTL()[0].getTriangles();
//            for (int dim = 0; dim < DIM; dim++)
//            {
//                OnePointVector(dim) = m_triangles->at(ret_index[0]).triangleCoord[0][dim] + shift_[dim] - pt(dim);
//                PickNormalVector(dim) = m_triangles->at(ret_index[0]).normal[dim];
//            }
//            PickTrianleID = ret_index[0];
//            PickGeomID = 0 /*TODO: fix hard code*/;
//
//
//            // scaling of vector
//            double scale = 0.0;
//            for (int dim = 0; dim < DIM; dim++)
//            {
//                scale += OnePointVector(dim) * PickNormalVector(dim);
//            }
//
//            for (int dim = 0; dim < DIM; dim++)
//            {
//                d[dim] = scale * PickNormalVector(dim);
//            }
//
//            GEOMETRY::Triangles m_triangle = imga_->getGeometries()[PickGeomID]->getSTL()[0].getTriangles()[PickTrianleID];
//            if (!CheckInside3DTriangle(pt, d, m_triangle, shift_))
//            {
//                ShortestDist2TriEdge(pt, m_triangle, shift_, d);
//                if (!CheckInside3DTriangle(pt, d, m_triangle, shift_)){
//                    ShortestDist2TriVertex(pt, m_triangle, shift_, d);
//                }
//            }
//
//
//#ifndef NDEBUG
//            std::ofstream fout("SurrogateGPandDistanceFunc.txt",  std::ios::out | std::ios::app);
//            fout << pt.x() << " " << pt.y() << " " << pt.z() << " " << d[0] << " " << d[1] << " " << d[2]
//                 << pt.x() + d[0] << " " << pt.y() + d[1] << " " << pt.z() + d[2] << "\n";
//            fout.close();
//#endif
//
//            break;
//        }
//    }
//
//    break;
//  }
//
//  /// returning zero distance when it is not specific SBM geometries
//  default:
//  {
//    d[0] = 0.0;
//    d[1] = 0.0;
//    d[2] = 0.0;
//    break;
//  }
//  }
//
//#endif
//}

/**
 * this is for sbm without lambda -> the gauss points on surrogate boundary will always treat as Inactive
 * @param normal
 * @param d
 */
void SBMCalc::NormalofGeo(ZEROPTV &normal, const double (&d)[DIM])
{
    assert(idata_->RatioGPSBM == 1);

  double R2 = 0;

  for (int dim = 0; dim < DIM; dim++)
  {
    R2 += pow(d[dim], 2);
  }
  for (int dim = 0; dim < DIM; dim++)
  {
    normal(dim) = -d[dim] / sqrt(R2);
  }

}



void SBMCalc::GetBC(const double (&d_)[DIM], double *BCValue, BCTypes &BCType) {
    const ZEROPTV pt = fe_.position();

    double x_true = pt.x() + d_[0];
    double y_true = pt.y() + d_[1];
#if (DIM == 3)
    double z_true = pt.z() + d_[2];
#endif
#if (DIM == 2)
    //// throw error if the dimension is 2
    throw std::runtime_error("Not implemented for 2D!");
#endif
    switch (idata_->SbmGeo) {

        case LEInputData::SBMGeo::BUNNY:
        {
            /// we have three cases if z_true is below some value z1 we fix it
            /// if z_true is above some value z2 we apply the force
            /// else we apply the no traction BC
            /// we have three cases if z_true is below some value z1 we fix it
            /// if z_true is above some value z2 we apply the force
            /// else we apply the no traction BC
            BCValue[0] = sin(M_PI * x_true) * cos(M_PI * y_true) / 50.0;  // Smooth oscillatory X-flow
            BCValue[1] = cos(M_PI * x_true) * sin(M_PI * z_true) / 50.0;  // Cross-flow in Y
            BCValue[2] = exp(-y_true * y_true) * sin(M_PI * z_true) / 100.0;  // Decaying vertical motion
            BCType = BCTypes::DIRICHLET;
            break;
        }
        case LEInputData::SBMGeo::PLANE: {
            /// if the y_value is less than -0.11 we fix the displacement
            if (y_true < -0.11) {
                BCValue[0] = 0.0;
                BCValue[1] = 0.0;
                BCValue[2] = 0.0;
                BCType = BCTypes::DIRICHLET;
                break;
            }
            //// above it we apply the neumann BC
            else
            {
                BCValue[0] = 1000*cos(M_PI * x_true) * sin(M_PI * y_true);
                BCValue[1] = 1000*sin(M_PI * x_true) * cos(M_PI * y_true);
                BCValue[2] = -1000;
                BCType = BCTypes::NEUMANN;
                break;
            }
        }
        case LEInputData::SBMGeo::PLANT:
        {
            BCValue[0] = sin(M_PI*x_true)*cos(M_PI*y_true)/10.0;
            BCValue[1] = cos(M_PI*x_true)*sin(M_PI*y_true)/10.0;
            BCType = BCTypes::DIRICHLET;
            break;
        }
        case LEInputData::SBMGeo::GYROID:
        {
            BCValue[0] = sin(M_PI*x_true)*cos(M_PI*y_true)/10.0;
            BCValue[1] = cos(M_PI*x_true)*sin(M_PI*y_true)/10.0;
            BCType = BCTypes::DIRICHLET;
            break;
        }
        case LEInputData::SBMGeo::EIFFEL:
        {
            // Simulated displacement pattern for structural deformation
            BCValue[0] = 0.1 * sin(M_PI * x_true) * cos(M_PI * y_true);  // Horizontal displacement
            BCValue[1] = 0.05 * sin(M_PI * y_true) * sin(M_PI * z_true); // Vertical wave displacement
            BCValue[2] = 0.0;  // No displacement in Z (fixed in the depth direction)

            BCType = BCTypes::DIRICHLET;
            break;
        }
        case LEInputData::SBMGeo::SPHERE:
        {
            BCValue[0] = sin(M_PI * (x_true )) * cos(M_PI * (y_true)) * sin(M_PI * (z_true)) / 10.0;
            BCValue[1] = cos(M_PI * (x_true )) * sin(M_PI * (y_true)) * sin(M_PI * (z_true)) / 10.0;
            BCValue[2] = sin(M_PI * (x_true )) * sin(M_PI * (y_true)) * cos(M_PI * (z_true)) / 20.0;

            BCType = BCTypes::DIRICHLET;

            break;

        }
        default:
        {
            break;
        }




    }
}

#endif // LE_KT_SBMCALC_H
