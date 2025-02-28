//
// Created by chenghau on 8/18/22.
//

#ifndef DENDRITEKT_SBMCalc_H
#define DENDRITEKT_SBMCalc_H
#include <IMGA/IMGA.h>
#include "nanoflann.hpp"

// construct a kd-tree index:
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
        PointCloud<double>, 3 /* dim */>;

using my_kd_tree_ptr = std::unique_ptr<my_kd_tree_t>;


class SBMCalc
{
private:
    double d[DIM];
    double DirichletBCValue;
    TALYFEMLIB::FEMElm fe_ = TALYFEMLIB::FEMElm(nullptr);
    const IMGA *imga_;
    std::vector<my_kd_tree_t*> kd_trees_;

    NSInputData *idata_;
    TALYFEMLIB::ZEROPTV shift_;

    TALYFEMLIB::ZEROPTV normal;
    const double min_domain = 0.0;
    const double max_domain = 1.0 - min_domain;
    const double min_domain_RotBox = -0.5;
    const double max_domain_RotBox = 0.5;
    const int CarvedOutGeomID_;

    void CalDist(double A, double B, double C, double x1, double y1, double &d);
    void CalVec(double A, double B, double C, double x1, double y1, double &d1, double &d2);
    DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe);

    /**
     * @brief function to check whether the shited point on true boundary is within the corresponding triangle
     * @param pt the carved-out based gauss point position
     * @param d the distance function between true boundary and surrogate boundary
     * @param m_triangle the corresponding triangle to check whether the point is inside
     * @param shift the initial displacement of the stl geometry
     * @return [bool] whether the point is on the plane of the triangle
     */
    bool CheckInside3DTriangle(const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift);

    /**
     * @brief calculate the shortest distance from octree-based Gauss Points to triangle edges
     * @param pt the carved-out based gauss point position
     * @param m_triangle the corresponding triangle to find the distance
     * @param shift the initial displacement of the stl geometry
     * @param d [out] the distance function from octree-based Gauss Points to triangle edges
     */
    void ShortestDist2TriEdge(const TALYFEMLIB::ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift, double (&d)[DIM]);

    void ShortestDist2TriVertex(const TALYFEMLIB::ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift, double (&d)[DIM]);

    bool isPointOnLineSegment(const TALYFEMLIB::ZEROPTV p1, const TALYFEMLIB::ZEROPTV p2, const TALYFEMLIB::ZEROPTV p) {
        double totalDist = std::sqrt(std::pow(p2.x() - p1.x(), 2) + std::pow(p2.y() - p1.y(), 2));
        double dist1 = std::sqrt(std::pow(p.x() - p1.x(), 2) + std::pow(p.y() - p1.y(), 2));
        double dist2 = std::sqrt(std::pow(p2.x() - p.x(), 2) + std::pow(p2.y() - p.y(), 2));
//        if (std::abs(dist1 + dist2 - totalDist) > 1e-9)
//        {
//            std::cout << "std::abs(dist1 + dist2 - totalDist) = " << std::abs(dist1 + dist2 - totalDist) << "\n";
//        }
        return std::abs(dist1 + dist2 - totalDist) < 1e-9;
    };

    int findIndexOfSmallest(const std::vector<double>& vec) {
        if(vec.empty())
            throw std::invalid_argument("Cannot get the minimum element from an empty vector.");

        auto minElementPtr = std::min_element(vec.begin(), vec.end());
        return std::distance(vec.begin(), minElementPtr);
    };

    ZEROPTV findNormalVector(const ZEROPTV& p1, const ZEROPTV& p2) {
        // Calculate the differences
        double dx = p2.x() - p1.x();
        double dy = p2.y() - p1.y();

        // The normal vector components
        double normalX = dy;
        double normalY = -dx;

        ZEROPTV normal{normalX, normalY, 0.0};
        normal.SafeNormalize();
        return normal;
    };

    template<class T>
    inline const T& clamp( const T& v, const T& lo, const T& hi )
    {
        assert( !(hi < lo) );
        return (v < lo) ? lo : (hi < v) ? hi : v;
    };

public:
    /**
     * @brief constructor
     * @param fe the element we use to access gauss points
     * @param idata input Data
     * @param imga imga context, and we use this to access geometry
     * @param CarvedOutGeomID geometry ID of carved-out geometry
     */
    SBMCalc(const TALYFEMLIB::FEMElm &fe, NSInputData *idata, const IMGA *imga, const int CarvedOutGeomID);
    SBMCalc(NSInputData *idata, const IMGA *imga, const int CarvedOutGeomID);
    SBMCalc(const TALYFEMLIB::FEMElm &fe, NSInputData *idata, const IMGA *imga, std::vector<my_kd_tree_t *> kd_trees, const int CarvedOutGeomID);
    SBMCalc(NSInputData *idata, const IMGA *imga, std::vector<my_kd_tree_t *> kd_trees, const int CarvedOutGeomID);

    /**
     * @brief calculate the distance function for different kinds of geometries
     * @param d [out] distance function
     */
    void Dist2Geo(double (&d)[DIM]);

    void PtDist2Geo(const TALYFEMLIB::ZEROPTV &pt, double (&d)[DIM]);

    /**
     * @brief calculate the true normal of the SBM geometry
     * @param normal true normal of SBM geometry
     * @param d distance function
     */
    void NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM],bool OuterBoundary);
    void NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM]);
    void NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const double (&d)[DIM]);
};

SBMCalc::SBMCalc(const TALYFEMLIB::FEMElm &fe, NSInputData *idata, const IMGA *imga, std::vector<my_kd_tree_t *> kd_trees,
                 const int CarvedOutGeomID)
        : fe_(fe), imga_(imga), shift_(idata->carved_out_geoms_def[CarvedOutGeomID].InitialDisplacement), kd_trees_(kd_trees),CarvedOutGeomID_(CarvedOutGeomID)
{
    idata_ = idata;
}

SBMCalc::SBMCalc(NSInputData *idata, const IMGA *imga, std::vector<my_kd_tree_t *> kd_trees,
                 const int CarvedOutGeomID)
        : imga_(imga), shift_(idata->carved_out_geoms_def[CarvedOutGeomID].InitialDisplacement), kd_trees_(kd_trees),CarvedOutGeomID_(CarvedOutGeomID)
{
    idata_ = idata;
}


SBMCalc::SBMCalc(const TALYFEMLIB::FEMElm &fe, NSInputData *idata, const IMGA *imga
        , const int CarvedOutGeomID)
        : fe_(fe), imga_(imga), shift_(idata->carved_out_geoms_def[CarvedOutGeomID].InitialDisplacement),CarvedOutGeomID_(CarvedOutGeomID)
{
    idata_ = idata;
}

SBMCalc::SBMCalc(NSInputData *idata, const IMGA *imga
        , const int CarvedOutGeomID)
        : imga_(imga), shift_(idata->carved_out_geoms_def[CarvedOutGeomID].InitialDisplacement),CarvedOutGeomID_(CarvedOutGeomID)
{
    idata_ = idata;
}

bool SBMCalc::CheckInside3DTriangle(const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM], const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift)
{

    const TALYFEMLIB::ZEROPTV ptMove{pt(0) + d[0] - shift(0), pt(1) + d[1] - shift(1), pt(2) + d[2] - shift(2)};

    const TALYFEMLIB::ZEROPTV ptA{m_triangle.triangleCoord[0][0], m_triangle.triangleCoord[0][1], m_triangle.triangleCoord[0][2]};
    const TALYFEMLIB::ZEROPTV ptB{m_triangle.triangleCoord[1][0], m_triangle.triangleCoord[1][1], m_triangle.triangleCoord[1][2]};
    const TALYFEMLIB::ZEROPTV ptC{m_triangle.triangleCoord[2][0], m_triangle.triangleCoord[2][1], m_triangle.triangleCoord[2][2]};

    TALYFEMLIB::ZEROPTV u;
    TALYFEMLIB::ZEROPTV v;
    TALYFEMLIB::ZEROPTV w;

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

void SBMCalc::ShortestDist2TriEdge(const TALYFEMLIB::ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift, double (&d)[DIM])
{
    std::vector<TALYFEMLIB::ZEROPTV> ShortestVector2Line(3);
    std::vector<TALYFEMLIB::ZEROPTV> ptri(3);
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

void SBMCalc::ShortestDist2TriVertex(const TALYFEMLIB::ZEROPTV &pt, const GEOMETRY::Triangles &m_triangle, const TALYFEMLIB::ZEROPTV &shift, double (&d)[DIM])
{
    std::vector<TALYFEMLIB::ZEROPTV> ShortestVector2Vertex(3);
    std::vector<TALYFEMLIB::ZEROPTV> ptri(3);
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

void SBMCalc::Dist2Geo(double (&d)[DIM])
{
    const TALYFEMLIB::ZEROPTV pt = fe_.position();
    PtDist2Geo(pt, d);
}

void SBMCalc::NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM])
{
    /*
     * [update] it also consider NS, the normal is pointing inward if the case is circle
     */
    double scale = 0;
    /// doing below because of the lambda considerations
    if (imga_->ifInside(pt.data()))
    { // IN means INactive
        scale = -1;
    }
    else
    { // Active
        scale = 1;
#ifndef NDEBUG
//        std::cout << "[this is for debug] active GP\n";
#endif
    }

    double R2 = 0;

    for (int dim = 0; dim < DIM; dim++)
    {
        R2 += pow(d[dim], 2);
    }
    for (int dim = 0; dim < DIM; dim++)
    {
        normal(dim) = scale * d[dim] / sqrt(R2);
    }
}

/*
 * this is for sbm without lambda -> the gauss points on surrogate boundary will always treat as Inactive
 */
void SBMCalc::NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const double (&d)[DIM])
{

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

void SBMCalc::CalDist(double A, double B, double C, double x1, double y1, double &d)
{
    d = fabs(A * x1 + B * y1 + C) / sqrt(A * A + B * B);
}

DENDRITE_REAL SBMCalc::ElementSize(const TALYFEMLIB::FEMElm &fe)
{
    return pow((pow(2, DIM) * fe.volume_jacc()), (double)1 / DIM);
}

void SBMCalc::CalVec(double A, double B, double C, double x1, double y1, double &d1, double &d2)
{
    d1 = -A * (A * x1 + B * y1 + C) / (A * A + B * B);
    d2 = -B * (A * x1 + B * y1 + C) / (A * A + B * B);
}

void SBMCalc::PtDist2Geo(const TALYFEMLIB::ZEROPTV &pt, double (&d)[DIM])
{
    double x = pt.x();
    double y = pt.y();

#if (DIM == 2)

    size_t num_results = 1;
            std::vector<uint32_t> ret_index(num_results);
            std::vector<double> out_dist_sqr(num_results);

            const double query_pt[3] = {x - shift_[0], y - shift_[1], 0.0}; // shift_

            num_results = kd_trees_[CarvedOutGeomID_]->knnSearch(
                    &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

            const std::vector<GEOMETRY::Lines> *m_lines = &imga_->getGeometries()[CarvedOutGeomID_]->getMSH()->getLines();

            TALYFEMLIB::ZEROPTV OnePointVector,OnePointVectorOtherEnd;
//            d[0] = m_lines->at(ret_index[0]).lineCoord[0][0] + shift_[0] - x;
//            d[1] = m_lines->at(ret_index[0]).lineCoord[0][1] + shift_[1] - y;
            for (int dim =0;dim<DIM;dim++) {
                OnePointVector(dim) = m_lines->at(ret_index[0]).lineCoord[0][dim] + shift_[dim] - pt(dim);
                OnePointVectorOtherEnd(dim) = m_lines->at(ret_index[0]).lineCoord[1][dim] + shift_[dim] - pt(dim);
            }

//            TALYFEMLIB::ZEROPTV PickNormalVector= findNormalVector(TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[0]).lineCoord[0][0],
//                                                                             m_lines->at(ret_index[0]).lineCoord[0][1],0.0},
//                                                                  TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[0]).lineCoord[1][0],
//                                                                             m_lines->at(ret_index[0]).lineCoord[1][1],0.0});
//            TALYFEMLIB::ZEROPTV PickNormalVector2= findNormalVector(TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[1]).lineCoord[0][0],
//                                                                                        m_lines->at(ret_index[1]).lineCoord[0][1],0.0},
//                                                                    TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[1]).lineCoord[1][0],
//                                                                                        m_lines->at(ret_index[1]).lineCoord[1][1],0.0});

            ZEROPTV PickNormalVector= ZEROPTV{m_lines->at(ret_index[0]).normal[0],m_lines->at(ret_index[0]).normal[1],0.0};

            // scaling of vector
            double scale = 0.0;
            for (int dim = 0; dim < DIM; dim++)
            {
                scale += OnePointVector(dim) * PickNormalVector(dim);
            }

            TALYFEMLIB::ZEROPTV DistanceVector = PickNormalVector * scale;
            for (int dim = 0; dim < DIM; dim++)
            {
                d[dim] = DistanceVector(dim);
            }

            int PickClosestNumber  =  0;

            std::vector<double> Distance2EndPoint = {OnePointVector.norm(),
                                                     OnePointVectorOtherEnd.norm()};

            if(!isPointOnLineSegment(TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[PickClosestNumber]).lineCoord[0][0] + shift_[0],
                                             m_lines->at(ret_index[PickClosestNumber]).lineCoord[0][1] + shift_[1],0.0},
                                     TALYFEMLIB::ZEROPTV{m_lines->at(ret_index[PickClosestNumber]).lineCoord[1][0] + shift_[0],
                                             m_lines->at(ret_index[PickClosestNumber]).lineCoord[1][1] + shift_[1],0.0},
                                     pt+TALYFEMLIB::ZEROPTV{d[0],d[1],0.0})){

//                std::cout << "projection point is not on the line segment\n";

                // distance to line endpoint

                int IndexOfSmallest = findIndexOfSmallest(Distance2EndPoint);
//                std::cout << " --------------------------------\n";
//                std::cout << " IndexOfSmallest = " << IndexOfSmallest << "\n";
//                std::cout << " --------------------------------\n";

                if (IndexOfSmallest == 0) {
                    for (int dim = 0; dim < DIM; dim++) {
                        d[dim] = OnePointVector(dim);
                    }
                }
                else if (IndexOfSmallest == 1){
                    for (int dim = 0; dim < DIM; dim++) {
                        d[dim] = OnePointVectorOtherEnd(dim);
                    }
                }
            }



//            std::ofstream fout("DistanceFunc.txt",std::ios::app);
//            fout << x << " " << y << " ";
//            fout << d[0] << " " << d[1] << " ";
//            fout << x+ d[0] << " " << y+d[1] << "\n";



#endif

#if (DIM == 3)

    /// beacause remember to implement how to calculate distance function if you implement new SBMGeo

    double z = pt.z();



    double MinDist = 100.0;
    TALYFEMLIB::ZEROPTV OnePointVector;
    TALYFEMLIB::ZEROPTV PickNormalVector;
    int PickGeomID = 0;
    int PickTrianleID = 0;


    size_t num_results = 1;
    std::vector<uint32_t> ret_index(num_results);
    std::vector<double> out_dist_sqr(num_results);



    const double query_pt[3] = {x - shift_[0], y - shift_[1], z - shift_[2]}; // shift_

    num_results = kd_trees_[CarvedOutGeomID_]->knnSearch(
            &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

    const std::vector<GEOMETRY::Triangles> *m_triangles = &imga_->getGeometries()[CarvedOutGeomID_]->getSTL()[0].getTriangles();
    for (int dim = 0; dim < DIM; dim++)
    {
        OnePointVector(dim) = m_triangles->at(ret_index[0]).triangleCoord[0][dim] + shift_[dim] - pt(dim);
        PickNormalVector(dim) = m_triangles->at(ret_index[0]).normal[dim];
    }
    PickTrianleID = ret_index[0];
    PickGeomID = CarvedOutGeomID_;


    // scaling of vector
    double scale = 0.0;
    for (int dim = 0; dim < DIM; dim++)
    {
        scale += OnePointVector(dim) * PickNormalVector(dim);
    }

    for (int dim = 0; dim < DIM; dim++)
    {
        d[dim] = scale * PickNormalVector(dim);
    }

    GEOMETRY::Triangles m_triangle = imga_->getGeometries()[PickGeomID]->getSTL()[0].getTriangles()[PickTrianleID];
    if (!CheckInside3DTriangle(pt, d, m_triangle, shift_))
    {
        ShortestDist2TriEdge(pt, m_triangle, shift_, d);
        if (!CheckInside3DTriangle(pt, d, m_triangle, shift_)){
            ShortestDist2TriVertex(pt, m_triangle, shift_, d);
        }
    }
#endif
}

void SBMCalc::NormalofGeo(TALYFEMLIB::ZEROPTV &normal, const TALYFEMLIB::ZEROPTV &pt, const double (&d)[DIM], bool OuterBoundary) {
    /*
     * this consider NSHT -> flow passing through cylinder cases
     */
    double scale = 0;
    if (imga_->ifInside(pt.data()))
    { // IN means INactive
        scale = -1;
    }
    else
    { // Active
        scale = 1;
    }

    /*
     * for NS, most cases the geom is not the outer boundary and we should flip it
     */
    int flip = (OuterBoundary)? 1: -1;

    double R2 = 0;

    for (int dim = 0; dim < DIM; dim++)
    {
        R2 += pow(d[dim], 2);
    }
    for (int dim = 0; dim < DIM; dim++)
    {
        normal(dim) = flip* scale * d[dim] / sqrt(R2);
    }
}

#endif // DENDRITEKT_SBMCalc_H
