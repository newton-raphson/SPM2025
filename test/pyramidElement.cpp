//
// Created by maksbh on 3/25/23.
//
#include <DendriteUtils.h>
#include <talyfem/grid/elem_types/elemPyramid.h>
#include <talyfem/grid/elem_types/elem2dtriangle.h>

int main(int argc, char *argv[]) {
    TALYFEMLIB::ELEM *taly_elem;
    TALYFEMLIB::GRID grid;

    static constexpr int numNodes = 5;
    int node_id_array[3];
    grid.redimArrays(numNodes, 1);
    for (int i = 0; i < numNodes; i++) {
        grid.node_array_[i] = new TALYFEMLIB::NODE();
        node_id_array[i] = i;
    }

    taly_elem = new TALYFEMLIB::ELEMPyramid();
    grid.elm_array_[0] = taly_elem;


    taly_elem->redim(numNodes, node_id_array);
//    https://www.dealii.org/current/doxygen/deal.II/reference__cell_8h_source.html

    static constexpr int dim = 3;
    static const Point<dim> coords[5] = {Point<dim>{-1.0, -1.0, 0.0},
                                           Point<dim>{+1.0, -1.0, 0.0},
                                           Point<dim>{-1.0, +1.0, 0.0},
                                           Point<dim>{+1.0, +1.0, 0.0},
                                           Point<dim>{+1.0, +1.0, 1.0}};
    for(int i = 0; i < numNodes; i++){
        grid.node_array_[i]->setCoor(coords[i].x(), coords[i].y(), coords[i].z());
    }

    TALYFEMLIB::FEMElm fe(&grid, TALYFEMLIB::BASIS_LINEAR | TALYFEMLIB::BASIS_FIRST_DERIVATIVE);
    fe.refill(0, 0);
    double val = 0.0;
    while (fe.next_itg_pt()){
        val += fe.detJxW();
    }
    std::cout << val << "\n";

    delete taly_elem;




}