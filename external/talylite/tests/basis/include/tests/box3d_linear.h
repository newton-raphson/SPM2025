/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once

#include <basistest.h>

// Yes, this should be split into a .cpp/.h pair.
// But this is only ever #included in one file (main.cpp) and I'm lazy.

class Box3DLinearBasisTest : public BasisTest
{
public:
  std::string name() override { return "3D Box (Linear)"; }
  ElemType elm_type() override { return kElem3dHexahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[8] = {
      ZEROPTV(-7.0/5,-4.0/5,-1.0/5), ZEROPTV(2.0/5,-4.0/5,-1.0/5), ZEROPTV(2.0/5,3.0/5,-1.0/5), ZEROPTV(-7.0/5,3.0/5,-1.0/5),
      ZEROPTV(-7.0/5,-4.0/5,2.0/5), ZEROPTV(2.0/5,-4.0/5,2.0/5), ZEROPTV(2.0/5,3.0/5,2.0/5), ZEROPTV(-7.0/5,3.0/5,2.0/5)
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[8] = {
      ZEROPTV(-sqrt(3.0)/3, -sqrt(3.0)/3, -sqrt(3.0)/3),
      ZEROPTV(sqrt(3.0)/3, -sqrt(3.0)/3, -sqrt(3.0)/3),
      ZEROPTV(-sqrt(3.0)/3,  sqrt(3.0)/3, -sqrt(3.0)/3),
      ZEROPTV(sqrt(3.0)/3,  sqrt(3.0)/3, -sqrt(3.0)/3),
      ZEROPTV(-sqrt(3.0)/3, -sqrt(3.0)/3,  sqrt(3.0)/3),
      ZEROPTV(sqrt(3.0)/3, -sqrt(3.0)/3,  sqrt(3.0)/3),
      ZEROPTV(-sqrt(3.0)/3,  sqrt(3.0)/3,  sqrt(3.0)/3),
      ZEROPTV(sqrt(3.0)/3,  sqrt(3.0)/3,  sqrt(3.0)/3)
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1;
  }

  double N(int itg_pt, int bf) override {
    static const double n[8][8] = {
      {
        sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
      },
      {
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, 
      },
      {
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
      },
      {
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
      },
      {
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
      },
      {
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
      },
      {
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, 
      },
      {
        sqrt(3.0)*(-5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(-1.0/3.6E1)+1.0/1.2E1, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
        sqrt(3.0)*(5.0/3.6E1)+1.0/4.0, sqrt(3.0)*(1.0/3.6E1)+1.0/1.2E1, 
      }
    };
    return n[itg_pt][bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const double dnde[8][8][3] = {
      {
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, -1.0/1.2E1, -1.0/1.2E1 },
        { 1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, -1.0/1.2E1 },
        { -1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, 1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, 1.0/1.2E1, 1.0/1.2E1 }
      },
      {
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, -1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, -1.0/1.2E1 },
        { -1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, 1.0/1.2E1 },
        { 1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, 1.0/1.2E1, 1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 }
      },
      {
        { -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, -1.0/1.2E1 },
        { 1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, 1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, -1.0/1.2E1, 1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, 1.0/1.2E1 },
        { -1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 }
      },
      {
        { -1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, -1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, 1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, -1.0/1.2E1, 1.0/1.2E1 },
        { 1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, 1.0/1.2E1 }
      },
      {
        { -1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, -1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, 1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, -1.0/1.2E1, 1.0/1.2E1 },
        { 1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, 1.0/1.2E1 }
      },
      {
        { -1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, -1.0/1.2E1 },
        { 1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, 1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, -1.0/1.2E1, 1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, 1.0/1.2E1 },
        { -1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 }
      },
      {
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, -1.0/1.2E1, -1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, -1.0/1.2E1 },
        { -1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, 1.0/1.2E1 },
        { 1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, 1.0/1.2E1, 1.0/1.2E1 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 }
      },
      {
        { sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0, sqrt(3.0)*(1.0/1.2E1)-1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, -1.0/1.2E1, -1.0/1.2E1 },
        { 1.0/1.2E1, 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0 },
        { -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0, -1.0/1.2E1 },
        { -1.0/1.2E1, -1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)+1.0/6.0 },
        { 1.0/1.2E1, sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, 1.0/1.2E1 },
        { sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0, sqrt(3.0)*(1.0/1.2E1)+1.0/6.0 },
        { sqrt(3.0)*(-1.0/1.2E1)-1.0/6.0, 1.0/1.2E1, 1.0/1.2E1 }
      },

    };
    return dnde[itg_pt][bf][axis];
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[3][3] = {
      { 9.0/1.0E1, 0.0, 0.0 },
      { 0.0, 7.0/1.0E1, 0.0 },
      { 0.0, 0.0, 3.0/1.0E1 }
    };
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[3][3] = {
      { 2.1E1/1.0E2, 0.0, 0.0 },
      { 0.0, 2.7E1/1.0E2, 0.0 },
      { 0.0, 0.0, 6.3E1/1.0E2 }
    };
    return cof_values[i][j];
  }

  double jacobian_det(int itg_pt) override {
    return 189.0/1000.0;
  }

  double dN(int itg_pt, int i, int axis) override {
    static const double dn[8][8][3] = {
      {
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, -5.0/4.2E1, -5.0/1.8E1 },
        { 5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, -5.0/1.8E1 },
        { -5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, 5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, 5.0/4.2E1, 5.0/1.8E1 }
      },
      {
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, -5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, -5.0/1.8E1 },
        { -5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, 5.0/1.8E1 },
        { 5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, 5.0/4.2E1, 5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 }
      },
      {
        { -5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, -5.0/1.8E1 },
        { 5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, 5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, -5.0/4.2E1, 5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, 5.0/1.8E1 },
        { -5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 }
      },
      {
        { -5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, -5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, 5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, -5.0/4.2E1, 5.0/1.8E1 },
        { 5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, 5.0/1.8E1 }
      },
      {
        { -5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, -5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, 5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, -5.0/4.2E1, 5.0/1.8E1 },
        { 5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, 5.0/1.8E1 }
      },
      {
        { -5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, -5.0/1.8E1 },
        { 5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, 5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, -5.0/4.2E1, 5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, 5.0/1.8E1 },
        { -5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 }
      },
      {
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, -5.0/4.2E1, -5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, -5.0/1.8E1 },
        { -5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, 5.0/1.8E1 },
        { 5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, 5.0/4.2E1, 5.0/1.8E1 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 }
      },
      {
        { sqrt(3.0)*(5.0/5.4E1)-5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)-5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)-5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)+5.0/2.7E1, -5.0/4.2E1, -5.0/1.8E1 },
        { 5.0/5.4E1, 5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)-5.0/9.0 },
        { -5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)+5.0/2.1E1, -5.0/1.8E1 },
        { -5.0/5.4E1, -5.0/4.2E1, sqrt(3.0)*(-5.0/1.8E1)+5.0/9.0 },
        { 5.0/5.4E1, sqrt(3.0)*(-5.0/4.2E1)-5.0/2.1E1, 5.0/1.8E1 },
        { sqrt(3.0)*(5.0/5.4E1)+5.0/2.7E1, sqrt(3.0)*(5.0/4.2E1)+5.0/2.1E1, sqrt(3.0)*(5.0/1.8E1)+5.0/9.0 },
        { sqrt(3.0)*(-5.0/5.4E1)-5.0/2.7E1, 5.0/4.2E1, 5.0/1.8E1 }
      }
    };
    return dn[itg_pt][i][axis];
  }
};

