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



#include "../elem.h"
#include "../elem_common.h"
#include "../../basis/basis.h"

namespace TALYFEMLIB {

/**
 * 2D triangle element.
 */
    class ELEMPyramid : public ELEM {
    public:
        ElemType elmType() const override { return kElemPyramid; }
        int GetNodesPerSurface() const override;
        const int* GetSurfaceCheckArray() const override;
        int GetSurfaceCount() const override { return 0; }
        ZEROPTV CalculateNormal(const GRID* grid, int surface_id) const override;
        int nsd() const override { return 3; }
        void Validate(const GRID* grid) override{};
        kBasisFunction basis_function() const override;
    };

}  // namespace TALYFEMLIB
