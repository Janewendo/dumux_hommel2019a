// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief The Verma-Pruess relationship for the calculation of a porosity-dependent permeability.
 */
#ifndef DUMUX_PERMEABILITY_VERMA_PRUESS_HH
#define DUMUX_PERMEABILITY_VERMA_PRUESS_HH

#include <cmath>
#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The Verma-Pruess relationship for the calculation of a porosity-dependent permeability.
 *        When the porosity is implemented as solution-independent, using this relationship for the
 *        permeability leads to unnecessary overhead.
 */
template<class PermeabilityType>
class PermeabilityVermaPruess
{
public:
    /*!
     * \brief calculates the permeability for a given sub-control volume
     * \param refPerm Reference permeability before porosity changes
     * \param refPoro The poro corresponding to the reference permeability
     * \param critPoro The poro corresponding to a zero permeability
     * \param poro The porosity for which permeability is to be evaluated
     */
    template<class Scalar>
    PermeabilityType evaluatePermeability(PermeabilityType refPerm, Scalar refPoro, Scalar critPoro, Scalar poro) const
    {
        using std::pow;
        auto factor = pow((poro-critPoro)/(refPoro-critPoro), 3);
        refPerm *= factor;
        return refPerm;
    }
};

} // namespace Dumux


#endif
