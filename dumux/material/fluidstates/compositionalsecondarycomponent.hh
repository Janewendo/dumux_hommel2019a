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
 *
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_COMPOSITIONAL_SEC_COMP_FLUID_STATE_HH
#define DUMUX_COMPOSITIONAL_SEC_COMP_FLUID_STATE_HH

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <cassert>

#include <dune/common/exceptions.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{
/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 *        The difference to the parent CompositionalFluidState is
 *        that secondary component mole fractions can be set, e.g.
 *        calculated from the dissociation of a regular component
 *        into several sub- or secondary components.
 *        E.g. total inorganic carbon into CO_2, HCO_3^-, CO_3^{2-}.
 */
template <class ScalarType, class FluidSystem>
class CompositionalSecCompFluidState : public CompositionalFluidState<ScalarType, FluidSystem>
{
    using ParentType = CompositionalFluidState<ScalarType, FluidSystem>;
public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int numSecComponents = FluidSystem::numSecComponents;

    //! export the scalar type
    using Scalar = ScalarType;

    //! default constructor
    CompositionalSecCompFluidState() = default;

    //! copy constructor from arbitrary fluid state
    template <class FluidState, typename std::enable_if_t<!std::is_same<FluidState, CompositionalSecCompFluidState>::value, int> = 0>
    CompositionalSecCompFluidState(const FluidState &fs)
    { assign(fs); }

    // copy and move constructor / assignment operator
    CompositionalSecCompFluidState(const CompositionalSecCompFluidState &fs) = default;
    CompositionalSecCompFluidState(CompositionalSecCompFluidState &&fs) = default;
    CompositionalSecCompFluidState& operator=(const CompositionalSecCompFluidState &fs) = default;
    CompositionalSecCompFluidState& operator=(CompositionalSecCompFluidState &&fs) = default;

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/

    /*!
     * \copydoc CompositionalFluidState::moleFraction
     *
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {   if(compIdx<numComponents)
            return ParentType::moleFraction_[phaseIdx][compIdx];
        else
            return
                moleFractionSecComp_[phaseIdx][compIdx-numComponents];
    }

    /*!
     * \copydoc CompositionalFluidState::massFraction
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        // calculate the mass fractions:
        // for "mass" models this is just a back calculation
        if(compIdx<numComponents)
            {
                return ParentType::sumMoleFractions_[phaseIdx]
                    * moleFraction(phaseIdx, compIdx)
                    * FluidSystem::molarMass(compIdx)
                    / ParentType::averageMolarMass_[phaseIdx];
            }
        else
            {
            return ParentType::sumMoleFractions_[phaseIdx]
                * moleFractionSecComp_[phaseIdx][compIdx-numComponents]
                * FluidSystem::molarMass(compIdx)
                / ParentType::averageMolarMass_[phaseIdx];
            }
    }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     *
     * \note If the other fluid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class FluidState>
    void assign(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            ParentType::averageMolarMass_[phaseIdx] = 0;
            ParentType::sumMoleFractions_[phaseIdx] = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                ParentType::moleFraction_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
                ParentType::fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
                ParentType::averageMolarMass_[phaseIdx] += ParentType::moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
                ParentType::sumMoleFractions_[phaseIdx] += ParentType::moleFraction_[phaseIdx][compIdx];
            }
            for (int compIdx = 0; compIdx < numSecComponents; ++compIdx) {
                moleFractionSecComp_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx + numComponents);
            }
            ParentType::averageMolarMass_[phaseIdx] = fs.averageMolarMass(phaseIdx);
            ParentType::pressure_[phaseIdx] = fs.pressure(phaseIdx);
            ParentType::saturation_[phaseIdx] = fs.saturation(phaseIdx);
            ParentType::density_[phaseIdx] = fs.density(phaseIdx);
            ParentType::enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            ParentType::viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
        }
        ParentType::wPhaseIdx_ = fs.wettingPhase();
    }

    /*!
     * \brief Set the mole fraction of a secondary component in a phase []
     */
    void setMoleFractionSecComp(int phaseIdx, int compIdx, Scalar value)
    {
        moleFractionSecComp_[phaseIdx][compIdx-numComponents] = value;
    }

protected:
    //! zero-initialize all data members with braces syntax
    Scalar moleFractionSecComp_[numPhases][numComponents] = {};
};

} // end namespace Dumux

#endif
