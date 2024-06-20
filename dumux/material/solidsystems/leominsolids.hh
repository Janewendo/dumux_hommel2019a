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
 * \ingroup SolidSystems
 * \brief @copybrief Dumux::SolidSystems::InertSolidPhase
 */
#ifndef DUMUX_SOLIDSYSTEMS_EICP_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_EICP_SOLID_PHASE_HH

#include <string>
#include <dune/common/exceptions.hh>

#include <dumux/material/components/urease.hh>
#include <dumux/material/components/calcite.hh>
#include <dumux/material/components/granite.hh>

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of a single inert solid component and two reactive solid components
 * \note a solid is considered inert if it can't dissolve in a liquid and
 *       and can't increase its mass by precipitation from a fluid phase.
 * \note inert components have to come after all non-inert components
 */
template <class Scalar>
class EnzymeMinSolidPhase
{
public:
    using Jbme = Components::Urease<Scalar>;
    using Calcite = Components::Calcite<Scalar>;
    using Granite = Components::Granite<Scalar>;

    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 3;
    static constexpr int numInertComponents = 1;
    static constexpr int JbmeIdx = 0;
    static constexpr int CalciteIdx = 1;
    static constexpr int GraniteIdx = 2;

    /*!
     * \brief Return the human readable name of a solid phase
     *
     * \param compIdx The index of the solid phase to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case JbmeIdx: return Jbme::name();
            case CalciteIdx: return Calcite::name();
            case GraniteIdx: return Granite::name();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief A human readable name for the solid system.
     */
    static std::string name()
    { return "EnzymeMinSolidPhase"; }

    /*!
     * \brief Returns whether the phase is incompressible
     */
    static constexpr bool isCompressible(int compIdx)
    { return false; }

    /*!
     * \brief Returns whether all components are inert (don't react)
     */
    static constexpr bool isInert()
    { return (numComponents == numInertComponents); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    const static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case JbmeIdx: return Jbme::molarMass();
            case CalciteIdx: return Calcite::molarMass();
            case GraniteIdx: return Granite::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
        const Scalar rho1 = Jbme::solidDensity(solidState.temperature());
        const Scalar rho2 = Calcite::solidDensity(solidState.temperature());
        const Scalar rho3 = Granite::solidDensity(solidState.temperature());
        const Scalar volFrac1 = solidState.volumeFraction(JbmeIdx);
        const Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        const Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (rho1*volFrac1+
                rho2*volFrac2+
                rho3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    const static Scalar density(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case JbmeIdx: return Jbme::solidDensity(solidState.temperature());
            case CalciteIdx: return Calcite::solidDensity(solidState.temperature());
            case GraniteIdx: return Granite::solidDensity(solidState.temperature());
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief The molar density of the solid phase at a given pressure and temperature.
     */
    template <class SolidState>
    const static Scalar molarDensity(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case JbmeIdx: return Jbme::solidDensity(solidState.temperature())/Jbme::molarMass();
            case CalciteIdx: return Calcite::solidDensity(solidState.temperature())/Calcite::molarMass();
            case GraniteIdx: return Granite::solidDensity(solidState.temperature())/Granite::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        const Scalar lambda1 = Jbme::solidThermalConductivity(solidState.temperature());
        const Scalar lambda2 = Calcite::solidThermalConductivity(solidState.temperature());
        const Scalar lambda3 = Granite::solidThermalConductivity(solidState.temperature());
        const Scalar volFrac1 = solidState.volumeFraction(JbmeIdx);
        const Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        const Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (lambda1*volFrac1+
                lambda2*volFrac2+
                lambda3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }

    /*!
     * \brief Specific isobaric heat capacity of the pure solids \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        const Scalar c1 = Jbme::solidHeatCapacity(solidState.temperature());
        const Scalar c2 = Calcite::solidHeatCapacity(solidState.temperature());
        const Scalar c3 = Granite::solidHeatCapacity(solidState.temperature());
        const Scalar volFrac1 = solidState.volumeFraction(JbmeIdx);
        const Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        const Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (c1*volFrac1+
                c2*volFrac2+
                c3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
