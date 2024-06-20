// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Components
 * \brief A class for the Fe(OH)2 mineral phase properties
 */
#ifndef DUMUX_FERROHYDRITE_HH
#define DUMUX_FERROHYDRITE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

#include <dumux/material/components/hydroxideion.hh>
#include <dumux/material/components/iron2ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Fe(OH)2 mineral phase properties
 */
template <class Scalar>
class Ferrohydrite
: public Components::Base<Scalar, Ferrohydrite<Scalar> >
, public Components::Solid<Scalar, Ferrohydrite<Scalar> >
{

public:
    using HydroxideIon = Components::HydroxideIon<Scalar>;
    using Iron2Ion = Components::Iron2Ion<Scalar>;
    /*!
     * \brief A human readable name for Ferrohydrite.
     */
    static std::string name()
    { return "Ferrohydrite"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Ferrohydrite.
     */
    static constexpr Scalar molarMass()
    { return HydroxideIon::molarMass() + Iron2Ion::molarMass(); } // kg/mol

    /*!
     * \brief Returns true if the solid phase is assumed to be compressible
     */
    static constexpr bool solidIsCompressible()
    { return false; }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static constexpr Scalar solidDensity(Scalar temperature)
    { return 3.4e3; }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    { return 0.5; }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    { return 420; }
};

} // end namespace Components
} // end namespace Dumux

#endif

