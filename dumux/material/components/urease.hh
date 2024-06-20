/*
 * urease.hh
 *
 *  Created on: 20.04.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the urease fluid properties
 */
#ifndef DUMUX_UREASE_HH
#define DUMUX_UREASE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Urease fluid properties
 */
template <class Scalar>
class Urease
: public Components::Base<Scalar, Urease<Scalar> >
, public Components::Solid<Scalar, Urease<Scalar> >
{
public:
    /*!
     * \brief A human readable name for Urease.
     */
    static std::string name()
    { return "Urease"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ urease
     */
    static Scalar molarMass()  //
    { return 1.0; } //545; //} // kg/mol 545 kDa for jack bean urease source: wikipedia, is that too high?

    /*!
     * \brief The density in [kg/m³] of attached urease
     */
    static Scalar solidDensity(Scalar temperature)  //
    { return 1100.0;} // TODO assumption check for the right density!

};

} // end namespace Components
} // end namespace Dumux

#endif
