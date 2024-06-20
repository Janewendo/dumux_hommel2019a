/*
 * HCO3.hh
 *
 *  Created on: 14.03.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the HCO3 fluid properties
 */
#ifndef DUMUX_HCO3_HH
#define DUMUX_HCO3_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the HCO3 fluid properties
 */
template <class Scalar>
class BicarbonateIon
: public Components::Base<Scalar, BicarbonateIon<Scalar> >
, public Components::Ion<Scalar, BicarbonateIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for HCO3.
    */
    static std::string name()
    { return "HCO3-"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of HCO3.
    */
    static Scalar molarMass()
    { return 61.01714e-3; } // kg/mol

   /*!
    * \brief The charge of HCO3.
    */
    static constexpr int charge()
    { return -1; }

};

} // end namespace Components
} // end namespace Dumux

#endif

