/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the Na2+ fluid properties
 */
#ifndef DUMUX_NA_HH
#define DUMUX_NA_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Na (Sodium Ion) fluid properties
 */
template <class Scalar>
class SodiumIon
: public Components::Base<Scalar, SodiumIon<Scalar> >
, public Components::Ion<Scalar, SodiumIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Na.
    */
    static std::string name()
    { return "Na+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Na.
    */
    static Scalar molarMass()
    { return 22.9898e-3; } // kgNa/molNa

   /*!
    * \brief The charge of the Na ion.
    */
    static constexpr int charge()
    {
        return +1;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

