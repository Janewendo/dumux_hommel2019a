/*
 * NH3.hh
 *
 *  Created on: 18.5.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the NH3 component properties
 */
#ifndef DUMUX_NH3_HH
#define DUMUX_NH3_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the NH3 component properties
 */
template <class Scalar>
class Ammonia
: public Components::Base<Scalar, Ammonia<Scalar> >
{
public:
   /*!
    * \brief A human readable name for the NH3.
    */
    static std::string name()
    { return "NH3"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of NH3.
    */
    static Scalar molarMass()
    { return 0.017031; } // kg/mol

   /*!
    * \brief The charge of NH3.
    */
    static Scalar charge()
    { return 0.0; }
};

} // end namespace Components
} // end namespace Dumux

#endif

