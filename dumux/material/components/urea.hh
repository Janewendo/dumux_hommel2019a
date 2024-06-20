/*
 * urea.hh
 *
 *  Created on: 20.04.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the urea fluid properties
 */
#ifndef DUMUX_UREA_HH
#define DUMUX_UREA_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the urea fluid properties
 */
template <class Scalar>
class Urea
: public Components::Base<Scalar, Urea<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Urea.
    */
    static std::string name()
    { return "Urea"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Urea.
    */
    static Scalar molarMass()
    { return 0.0606; } // kg/mol

   /*!
    * \brief The charge of Urea.
    */
    static Scalar charge()
    {
        return 0.0;
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
