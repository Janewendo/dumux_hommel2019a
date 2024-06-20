/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the Cl- fluid properties
 */
#ifndef DUMUX_CL_HH
#define DUMUX_CL_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Na (Natrium ion) fluid properties
 */
template <class Scalar>
class ChlorideIon
: public Components::Base<Scalar, ChlorideIon<Scalar> >
, public Components::Ion<Scalar, ChlorideIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Cl.
    */
    static std::string name()
    { return "Cl-"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Cl.
    */
    static Scalar molarMass()
    { return 35.453e-3; } // kgNa/molNa
   /*!
    * \brief The charge of the Cl ion.
    */
    static constexpr int charge()
    {
        return -1;
    }

};

} // end namespace Components
} // end namespace Dumux

#endif

