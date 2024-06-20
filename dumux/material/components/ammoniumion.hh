/*
 * NH4.hh
 *
 *  Created on: 20.04.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the NH4 fluid properties
 */
#ifndef DUMUX_NH4_HH
#define DUMUX_NH4_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the NH4 fluid properties
 */
template <class Scalar>
class AmmoniumIon
: public Components::Base<Scalar,  AmmoniumIon<Scalar> >
, public Components::Ion<Scalar,  AmmoniumIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for NH4.
    */
    static std::string name()
    { return "NH4+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of NH4.
    */
    static Scalar molarMass()
    { return 0.018039; } // kg/mol

   /*!
    * \brief The charge of H.
    */
    static constexpr int charge()
    { return +1; }

};

} // end namespace Components
} // end namespace Dumux

#endif

