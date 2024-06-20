/*
 * hPlus.h
 *
 *  Created on: 14.03.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the H+ fluid properties
 */
#ifndef DUMUX_HYDRON_HH
#define DUMUX_HYDRON_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the H+ (hydron or proton) ion properties
 */
template <class Scalar>
class HydronIon
: public Components::Base<Scalar, HydronIon<Scalar> >
, public Components::Ion<Scalar, HydronIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for H+.
    */
    static const char *name()
    { return "H+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of H+.
    */
    static Scalar molarMass()
    { return 1.00794e-3; }

   /*!
    * \brief The charge of H.
    */
    static constexpr int charge()
    { return +1; }

};

} // end namespace Components
} // end namespace Dumux

#endif
