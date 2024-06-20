/*
 * OH.hh
 *
 *  Created on: 14.03.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the OH- fluid properties
 */
#ifndef DUMUX_HYDROXIDE_HH
#define DUMUX_HYDROXIDE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Hydroxide ion properties
 */
template <class Scalar>
class HydroxideIon
: public Components::Base<Scalar, HydroxideIon<Scalar> >
, public Components::Ion<Scalar, HydroxideIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the OH.
     */
    static std::string name()
    { return "OH-"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of OH.
     */
    static Scalar molarMass()
    { return 16.9994; } // kg/mol

   /*!
    * \brief The charge of OH.
    */
    static constexpr int charge()
    { return -1; }

};

} // end namespace Components
} // end namespace Dumux

#endif

