/*
 * Fe.hh
 *
 *  Created on: 16.06.2024
 *      Author: Jianwen
 */

/*!
 * \file
 *
 * \brief A class for the Fe2+ fluid properties
 */
#ifndef DUMUX_Fe2_HH
#define DUMUX_Fe2_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Fe2 (Iron ion) fluid properties
 */
template <class Scalar>
class Iron2Ion
: public Components::Base<Scalar, Iron2Ion<Scalar> >
, public Components::Ion<Scalar, Iron2Ion<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Fe2.
    */
    static std::string name()
    { return "Fe2+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Fe2.
    */
    static Scalar molarMass()
    { return 55.845e-3; } // kgFe/molFe
   /*!
    * \brief The charge of the Fe2 ion.
    */
    static constexpr int charge()
    {
        return +2;
    }

};

} // end namespace Components
} // end namespace Dumux

#endif


