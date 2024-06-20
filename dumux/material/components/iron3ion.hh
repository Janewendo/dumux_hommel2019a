/*
 * Fe.hh
 *
 *  Created on: 16.06.2024
 *      Author: Jianwen
 */

/*!
 * \file
 *
 * \brief A class for the Fe3+ fluid properties
 */
#ifndef DUMUX_Fe3_HH
#define DUMUX_Fe3_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Fe3 (Iron ion) fluid properties
 */
template <class Scalar>
class Iron3Ion
: public Components::Base<Scalar, Iron3Ion<Scalar> >
, public Components::Ion<Scalar, Iron3Ion<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Fe3.
    */
    static std::string name()
    { return "Fe3+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Fe3.
    */
    static Scalar molarMass()
    { return 55.845e-3; } // kgFe/molFe
   /*!
    * \brief The charge of the Fe3 ion.
    */
    static constexpr int charge()
    {
        return +3;
    }

};

} // end namespace Components
} // end namespace Dumux

#endif


