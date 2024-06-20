// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_ICP_SPATIAL_PARAMS_HH
#define DUMUX_ICP_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilityvermapruess.hh>


namespace Dumux{
//forward declaration
template<class TypeTag>
class ICPSpatialParams;

/*!
* \ingroup 2PICPModel
* \brief Definition of the spatial parameters for the two-phase
* induced calcium carbonate precipitation problem.
*/
// template<class TypeTag>
// class ICPSpatialParams : public FVSpatialParams<TypeTag>
template<class TypeTag>
class ICPSpatialParams
: public FVSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                         GetPropType<TypeTag, Properties::Scalar>,
                         ICPSpatialParams<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ICPSpatialParams<TypeTag>>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using CoordScalar = typename GridView::ctype;
    enum { dimWorld=GridView::dimensionworld };

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

    ICPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {

        try
        {
        referencePorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity", 0.4);
        referencePermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability", 2.e-10);
        critPorosity_ = getParam<Scalar>("SpatialParams.CritPorosity", 0.0);

        }
        catch (Dumux::ParameterException &e)
        {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }

        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(0.05);

        // parameters for the Brooks-Corey law
        materialParams_.setPe(1e4);
        materialParams_.setLambda(2.0);
/*
                //thermalconductivity
                lambdaSolid_ = 14.7; //[W/(m*K)] Acosta et al. [2006]*/
    }

    /*!
     *  \brief Define the reference porosity \f$[-]\f$ distribution.
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }


    /*!
     *  \brief Define the volume fraction of the inert component
     *
     *  \param globalPos The global position in the domain
     *  \param compIdx The index of the inert solid component
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 1.0-referencePorosity_; }

    /*!
     *  \brief Return the actual recent porosity \f$[-]\f$ accounting for
     *  clogging caused by mineralization
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol, referencePorosity_); }


    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *
     *  Solution dependent permeability function
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        const auto poro = porosity(element, scv, elemSol);
        return permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, critPorosity_, poro);
    }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    MaterialLawParams materialParams_;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PorosityLaw = PorosityPrecipitation<Scalar, ModelTraits::numFluidComponents(), SolidSystem::numComponents - SolidSystem::numInertComponents>;
    PorosityLaw poroLaw_;
    PermeabilityVermaPruess<PermeabilityType> permLaw_;
    Scalar referencePorosity_;
    Scalar critPorosity_;
    PermeabilityType referencePermeability_ = 0.0;
};

}

#endif
