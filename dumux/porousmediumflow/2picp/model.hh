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
/*!
 * \file
 *
 * \ingroup TwoPICPModel, Mineralization
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase induced calcium carbonate precipitation model.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 * Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 * The primary variable of the solid phases is the volume fraction \f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
 * Additionally, secondary components are used, which have no mass balance, but are calculated using
 * chemical equilibrium conditions (electroneutrality of the solution, laws of mass action) from the
 * respective integral component which has a mass balance solved for.
 * E.g. total inorganic carbon with the secondary components CO2 or H2CO3, HCO3, CO3.
 *
 * The source an sink terms link the mass balances of the n-transported component to the solid phases.
 * The porosity \f$\phi\f$ is updated according to the reduction of the initial (or solid-phase-free porous medium) porosity \f$\phi_0\f$
 * by the accumulated volume fractions of the solid phases:
 * \f$ \phi = \phi_0 - \sum (\phi_\lambda)\f$
 * Additionally, the permeability is updated depending on the current porosity.
 */
#ifndef DUMUX_2PICP_MODEL_HH
#define DUMUX_2PICP_MODEL_HH

#include <dumux/porousmediumflow/2pnc/model.hh>

#include <dumux/material/solidstates/compositionalsolidstate.hh>
#include <dumux/material/fluidstates/compositionalsecondarycomponent.hh>

#include <dumux/porousmediumflow/mineralization/model.hh>
#include <dumux/porousmediumflow/mineralization/localresidual.hh>
#include <dumux/porousmediumflow/mineralization/volumevariables.hh>
#include <dumux/porousmediumflow/mineralization/iofields.hh>

#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

    /*!
 * \ingroup TwoPICPModel
 * \brief Specifies a number properties of two-phase n-component models.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use molar or mass balances
 * \tparam setMoleFractionForFP whether to set mole fractions for first or second phase
 * \tparam formulation the formulation.
 * \tparam nSecComp the number of secondary components to be considered, which do not have their own mass balance, but are calculated from equilibrium chemistry.
 */
template<int nComp, bool useMol, bool setMoleFractionForFP, TwoPFormulation formulation, int nSecComp, int repCompEqIdx = nComp>
struct TwoPICPModelTraits : public TwoPNCModelTraits <nComp, useMol, setMoleFractionForFP, formulation, repCompEqIdx>
{
    //specific stuff
    static constexpr int numSecComponents() { return nSecComp; }
    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int numPhases() { return 2; }
};


/*!
 * \ingroup TwoPICPModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class CO2Tab, class Chem>
struct TwoPICPVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    //specific stuff
    using CO2Tables = CO2Tab;
    using Chemistry = Chem;
};


namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct TwoPICP { using InheritsFrom = std::tuple<TwoPNC>; };
struct TwoPICPNI { using InheritsFrom = std::tuple<TwoPICP>; };
} // end namespace TTag
NEW_PROP_TAG(CO2Tables); //!< The CO2 Tables that are used
// template<class TypeTag, class MyTypeTag>
// struct CO2Tables { using type = UndefinedProperty; }; //!< The CO2 Tables that are used

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 2picp model
//////////////////////////////////////////////////////////////////
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::TwoPICP>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! use the mineralization volume variables together with the 2picp vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPICP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using Chemistry = GetPropType<TypeTag, Properties::Chemistry>;

    using Traits = TwoPICPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, CO2Tables, Chemistry>;
    using NonMinVolVars = TwoPICPVolumeVariables<Traits>;
public:
    using type = MineralizationVolumeVariables<Traits, NonMinVolVars>;
};

//! Set the base model traits
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::TwoPICP>
{
private:
    //! we use the number of components specified by the fluid system here
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numPhases == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
public:
    using type = TwoPICPModelTraits<FluidSystem::numComponents,
                                                      GET_PROP_VALUE(TypeTag, UseMoles),
                                                      GET_PROP_VALUE(TypeTag, SetMoleFractionsForFirstPhase),
                                                      GET_PROP_VALUE(TypeTag, Formulation),
                                                      FluidSystem::numSecComponents,
                                                      GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx)>;
};
// template<class TypeTag>
// struct ModelTraits<TypeTag, TTag::TwoPNC> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPICP>
{
private:
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using NonMineralizationTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
//     using NonMineralizationTraits = TwoPICPModelTraits<FluidSystem::numComponents,
//                                                       GET_PROP_VALUE(TypeTag, UseMoles),
//                                                       GET_PROP_VALUE(TypeTag, SetMoleFractionsForFirstPhase),
//                                                       GET_PROP_VALUE(TypeTag, Formulation),
//                                                       FluidSystem::numSecComponents,
//                                                       GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx)>;
public:
    using type = MineralizationModelTraits<NonMineralizationTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;

};

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPICP> { using type = MineralizationIOFields<TwoPICPIOFields>; };

// use the mineralization local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPICP> { using type = MineralizationLocalResidual<TypeTag>; };

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::TwoPICP> { static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::numComponents; }; //!< Per default, no component mass balance is replaced

//! Default formulation is pw-Sn, overwrite if necessary
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPICP>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct SetMoleFractionsForFirstPhase<TypeTag, TTag::TwoPICP> { static constexpr bool value = true; };  //!< Set the primary variables mole fractions for the wetting or non-wetting phase
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPICP> { static constexpr bool value = true; };                         //!< Use mole fractions in the balance equations by default

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::TwoPICP> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! This model uses the compositional fluid state
template<class TypeTag>
struct FluidState<TypeTag, TTag::TwoPICP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalSecCompFluidState<Scalar, FluidSystem>;
};

//! This model uses a compositional solid state
template<class TypeTag>
struct SolidState<TypeTag, TTag::TwoPICP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = CompositionalSolidState<Scalar, SolidSystem>;
};

//////////////////////////////////////////////////////////////////
// Properties for the non-isothermal 2picp model
//////////////////////////////////////////////////////////////////

//! Set the non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPICPNI>
{
private:
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using NonMineralizationTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
//     using NonMineralizationTraits = TwoPICPModelTraits<FluidSystem::numComponents,
//                                                     GET_PROP_VALUE(TypeTag, UseMoles),
//                                                     GET_PROP_VALUE(TypeTag, SetMoleFractionsForFirstPhase),
//                                                     GET_PROP_VALUE(TypeTag, Formulation),
//                                                     FluidSystem::numSecComponents,
//                                                     GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx)>;
    using IsothermalTraits = MineralizationModelTraits<NonMineralizationTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;
public:
    // the mineralization traits, based on 2picp traits, are the isothermal traits
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPICPNI> { using type = EnergyIOFields<TwoPICPIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPICPNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};
} // end namespace Properties
} // end namespace Dumux

#endif
