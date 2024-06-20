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
 * \ingroup TwoPICPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase induced calcium carbonate precipitation model.
 */
#ifndef DUMUX_2PICP_VOLUME_VARIABLES_HH
#define DUMUX_2PICP_VOLUME_VARIABLES_HH

// #include <iostream>
// #include <vector>

// #include <dumux/common/math.hh>
// #include <dumux/common/properties.hh>
// #include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/2pnc/volumevariables.hh>
#include <dumux/material/fluidstates/compositionalsecondarycomponent.hh>
// #include <dumux/material/constraintsolvers/miscible2pnccomposition.hh>
// #include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

// #include "indices.hh" // for formulation

#include <dumux/material/binarycoefficients/brine_co2.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPICPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase induced calcium carbonate precipitation model.
 */
template <class Traits>
class TwoPICPVolumeVariables
: public TwoPNCVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPICPVolumeVariables<Traits> >
{
    using ParentType = TwoPNCVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPICPVolumeVariables<Traits> >;
//     using Implementation = GetPropType<TypeTag, Properties::VolumeVariables>;
//     using GridView = GetPropType<TypeTag, Properties::GridView>;
//     using Problem = GetPropType<TypeTag, Properties::Problem>;
//     using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
//     using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
//     using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using FS =  typename Traits::FluidSystem;
    using ModelTraits = typename Traits::ModelTraits;

    using CO2Tables = typename Traits::CO2Tables;
    using Chemistry = typename Traits::Chemistry;
    using Brine_CO2 = Dumux::BinaryCoeff::Brine_CO2<Scalar, CO2Tables, true>;

    static constexpr int numFluidComps = ParentType::numFluidComponents();

    enum
    {
//         numPhases = ModelTraits::numPhases(),
        numMajorComponents = ModelTraits::numPhases(),
        numComponents = ModelTraits::numFluidComponents(),
        numSecComponents = ModelTraits::numSecComponents(),

        // phase indices
        wPhaseIdx = FS::wPhaseIdx,
        nPhaseIdx = FS::nPhaseIdx,
        phase0Idx = wPhaseIdx,
        phase1Idx = nPhaseIdx,

        // component indices
        wCompIdx = FS::wCompIdx,
        nCompIdx = FS::nCompIdx,
        comp0Idx = wCompIdx,
        comp1Idx = nCompIdx,
        NaIdx = FS::NaIdx,
        CaIdx = FS::CaIdx,
        ClIdx = FS::ClIdx,
        CO2Idx = FS::CO2Idx,

        // phase presence enums
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases,
        nPhaseOnly = secondPhaseOnly,
        wPhaseOnly = firstPhaseOnly,

        // primary variable indices
        pressureIdx = ModelTraits::Indices::pressureIdx,
        switchIdx = ModelTraits::Indices::switchIdx
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();
    static constexpr bool setFirstPhaseMoleFractions = ModelTraits::setMoleFractionsForFirstPhase();

    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FS>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FS>;

public:
    //! return number of secondary components considered by the model
    static constexpr int numSecFluidComponents() { return Traits::ModelTraits::numSecComponents(); }
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! return whether moles or masses are balanced
    static constexpr bool useMoles() { return Traits::ModelTraits::useMoles(); }
    //! return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }
    // check for permissive specifications
    static_assert(useMoles(), "use moles has to be set true in the 2pnc model");
    static_assert(ModelTraits::numPhases() == 2, "NumPhases set in the model is not two!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        //we need porosity and permeability in completeFluidState for Leverett scaling of capillary pressure!
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting-phase saturation as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));

            mobility_[phaseIdx] = kr / fluidState_.viscosity(phaseIdx);

            int compIIdx = phaseIdx;
            for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                // binary diffusion coefficients
                if(compIIdx!= compJIdx)
                {
                    setDiffusionCoefficient_(phaseIdx, compJIdx,
                                             FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                     paramCache,
                                                                                     phaseIdx,
                                                                                     compIIdx,
                                                                                     compJIdx));
                }
            }
        }

        // calculate the remaining quantities
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
    }

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // set the saturations
        if (phasePresence == secondPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 0.0);
            fluidState.setSaturation(phase1Idx, 1.0);
        }
        else if (phasePresence == firstPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 1.0);
            fluidState.setSaturation(phase1Idx, 0.0);
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == TwoPFormulation::p0s1)
            {
                fluidState.setSaturation(phase1Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase0Idx, 1 - priVars[switchIdx]);
            }
            else
            {
                fluidState.setSaturation(phase0Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase1Idx, 1 - priVars[switchIdx]);
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        // set pressures of the fluid phases
        refPC_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
        refPorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity");
        refPermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability");
        pc_ = refPC_ * sqrt((refPermeability_*porosity_)/
                               (permeability_*refPorosity_));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase1Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] + pc_
                                                                       : priVars[pressureIdx] - pc_);
        }
        else
        {
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase0Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] - pc_
                                                                       : priVars[pressureIdx] + pc_);
        }

        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        Dune::FieldVector<Scalar, numComponents + numSecComponents> moleFrac(0.0);

        if (phasePresence == bothPhases)
        {
            // both phases are present, phase composition results from
            // the first <-> second phase equilibrium. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the MiscibleMultiPhaseComposition constraint solver

            const int knownPhaseIdx = setFirstPhaseMoleFractions ? phase0Idx : phase1Idx;

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(knownPhaseIdx, compIdx, priVars[compIdx]);
//                 std::cout <<"both phases, priVars["<< FluidSystem::componentName(compIdx)<<"]"<<priVars[compIdx]<< std::endl;
            }
            //TODO: Only CO2 and water are present in the non-wetting phase, which is why we do not use the constraintsolver
            MiscibleMultiPhaseComposition::solve(fluidState,
                                                 paramCache,
                                                 knownPhaseIdx);

            // set the fluid state for secondary components to zero for now. will be calculated later in the chemistry
            for (int compIdx=numComponents; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFractionSecComp(wPhaseIdx, compIdx, 0);
                fluidState.setMoleFractionSecComp(nPhaseIdx, compIdx, 0);
            }
        }
        //we don't really care about the unrealistic nPhaseOnly-case, most comonents are only in the water phase
        else if (phasePresence == nPhaseOnly)
        {
            moleFrac[comp0Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp0Idx];

            for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp1Idx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setMoleFraction(phase1Idx, compIdx, moleFrac[compIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase1Idx);

            // set the fluid state for secondary components (only in water phase) to zero in the nPhaseOnly-case
            for (int compIdx=numComponents; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFractionSecComp(wPhaseIdx, compIdx, 0);
                fluidState.setMoleFractionSecComp(nPhaseIdx, compIdx, 0);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            // extract _mass_ fractions in the non-wetting phase

            moleFrac[comp1Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp1Idx];
            for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
//                 std::cout <<"w phase only, priVars["<< FluidSystem::componentName(compIdx)<<"]"<<priVars[compIdx]<< std::endl;

                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp0Idx] = 1 - sumMoleFracOtherComponents;

            // convert mass to mole fractions and set the fluid state
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setMoleFraction(phase0Idx, compIdx, moleFrac[compIdx]);
            // set the fluid state for secondary components to zero for now. will be calculated later in the chemistry
            for (int compIdx=numComponents; compIdx<numComponents + numSecComponents; ++compIdx)
            {
                moleFrac[compIdx] = 0;
            }

//             // calculate the composition of the remaining phases (as
//             // well as the densities of all phases). this is the job
//             // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase0Idx);

            //TODO:later use Brine_CO2 again
            Scalar XwSalinity= 0.0;
            for (int compIdx = NaIdx; compIdx<= CaIdx ; compIdx++)  //salinity = XlNa + XlCl + XlCa
            {
                if(fluidState.massFraction(wPhaseIdx, compIdx)>0)
                {
                    XwSalinity+= fluidState.massFraction(wPhaseIdx, compIdx);
                }
            }
                        Scalar xnH2O;
                        Scalar xwCO2;

             Brine_CO2::calculateMoleFractions(fluidState.temperature(),
                                        fluidState.pressure(nPhaseIdx),
                                        XwSalinity,
                                        /*knownPhaseIdx=*/-1,
                                        xwCO2,
                                        xnH2O);
            // normalize the phase compositions
            xwCO2 = std::max(0.0, std::min(1.0, xwCO2));
            xnH2O = std::max(0.0, std::min(1.0, xnH2O));

            // set the fluid state
            for (int compIdx=0; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
                fluidState.setMoleFraction(nPhaseIdx, compIdx, 0);
            }
            //only CO2 and water are present in the non-wetting Phase
            fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnH2O);
            fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1-xnH2O);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
            fluidState.setViscosity(phaseIdx, mu);
            fluidState.setEnthalpy(phaseIdx, h);
        }

        //calculate the actual equilibrium aqueous phase composition including secondary components
//         if (phasePresence == bothPhases)
//         {
            Chemistry chemistry;
            chemistry.calculateEquilibriumChemistry(fluidState, phasePresence, moleFrac);
//         }
            for (int compIdx=numComponents; compIdx<numComponents + numSecComponents; ++compIdx)
            {
                fluidState.setMoleFractionSecComp(phase0Idx, compIdx, moleFrac[compIdx]);
                fluidState.setMoleFractionSecComp(phase1Idx, compIdx, 0);
            }
            //all inorganic carbon in the gas phase is CO2:
            fluidState.setMoleFractionSecComp(phase1Idx, CO2Idx, fluidState.moleFraction(phase1Idx, nCompIdx));
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.density(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the kinematic viscosity of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.viscosity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.molarDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(FluidSystem::nPhaseIdx) - fluidState_.pressure(FluidSystem::wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }


    /*!
     * \brief Returns the diffusion coefficient
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient called for phaseIdx = compIdx");
    }

    /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     */
     Scalar molarity(int phaseIdx, int compIdx) const // [moles/m^3]
    { return fluidState_.molarity(phaseIdx, compIdx);}

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     { return fluidState_.massFraction(phaseIdx, compIdx); }

     /*!
      * \brief Returns the mole fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     { return fluidState_.moleFraction(phaseIdx, compIdx); }

//      int numSecFluidComponents()
//      { return numSecFluidComps;}

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient for phaseIdx = compIdx doesn't exist");
    }

    Scalar refPC_;                     //!< The reference capillary pressure
    Scalar refPorosity_;               //!< The reference effective porosity within the control volume
    PermeabilityType refPermeability_; //!> The reference effective permeability within the control
    Scalar pc_;                        //!< The capillary pressure
    Scalar porosity_;                  //!< Effective porosity within the control volume
    PermeabilityType permeability_;    //!> Effective permeability within the control volume
    Scalar mobility_[ModelTraits::numPhases()]; //!< Effective mobility within the control volume
    std::array<std::array<Scalar, numComponents-1>, ModelTraits::numPhases()> diffCoefficient_;
};

} // end namespace Dumux

#endif
