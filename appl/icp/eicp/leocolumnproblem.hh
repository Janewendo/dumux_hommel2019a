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

#ifndef DUMUX_LEO_COLUMN_PROBLEM_HH
#define DUMUX_LEO_COLUMN_PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/leomin.hh>
#include <dumux/material/solidsystems/leominsolids.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2picp/model.hh>

#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/chemistry/biogeochemistry/leocarbonicacid.hh>

#include <appl/icp/icpspatialparams.hh>
#include <appl/icp/co2tableslaboratoryhightemp.hh>

#include "dumux/linear/seqsolverbackend.hh"

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class LEOColumnProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
#if NONISOTHERMAL
struct LEOColumnTypeTag { using InheritsFrom = std::tuple<TwoPICPNI>; };
struct LEOColumnCCTpfaTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, CCTpfaModel>; };
struct LEOColumnBoxTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, BoxModel>; };
#else  LEO
struct LEOColumnTypeTag { using InheritsFrom = std::tuple<TwoPICP>; };
struct LEOColumnCCTpfaTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, CCTpfaModel>; };
struct LEOColumnBoxTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, BoxModel>; };
#endif
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::YaspGrid<1>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LEOColumnTypeTag> { using type = LEOColumnProblem<TypeTag>; };

//Set the CO2 tables used.
SET_TYPE_PROP(LEOColumnTypeTag, CO2Tables, Dumux::ICP::CO2Tables);

// set the fluidSystem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LEOColumnTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::LeoMinFluid<Scalar, CO2Tables, H2OTabulated>;
};

// set the solidSystem
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::LEOColumnTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::LeoMinSolidPhase<Scalar>;
};


//Set the problem chemistry
template<class TypeTag>
struct Chemistry<TypeTag, TTag::LEOColumnTypeTag>
{
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = Dumux::LeoMinCarbonicAcid<TypeTag, CO2Tables, ModelTraits>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LEOColumnTypeTag> { using type = ICPSpatialParams<TypeTag>; };

template<class TypeTag>
struct Formulation<TypeTag, TTag::LEOColumnTypeTag>
{ static constexpr auto value = TwoPFormulation::p0s1; };

}

/*!
 * \ingroup TwoPNCSecCompMinModel
 * \ingroup ImplicitTestProblems
 * \brief Problem for enzyme-induced calcium carbonate precipitation
 *  */
template <class TypeTag>
class LEOColumnProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
        xwNaIdx = FluidSystem::NaIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwFe2Idx = FluidSystem::Fe2Idx,
        xwUreaIdx = FluidSystem::UreaIdx,
        xwTNHIdx = FluidSystem::TNHIdx,
        xwUreaseIdx = FluidSystem::UreaseIdx,
        phiImmUreaseIdx = numComponents,
        phiCalciteIdx = numComponents +1,
        phiFerrohydriteIdx = numComponents +2,
		
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

        //Indices of the components
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
        NaIdx = FluidSystem::NaIdx,
        ClIdx = FluidSystem::ClIdx,
        CaIdx = FluidSystem::CaIdx,
        Fe2Idx = FluidSystem::Fe2Idx,
        UreaIdx = FluidSystem::UreaIdx,
        TNHIdx = FluidSystem::TNHIdx,
        UreaseIdx = FluidSystem::UreaseIdx,

        NH4Idx = FluidSystem::NH4Idx,
        CO3Idx = FluidSystem::CO3Idx,
        HCO3Idx = FluidSystem::HCO3Idx,
        CO2Idx = FluidSystem::CO2Idx,
        HIdx = FluidSystem::HIdx,
        OHIdx = FluidSystem::OHIdx,

        //Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase State
        wPhaseOnly = Indices::firstPhaseOnly,
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Chemistry = GetPropType<TypeTag, Properties::Chemistry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

public:
    LEOColumnProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);

        name_  = getParam<std::string>("Problem.Name");

        //initial values
        densityW_ = getParam<Scalar>("Initial.initDensityW");
        initPressure_ = getParam<Scalar>("Initial.initPressure");

        initxwTC_ = getParam<Scalar>("Initial.initxwTC");
        initxwNa_ = getParam<Scalar>("Initial.initxwNa");
        initxwCl_ = getParam<Scalar>("Initial.initxwCl");
        initxwCa_ = getParam<Scalar>("Initial.initxwCa");
        initxwFe2_ = getParam<Scalar>("Initial.initxwFe2");
        initxwUrea_ = getParam<Scalar>("Initial.initxwUrea");
        initxwTNH_ = getParam<Scalar>("Initial.initxwTNH");
        initxwEnzymeSource_ = getParam<Scalar>("Initial.initxwEnzymeSource");
        initCalcite_ = getParam<Scalar>("Initial.initCalcite");
        initFerrohydrite_ = getParam<Scalar>("Initial.initFerrohydrite");
        initImmUrease_ = getParam<Scalar>("Initial.initImmUrease");
        initTemperature_ = getParam<Scalar>("Initial.initTemperature");

        xwNaCorr_ = getParam<Scalar>("Initial.xwNaCorr");
        xwClCorr_ = getParam<Scalar>("Initial.xwClCorr");

        //injection values
        injQ_ = getParam<Scalar>("Injection.injVolumeflux");

        injTC_ = getParam<Scalar>("Injection.injTC");
        injNa_ = getParam<Scalar>("Injection.injNa");
        injCa_ = getParam<Scalar>("Injection.injCa");
        injFe2_ = getParam<Scalar>("Injection.injFe2");
        injUrea_ = getParam<Scalar>("Injection.injUrea");
        injTNH_ = getParam<Scalar>("Injection.injTNH");
        injEnzymeSource_= getParam<Scalar>("Injection.injEnzymeSource");

        injNaCorr_ = getParam<Scalar>("Injection.injNaCorr");
        injTemperature_ = getParam<Scalar>("Injection.injTemperature");
        injPressure_ = getParam<Scalar>("Injection.injPressure");

        numInjections_ = getParam<int>("Injection.numInjections");
        injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(gridGeometry->gridView().size(codim));
        calcium_.resize(gridGeometry->gridView().size(codim));
        urea_.resize(gridGeometry->gridView().size(codim));

        std::ifstream injectionData;
        std::string row;
        injectionData.open( injectionParameters_); // open the Injection data file
        if (not injectionData.is_open())
        {
            std::cerr << "\n\t -> Could not open file '"
                    << injectionParameters_
                    << "'. <- \n\n\n\n";
            exit(1) ;
        }
        int tempType = 0;

        // print file to make sure it is the right file
        std::cout << "Read file: " << injectionParameters_ << " ..." << std::endl << std::endl;
        while(!injectionData.eof())
        {
            getline(injectionData, row);
            std::cout << row << std::endl;
        }
        injectionData.close();

        //read data from file
        injectionData.open(injectionParameters_);

        while(!injectionData.eof())
        {
            getline(injectionData, row);

            if(row == "InjectionTypes")
            {
                getline(injectionData, row);
                while(row != "#")
                {
                    if (row != "#")
                        {
                        std::istringstream ist(row);
                        ist >> tempType;
                        injType_.push_back(tempType);
//                         std::cout << "size of injType: "<<injType_.size() << std::endl;
                        }
                    getline(injectionData, row);
                }
            }
        }

        injectionData.close();

//      check the injection data against the number of injections specified in the parameter file
        if (injType_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection types specified in the injection data file do not match!"
                    <<"\n numInjections from parameter file = "<<numInjections_
                    <<"\n numInjTypes from injection data file = "<<injType_.size()
                    <<"\n Abort!\n";
            exit(1) ;
        }

#if NONISOTHERMAL
        FluidSystem::init(/*startTemp=*/295.15, /*endTemp=*/445.15, /*tempSteps=*/151,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
#else
        FluidSystem::init(/*startTemp=*/initTemperature_ -5.0, /*endTemp=*/initTemperature_ +5.0, /*tempSteps=*/5,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
#endif
    }

    void setTime( Scalar time )
    {
        time_ = time;
    }

    void setTimeStepSize( Scalar timeStepSize )
    {
        timeStepSize_ = timeStepSize;
    }

    void setEpisodeIdx( int epiIdx )
    {
        episodeIdx_ = epiIdx;
    }

    int injectionType(int episodeIdx)
    {
        return injType_[episodeIdx];
    }

   /*!
    * \name Problem parameters
    */


   /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
//    const char *name() const
    const std::string name() const
    { return name_; }

 #if !NONISOTHERMAL
   /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 25 degrees Celsius.
    */
    Scalar temperature() const
    {
        return initTemperature_; //
    };
#endif

    // \}

   /*!
    * \name Boundary conditions
    */
    // \{

    /*!
    * \brief Specifies which kind of boundary condition should be
    *        used for which equation on a given boundary segment.
    */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        Scalar zmax = this->gridGeometry().bBoxMax()[dim - 1];
        bcTypes.setAllNeumann();
        if(globalPos[dim - 1] > zmax - eps_)
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

   /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        boundary segment.
    */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

   /*!
    * \brief Evaluate the initial value for a control volume.
    *
    * \param globalPos The global position
    *
    * For this method, the \a values parameter stores primary
    * variables.
    */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return initial_(globalPos);
    }

   /*!
    * \brief Evaluate the boundary conditions for a Neumann
    *        boundary segment.
    *
    * For this method, the \a values parameter stores the mass flux
    * in normal direction of each component. Negative values mean
    * influx.
    *
    * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
    */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        Scalar waterFlux = injQ_/4.640675e-4; //[m/s]
        //4.640675e-4 m^2 = 0.957" diameter column cross-sectional area.

        int injProcess = injType_[episodeIdx_];

//         negative values for injection
        if(globalPos[1]<= eps_)
        {
           //basic rinse injection (injProcess == -1 )
           values[conti0EqIdx + wCompIdx] = -waterFlux * 996/FluidSystem::molarMass(wCompIdx);
           values[conti0EqIdx + nCompIdx] = -waterFlux * injTC_*996 /FluidSystem::molarMass(nCompIdx);
           values[conti0EqIdx + xwCaIdx] = 0;
           values[conti0EqIdx + xwUreaIdx] = 0;
           values[conti0EqIdx + xwUreaseIdx] = 0;
           values[conti0EqIdx + xwTNHIdx] = 0;
           values[conti0EqIdx + phiCalciteIdx] = 0;
           values[conti0EqIdx + xwFe2Idx] = 0;
           values[conti0EqIdx + phiFerrohydriteIdx] = 0;
           values[conti0EqIdx + phiImmUreaseIdx] = 0;
           values[conti0EqIdx + xwNaIdx] = -waterFlux * (injNa_+injNaCorr_) /FluidSystem::molarMass(NaIdx);
           values[conti0EqIdx + xwClIdx] = -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
#if NONISOTHERMAL
            values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
                                    injTemperature_, injPressure_,
                                   (injNa_
                                   +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
#endif

            if (injProcess == -1)       // rinse,
            {
                //only NH4Cl
                values[conti0EqIdx + xwTNHIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
                values[conti0EqIdx + xwClIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
#if NONISOTHERMAL
                values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
                                    injTemperature_, injPressure_,
                                    (injTNH_/FluidSystem::molarMass(TNHIdx)*FluidSystem::molarMass(ClIdx)
                                    +injNa_
                                    +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
#endif
            }

           else if (injProcess <-89)       // no injection
           {
            values = 0.0; //mol/m²/s
           }

           else if (injProcess == 1)              //ca-rich injection: ca and urea injected Na (pH) and Cl are different(CaCl2)
           {
               values[conti0EqIdx + wCompIdx] = - waterFlux * 0.8716 * densityW_ /FluidSystem::molarMass(wCompIdx);       //TODO 0.8716 check factor!!!
               values[conti0EqIdx + nCompIdx] = - waterFlux * injTC_ * densityW_ /FluidSystem::molarMass(nCompIdx);
               values[conti0EqIdx + xwCaIdx] = - waterFlux * injCa_/FluidSystem::molarMass(CaIdx);
               values[conti0EqIdx + xwFe2Idx] = - waterFlux * injFe2_/FluidSystem::molarMass(Fe2Idx);			   
               values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
               values[conti0EqIdx + xwClIdx] += - waterFlux * 2 * injCa_/FluidSystem::molarMass(CaIdx);

#if NONISOTHERMAL
            values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
                                    injTemperature_, injPressure_,
                                   (injCa_
                                   +2*injCa_ /FluidSystem::molarMass(CaIdx)*FluidSystem::molarMass(ClIdx)
                                   +injNa_
                                   +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
#endif
       }
              else if(injProcess == 2)              //injection of Urease and Jack Bean Meal
       {
                values[conti0EqIdx + xwUreaseIdx] = -waterFlux * injEnzymeSource_ /FluidSystem::molarMass(UreaseIdx);
       }
               else
               {
                      DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << injProcess);
               }
        }
        else
        {
               values = 0.0; //mol/m²/s
        }
        return values;
    }

   /*!
    * \name Volume terms
    */
    // \{

   /*!
    * \brief Evaluate the source term for all phases within a given
    *        sub-control-volume.
    *
    * This is the method for the case where the source term is
    * potentially solution dependent and requires some quantities that
    * are specific to the fully-implicit method.
    *
    * \param values The source and sink values for the conservation equations in units of
    *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
    * \param element The finite element
    * \param fvGeometry The finite-volume geometry
    * \param elemVolVars All volume variables for the element
    * \param scv The subcontrolvolume
    *
    * For this method, the \a values parameter stores the conserved quantity rate
    * generated or annihilate per volume unit. Positive values mean
    * that the conserved quantity is created, negative ones mean that it vanishes.
    * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
    */
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        Chemistry chemistry;
        source = chemistry.reactionSource(elemVolVars[scv],
                        timeStepSize_);
        return source;
    }

   /*!
    * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
    */

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }
    const std::vector<Scalar>& getCalcium()
    {
        return calcium_;
    }
    const std::vector<Scalar>& getUrea()
    {
        return urea_;
    }



    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = volVars.permeability();
                calcium_[dofIdxGlobal] = volVars.moleFraction(0,CaIdx)* volVars.molarDensity(0) * FluidSystem::molarMass(CaIdx);
                urea_[dofIdxGlobal] = volVars.moleFraction(0,UreaIdx)* volVars.molarDensity(0) * FluidSystem::molarMass(UreaIdx);
            }
        }
    }

    void setGridVariables(std::shared_ptr<GridVariables> gv)
    { gridVariables_ = gv; }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        priVars[pressureIdx] = initPressure_ ;
        priVars[switchIdx] = initxwTC_;
        priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
        priVars[xwCaIdx] = initxwCa_;
        priVars[xwUreaIdx] = initxwUrea_;
        priVars[xwTNHIdx] = initxwTNH_;
        priVars[xwFe2Idx] = initxwFe2_;
        priVars[xwUreaseIdx]= initxwEnzymeSource_;
        priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]
        priVars[phiFerrohydriteIdx] = initFerrohydrite_; // [m^3/m^3]
        priVars[phiImmUreaseIdx] = initImmUrease_; // [m^3/m^3]

#if NONISOTHERMAL
        priVars[temperatureIdx] = initTemperature_;
#endif
        return priVars;
    }
    /*!
        * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
        *
        * \param XwNaCl the XwNaCl [kg NaCl / kg solution]
        */
    static Scalar massTomoleFrac_(Scalar XwNaCl)
    {
        const Scalar Mw = FluidSystem::molarMass(wCompIdx);  // 18.015e-3; /* molecular weight of water [kg/mol] */
        const Scalar Ms = FluidSystem::molarMass(NaIdx) + FluidSystem::molarMass(ClIdx); // 58.44e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = XwNaCl;
        /* XwNaCl: conversion from mass fraction to mol fraction */
        const Scalar xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return xwNaCl;
    }

    static constexpr Scalar eps_ = 1e-6;

    Scalar initPressure_;
    Scalar densityW_;

    Scalar initxwTC_;
    Scalar initxwNa_;
    Scalar initxwCl_;
    Scalar initxwCa_;
    Scalar initxwUrea_;
    Scalar initxwTNH_;
    Scalar initxwFe2_;
    Scalar initxwEnzymeSource_;
    Scalar xwNaCorr_;
    Scalar xwClCorr_;

    Scalar initCalcite_;
    Scalar initFerrohydrite_;
    Scalar initImmUrease_;
    Scalar initImmJBM_;
    Scalar initPorosity_;
    Scalar initTemperature_;

    Scalar injQ_;

    Scalar injTC_;
    Scalar injNa_;
    Scalar injCa_;
    Scalar injFe2_;
    Scalar injUrea_;
    Scalar injTNH_;
    Scalar injNaCorr_;
    Scalar injEnzymeSource_;
    Scalar injTemperature_;
    Scalar injPressure_;

    Scalar ureaseInJBM_;

    int numInjections_;
    std::string injectionParameters_;

    std::vector<int> injType_;
    std::string name_;

    std::vector<Scalar> permeability_;
    std::vector<Scalar> calcium_;
    std::vector<Scalar> urea_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
    std::shared_ptr<GridVariables> gridVariables_;
};
} //end namespace

#endif

