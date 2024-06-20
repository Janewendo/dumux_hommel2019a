
#ifndef LEO_CARBONIC_ACID_HH_
#define LEO_CARBONIC_ACID_HH_

#include <dumux/common/exceptions.hh>
#include <dumux/material/fluidsystems/leomin.hh>

#include <cmath>
#include <iostream>
#include <dumux/common/math.hh>

namespace Dumux
{
/*!
 * \brief The equilibrium chemistry is calculated in this class. The function calculateEquilbriumChemistry is used to
 * control the Newton Solver "newton1D". The chemical functions and derivations are implemented in the private part of
 * class.
 */
template <class TypeTag, class CO2Tables, class ModelTraits>
class LeoMinCarbonicAcid
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Sources = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using ThisType = LeoMinCarbonicAcid<TypeTag, CO2Tables, ModelTraits>;
    using H2O = Components::H2O<Scalar>;

    enum
    {
        // phase presence enums
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases,
        nPhaseOnly = secondPhaseOnly,
        wPhaseOnly = firstPhaseOnly,
    };

public:

    LeoMinCarbonicAcid()
{

    pKaFactor_ = getParam<Scalar>("Geochem.pKaFactor");
    //what is the enzyme source, heat-killed cells or jack bean?
    useHeatKilledCells_ = getParam<bool>("Problem.UseHeatKilledCells");
    useJackBeans_ = getParam<bool>("Problem.UseJackBeans");

    // calcite parameters
    ac_        = getParam<Scalar>("CalciteCoefficients.ac");
    kdiss1_    = getParam<Scalar>("CalciteCoefficients.kdiss1");
    kdiss2_    = getParam<Scalar>("CalciteCoefficients.kdiss2");
    kprec_     = getParam<Scalar>("CalciteCoefficients.kprec");
    ndiss_     = getParam<Scalar>("CalciteCoefficients.ndiss");
    nprec_     = getParam<Scalar>("CalciteCoefficients.nprec");
    Asw0_      = getParam<Scalar>("CalciteCoefficients.Asw0");

    // ferrohydrite parameters
    fac_        = getParam<Scalar>("FerrohydriteCoefficients.fac");
    fkdiss1_    = getParam<Scalar>("FerrohydriteCoefficients.fkdiss1");
    fkdiss2_    = getParam<Scalar>("FerrohydriteCoefficients.fkdiss2");
    fkprec_     = getParam<Scalar>("FerrohydriteCoefficients.fkprec");
    fndiss_     = getParam<Scalar>("FerrohydriteCoefficients.fndiss");
    fnprec_     = getParam<Scalar>("FerrohydriteCoefficients.fnprec");
    fAsw0_      = getParam<Scalar>("FerrohydriteCoefficients.fAsw0");
	
    // thermal ureolysis parameters
    cu_        = getParam<Scalar>("UreolysisCoefficients.cu");
    cuT_       = getParam<Scalar>("UreolysisCoefficients.cuT");

    // urease parameters
    if(useHeatKilledCells_)
    {
    //     ureaseInEnzymeSource_ = getParam<Scalar>("UreaseCoefficients.ureaseInEnzymeSource");
        cia_       = getParam<Scalar>("UreaseCoefficientsHKC.cia");
        ciaT_      = getParam<Scalar>("UreaseCoefficientsHKC.ciaT");
        cureaseT_  = getParam<Scalar>("UreaseCoefficientsHKC.cureaseT");
        kurease_   = getParam<Scalar>("UreaseCoefficientsHKC.kurease");
        ciaPrec_   = getParam<Scalar>("UreaseCoefficientsHKC.ciaPrec");

        //attachment and detachment parameters
        ka_urease_ = getParam<Scalar>("UreaseCoefficientsHKC.ka_urease");
        kd_urease_ = getParam<Scalar>("UreaseCoefficientsHKC.kd_urease");
    }
    if(useJackBeans_)
    {
    //     ureaseInEnzymeSource_ = getParam<Scalar>("UreaseCoefficients.ureaseInEnzymeSource");
        cia_       = getParam<Scalar>("UreaseCoefficientsJB.cia");
        ciaT_      = getParam<Scalar>("UreaseCoefficientsJB.ciaT");
        cureaseT_  = getParam<Scalar>("UreaseCoefficientsJB.cureaseT");
        kurease_   = getParam<Scalar>("UreaseCoefficientsJB.kurease");
        ciaPrec_   = getParam<Scalar>("UreaseCoefficientsJB.ciaPrec");

        //attachment and detachment parameters
        ka_urease_ = getParam<Scalar>("UreaseCoefficientsJB.ka_urease");
        kd_urease_ = getParam<Scalar>("UreaseCoefficientsJB.kd_urease");
    }
}

    static const int wPhaseIdx    = FluidSystem::wPhaseIdx;
    static const int nPhaseIdx    = FluidSystem::nPhaseIdx;

    static const int wCompIdx     = FluidSystem::wCompIdx;
    static const int nCompIdx     = FluidSystem::nCompIdx;

    static const int H2OIdx       = FluidSystem::H2OIdx;
    static const int CTotIdx      = FluidSystem::TCIdx;
    static const int CaIdx        = FluidSystem::CaIdx;
    static const int NaIdx        = FluidSystem::NaIdx;
    static const int ClIdx        = FluidSystem::ClIdx;
    static const int HIdx         = FluidSystem::HIdx;
    static const int OHIdx        = FluidSystem::OHIdx;
    static const int CO2Idx       = FluidSystem::CO2Idx;
    static const int HCO3Idx      = FluidSystem::HCO3Idx;
    static const int CO3Idx       = FluidSystem::CO3Idx;
    static const int Fe2Idx       = FluidSystem::Fe2Idx;
	
    static const int UreaIdx      = FluidSystem::UreaIdx;
    static const int UreaseIdx    = FluidSystem::UreaseIdx;

    static const int TNHIdx       = FluidSystem::TNHIdx;
    static const int NH4Idx       = FluidSystem::NH4Idx;

    static const int numComponents      = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numMajorComponents;
    static const int numSecComponents   = FluidSystem::numSecComponents;
    static const int numTotComponents   = numComponents + numSecComponents;
    static const int numPhases          = FluidSystem::numPhases;

    static const int cPhaseIdx          = SolidSystem::CalciteIdx;
    static const int uPhaseIdx          = SolidSystem::JbmeIdx;
    static const int fPhaseIdx          = SolidSystem::FerrohydriteIdx;	
    static const int numSolidComponents = SolidSystem::numComponents;
    static const int numInertComponents = SolidSystem::numInertComponents;

    static const int phiCalciteIdx      = numComponents + cPhaseIdx;
    static const int phiImmUreaseIdx    = numComponents + uPhaseIdx;
    static const int phiFerrohydriteIdx    = numComponents + fPhaseIdx;

    typedef Dune::FieldVector<Scalar, 4> Vector;   // Ionic Strength with NH4/totalnh
    typedef Dune::FieldVector<Scalar, 2> SolVector;
    typedef Dune::FieldVector<Scalar, numTotComponents> CompVector;

    typedef CompositionalSecCompFluidState<Scalar, FluidSystem> FluidState;

    template <class FluidState>
    void calculateEquilibriumChemistry(const FluidState &fluidState, int phaseState, CompVector &variable)
    {

        gammaCO2_ = 1.0;
        h2o_ = 55.508; //molH2O/kgH2O
        pressure_ = fluidState.pressure(wPhaseIdx);
        temperature_ = fluidState.temperature();

        Scalar moleFracSalinity = variable[NaIdx] + variable[ClIdx] + variable[CaIdx];

        for (int i = 0; i < numComponents + numSecComponents; ++i)
            {
                if(variable[i]<0)
                    variable[i]=0;

                if(std::isnan(variable[i]))
                {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid component mole fraction " << variable); break;
                }
            }

        if (phaseState == bothPhases) //both Phases: solve an open system with co2 concentration constant
        {
            salinity_ = moleFracToMolality(moleFracSalinity, moleFracSalinity, variable[nCompIdx]);

            co2_ =  (salinity_ + h2o_)/(1/variable[nCompIdx] - 1);
            ca_ = moleFracToMolality(variable[CaIdx], moleFracSalinity, variable[nCompIdx]);
            na_ = moleFracToMolality(variable[NaIdx], moleFracSalinity, variable[nCompIdx]);
            cl_ = moleFracToMolality(variable[ClIdx], moleFracSalinity, variable[nCompIdx]);
            fe2_ = moleFracToMolality(variable[Fe2Idx], moleFracSalinity, variable[nCompIdx]);
            totalnh_ = moleFracToMolality(variable[TNHIdx], moleFracSalinity, variable[nCompIdx]);

            Scalar m = na_ + ca_;
            Scalar Temp = fluidState.temperature();
            /* Millero et al. 2007: The dissociation of carbonic acid */
            /* in NaCl solutions as a function of concentration and temperature */

            /*for pK1*/
            Scalar a0 = 31.3616;Scalar  a1 = 0.86644;Scalar  a2 = -0.33611;Scalar  a3 = 0.05888;
            Scalar  b0 = -1422.25317; Scalar c0 = -4.84141;

//             Scalar A = a0*sqrt(m) + a1*m+ a2*pow(m,1.5) + a3*m*m ;
            Scalar A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
//             Scalar B = b0 * pow (m,0.5); Scalar C = c0*pow(m,0.5);
            Scalar B = b0 * sqrt(m);
            Scalar C = c0*sqrt(m);

            Scalar dpK1 = A +  B/Temp + C*log(Temp);
            Scalar pK1 = - 402.56788 + 11656.46/Temp + 72.173*log(Temp) - 0.161325*Temp + 7.5526E-5*Temp*Temp;

            pK1 = pK1 + dpK1;

            /*for pK2*/
            a0 = 36.88545; a1 = 1.66599; a2 = -0.68730; a3 = 0.12070;
            b0 = -1669.55918; c0 = -5.83555;

            A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
            B = b0 * sqrt(m);
            C = c0*sqrt(m);

            Scalar dpK2 = A +  B/Temp + C*log(Temp);
            Scalar  pK2 = -122.4994 + 5811.18/Temp + 20.5263*log(Temp) - 0.0120897*Temp;

            pK2 = pK2 + dpK2;

            /*Bell et al. 2008: Ammonia/ammonium dissociation coefficient*/
            /*in seawater: A significant numerical correction*/
            Scalar I_f = 0.5*(na_ + 4.*ca_ + cl_);                         /*ionic strength of salt solution: here, equivalent to m, neglecting other ions*/
//             Scalar pKa = 10.0423 - 0.0315536*(Temp-273.15) + 0.14737*I_f;
//             Olofsson 1975: Thermodynamic quantities for the dissociation of the ammonium Ion and for the ionization of aqueous ammonia over a wide temperature range (up to 570K!)
            Scalar pKa = 2533/Temp - 0.5936* log(Temp) + 4.127;
            //using the ionic strength correction by Bell et al. 2008
            pKa += 0.14737*I_f;
            pKa *= pKaFactor_;

            apparentk1_ = pow (10.,-pK1);
            apparentk2_ = pow (10.,-pK2);
            apparentka_ = pow (10.,-pKa);

            initH_ = 1e-5; //Initial guess
            Scalar activityH = initH_;

            //Anozies apparent constants
            k1_ = apparentk1_;
            k2_ = apparentk2_;
            ka_ = apparentka_;
            kw_ = constW(pressure_, temperature_);

            Scalar tolAbs = 1e-20;
            Scalar tolRel = 1e-15;
            int maxIter = 30;

            //Do the Newton iteration and calculate the components molalities and update the mass fraction array and
            //the concentration of the changed primary variables
            if(newton1D(activityH, &ThisType::H_CO2, tolAbs, tolRel, maxIter) == false) //Alex' Newton
//            if(newton1D(activityH, tolAbs, maxIter) == false) //Anozies Newton
            {
                initH_ = 1e-5;
                activityH = initH_;
                Scalar a0 = 0.0;
                Scalar b0 = 1e-1;
                Scalar tol = 1e-15;
                if(bisection1D(activityH, &ThisType::H_CO2, a0, b0, tol) == false) //Alex' bisection
//                if(bisection1D(tol) == false) //Anozies bisection
                {
                    DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
                }
            }
            H_CO2(activityH); //update component molalities

            //update mole fractions in the variable vector for the open system

            cTot_ = co2_ + hco3_ + co3_; //calculate the molality of cTot from the other c-components

            Scalar moleFracCTot = molalityToMoleFrac(cTot_, moleFracSalinity, variable[nCompIdx]);

//             Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracWater);
            Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracSalinity, variable[nCompIdx]);
            Scalar urease = moleFracToMolality(variable[UreaseIdx], moleFracSalinity, variable[nCompIdx]);

                  Scalar totalMolality = h2o_ + cTot_ + na_ + cl_ + ca_ + totalnh_ + urea + urease + fe2_;
                  variable[CTotIdx] = cTot_/totalMolality; //calculate the mole fraction of cTot in terms of mol CO2 / mol solution
                  variable[CO2Idx] = co2_/totalMolality;
                  variable[HCO3Idx] = hco3_/totalMolality;
                  variable[CO3Idx] = co3_/totalMolality;
                  variable[NH4Idx] = nh4_/totalMolality;
                  variable[OHIdx] = oh_/totalMolality;
                  variable[HIdx] = h_/totalMolality;

        }

        else if (phaseState == wPhaseOnly) //wPhaseOnly: solve a closed system with cTot concentration constant
        {
            salinity_ = moleFracToMolality(moleFracSalinity,  moleFracSalinity, variable[nCompIdx]);
            cTot_ =  moleFracToMolality(variable[nCompIdx], moleFracSalinity, variable[nCompIdx]);
            ca_ = moleFracToMolality(variable[CaIdx], moleFracSalinity, variable[nCompIdx]);
            na_ = moleFracToMolality(variable[NaIdx], moleFracSalinity, variable[nCompIdx]);
            cl_ = moleFracToMolality(variable[ClIdx],  moleFracSalinity, variable[nCompIdx]);
            fe2_ = moleFracToMolality(variable[Fe2Idx], moleFracSalinity, variable[nCompIdx]);         
		    totalnh_ = moleFracToMolality(variable[TNHIdx], moleFracSalinity, variable[nCompIdx]);

            Scalar m = na_ + ca_;
            Scalar Temp = fluidState.temperature();
            /* Millero et al. 2007: The dissociation of carbonic acid */
            /* in NaCl solutions as a function of concentration and temperature */

            /*for pK1*/
            Scalar a0 = 31.3616;Scalar  a1 = 0.86644;Scalar  a2 = -0.33611;Scalar  a3 = 0.05888;
            Scalar  b0 = -1422.25317; Scalar c0 = -4.84141;

            Scalar A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
            Scalar B = b0 * sqrt(m);
            Scalar C = c0*sqrt(m);

            Scalar dpK1 = A +  B/Temp + C*log(Temp);
            Scalar pK1 = - 402.56788 + 11656.46/Temp + 72.173*log(Temp) - 0.161325*Temp + 7.5526E-5*Temp*Temp;

            pK1 = pK1 + dpK1;

            /*for pK2*/
            a0 = 36.88545; a1 = 1.66599; a2 = -0.68730; a3 = 0.12070;
            b0 = -1669.55918; c0 = -5.83555;

//             A = a0*pow(m,0.5) + a1*m+ a2*pow(m,1.5) + a3*m*m;
            A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
//             B = b0 * pow (m,0.5); C = c0*pow(m,0.5);
            B = b0 * sqrt(m);
            C = c0*sqrt(m);

            Scalar dpK2 = A +  B/Temp + C*log(Temp);
            Scalar  pK2 = -122.4994 + 5811.18/Temp + 20.5263*log(Temp) - 0.0120897*Temp;

            pK2 = pK2 + dpK2;

            /*Bell et al. 2008: Ammonia/ammonium dissociation coefficient*/
            /*in seawater: A significant numerical correction*/
            Scalar I_f = 0.5*(na_ + 4.*ca_ + cl_);                         /*ionic strength of salt solution: here, equivalent to m, neglecting other ions*/
//             Scalar pKa = 10.0423 - 0.0315536*(Temp-273.15) + 0.14737*I_f;
//             Olofsson 1975: Thermodynamic quantities for the dissociation of the ammonium Ion and for the ionization of aqueous ammonia over a wide temperature range (up to 570K!)
            Scalar pKa = 2533/Temp - 0.5936* log(Temp) + 4.127;
            //using the ionic strength correction by Bell et al. 2008
            pKa += 0.14737*I_f;
            pKa *= pKaFactor_;

            apparentk1_ = pow (10.,-pK1);
            apparentk2_ = pow (10.,-pK2);
            apparentka_ = pow (10.,-pKa);

            //Anozies apparent constants
            k1_ = apparentk1_;
            k2_ = apparentk2_;
            ka_ = apparentka_;
            kw_ = constW(pressure_, temperature_);

            //Parameters for the newton solver
//            Scalar CTotLow = 1e-4;
//            Scalar CTotHigh = 1e-2;
            Scalar tolAbs = 1e-20; //1e-11;old
            Scalar tolRel = 1e-20; //1e-11;old
            int maxIter = 40;
            initH_ = 1e-7;
            Scalar activityH = initH_;
            if(newton1D(activityH, &ThisType::H_Ctot, tolAbs, tolRel, maxIter) == true) //Alex' Newton
//          if(newton1D(activityH, tolAbs, maxIter) == true)  //Anozies Newton
            {
                //update all component molalities
                H_Ctot(activityH);
            }
            else //solve with the bisection method and hco3 as primary variable
            {
                Scalar a0 = 0;
                Scalar b0 = 1e-3;//1e-4; old
//                initH_ = 1e-5;
                Scalar tol = 1e-12;//1e-8; old
                Scalar activityH = b0;
                if(bisection1D(activityH, &ThisType::H_Ctot, a0, b0, tol) == true) //Alex' Bisection
//                if(bisection1D(tol) == true) //Anozies Bisection
                {
                    H_Ctot(activityH);//   CTot_HCO3(activityHCO3); //update all component molalities
                }
                else
                {
                    DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
                }

            }

            Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracSalinity, variable[nCompIdx]);
            Scalar urease = moleFracToMolality(variable[UreaseIdx], moleFracSalinity, variable[nCompIdx]);

            Scalar totalMolality = h2o_ + cTot_ + na_ + cl_ + ca_ + totalnh_ + urea + urease + fe2_;

            // calculate the secondary component mole fractions
            variable[CO2Idx] = co2_/totalMolality;
            variable[HCO3Idx] = hco3_/totalMolality;
            variable[CO3Idx] = co3_/totalMolality;
            variable[NH4Idx] = nh4_/totalMolality;
            variable[OHIdx] = oh_/totalMolality;
            variable[HIdx] = h_/totalMolality;
			
            Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_ + 2*fe2_;
            Scalar fmolfrac = 2*variable[CaIdx] + variable[NaIdx] + variable[NH4Idx] + variable[HIdx] - variable[ClIdx] - variable[HCO3Idx] - 2*variable[CO3Idx] - variable[OHIdx] + 2*variable[Fe2Idx];

           for (int i = 0; i < numComponents + numSecComponents; ++i)
           {
              if(std::isnan(variable[i]))
              {
                  std::cout<<"wPhaseOnly, moleFrac of  "<< FluidSystem::componentName(i) << "  in chemistry is: "<< variable[i]<< " all moleFracs are: \n"  << variable<<"\n" <<std::endl;
              }
           }

        }

        else if (phaseState == nPhaseOnly) //no secondary components in the gas phase, except CTot = CO2!
        {
            variable[CO2Idx] = variable[CTotIdx]; //all cTot in the gas is CO2
            variable[HCO3Idx] = 0.0;
            variable[CO3Idx] = 0.0;
            variable[NH4Idx] = 0.0;
            variable[OHIdx] = 0.0;
            variable[HIdx] = 0.0;
        }

        else
        {
            DUNE_THROW(Dune::InvalidStateException, "Invalid phaseState" );
        }
    }

    //Return equlibrium constant for chemical equation:
    //H2CO3 <--> H + HCO3
    static Scalar const1(const Scalar pw, const Scalar T)
    {
        return 5.01187e-7;//return(pow(10,-6.3));
    }

    //Return equlibrium constant for chemical equation:
    //HCO3 <--> H + CO3
    static Scalar const2(const Scalar pw, const Scalar T)
    {
        return 5.01187e-11;//return(pow(10,-10.3));
    }

//    Return equlibrium constant for dissolution reaction:
//    CaCO3(s) <--> Ca + CO3
    static Scalar solubilityProductCaCO(const Scalar pw, const Scalar T)
    {
        return 4.8e-9;
    }

//    Return equlibrium constant for dissolution reaction:
//    Fe(OH)2(s) <--> Fe + 2OH
    static Scalar solubilityProductFeOH2(const Scalar pw, const Scalar T)
    {
        return(pow(10,-4.89));
    }

    //Return equlibrium constant for chemical equation:
    // H2O <--> H + OH
    static Scalar constW(const Scalar pw, const Scalar T)
    {
        return 1e-14;
    }
    //Return equlibrium constant for chemical equation:
    // NH4 <--> H + NH3
    /*static*/ Scalar consta(const Scalar pw, const Scalar T)
    {
        return 5.12861e-10;//return(pow(10,-9.29)); //pow(10,-9.25)
    }


    static Scalar massFracToMolality(const Scalar massFracX, const Scalar molarMassX, const Scalar massFracSalinity,
            const Scalar massFracC)
    {
        Scalar molalityX = massFracX/molarMassX/(1- massFracSalinity - massFracC);
        return molalityX;
    }


    /*!
     * \brief Returns the mass fraction of a component x (kg x / kg solution) for a given
     * molality fraction (mol x / mol solution)
     * The salinity and the mole Fraction of CO2 are considered
     *
     */

    static Scalar molalityToMassFrac(Scalar molalityX, Scalar molarMassX, Scalar massFracSalinity, Scalar massFracCTot)
    {
        Scalar massFracX = molalityX * molarMassX * (1 - massFracSalinity - massFracCTot);
        return massFracX;
    }

    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        if(moleFracX<0)
            moleFracX=0;
        if(moleFracSalinity<0)
            moleFracSalinity=0;
        if(moleFracCTot<0)
            moleFracCTot=0;

        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

    static Scalar molalityToMoleFrac(Scalar molalityX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar moleFracX = molalityX * (1 - moleFracSalinity - moleFracCTot) * FluidSystem::molarMass(H2OIdx);
        return moleFracX;
    }

    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracWater)
    {
        if(moleFracX<0)
            moleFracX=0;
        if(moleFracWater<0)
            moleFracWater=0;
        Scalar molalityX = moleFracX / moleFracWater / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

    static Scalar molalityToMoleFrac(Scalar molalityX, Scalar moleFracWater)
    {
        Scalar moleFracX = molalityX * moleFracWater * FluidSystem::molarMass(H2OIdx);
        return moleFracX;
    }


    /*!
     * \brief The ionic strength of a substance is calculated only with the salinity until now!
     */

    static Scalar ionicStrength(Scalar molalitySalinity)
    {
        Scalar ionicStrength = 0.0;
        //Scalar Product
        for (int compIdx = 0; compIdx < 2; ++compIdx)
        {
            ionicStrength += molalitySalinity;
        }
        ionicStrength *= 0.5;

        return ionicStrength;
    }
    static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa, Scalar mNH4, Scalar mFe2 )
    {
        Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
        + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
        + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx)
        + mNH4  * FluidSystem::charge(NH4Idx) * FluidSystem::charge(NH4Idx)
		+ mFe2  * FluidSystem::charge(Fe2Idx) * FluidSystem::charge(Fe2Idx));
        return ionicStrength;
    }
    static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa, Scalar mNH4 )
    {
        Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
        + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
        + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx)
        + mNH4  * FluidSystem::charge(NH4Idx) * FluidSystem::charge(NH4Idx));

        return ionicStrength;
    }
    static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa )
    {
        Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
        + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
        + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx));

        return ionicStrength;
    }

    void ionicStrength()
    {
        ionicStrength_ = 0.0;
        //Scalar Product
        for (int compIdx = 0; compIdx < 4; ++compIdx)
        {
            ionicStrength_ += molality_[compIdx] * charge_[compIdx] * charge_[compIdx];
        }
        ionicStrength_ *= 0.5;
    }

    //Calculates the activity with a modified Debye-Hückel equation after Parkhurst (1990) for
    //ionic strengths up to 2.
    static Scalar activityCoefficient(Scalar ionicStrength, Scalar temperatureK, int compIdx)
    {
        if (ionicStrength<0)
        {
            ionicStrength = 0;
        }
        Scalar charge = FluidSystem::charge(compIdx);
        Scalar ai = FluidSystem::ai(compIdx);
        Scalar bi = FluidSystem::bi(compIdx);
        Scalar A = 0.5085;
        Scalar B = 0.3285e10;
        // The actual modified Debye Hückel equation
        Scalar logActivityCoefficient = -A*(charge*charge)*sqrt(ionicStrength)/(1 + B*ai*sqrt(ionicStrength))
                + bi*ionicStrength;

        return pow(10, logActivityCoefficient);
//          return 1.0;
    }

    static Scalar J(Scalar x)
    {

      Scalar c[5], res;

      /*Pitzer 1974, Thermodaynamics of Electrolytes V*/

      c[1]=4.581;  c[2]=0.7237;  c[3]=0.0120;  c[4]=0.528;

      res = x/(4. + c[1]*pow(x,-c[2])*exp(-c[3]*pow(x,c[4])));


      return(res);
    }
    static Scalar Jprime(Scalar x)
    {

      Scalar res, eps;

      eps = 1.E-3;

      res = (J(x+eps) - J(x))/eps;


      return(res);
    }

    static Scalar Appa_Ksp(Scalar mNa, Scalar mCa, Scalar mNH4, Scalar mHCO3, Scalar mCO3, Scalar mCl, Scalar temp)
    {

      Scalar f, B_cacl, C_cacl, B1_cacl, I, sqrt_I, gamma_Ca, gamma_CO3, Ksp;
      Scalar beta_cacl_0, beta_cacl_1, C_cacl_phi;
      Scalar beta_nacl_0, beta_nacl_1, C_nacl_phi;
      Scalar beta_nahco3_0, beta_nahco3_1, C_nahco3_phi;
      Scalar beta_naco3_0, beta_naco3_1, C_naco3_phi;
      Scalar psi_canacl, psi_co3nacl, theta_naca, theta_clco3;
      Scalar B1_nacl, C_nacl, B1_nahco3, C_nahco3, B1_naco3, C_naco3, B_naco3;
      Scalar A_phi, a[6], T,x_clco3,x_clcl, x_co3co3,x_cana,x_caca,x_nana;
      Scalar E_theta_cana, E_theta_clco3, E1_theta_cana, E1_theta_clco3;

      Scalar beta_nh4cl_0, beta_nh4cl_1, beta_nh4co3_0, beta_nh4co3_1, beta_nh4hco3_0, beta_nh4hco3_1; /*new*/
      Scalar B_nh4cl, B_nh4co3, B_nh4hco3, B1_nh4cl, B1_nh4co3, B1_nh4hco3, C_nh4cl, C_nh4co3, C_nh4cl_phi, C_nh4co3_phi; /*new*/

      I = 0.5*( mNa + 4.*mCa + mNH4 + mHCO3 + 4*mCO3 + mCl) + 1.E-20;
      sqrt_I = sqrt(I);

      T = temp;
      a[0]=-8.1765300E-1; a[1]=-8.6852760E-1; a[2]=1.9251000E+4; a[3]=5.2514840E-3; a[4]=-7.1493970E-6; a[5]=9.3385590E-12;

      A_phi = a[0] + a[1]/(T-222.) + a[2]/(T*T) + a[3]*T + a[4]*T*T + a[5]*T*T*T*T;
      /*MODELING AND NUMERICAL SIMULATION OF SALT TRANSPORT AND PHASE TRANSITIONS IN UNSATURATED POROUS BUILDING MATERIALS By Andreas Nicolai*/

      beta_cacl_0 = 0.3159;  beta_cacl_1 = 1.614; C_cacl_phi = -0.00034;
      beta_nacl_0 = 0.0765; beta_nacl_1 = 0.2664; C_nacl_phi = 0.00127;
      beta_nahco3_0 = 0.0277; beta_nahco3_1 = 0.0411; C_nahco3_phi = 0.0;
      beta_naco3_0 = 0.1898; beta_naco3_1 = 0.846; C_naco3_phi = -0.048;
      psi_canacl = -0.014; psi_co3nacl = 0.016;
      theta_naca = 0.07; theta_clco3 = -0.053;


      beta_nh4co3_0 = 0.1288; beta_nh4co3_1 = 1.433; C_nh4co3_phi = 0.0005; /*new*/
      beta_nh4hco3_0 = -0.038; beta_nh4hco3_1 = 0.07;                       /*new*/
      beta_nh4cl_0 = 0.0522; beta_nh4cl_1 = 0.1918; C_nh4cl_phi = 0.003;    /*new*/


      x_clco3 = 6.*(-1.)*(-2.)*A_phi*sqrt_I;
      x_clcl = 6.*(-1.)*(-1.)*A_phi*sqrt_I;
      x_co3co3 = 6.*(-2.)*(-2.)*A_phi*sqrt_I;
      x_cana = 6.*(+2.)*(+1.)*A_phi*sqrt_I;
      x_caca = 6.*(+2.)*(+2.)*A_phi*sqrt_I;
      x_nana = 6.*(+1.)*(+1.)*A_phi*sqrt_I;


      E_theta_cana = ((+2.)*(+1.)/(4.*I))*( J(x_cana) - 0.5*J(x_caca) - 0.5*J(x_nana) );
      E_theta_clco3 = ((-1.)*(-2.)/(4.*I))*( J(x_clco3) - 0.5*J(x_clcl) - 0.5*J(x_co3co3) );

      E1_theta_cana = -(E_theta_cana/I) + ((+2)*(+1)/(8*I*I))*( x_cana*Jprime(x_cana) - 0.5*x_caca*Jprime(x_caca) - 0.5*x_nana*Jprime(x_nana) );
      E1_theta_clco3 = -(E_theta_clco3/I) + ((-1)*(-2)/(8*I*I))*( x_clco3*Jprime(x_clco3) - 0.5*x_clcl*Jprime(x_clcl) - 0.5*x_co3co3*Jprime(x_co3co3) );

        f = -A_phi * ( sqrt_I/(1. + 1.2*sqrt_I) + (2./1.2)*log(1. + 1.2*sqrt_I) );
        B_cacl = beta_cacl_0 + (beta_cacl_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));
        B1_cacl = (beta_cacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_cacl = C_cacl_phi  / (2.*pow ( (2.*1.) , 0.5 ));
        C_cacl = C_cacl_phi  / (2.*sqrt(2.*1.));

        B1_nacl = (beta_nacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_nacl = C_nacl_phi  / (2.*pow ( (1.*1.) , 0.5 ));
        C_nacl = C_nacl_phi  / (2.*sqrt(1.*1.));

        B1_nahco3 = (beta_nahco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_nahco3 = C_nahco3_phi  / (2.*pow ( (1.*1.) , 0.5 ));
        C_nahco3 = C_nahco3_phi  / (2.*sqrt(1.*1.));

        B_naco3 = beta_naco3_0 + (beta_naco3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));
        B1_naco3 = (beta_naco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_naco3 = C_naco3_phi  / (2.*pow ( (1.*2.) , 0.5 ));
        C_naco3 = C_naco3_phi  / (2.*sqrt(1.*2.));


        B_nh4cl = beta_nh4cl_0 + (beta_nh4cl_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        B1_nh4cl = (beta_nh4cl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));    /*new*/
//         C_nh4cl = C_nh4cl_phi  / (2.*pow ( (2.*1.) , 0.5 ));       /*new*/
        C_nh4cl = C_nh4cl_phi  / (2.*sqrt(2.*1.));                                         /*new*/

        B_nh4co3 = beta_nh4co3_0 + (beta_nh4co3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        B1_nh4co3 = (beta_nh4co3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));      /*new*/
//         C_nh4co3 = C_nh4co3_phi  / (2.*pow ( (1.*2.) , 0.5 ));
        C_nh4co3 = C_nh4co3_phi  / (2.*sqrt(1.*2.));                                          /*new*/

        B_nh4hco3 = beta_nh4hco3_0 + (beta_nh4hco3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        B1_nh4hco3 = (beta_nh4hco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));      /*new*/



        gamma_Ca = exp (
                4.*f
                + mCl*(2.*B_cacl + mCl*C_cacl)
                + mNa*mCl*(4.*B1_nacl + 2.*C_nacl)
                + mNa*mHCO3*(4.*B1_nahco3 + 2.*C_nahco3)
                + mNa*mCO3*(4.*B1_naco3 + 2.*C_naco3)
                + mCa*mCl*(4.*B1_cacl + 2.*C_cacl)

                        + mNH4*mCl*(4.*B1_nh4cl + 2.*C_nh4cl)        /*new*/
                        + mNH4*mHCO3*(4.*B1_nh4hco3)  /*new*/
                        + mNH4*mCO3*(4.*B1_nh4co3 + 2.*C_nh4co3)     /*new*/


                + mNa*(2.*theta_naca + 2.*E_theta_cana + mCl*psi_canacl)
                + 4.*mNa*mCa*E1_theta_cana + 4.*mCl*mCO3*E1_theta_clco3
                );
        gamma_CO3 = exp (
                 4.*f
                 + mNa*(2.*B_naco3 + mNa*C_naco3)
                 + mNa*mCl*(4.*B1_nacl + 2.*C_nacl)
                 + mNa*mHCO3*(4.*B1_nahco3 + 2.*C_nahco3)
                 + mNa*mCO3*(4.*B1_naco3 + 2.*C_naco3)

                         + mNH4*(2.*B_nh4co3 + mNH4*C_nh4co3)        /*new*/
                         + mNH4*mCl*(4.*B1_nh4cl + 2.*C_nh4cl)        /*new*/
                         + mNH4*mHCO3*(4.*B1_nh4hco3)                /*new*/
                         + mNH4*mCO3*(4.*B1_nh4co3 + 2.*C_nh4co3)       /*new*/

                 + mCa*mCl*(4.*B1_cacl + 2.*C_cacl)
                 + mCl*(2.*theta_clco3 + 2.*E_theta_clco3 + mNa*psi_co3nacl)
                 + 4.*mNa*mCa*E1_theta_cana + 4.*mCl*mCO3*E1_theta_clco3
                 );

//        Ksp = pow (10.,-8.48)/(gamma_Ca*gamma_CO3);
        Ksp = 3.31131e-9/(gamma_Ca*gamma_CO3);
      return(Ksp);
    }

   Scalar pH(const VolumeVariables &volVars)
   {
      Scalar mH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HIdx), volVars.moleFracSalinity(), volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_H/kg_H2O]
      Scalar pH = -log10(mH);
         return pH;
   }
    Scalar Omega(const Scalar mNa,
            const Scalar mCa,
            const Scalar mNH4,
            const Scalar mHCO3,
            const Scalar mCO3,
            const Scalar mCl,
            const Scalar temperature)
    {
     Scalar Ksp = Appa_Ksp( mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
     Scalar Omega_ = mCa * mCO3 / Ksp;
     return Omega_;
    }
    Scalar OmegaApprox(const Scalar mCa,
            const Scalar mCO3)
    {
     Scalar Omega_ = mCa * mCO3 / pow (10.,-8.48); // = 3.3e-9= Ksp(Standard) = pow (10.,-8.48);
     return Omega_;
    }
    Scalar rdiss(const Scalar initialPorosity,
            const Scalar volFracCalcite,
            const Scalar mNa,
            const Scalar mCa,
            const Scalar mNH4,
            const Scalar mHCO3,
            const Scalar mCO3,
            const Scalar mCl,
            const Scalar temperature,
            const Scalar mH)
    {
        Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        if (Asw < 1e-8 || isnan(Asw))
        {
            std::cout<< "Asw = "<<Asw<<std::endl;
            Asw = 0;
            std::cout<< "Asw, corrected = "<<Asw<<std::endl;
        }
        Scalar Acw = ac_ * volFracCalcite;
        if (ac_ * volFracCalcite > Asw)
            Acw = Asw;
        Scalar rdiss_ = 0;
        Scalar Omega_ = Omega(mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
        if (Omega_ <= 1)
        {
              rdiss_ = (kdiss1_ * mH + kdiss2_) * Acw * pow((1 - Omega_),ndiss_); //[mol/dm³s]
              rdiss_ *= 1000; // rdiss [mol/m³s]
        }

        return rdiss_;
    }
    Scalar rprec(const Scalar initialPorosity,
            const Scalar volFracCalcite,
            const Scalar mNa,
            const Scalar mCa,
            const Scalar mNH4,
            const Scalar mHCO3,
            const Scalar mCO3,
            const Scalar mCl,
            const Scalar temperature)
    {
        Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
             if (Asw < 1e-8 || isnan(Asw))
             {
                 std::cout<< "Asw = "<<Asw<<std::endl;
                 Asw = 0;
                 std::cout<< "Asw, corrected = "<<Asw<<std::endl;
             }
             Scalar rprec_ = 0;
             Scalar Omega_ = Omega(mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
             if (Omega_ >= 1)
             {
                 rprec_ = kprec_ * Asw * pow(Omega_ - 1 , nprec_);//[mol/dm³s]
                 rprec_ *= 1000; // rprec [mol/m³s]
             }


        return rprec_;
    }

    Scalar Fe2OmegaApprox(const Scalar mFe2,
            const Scalar mOH)
    {
     Scalar Fe2Omega_ = mFe2 * mOH * mOH / pow (10.,-4.89); // = 3.3e-9= Ksp(Standard) = pow (10.,-8.48);
     return Fe2Omega_;
    }
    Scalar frdiss(const Scalar initialPorosity,
            const Scalar volFracCalcite,
            const Scalar volFracFerrohydrite,
            const Scalar mNa,
            const Scalar mCa,
            const Scalar mNH4,
            const Scalar mHCO3,
            const Scalar mCO3,
            const Scalar mCl,
            const Scalar mFe2,
            const Scalar mOH,
            const Scalar temperature,
            const Scalar mH)
    {
        Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        if (fAsw < 1e-8 || isnan(fAsw))
        {
            std::cout<< "fAsw = "<<fAsw<<std::endl;
            fAsw = 0;
            std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
        }
        Scalar fAcw = fac_ * volFracFerrohydrite;
        if (fac_ * volFracFerrohydrite > fAsw)
            fAcw = fAsw;
        Scalar frdiss_ = 0;
        Scalar Fe2OmegaApprox_ = Fe2Omega(mFe2,  mOH);
        if (Fe2OmegaApprox_ <= 1)
        {
              frdiss_ = (fkdiss1_ * mH + fkdiss2_) * fAcw * pow((1 - Fe2OmegaApprox_),fndiss_); //[mol/dm³s]
              frdiss_ *= 1000; // rdiss [mol/m³s]
        }

        return frdiss_;
    }
    Scalar frprec(const Scalar initialPorosity,
            const Scalar volFracCalcite,
            const Scalar volFracFerrohydrite,
            const Scalar mNa,
            const Scalar mCa,
            const Scalar mNH4,
            const Scalar mHCO3,
            const Scalar mCO3,
            const Scalar mCl,
            const Scalar mFe2,
            const Scalar mOH,
            const Scalar temperature)
    {
        Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
             if (fAsw < 1e-8 || isnan(fAsw))
             {
                 std::cout<< "fAsw = "<<fAsw<<std::endl;
                 fAsw = 0;
                 std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
             }
             Scalar frprec_ = 0;
             Scalar Fe2OmegaApprox_ = Fe2Omega(mFe2,  mOH);
             if (Fe2OmegaApprox_ >= 1)
             {
                 frprec_ = fkprec_ * fAsw * pow(Fe2OmegaApprox_ - 1 , fnprec_);//[mol/dm³s]
                 frprec_ *= 1000; // rprec [mol/m³s]
             }


        return frprec_;
    }

   Sources reactionSource(const VolumeVariables &volVars,
            const Scalar dt)
    {
        Sources q(0.0);
    //  //define and compute some parameters for simplicity:
        Scalar temperature = volVars.temperature(); //temperature
        Scalar temperatureC = volVars.temperature() - 273.15; //temperature in °C
        Scalar porosity = volVars.porosity();
        Scalar initialPorosity = 1.0;
        for (int i=numSolidComponents-numInertComponents; i<numSolidComponents ; ++i)
        {
            initialPorosity   -= volVars.solidVolumeFraction(i);
        }
        Scalar Sw   =  volVars.saturation(wPhaseIdx);
        Scalar xWWater = volVars.moleFraction(wPhaseIdx,wCompIdx);
        Scalar volFracCalcite = volVars.solidVolumeFraction(cPhaseIdx);
        if (volFracCalcite < 0)
        volFracCalcite = 0;
	
        Scalar volFracFerrohydrite = volVars.solidVolumeFraction(fPhaseIdx);
        if (volFracFerrohydrite < 0)
        volFracFerrohydrite = 0;
	
        Scalar massImmUrease = volVars.solidVolumeFraction(uPhaseIdx)*volVars.solidComponentDensity(uPhaseIdx);
        if (massImmUrease < 0)
        massImmUrease = 0;

        Scalar cUrea = volVars.moleFraction(wPhaseIdx, UreaIdx) * volVars.molarDensity(wPhaseIdx) * FluidSystem::molarMass(UreaIdx);
        Scalar cUrease = volVars.moleFraction(wPhaseIdx, UreaseIdx) * volVars.molarDensity(wPhaseIdx) * FluidSystem::molarMass(UreaseIdx);

        Scalar xlSalinity = volVars.moleFraction(wPhaseIdx,NaIdx)
                            + volVars.moleFraction(wPhaseIdx,CaIdx)
                            + volVars.moleFraction(wPhaseIdx,ClIdx);
        Scalar mH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_H/kg_H2O]
        Scalar mNa = moleFracToMolality(volVars.moleFraction(wPhaseIdx,NaIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_sodium/kg_H2O]
        if (mNa < 0)
            mNa = 0;
        Scalar mCl = moleFracToMolality(volVars.moleFraction(wPhaseIdx,ClIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_chloride/kg_H2O]
        if (mCl < 0)
            mCl = 0;
        Scalar mUrea = moleFracToMolality(volVars.moleFraction(wPhaseIdx,UreaIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_urea/kg_H2O]
        if (mUrea < 0)
            mUrea = 0;
        Scalar mNH4 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,NH4Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_NH4/kg_H2O]
        if (mNH4 < 0)
            mNH4 = 0;
        Scalar mCa = moleFracToMolality(volVars.moleFraction(wPhaseIdx,CaIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_calcium/kg_H2O]
        if (mCa < 0)
            mCa = 0;
        Scalar mCO3 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,CO3Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_CO3/kg_H2O]
        if (mCO3 < 0)
            mCO3 = 0;
        Scalar mHCO3 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HCO3Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_HCO3/kg_H2O]
        if (mHCO3 < 0)
            mHCO3 = 0;
        Scalar mFe2 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,Fe2Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_HCO3/kg_H2O]
        if (mFe2 < 0)
            mFe2 = 0;
        Scalar mOH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,OHIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_HCO3/kg_H2O]
        if (mOH < 0)
            mOH = 0;
        // compute dissolution and precipitation rate of ferrohydrite
        // Scalar Ksp = Appa_Ksp( mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, volVars.temperature());
        // Scalar Omega = mCa * mCO3 / Ksp;
        Scalar Fe2OmegaApprox_ = mFe2 * mOH * mOH / pow (10.,-4.89);
        Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        if (fAsw < 1e-8 || std::isnan(fAsw))
        {
        std::cout<< "fAsw = "<<fAsw<<std::endl;
        fAsw = 0;
        std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
        }
        Scalar fAcw = fac_ * volFracFerrohydrite;
        if (fac_ * volFracFerrohydrite > fAsw)
            fAcw = fAsw;

        Scalar frdiss = 0;
        Scalar frprec = 0;
        if (Fe2OmegaApprox_ >= 1)
        {
        frdiss = 0;
        frprec = fkprec_ * fAsw * pow(Fe2OmegaApprox_ - 1 , fnprec_);//[mol/dm³s]
        frprec *= 1000; // rprec [mol/m³s]
        }
        else
        {
//             rdiss = (kdiss1_ * mH + kdiss2_) * Acw * pow((1 - Omega),ndiss_); //[mol/dm³s]
//             rdiss *= 1000; // rdiss [mol/m³s]
            frprec = 0;
        }
        if(frprec >
            volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        {
            frprec =  volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        }

        // compute dissolution and precipitation rate of calcite
        Scalar Ksp = Appa_Ksp( mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, volVars.temperature());
        Scalar Omega = mCa * mCO3 / Ksp;
        Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        if (Asw < 1e-8 || std::isnan(Asw))
        {
        std::cout<< "Asw = "<<Asw<<std::endl;
        Asw = 0;
        std::cout<< "Asw, corrected = "<<Asw<<std::endl;
        }
        Scalar Acw = ac_ * volFracCalcite;
        if (ac_ * volFracCalcite > Asw)
            Acw = Asw;

        Scalar rdiss = 0;
        Scalar rprec = 0;
        if (Omega >= 1)
        {
        rdiss = 0;
        rprec = kprec_ * Asw * pow(Omega - 1 , nprec_);//[mol/dm³s]
        rprec *= 1000; // rprec [mol/m³s]
        }
        else
        {
//             rdiss = (kdiss1_ * mH + kdiss2_) * Acw * pow((1 - Omega),ndiss_); //[mol/dm³s]
//             rdiss *= 1000; // rdiss [mol/m³s]
            rprec = 0;
        }
        if(rprec >
            volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        {
            rprec =  volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        }

        // compute first-order rate of ureolysis:
        // calculate the temperature-dependent rate coefficient
        Scalar kureaseT = kurease_ * exp(cureaseT_ / temperature);

        Scalar rurea_urease = kureaseT * cUrease * porosity * Sw * cUrea;
//         Scalar rurea_immUrease = kureaseT * massImmUrease * ureaseInEnzymeSource_ * cUrea;
        Scalar rurea_immUrease = kureaseT * massImmUrease * cUrea;

        //[mol_urea/m³s]
        rurea_urease /=FluidSystem::molarMass(UreaIdx);
        rurea_immUrease /=FluidSystem::molarMass(UreaIdx);

        //compute rate of enzyme inactivation
        Scalar kia = cia_ * exp(ciaT_ / temperature); //[1/s]
        //kg_urease/m³s]
        Scalar ria_urease = kia * cUrease * porosity * Sw;
        Scalar ria_immUrease = (kia  +
            pow(rprec * SolidSystem::molarMass(cPhaseIdx) /
            (volVars.solidComponentDensity(cPhaseIdx) * (initialPorosity - volFracCalcite)), ciaPrec_)
                               )
            * massImmUrease;


        // compute rate of ureolysis due to temperature:
        //using actual parameters fitted to experimental data:
        // calculate the temperature-dependent rate coefficient October 2019 data fitted
        Scalar ku = cu_ * exp(cuT_ / temperature - 0.5*mCa);
        Scalar rurea_temp = ku * cUrea/FluidSystem::molarMass(UreaIdx);

        //compute the combined rate
        Scalar rurea = rurea_urease + rurea_immUrease + rurea_temp;

        // compute attachment rates:
        Scalar ra_urease = ka_urease_ * cUrease * porosity * Sw;          //[kg/m³s]

        // compute detachment rates:
        Scalar rd_urease = kd_urease_ * massImmUrease;                      //[kg/m³s]

        // rdiss+rprec[mol/m³s]
        // rurea[mol/m³s]
        // q[kg/m³s]
        q[wCompIdx] += 0;
        q[nCompIdx] += rurea - rprec + rdiss;
        q[NaIdx] += 0;
        q[ClIdx] += 0;
        q[CaIdx] += - rprec + rdiss;
        q[Fe2Idx] += - frprec + frdiss;
        q[UreaIdx] += - rurea;
        q[TNHIdx] += 2 * rurea;
        q[UreaseIdx] += (- ra_urease + rd_urease - ria_urease)/FluidSystem::molarMass(UreaseIdx);
//      q[JBMIdx] += (- ra_iJBM + rd_iJBM)/FluidSystem::molarMass(JBMIdx);
        q[phiCalciteIdx] += + rprec - rdiss;
        q[phiFerrohydriteIdx] += + frprec - frdiss;
        q[phiImmUreaseIdx] += (ra_urease - rd_urease - ria_immUrease)/SolidSystem::molarMass(uPhaseIdx);

        return q;
    }

private:
    //Newton Solver which returns true if convergence is reached and false if not.
    // x(i+1) = x(i) - f(i)/df(i) = x(i) + h(i)

    bool newton1D(Scalar &xVar, const Scalar tolAbs, const int maxIter)
    {
         /*NEWTON*/
        bool converge = false;
          Scalar eps = 1.e-3;
          Scalar eps2 = 1.e-10;
          Scalar b=0;
          Scalar c=100;
          Scalar r;
          Scalar pHc = - log(xVar);
          Scalar pHb = pHc+eps;
          Scalar Hb,Hc;
          Scalar CO3l,CO3r,CO3b, CO3c;
          Scalar NH3l,NH3r,NH3b, NH3c;
          Scalar error =100;
          iter_ = 0;
          int i = 0;
//        *HCO3 = *NH4 = CO3c = 0.;
          Scalar oh,hco3,nh4,co3,co2;
          oh=hco3=nh4=co3=co2=0;
          while (absolute(c) > tolAbs)
            {
              Hb = pow(10.,-pHb);
              CO3l = 0.; CO3r = cTot_;
              CO3b = (CO3l+CO3r)/2.;
              while (absolute(error)>1.E-11)
            {
              CO3r = CO3b + eps2;
              r = cTot_ - (Hb * CO3r / k2_)  - (Hb * (Hb*CO3r/k2_) / k1_) - CO3r;
              error = cTot_ - (Hb * CO3b / k2_)  - (Hb * (Hb*CO3b/k2_) / k1_) - CO3b;
              CO3b = CO3b - (eps2*error)/(r-error);
              i++;  if (i>1.E2) break;
            }

              hco3 = Hb * CO3b / k2_;
              co2 = Hb * hco3 / k1_;

              NH3l = 0.;
              NH3r = totalnh_;
              error =100.;
              i = 0;
              NH3b = (NH3l+NH3r)/2.;
              while (absolute(error)>1.E-11)
            {
              NH3r = NH3b + eps2;
              r = totalnh_ - Hb * NH3r / ka_ - NH3r;
              error = totalnh_ - Hb * NH3b / ka_ - NH3b;
              NH3b = NH3b - (eps2*error)/(r-error);
              i++;  if (i>1.E2) break;
            }
              nh4 = Hb * NH3b / ka_;

              oh = kw_ / Hb;

              b = - Hb + 2*CO3b + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_ - 2.*fe2_;

              pHc = pHc - (eps*c)/(b-c);

              Hc = pow(10.,-pHc);
              CO3l = 0.;
              CO3r = cTot_;
              error =100.; i = 0;
              CO3c = (CO3l+CO3r)/2.;
              while (absolute(error)>1.E-11)
            {
              CO3r = CO3c + eps2;
              r = cTot_ - (Hc * CO3r / k2_)  - (Hc * (Hc*CO3r/k2_) / k1_) - CO3r;
              error = cTot_ - (Hc * CO3c / k2_)  - (Hc * (Hc*CO3c/k2_) / k1_) - CO3c;
              CO3c = CO3c - (eps2*error)/(r-error);
              i++; if (i>1.E2) break;
            }
              hco3 = Hc * CO3c / k2_;
              co2 = Hc * hco3 / k1_;

              NH3l = 0.;
              NH3r = totalnh_;
              error =100.;
              i = 0;
              NH3c = (NH3l+NH3r)/2.;
              while (absolute(error)>1.E-11)
            {
              NH3r = NH3c + eps2;
              r = totalnh_ - Hc * NH3r / ka_ - NH3r;
              error = totalnh_ - Hc * NH3c / ka_ - NH3c;
              NH3c = NH3c - (eps2*error)/(r-error);
              i++; if (i>1.E2) break;
            }
              nh4 = Hc * NH3c / ka_;

              oh = kw_ / Hc;
              c = - Hc + 2*CO3c + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;

              pHb = pHc+eps;
              iter_+=1;
              if (iter_>maxIter || isnan(error) || isnan(c))
            {
              /*sprintf(buf, "Bisection pH: %4.2f \n", pHc);
              UserWrite(buf);*/
              break;
            }
            }
          h_ = Hc;
          oh_ = kw_ / Hc;
          nh4_ = Hc * NH3c / ka_;
          co3_ = CO3c;
          hco3_ = Hc * CO3c / k2_;
          co2_ = Hc * hco3 / k1_;
//        (this->*funcPtr)(xVar);
//        Scalar h = -fdf_[0]/fdf_[1]; // h = x(i) - x(i-1)
//        Scalar hLast = h*0.5; //initial Step
//        iter_ = 0;
//        bool converge = false;
//        if (std::isnan(h))
//        {
//            return converge = false;
//        }
//
//        while(absolute(h) > tolAbs || absolute(h/hLast)  > 1 + tolRel)
//        {
//            if(iter_ > maxIter){break;}
//
//            if(iter_ > 0)
//            {
//                (this->*funcPtr)(xVar);
//                hLast = h;
//                h = -fdf_[0]/fdf_[1];
//            }
//            if (std::isnan(h))
//            {
//                return converge = false;
//            }
//
//            xVar = xVar + h;
//            iter_ = iter_ + 1;
//        }
        if(Hc < 0.0) {return converge = false;}
        if(iter_ <= maxIter) {converge = true;}
        return converge;

    }

    bool newton1D(Scalar &xVar, void (ThisType::*funcPtr)(Scalar), const Scalar tolAbs, const Scalar tolRel, const int maxIter)
    {
        if (!Valgrind::CheckDefined(xVar))
        {
            std::cout << "----!Valgrind::CheckDefined(xVar) in chemistry \n";
            DUNE_THROW(Dune::InvalidStateException, "xVar is not defined.");
        }
        (this->*funcPtr)(xVar);

        Scalar h = -fdf_[0]/fdf_[1]; // h = x(i) - x(i-1)
        Scalar hLast = h*0.5; //initial Step
        iter_ = 0;
        bool converge = false;
        if (std::isnan(h))
        {
            return converge = false;
        }

        while(absolute(h) > tolAbs || absolute(h/hLast)  > 1 + tolRel)
        {
            if(iter_ > maxIter){break;}

            if(iter_ > 0)
            {
                (this->*funcPtr)(xVar);
                hLast = h;
                h = -fdf_[0]/fdf_[1];
            }
            if (std::isnan(h))
            {
                return converge = false;
            }

            xVar = xVar + h;
            iter_ = iter_ + 1;
        }
        if(xVar < 0.0) {return converge = false;}
        if(iter_ <= maxIter) {converge = true; newtonOrBisection_ = true; }
        return converge;

    }

    //Bisection Method Solver returns true if convergence is reached and false if not.
    //xVar is the variable for which the system is solved
    //funcPtr is the pointer to the function which is to be solved
    //a0 is the lower starting value, b0 is the upper starting value. The root must be inside the interval [a0, b0]
    //tol is the stopping critium a-b
    bool bisection1D(Scalar &xVar, void (ThisType::*funcPtr)(Scalar), const Scalar a0, const Scalar b0, const Scalar tol)
    {
        Scalar iterNo = 0;
        int maxIter = 200;
        bool converge = false;
        int sfb, sfx;

        Scalar a = a0;
        Scalar b = b0;
        (this->*funcPtr)(b);
        sfb = sign(fdf_[0]);

        while(b-a > tol)
        {
            if(iterNo > maxIter)
            {
                return converge;
            }
            xVar = (b + a)/2;
            (this->*funcPtr)(xVar);
            sfx = sign(fdf_[0]);
            iterNo = iterNo + 1;
            if (sfx == 0)
                break;
            else
                {
                    if(sfx == sfb)
                    {
                        b = xVar;
                    }
                    else
                    {
                        a = xVar;
                    }
                }
        }
        newtonOrBisection_ = false;
        converge = true;
        return converge;
    }

    bool bisection1D(const Scalar tol)
    {
        bool converge = false;
      Scalar eps = 1.e-3;
      Scalar eps2 = 1.e-10;
        Scalar pHc = 7.;
        Scalar pHa = -1.;
        Scalar pHb = 15.;
        Scalar Ha,Hb,Hc;
        Scalar CO3r,CO3l,CO3a,CO3b,CO3c;
        Scalar NH3l,NH3r,NH3a,NH3b,NH3c;
        Scalar c=100.;
        Scalar a,b;
        Scalar error=100;
        Scalar r;
        iter_ = 0;
        int i = 0;
        Scalar oh,hco3,nh4,co3,co2;
        oh=hco3=nh4=co3=co2=0;
        while (absolute(c) > tol)
    {
      Ha =pow(10.,-pHa);
      CO3l = 0.;
      CO3r = cTot_;
      error =100.;
      i = 0;
      CO3a = (CO3l+CO3r)/2.;
      while (absolute(error)>1.E-10)
        {
          CO3r = CO3a + eps2;
          r = cTot_ - (Ha * CO3r / k2_)  - (Ha * (Ha*CO3r/k2_) / k1_) - CO3r;
          error = cTot_ - (Ha * CO3a / k2_)  - (Ha * (Ha*CO3a/k2_) / k1_) - CO3a;
          CO3a = CO3a - (eps2*error)/(r-error);
          i++;  if (i>1.E2) break;
        }

      hco3 = Ha * CO3a / k2_;
      co2 = Ha * hco3 / k1_;

      NH3l = 0.;
      NH3r = totalnh_;
      error =100.;
      i=0;
      NH3a = (NH3l+NH3r)/2.;
      while (absolute(error)>1.E-10)
        {
          NH3r = NH3a + eps2;
          r = totalnh_ - Ha * NH3r / ka_ - NH3r;
          error = totalnh_ - Ha * NH3a / ka_ - NH3a;
          NH3a = NH3a - (eps2*error)/(r-error);
          i++;  if (i>1.E2) break;
        }
      nh4 = Ha * NH3a / ka_;

      oh = kw_ / Ha;
      a = - Ha + 2*CO3a + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;

      Hb = pow(10.,-pHb);
      CO3l = 0.;
      CO3r = cTot_;
      error =100.;
      i = 0;
      CO3b = (CO3l+CO3r)/2.;
      while (absolute(error)>1.E-11)
        {
          CO3r = CO3b + eps2;
          r = cTot_ - (Hb * CO3r / k2_)  - (Hb * (Hb*CO3r/k2_) / k1_) - CO3r;
          error = cTot_ - (Hb * CO3b / k2_)  - (Hb * (Hb*CO3b/k2_) / k1_) - CO3b;
          CO3b = CO3b - (eps2*error)/(r-error);
          i++;  if (i>1.E2) break;
        }

      hco3 = Hb * CO3b / k2_;
      co2 = Hb * hco3 / k1_;

      NH3l = 0.;
      NH3r = totalnh_;
      error =100.;
      i = 0;
      NH3b = (NH3l+NH3r)/2.;
      while (absolute(error)>1.E-11)
        {
          NH3r = NH3b + eps2;
          r = totalnh_ - Hb * NH3r / ka_ - NH3r;
          error = totalnh_ - Hb * NH3b / ka_ - NH3b;
          NH3b = NH3b - (eps2*error)/(r-error);
          i++;  if (i>1.E2) break;
        }
      nh4 = Hb * NH3b / ka_;

      oh = kw_ / Hb;
      b = - Hb + 2*CO3b + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;

      pHc = (pHa + pHb)/2.;

      Hc = pow(10.,-pHc);
      CO3l = 0.;
      CO3r = cTot_;
      error =100.;
      i = 0;
      CO3c = (CO3l+CO3r)/2.;
      while (absolute(error)>1.E-10)
        {
          CO3r = CO3c + eps2;
          r = cTot_ - (Hc * CO3r / k2_)  - (Hc * (Hc*CO3r/k2_) / k1_) - CO3r;
          error = cTot_ - (Hc * CO3c / k2_)  - (Hc * (Hc*CO3c/k2_) / k1_) - CO3c;
          CO3c = CO3c - (eps2*error)/(r-error);
          i++; if (i>1.E2) break;
        }

      hco3 = Hc * CO3c / k2_;
      co2 = Hc * hco3 / k1_;

        NH3l = 0.;
        NH3r = totalnh_;
        error =100.;
        i = 0;
        NH3c = (NH3l+NH3r)/2.;
        while (absolute(error)>1.E-11)
    {
      NH3r = NH3c + eps2;
      r = totalnh_ - Hc * NH3r / ka_ - NH3r;
      error = totalnh_ - Hc * NH3c / ka_ - NH3c;
      NH3c = NH3c - (eps2*error)/(r-error);
      i++; if (i>1.E2) break;
    }
      nh4 = Hc * NH3c / ka_;

      oh = kw_ / Hc;
      c = - Hc + 2*CO3c + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;

      if (a*c<0.) pHb = pHc;
      else pHa = pHc;
      iter_+=1;

    }

      h_ = Hc;
      oh_ = kw_ / Hc;
      nh4_ = Hc * NH3c / ka_;
      co3_ = CO3c;
      hco3_ = Hc * CO3c / k2_;
      co2_ = Hc * hco3 / k1_;

        converge = true;
        return converge;
    }

    //Function solves electro neutrality equation f and derivative df/dH for H with constant CO2
    void H_CO2(Scalar activityH)
    {

        h_ = activityH;
        oh_ = kw_/h_;
        hco3_ = k1_*co2_/h_;
//         co3_ = k1_*k2_*co2_/pow(h_, 2.);
        co3_ = k1_*k2_*co2_/(h_*h_);
        nh4_ = totalnh_/(1+ka_/h_);

        //Solve the function
        Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = h_ + eps*h_; // x + dx
        Scalar xLeft = h_ - eps*h_; // x - dx
        Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - k1_*co2_/xRight - 2*k1_*k2_*co2_/(xRight*xRight) - cl_ + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft -  k1_*co2_/xLeft - 2*k1_*k2_*co2_/(xLeft*xLeft) - cl_ + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }


    void H_Ctot(Scalar activityH)
        {

        h_ = activityH;
        oh_ = kw_/h_;
        hco3_ = cTot_/(h_/k1_ + 1 + k2_/h_);
        co3_ = cTot_/((h_*h_)/k1_/k2_ + h_/k2_ + 1);
        co2_ = cTot_-co3_-hco3_;//h_*hco3_/k1_;
        nh4_ = totalnh_/(1+ka_/h_);

        //Solve the function
        Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = h_ + eps*h_; // x + dx
        Scalar xLeft = h_ - eps*h_; // x - dx
        Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cTot_/(xRight/k1_ + 1 + k2_/xRight) - 2*cTot_/((xRight*xRight)/k1_/k2_ + xRight/k2_ + 1) - cl_ + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_  - kw_/xLeft - cTot_/(xLeft/k1_ + 1 + k2_/xLeft) - 2*cTot_/((xLeft*xLeft)/k1_/k2_ + xLeft/k2_ + 1) - cl_ + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
    }

    //Value of numerical derivative at xVar
    /*static*/ Scalar equationNumDeri(Scalar xVar)
    {
        Scalar eps = 1e-8;
        Scalar xRight = xVar + eps*xVar; // x + dx
        Scalar xLeft = xVar - eps*xVar; // x - dx
        Scalar fRight = equationValue(xRight); // f(x+dx)
        Scalar fLeft = equationValue(xLeft); // f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/xVar; // {f(x+dx) - f(x-dx)}/2dx
        return df;
    }



    Scalar absolute(Scalar x)
    {
        if(x<0.0)
        {
            return x*(-1);
        }
        else return x;
    }

    Scalar sign(Scalar x)
    {
        if(x > 0.0)
        {
           return 1;
        }
        else if (x < 0.0)
        {
           return -1;
        }
        else
        {
            return 0.0;
        }
    }



    int iter_; //Number of iterations the Newton solver needs until convergence
    Scalar pressure_;
    Scalar temperature_;
    Scalar salinity_;
    Scalar h2o_;
    Scalar co2_;
    Scalar hco3_;
    Scalar co3_;
    Scalar oh_;
    Scalar h_;
    Scalar ca_;
    Scalar na_;
    Scalar cl_;
    Scalar fe2_;
    Scalar totalnh_;
    Scalar nh4_;
    Scalar initH_;
    Scalar ionicStrength_;
    Scalar cTot_;
    Scalar gammaH_;
    Scalar gammaCO2_;
    Scalar gammaCa_;
    Scalar gammaOH_;
    Scalar gammaHCO3_;
    Scalar gammaCO3_;
    Scalar gammaNH3_;
    Scalar gammaNH4_;
    SolVector fdf_; //Solution vector for the newtons solver every equation f solved by the newton solver for an unknown x
    // has to store f(x) in fdf_[0] and df/dx in fdf[1]
    Vector molality_;
    Vector charge_;
    Scalar x_;
    Scalar y_;
    Scalar k1_;
    Scalar k2_;
    Scalar kw_;
    Scalar ka_;
    Scalar apparentk1_;
    Scalar apparentk2_;
    Scalar apparentka_;
    bool newtonOrBisection_;

    static constexpr Scalar KpHb_ = 0;//9.14e-8;//[mol/kgH2O] Kim et al. 2000 //Not implemented by Anozie!!

    //attachment and detachment parameters
        Scalar ka_urease_;
        Scalar kd_urease_;
        Scalar ka_iJBM_;
        Scalar kd_iJBM_;

    // calcite parameters
        Scalar ac_;
        Scalar kdiss1_;
        Scalar kdiss2_;
        Scalar kprec_;
        Scalar ndiss_;
        Scalar nprec_;
        Scalar Asw0_;

    // ferrohydrite parameters
        Scalar fac_;
        Scalar fkdiss1_;
        Scalar fkdiss2_;
        Scalar fkprec_;
        Scalar fndiss_;
        Scalar fnprec_;
        Scalar fAsw0_;

    // urease parameters
        bool useHeatKilledCells_;
        bool useJackBeans_;
        Scalar kub_;
        Scalar kurease_;
        Scalar nub_;
        Scalar Keu1_;
        Scalar Keu2_;
        Scalar cia_;
        Scalar ciaT_;
        Scalar cureaseT_;
//         Scalar sorptionCoeffUrease_;
//         Scalar maxSorptionUrease_;
//         Scalar pmDensity_;
        Scalar cu_,cuT_;
        Scalar ciaPrec_;
        Scalar pKaFactor_;

public:

    // calcite parameters
       Scalar ac()    {       return ac_; }
       Scalar kdiss1()    {    return kdiss1_; }
       Scalar kdiss2()    {    return kdiss2_; }
       Scalar kprec()    {       return kprec_; }
       Scalar ndiss()    {       return ndiss_; }
       Scalar nprec()    {       return nprec_; }
       Scalar Asw0()    {       return Asw0_; }

    // ferrohydrite parameters
       Scalar fac()    {       return fac_; }
       Scalar fkdiss1()    {    return fkdiss1_; }
       Scalar fkdiss2()    {    return fkdiss2_; }
       Scalar fkprec()    {       return fkprec_; }
       Scalar fndiss()    {       return fndiss_; }
       Scalar fnprec()    {       return fnprec_; }
       Scalar fAsw0()    {       return fAsw0_; }

    // urease parameters
        Scalar kub()    {       return kub_; }
        Scalar kurease()    {   return kurease_; }
        Scalar nub()    {       return nub_; }
        Scalar Keu1()    {       return Keu1_; }
        Scalar Keu2()    {       return Keu2_; }
        Scalar cia()    {       return cia_; }
        Scalar ciaT()    {       return ciaT_; }
        Scalar cureaseT()    {       return cureaseT_; }


public:
    Scalar kprec() const
    {   return kprec_;}
    Scalar kub() const
    {   return kub_;}
    Scalar kurease() const
    {   return kurease_;}
    Scalar nprec() const
    {   return nprec_;}
    Scalar Asw0() const
    {   return Asw0_;}
    Scalar fnprec() const
    {   return fnprec_;}
    Scalar fAsw0() const
    {   return fAsw0_;}
    Scalar Keu1() const
    {   return Keu1_;}
    Scalar Keu2() const
    {   return Keu2_;}
    Scalar cia() const
    {   return cia_;}
    Scalar ciaT() const
    {   return ciaT_;}
    Scalar cureaseT() const
    {   return cureaseT_;}

  /*!
   * \brief Returns the mole fraction of NaCl \f$\mathrm{[mol \ NaCl / mol \ solution]}\f$  for a given mole fraction
   *
   * \param salinity the salinity \f$\mathrm{[kg \ NaCl / kg \ solution]}\f$
   */
  static Scalar salinityToMolFrac_(Scalar salinity) {

        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = salinity;
        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return x_NaCl;
    }
};

} // end namespace

#endif
