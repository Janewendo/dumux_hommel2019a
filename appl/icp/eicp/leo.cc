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
 * \brief Test for the the fully implicit two-phase enzymatically induced calcium carbonate
 * precipitation model.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "leocolumnproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>


/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-ParameterFile Parameter file (Input file) \n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
        const auto fileName = getParam<std::string>("Restart.File");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem, SolidSystem>();
        loadSolution(x, fileName, pvName, *gridGeometry);
    }
    else
        problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    //add specific output
    vtkWriter.addField(problem->getPermeability(), "Permeability");
    vtkWriter.addField(problem->getCalcium(), "CCalcium");
    vtkWriter.addField(problem->getUrea(), "CUrea");
    // update the output fields before write
    problem->updateVtkOutput(x);
    vtkWriter.write(restartTime);

    // instantiate time loop
    int episodeIdx_ = 0;
    // set the initial episodeIdx in the problem to zero
    problem->setEpisodeIdx(episodeIdx_);
    // set the time loop with episodes (check points)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // do it with episodes if an injection parameter file is specified
    if (hasParam("Injection.InjectionParamFile"))
    {
        // first read the times of the checkpoints from the injection file
        std::vector<Scalar> epiEnd_;
        std::string injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");
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
        Scalar tempTime = 0;

        // print file to make sure it is the right file
        std::cout << "Read file: " << injectionParameters_ << " ..." << std::endl;
        while (!injectionData.eof())
        {
            getline(injectionData, row);
            std::cout << row << std::endl;
        }
        injectionData.close();

        // read data from file
        injectionData.open(injectionParameters_);

        while (!injectionData.eof())
        {
            getline(injectionData, row);
            if (row == "EndTimes")
            {
                getline(injectionData, row);
                while (row != "#")
                {
                    if (row != "#")
                    {
                        std::istringstream ist(row);
                        ist >> tempTime;
                        // injectionParameters_ contains the time in seconds, MICP columns in hours!
                        epiEnd_.push_back(tempTime);
                    }
                    getline(injectionData, row);
                }
            }
        }
        injectionData.close();
        // check the injection data against the number of injections specified in the parameter file
        int numInjections_ = getParam<Scalar>("Injection.numInjections");
        if (epiEnd_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection end times specified in the injection data file do not match!"
                        <<"\n numInjections from parameter file = "<<numInjections_
                        <<"\n numEpisodes from injection data file = "<<epiEnd_.size()
                        <<"\n Abort!\n";
            exit(1) ;
        }
        // set the time loop with episodes of various lengths as specified in the injection file
        // set the episode ends /check points:
        for (int epiIdx = 0; epiIdx < epiEnd_.size(); ++epiIdx)
        {
            timeLoop->setCheckPoint(epiEnd_[epiIdx]);
        }
        // set the initial episodeIdx in the problem to zero
        problem->setEpisodeIdx(episodeIdx_);
        std::cout<<"test"<<std::endl;
    }

    // fixed episode lengths for convenience in output only if this is specified
    else if (hasParam("TimeLoop.EpisodeLength"))
    {
        // set the time loop with periodic episodes of constant length
        timeLoop->setPeriodicCheckPoint(getParam<Scalar>("TimeLoop.EpisodeLength"));
    }

    // no episodes
    else
    {
        // set the time loop with one big episodes
        timeLoop->setCheckPoint(tEnd);
    }

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop);

    // the linear solver
//     using LinearSolver = AMGBackend<TypeTag>;
//     auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
    using LinearSolver = Dumux::ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // set time for problem for implicit Euler scheme
        problem->setTime( timeLoop->time() + timeLoop->timeStepSize() );
        problem->setTimeStepSize( timeLoop->timeStepSize() );

        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        //only write output when desired
        if (hasParam("Problem.Inverse") && getParam<bool>("Problem.Inverse") == true)
        {
            //measurements are at the end of episodes -9x, only output when we have measurements!
            if (timeLoop->isCheckPoint())
                if (problem->injectionType(episodeIdx_)<-89
                    && problem->injectionType(episodeIdx_)>-100)
                {
                    // update the output fields before write
                    problem->updateVtkOutput(x);

                    // write vtk output
                    vtkWriter.write(timeLoop->time());
                }
        }
        else
        {
            if (timeLoop->timeStepIndex() % 10 == 0
                || timeLoop->isCheckPoint()
                )
            {
                // update the output fields before write
                problem->updateVtkOutput(x);

                // write vtk output
                vtkWriter.write(timeLoop->time());
            }
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // report and count episodes and tell the problem when advancing to the next episode
        if (hasParam("Injection.InjectionParamFile") || hasParam("TimeLoop.EpisodeLength"))
            if(timeLoop->isCheckPoint())
            {
                std::cout << "\n Episode " << episodeIdx_ << " done \n" << std::endl;

                episodeIdx_++;
                problem->setEpisodeIdx(episodeIdx_);
            }

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
