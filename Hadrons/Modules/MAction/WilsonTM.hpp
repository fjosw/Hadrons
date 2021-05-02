/*
 * WilsonTM.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MAction_WilsonTM_hpp_
#define Hadrons_MAction_WilsonTM_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            TWilsonTM quark action                            *
 * Note that the definitions of the twisted-mass in GRID and openQCD            *
 * differ by a minus sign.                                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class WilsonTMPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonTMPar,
                                    std::string, gauge,
                                    double     , mass,
                                    double     , mu,
                                    std::string, boundary,
                                    std::string, string,
                                    std::string, twist);
};

template <typename FImpl>
class TWilsonTM: public Module<WilsonTMPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWilsonTM(const std::string name);
    // destructor
    virtual ~TWilsonTM(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonTM, TWilsonTM<FIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(WilsonTMF, TWilsonTM<FIMPLF>, MAction);
#endif

/******************************************************************************
 *                     TWilsonTM template implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWilsonTM<FImpl>::TWilsonTM(const std::string name)
: Module<WilsonTMPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWilsonTM<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWilsonTM<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilsonTM<FImpl>::setup(void)
{
    LOG(Message) << "Setting up Wilson twisted-mass fermion matrix with m= " << par().mass
                 << " and mu= " << par().mu
                 << " using gauge field '" << par().gauge << "'" << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename WilsonTMFermion<FImpl>::ImplParams implParams;
    if (!par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L << std::endl;
    if (implParams.boundary_phases.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }
    envCreateDerived(FMat, WilsonTMFermion<FImpl>, getName(), 1, U, grid, gridRb,
                     par().mass, par().mu, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilsonTM<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WilsonTM_hpp_
