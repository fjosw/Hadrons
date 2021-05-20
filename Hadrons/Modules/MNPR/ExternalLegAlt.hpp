/*
 * ExternalLegAlt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
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

#ifndef Hadrons_ExternalLegAlt_hpp_
#define Hadrons_ExternalLegAlt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/Eigen/LU>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TExternalLegAlt                                       *
        Computes the propagator from a momentum source with
        appropriate phase corrections. Output is a SpinColourMatrix.
******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class ExternalLegAltPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExternalLegAltPar,
                                    std::string,    qIn, 
                                    std::string,    pIn,
                                    std::string,    output);
};

template <typename FImpl>
class TExternalLegAlt: public Module<ExternalLegAltPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string,  pIn);
    };
    typedef Correlator<Metadata, SpinColourMatrix> Result;
public:
    // constructor
    TExternalLegAlt(const std::string name);
    // destructor
    virtual ~TExternalLegAlt(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExternalLegAlt, ARG(TExternalLegAlt<FIMPL>), MNPR);

/******************************************************************************
 *                           TExternalLegAlt implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExternalLegAlt<FImpl>::TExternalLegAlt(const std::string name)
: Module<ExternalLegAltPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExternalLegAlt<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().qIn};
    
    return input;
}

template <typename FImpl>
std::vector<std::string> TExternalLegAlt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExternalLegAlt<FImpl>::execute(void)
{
    LOG(Message) << "Computing propagator '" << getName() << "' using"
                 << " momentum '" << par().pIn << "'"
                 << std::endl;
    BinaryWriter                    writer(par().output);
    auto                            &qIn    = envGet(PropagatorField, par().qIn);
    LatticeSpinColourMatrix         qIn_phased(env().getGrid());
    std::vector<int>                pIn  = strToVec<int>(par().pIn);
    Coordinate                      latt_size = GridDefaultLatt();  
    LatticeComplex                  pdotxin(env().getGrid()), coor(env().getGrid());
    LOG(Message) << "Propagators set up " << std::endl;
    Gamma                           g5(Gamma::Algebra::Gamma5);
    Complex                         Ci(0.0,1.0);
    Result                          r;

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    pdotxin=Zero();
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
        Real TwoPiL =  M_PI * 2.0 / latt_size[mu];
        LatticeCoordinate(coor,mu);
        pdotxin = pdotxin + (TwoPiL * pIn[mu]) * coor;
    }
    qIn_phased = qIn * exp(-Ci * pdotxin); // phase corrections

    SpinColourMatrix qIn_mom = (1.0 / volume) * sum(qIn_phased); // Divide by volume?
    LOG(Message) << "summed over lattice" << std::endl;

    r.info.pIn  = par().pIn;
    r.corr.push_back( qIn_mom );

    saveResult(par().output, "ExternalLegAlt", r);
    LOG(Message) << "Complete. Writing results to " << par().output << std:: endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ExternalLegAlt_hpp_
