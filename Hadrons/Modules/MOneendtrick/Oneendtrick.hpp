/*
 * Oneendtrick.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef Hadrons_MOneendtrick_Oneendtrick_hpp_
#define Hadrons_MOneendtrick_Oneendtrick_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                Oneendtrick                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MOneendtrick)

class OneendtrickPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(OneendtrickPar,
                                    std::string, solver,
                                    unsigned int, t,
                                    std::string, output);
};

template <typename FImpl>
class TOneendtrick: public Module<OneendtrickPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TOneendtrick(const std::string name);
    // destructor
    virtual ~TOneendtrick(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void prepareZ2source(FermionField &src);
    void solveFermion(FermionField &result, FermionField &propPhysical,
                         const FermionField &source);
private:
    bool        hasT_{false};
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(Oneendtrick, TOneendtrick<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZOneendtrick, TOneendtrick<ZFIMPL>, MFermion);

/******************************************************************************
 *                      TOneendtrick implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TOneendtrick<FImpl>::TOneendtrick(const std::string name)
: Module<OneendtrickPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TOneendtrick<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().solver, par().output};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TOneendtrick<FImpl>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return out;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TMeson<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOneendtrick<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);

    envTmpLat(LatticeFermion, "eta"); // Fermion Field or lattice fermion?
    
    envTmpLat(FermionField, "tmp"); 
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "source", Ls_);
        envTmpLat(FermionField, "sol", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "source");
        envTmpLat(FermionField, "sol");
    }
    // if (envHasType(PropagatorField, par().source))
    // {
    //     envCreateLat(PropagatorField, getName());
    //     if (Ls_ > 1)
    //     {
    //         envCreateLat(PropagatorField, getName() + "_5d", Ls_);
    //     }
    // }
    // else
    // {
    //     HADRONS_ERROR_REF(ObjectType, "object '" + par().source 
    //                       + "' has an incompatible type ("
    //                       + env().getObjectType(par().source)
    //                       + ")", env().getObjectAddress(par().source))
    // }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOneendtrick<FImpl>::prepareZ2source(FermionField &src)
{
    //auto    &src = envGet(FermionField, getName());
    auto    &t_tmp   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex shift(1., 1.);

    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }

    envGetTmp(LatticeFermion, eta);
    bernoulli(rng4d(), eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    eta = where((t_tmp >= par().t) and (t_tmp <= par().t), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}

template <typename FImpl>
void TOneendtrick<FImpl>::prepareU1source(FermionField &src)
{
    auto    &t_tmp   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex                         Ci(0.0,1.0);
    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }

    envGetTmp(LatticeFermion, eta);
    random(rng4d(), eta); // Uniform complex random number
    eta = exp(Ci*2*M_PI*real(eta));
    eta = where((t_tmp >= par().t) and (t_tmp <= par().t), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}

template <typename FImpl>
void TOneendtrick<FImpl>::solveFermion(PropagatorField &prop, 
                                        PropagatorField &propPhysical,
                                        const PropagatorField &fullSrc)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    LOG(Message) << "Inverting using solver '" << par().solver << "'" 
                 << std::endl;
    LOG(Message) << "Import source" << std::endl;
    if (!env().isObject5d(par().source))
    {
        if (Ls_ == 1)
        {
           //PropToFerm<FImpl>(source, fullSrc, s, c);
        }
        else
        {
            //PropToFerm<FImpl>(tmp, fullSrc, s, c);
            mat.ImportPhysicalFermionSource(tmp, source);
        }
    }
    // source conversion for 5D sources
    else
    {
        if (Ls_ != env().getObjectLs(par().source))
        {
            HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
        }
        //else
        //{
        //    PropToFerm<FImpl>(source, fullSrc, s, c);
        //}
    }
    sol = Zero();
    LOG(Message) << "Solve" << std::endl;
    solver(sol, source);
    LOG(Message) << "Export solution" << std::endl;
    //FermToProp<FImpl>(prop, sol, s, c);
    // create 4D propagators from 5D one if necessary
    if (Ls_ > 1)
    {
        mat.ExportPhysicalFermionSolution(sol, tmp);
        //FermToProp<FImpl>(propPhysical, tmp, s, c);
    }
    //}
}

template <typename FImpl>
void TOneendtrick<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");

    auto &prop         = envGet(PropagatorField, propName);
    auto &propPhysical = envGet(PropagatorField, getName());
    auto &fullSrc      = envGet(PropagatorField, par().source);

    LOG(Message) << "Using source '" << par().source << "'" << std::endl;
    solveFermion(prop, propPhysical, fullSrc);
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_Oneendtrick_hpp_
