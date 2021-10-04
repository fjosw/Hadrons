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
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TOneendtrick(const std::string name);
    // destructor
    virtual ~TOneendtrick(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void prepareU1source(FermionField &src);
    void solveFermion(FermionField &solution,
                    const FermionField &source);
private:
    bool        hasT_{false};
    std::string tName_;
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(Oneendtrick, TOneendtrick<FIMPL>, MOneendtrick);

/******************************************************************************
 *                      TOneendtrick implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TOneendtrick<FImpl>::TOneendtrick(const std::string name)
: Module<OneendtrickPar>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TOneendtrick<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TOneendtrick<FImpl>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

template <typename FImpl>
std::vector<std::string> TOneendtrick<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOneendtrick<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);

    envTmpLat(FermionField, "eta");
    envTmpLat(FermionField, "chi");
    envTmpLat(FermionField, "psi");

    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(FermionField, "rng_field");

    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "tmp_source", Ls_);
        envTmpLat(FermionField, "tmp_solution", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "tmp_source");
        envTmpLat(FermionField, "tmp_solution");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOneendtrick<FImpl>::prepareU1source(FermionField &src)
{
    auto    &t_tmp   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex                         Ci(0.0,1.0);

    if (!hasT_)
    {
        LatticeCoordinate(t_tmp, Tp);
        hasT_ = true;
    }

    envGetTmp(FermionField, rng_field);

    LOG(Message) << "Preparing U1 source" << std::endl;
    random(rng4d(), rng_field); // Uniform complex random number
    rng_field = exp(Ci*2.0*M_PI*real(rng_field));
    rng_field = where((t_tmp >= par().t) and (t_tmp <= par().t), rng_field, 0.*rng_field);
    src = rng_field;
}

template <typename FImpl>
void TOneendtrick<FImpl>::solveFermion(FermionField &solution,
                                        const FermionField &source)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();

    envGetTmp(FermionField, tmp_source);
    envGetTmp(FermionField, tmp_solution);


    LOG(Message) << "Inverting using solver '" << par().solver << "'" 
                 << std::endl;
    if (Ls_ > 1)
    {
        mat.ImportPhysicalFermionSource(source, tmp_source);
    }
    else
    {
        tmp_source = source;
    }

    tmp_solution = Zero();
    solver(tmp_solution, tmp_source);

    // create 4D FermionField from 5D one if necessary
    if (Ls_ > 1)
    {
        mat.ExportPhysicalFermionSolution(tmp_solution, solution);
    }
    else
    {
        solution = tmp_solution;
    }
}

template <typename FImpl>
void TOneendtrick<FImpl>::execute(void)
{
    int                             nt = env().getDim(Tp);
    Gamma                           g5(Gamma::Algebra::Gamma5);
    std::vector<Gamma::Algebra>     gammaList;
    std::vector<ComplexD>           res_vector;
    std::vector<Result>             result;

    const std::array<const Gamma, 2> grelevant = {{
      Gamma(Gamma::Algebra::Gamma5),      
      Gamma(Gamma::Algebra::GammaTGamma5)}};

    for (auto &G: grelevant)
    {
        gammaList.push_back(G.g);
    }

    result.resize(gammaList.size());
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma_snk = gammaList[i];
        result[i].gamma_src = gammaList[0]; // Gamma5
        result[i].corr.resize(nt);
    }

    LOG(Message) << "Computing one end trick '" << getName() << "'"
                 << std::endl;
    
    envGetTmp(FermionField, eta);
    envGetTmp(FermionField, chi);
    envGetTmp(FermionField, psi);

    prepareU1source(eta);    

    solveFermion(chi, eta);

    LOG(Message) << "Computing contractions" << std::endl;

    for (unsigned int i = 0; i < result.size(); ++i)
    {
        Gamma       gSnk(gammaList[i]);

        psi = gSnk*g5*chi;

        sliceInnerProductVector(res_vector, psi, chi, Tp);

        for (unsigned int t = 0; t < nt; ++t)
        {
            result[i].corr[t] = res_vector[t];
        }
    }

    saveResult(par().output, "Oneendtrick", result);
    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MOneendtrick_Oneendtrick_hpp_
