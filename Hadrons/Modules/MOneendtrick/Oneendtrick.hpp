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
    typedef Correlator<Metadata, SpinColourMatrix> Result;
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
    // void prepareZ2source(FermionField &src);
    void prepareU1source(FermionField &src);
    void solveFermion(FermionField &result, FermionField &propPhysical,
                         const FermionField &source);
private:
    bool        hasT_{false};
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(Oneendtrick, TOneendtrick<FIMPL>, MFermion);
// MODULE_REGISTER_TMP(ZOneendtrick, TOneendtrick<ZFIMPL>, MFermion);

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

    envTmpLat(FermionField, "eta");
    envTmpLat(FermionField, "chi");

    envTmpLat(LatticeComplex, "c");

    envTmpLat(FermionField, "rng_field"); // Fermion Field or lattice fermion?

    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "tmp_source", Ls_);
        envTmpLat(FermionField, "tmp_sol", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "tmp_source");
        envTmpLat(FermionField, "tmp_sol");
    }
}

// execution ///////////////////////////////////////////////////////////////////
// template <typename FImpl>
// void TOneendtrick<FImpl>::prepareZ2source(FermionField &src)
// {
//     //auto    &src = envGet(FermionField, getName());
//     auto    &t_tmp   = envGet(Lattice<iScalar<vInteger>>, tName_);
//     Complex shift(1., 1.);

//     if (!hasT_)
//     {
//         LatticeCoordinate(t, Tp);
//         hasT_ = true;
//     }

//     envGetTmp(LatticeFermion, eta);
//     bernoulli(rng4d(), eta);
//     eta = (2.*eta - shift)*(1./::sqrt(2.));
//     eta = where((t_tmp >= par().t) and (t_tmp <= par().t), eta, 0.*eta);
//     src = 1.;
//     src = src*eta;
// }

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

    envGetTmp(FermionField, rng_field);
    random(rng4d(), rng_field); // Uniform complex random number
    rng_field = exp(Ci*2*M_PI*real(rng_field));
    rng_field = where((t_tmp >= par().t) and (t_tmp <= par().t), rng_field, 0.*rng_field);
    src = 1.;
    src = src*rng_field;
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
    // LOG(Message) << "Import source" << std::endl;

    if (Ls_ > 1)
    {
        mat.ImportPhysicalFermionSource(source, tmp_source);
    }
    else
    {
        tmp_source = source;
    }

    tmp_solution = Zero();
    LOG(Message) << "Solve" << std::endl;
    solver(tmp_solution, tmp_source);
    LOG(Message) << "Export solution" << std::endl;

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

#define mesonConnected(q1, q2, gSnk, gSrc) \
(g5*(gSnk))*(q1)*(adj(gSrc)*g5)*adj(q2)

template <typename FImpl>
void TOneendtrick<FImpl>::execute(void)
{
    int                    nt = env().getDim(Tp);
    std::vector<TComplex>  buf;
    std::vector<Gamma>     gammaList;
    std::vector<Result>    result;

    for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
    {
        gammaList.push_back((Gamma::Algebra)i);
    }

    result.resize(gammaList.size());
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma_snk = gammaList[i];
        result[i].gamma_src = g5.G;
        result[i].corr.resize(nt);
    }

    LOG(Message) << "Computing one end trick '" << getName() << "'"
                 << std::endl;
    
    envGetTmp(FermionField, eta);
    envGetTmp(FermionField, chi);

    prepareU1source(eta);    

    solveFermion(chi, eta);

    envGetTmp(LatticeComplex, c);
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        Gamma       gSnk(gammaList[i]);

        c = trace(mesonConnected(q1, q2, gSnk, g5));
        sliceSum(c, buf, Tp);

        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[i].corr[t] = TensorRemove(buf[t]);
        }
    }

    // for (auto &G: Gamma::gall)
    // {

    //     {
    //         SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);
            
    //         c   = trace(mesonConnected(q1, q2, gSnk, gSrc));
    //         buf = sink(c);
    //     }
    //     for (unsigned int t = 0; t < buf.size(); ++t)
    //     {
    //         result[i].corr[t] = TensorRemove(buf[t]);
    //     }

    //     r.corr.push_back( trace(mesonConnected(chi, chi, G, g5)); );
    //     result.push_back(r);
    //     r.corr.erase(r.corr.begin());
    // }

    //////////////////////////////////////////////////
    saveResult(par().output, "Oneendtrick", result);
    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_Oneendtrick_hpp_
