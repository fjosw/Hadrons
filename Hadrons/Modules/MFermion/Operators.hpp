/*
 * Operators.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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
#ifndef Hadrons_MFermion_Operators_hpp_
#define Hadrons_MFermion_Operators_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Operators                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class OperatorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(OperatorsPar,
                                    std::string, action);
};

template <typename FImpl>
class TOperators: public Module<OperatorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TOperators(const std::string name);
    // destructor
    virtual ~TOperators(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Operators, TOperators<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZOperators, TOperators<ZFIMPL>, MFermion);

/******************************************************************************
 *                 TOperators implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TOperators<FImpl>::TOperators(const std::string name)
: Module<OperatorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TOperators<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TOperators<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_herm", getName() + "_schur"};
    
    return out;
}

template <typename FImpl>
DependencyMap TOperators<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName() + "_herm"});
    dep.insert({par().action, getName() + "_schur"});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOperators<FImpl>::setup(void)
{
    unsigned int Ls   = env().getObjectLs(par().action);
    auto         &mat = envGet(FMat, par().action);
    envCreateDerived(FBaseOp, FHermOp, getName() + "_herm", Ls, mat);
    envCreateDerived(FBaseOp, FSchurOp, getName() + "_schur", Ls, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOperators<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_Operators_hpp_
