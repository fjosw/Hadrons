/*
 * NPRUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Felix Erben <ferben@ed.ac.uk>
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

#ifndef Hadrons_MNPR_NPRUtils_hpp_
#define Hadrons_MNPR_NPRUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MNPR)

template <typename FImpl>
class NPRUtils
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    static void sumPropagator(SpinColourMatrix &res, PropagatorField &prop);
    static void sumFourQuark(SpinColourSpinColourMatrix &res, SpinColourSpinColourMatrixField &prop);
    static void tensorProd(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b);
    static void tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b);
    // covariant derivative
    static void dslash(PropagatorField &in, const PropagatorField &out,
        const GaugeField &Umu);
    static void phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut);
    static void dot(ComplexField &pDotX, std::vector<Real> p);
};

// Decompose the sum for a PropagatorField to be able to run the code on gpus
template <typename FImpl>
void NPRUtils<FImpl>::sumPropagator(SpinColourMatrix &res, PropagatorField &prop)
{
    for(int si=0; si < Ns; ++si)
    {
        for(int sj=0; sj < Ns; ++sj)
        {
            auto pjs = peekSpin(prop, si, sj);
            for (int ci=0; ci < Nc; ++ci)
            {
                for (int cj=0; cj < Nc; ++cj)
                {
                    const ComplexD val = sum(peekColour(pjs, ci, cj));
                    res()(si,sj)(ci,cj) = val;
                }
            }
        }
    }
}

// template <typename FImpl>
// void NPRUtils<FImpl>::sumFourQuark(SpinColourSpinColourMatrix &res, SpinColourSpinColourMatrixField &prop)
// {
//     for(int si=0; si < Ns; ++si)
//     {
//         for(int sj=0; sj < Ns; ++sj)
//         {
//             auto pjs = peekSpin(prop, si, sj);
//             for (int ci=0; ci < Nc; ++ci)
//             {
//                 for (int cj=0; cj < Nc; ++cj)
//                 {
//                     SpinColourMatrix sub_sum;
//                     NPRUtils<FImpl>::sumPropagator(sub_sum, peekColour(pjs, ci, cj););
//                     res()(si,sj)(ci,cj) = sub_sum();
//                 }
//             }
//         }
//     }
// }

// Decompose the sum for a SpinColourSpinColourMatrixField to be able to run the code on gpus
template <typename FImpl>
void NPRUtils<FImpl>::sumFourQuark(SpinColourSpinColourMatrix &res, SpinColourSpinColourMatrixField &prop)
{
    for(int si=0; si < Ns; ++si)
    {
        for(int sj=0; sj < Ns; ++sj)
        {
            auto pa = PeekIndex<1>(prop, si, sj);
            for (int ci=0; ci < Nc; ++ci)
            {
                for (int cj=0; cj < Nc; ++cj)
                {
                    auto pb = PeekIndex<2>(pa, ci, cj);
                    for(int sk=0; sk < Ns; ++sk)
                    {
                        for(int sl=0; sl < Ns; ++sl)
                        {
                            auto pc = PeekIndex<3>(pb, sk, sl);
                            for (int ck=0; ck < Nc; ++ck)
                            {
                                for (int cl=0; cl < Nc; ++cl)
                                {
                                    const ComplexD val = sum(PeekIndex<4>(pc, ck, cl));
                                    res()(si,sj)(ci,cj)(sk,sl)(ck,cl) = val;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Tensor product of two PropagatorFields (Lattice Spin Colour Matrices in many FImpls)
template <typename FImpl>
void NPRUtils<FImpl>::tensorProd(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b)
{
    autoView(lret_v, lret, AcceleratorWrite);
    autoView(a_v, a, AcceleratorRead);
    autoView(b_v, b, AcceleratorRead);

    accelerator_for( site, lret_v.size(), a.Grid()->Nsimd(), {
        vTComplex left;
        for(int si=0; si < Ns; ++si)
        {
            for(int sj=0; sj < Ns; ++sj)
            {
                for (int ci=0; ci < Nc; ++ci)
                {
                    for (int cj=0; cj < Nc; ++cj)
                    {
                        left()()() = a_v[site]()(si,sj)(ci,cj);
                        lret_v[site]()(si,sj)(ci,cj)=left()*b_v[site]();
                    }
                }
            }
        }
    });
}

// Tensor product on a single site only
template <typename FImpl>
void NPRUtils<FImpl>::tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b)
{
    for(int si=0; si < Ns; ++si)
    {
    for(int sj=0; sj < Ns; ++sj)
    {
        for (int ci=0; ci < Nc; ++ci)
	{
        for (int cj=0; cj < Nc; ++cj)
	{
            const ComplexD val = TensorRemove(a()(si,sj)(ci,cj));
            lret()(si,sj)(ci,cj) = val * b();
        }}
    }}
}

// Computes gamma^mu D_mu for the given input field. Currently uses the
// symmetric derivative, though this could change in the future.
template <typename FImpl>
void NPRUtils<FImpl>::dslash(PropagatorField &out, const PropagatorField &in,
        const GaugeField &Umu)
{
    assert(&out != &in);
    out = Zero();
    PropagatorField tmp(Umu.Grid());
    typename FImpl::GaugeLinkField U(Umu.Grid());
    for (int mu = 0; mu < Nd; mu++) 
    {
        // Overall formula:
        // tmp(x) = U_\mu(x) in(x + \hat{\mu}) - U_\mu^\dag(x - \hat{\mu}) in(x - \hat{\mu})
        U = peekLorentz(Umu, mu);
        tmp = FImpl::CovShiftForward(U, mu, in);
        tmp = tmp - FImpl::CovShiftBackward(U, mu, in);

        Gamma gamma_mu = Gamma::gmu[mu];
        out += gamma_mu * tmp;
    }
    out = 0.5 * out;
}


//// Compute phases for phasing propagators
// bilinearPhase = exp(-i (pIn - pOut) \cdot x)
template <typename FImpl>
void NPRUtils<FImpl>::phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut)
{
    bilinearPhase = Zero();
    ComplexField coordinate(bilinearPhase.Grid());
    Coordinate                  latt_size = GridDefaultLatt(); 
    for (int mu = 0; mu < Nd; mu++) 
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;

        bilinearPhase += coordinate * (pIn[mu] - pOut[mu]);
    }
    Complex Ci = Complex(0.0, 1.0);
    bilinearPhase = exp(-Ci * bilinearPhase);
}


// pDotX = p \cdot x
template <typename FImpl>
void NPRUtils<FImpl>::dot(ComplexField &pDotX, std::vector<Real> p)
{
    ComplexField coordinate(pDotX.Grid());
    Coordinate                  latt_size = GridDefaultLatt(); 
    pDotX = Zero();
    for (int mu = 0; mu < Nd; mu++) 
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;
        pDotX += coordinate * p[mu];
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
