// DW model class
// Created:     2/2/2023
// Modified:    3/2/2023
// Author:      Matthew Ellis
//

#ifndef __DW_MODEL_H__
#define __DW_MODEL_H__

#include <cmath>
#include <iostream>
#include <string>

#include "../hdr/utility.h"
#include "../hdr/calculate.h"

#include "../hdr/input_map.h"

class DW_model_t
{
public:
    double _Aex;
    double _Ms;
    double _alpha;
    double3 _L;
    double _S;
    double _gamma_p;

    double3 _N;

    double _A2;
    double _A4;

    double _DW_width;

    double _beta;
    double _P;
    double _j_dens;

    double _H_const;
    double _B_const;

    double _H_app;
    double _B_app;
    double _freq;
    double _omega;

    double _Temperature;

    double _B_k;
    double _ap;
    double _bp;


    DW_model_t()
        : _Aex(0), _Ms(0), _alpha(0), _L(0), _S(0), _N(0), _gamma_p(0), _A2(0), _A4(0), _DW_width(0), _beta(0),
        _P(0), _j_dens(0), _H_const(0), _B_const(0), _H_app(0), _B_app(0), _freq(0), _omega(0), _Temperature(0),
        _B_k(0), _ap(0), _bp(0)
        {};

    DW_model_t(double A, double Ms, double alpha, double Lx, double Ly, double Lz)
    {
        _Aex = A;
        _Ms = Ms;
        _alpha = alpha;
        _L = double3(Lx, Ly, Lz);
        _S = _L[1] * _L[2];
    }

    ~DW_model_t(){};

    void initialise();

    double K_eff(const double phi);

    double DW_width(const double phi);

    double B_q(const double q, const double phi, const double t);

    double B_phi(const double q, const double phi, const double t);

    void gradient(double &dq, double &dphi, const double q, const double phi, const double t);

};

double DW_model_t::K_eff(const double phi)
{
    return ::util::mu_0 * sqr(_Ms * sin(phi)) + _Ms * _B_k;
}

double DW_model_t::DW_width(const double phi)
{
    return 1e9 * Pi * sqrt( _Aex / K_eff(phi) );
}

void DW_model_t::initialise()
{
    // Input is read and stored in a map class
    input_map_t inputs;
    if (inputs.read_file("input") != 0)
    {
        std::cerr << "Error reading input file, exiting." << std::endl;
        exit(-1);
    }

    _Ms = inputs.get<double>("Ms");

    _Aex = inputs.get<double>("Aex");
    _alpha = inputs.get<double>("alpha");

    _L[0] = inputs.get<double>("L");
    _L[1] = inputs.get<double>("Ly");
    _L[2] = inputs.get<double>("Lz");

    _S = _L[1] * _L[2];

    _Temperature = inputs.get<double>("Temperature");
    _freq = inputs.get<double>("f");

    // Convert to per ns
    _omega = 2 * Pi * _freq * 1e-9;

    _P = inputs.get<double>("P");
    _j_dens = inputs.get<double>("j");
    _beta = inputs.get<double>("beta");

    _H_const = inputs.get<double>("H_const", 0.0);
    _B_const = ::util::mu_0 * _H_const;

    _H_app = inputs.get<double>("H");
    _B_app = ::util::mu_0 * _H_app;

    _A2 = inputs.get<double>("a2", -1.29e-7);
    _A4 = inputs.get<double>("a4", 1.632e8);

    _gamma_p = ::util::gyro_ns / (1.0 + sqr(_alpha));

    calculate::compute_demag_factor(_N, _L);
    std::cout << "Nx = " << _N[0] << ", Ny = " << _N[1] << ", Nz = " << _N[2] << std::endl;

    _B_k = util::mu_0 * _Ms * (_N[2] - _N[1]);

    _ap = (_A2 / (_Ms * _S))*1e-9;
    _bp = (2 * _A4 / (_Ms * _S))*1e-27;

    std::cout << "ap = " << _ap << std::endl;
    std::cout << "bp = " << _bp << std::endl;


}

double DW_model_t::B_q(const double q, const double phi, const double t)
{
    double B_q = _B_app * sin(_omega*t)  - _ap * q - _bp * cube(q);
    return B_q;    
}

double DW_model_t::B_phi(const double q, const double phi, const double t)
{
    double B_phi = - 0.5*_B_k *sin(2*phi);
    return B_phi;
}

void DW_model_t::gradient(double &dq, double &dphi, const double q, const double phi, const double t)
{
    const double Bq = B_q(q, phi, t);
    const double Bphi = B_phi(q, phi, t);
    const double D = DW_width(phi);

    dq = D * _gamma_p * (_alpha * Bq - Bphi);
    dphi = _gamma_p * (Bq + _alpha * Bphi);
}

#endif
