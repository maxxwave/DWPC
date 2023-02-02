// DW model class
// Created:     2/2/2023
// Modified:    2/2/2023
// Author:      Matthew Ellis
//

#ifndef __DW_MODEL_H__
#define __DW_MODEL_H__


class DW_model_t {
    public:

        double _Aex;
        double _Ms;
        double _alpha;
        double _Lx;
        double _Ly;
        double _Lz;
        double _S;

        DW_model_t ()
            : _Aex(0), _Ms(0), _alpha(0), _Lx(0), _Ly(0), _Lz(0), _S(0)
        {};

        DW_model_t( double A, double Ms, double alpha, double Lx, double Ly, double Lz)
        {
            _Aex = A;
            _Ms = Ms;
            _alpha = alpha;
            _Lx = Lx;
            _Ly = Ly;
            _Lz = Lz;
            S = Ly*Lz;
        }

        ~DW_model_t() {};



};



#endif
