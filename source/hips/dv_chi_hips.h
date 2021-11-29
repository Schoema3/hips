/**
 * @file dv_chi_hips.h
 * Header file for class dv_chi_hips
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_chi_hips of parent dv object.
 */

class dv_chi_hips : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        vector<double> tLC;         ///< time of last change of parcels

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

    virtual void setVar(const double vtilde, const double vhat, const int ipt=-1);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_chi_hips(){}
        dv_chi_hips(domain      *line,
                     const string s,
                     const bool   Lt,
                     const bool   Lo=true);

        virtual ~dv_chi_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

