/**
 * @file domaincase_hips_soot.h
 * Header file for class domaincase_hips_soot
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_hips_soot of parent domaincase object.
 *
 *  @author Isaac Wheeler
 */

class domaincase_hips_soot : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_hips_soot(){}
        virtual void init(domain *p_domn);
        ~domaincase_hips_soot(){}

};


////////////////////////////////////////////////////////////////////////////////

