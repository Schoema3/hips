/**
 * @file randomGenerator.h
 * Header file for classes randomGenerator
 */

///////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "MersenneTwister.h"
<<<<<<< HEAD


=======
#include "processor.h"
extern processor proc;
>>>>>>> Edit_hips
/** A random number generator class. This sets up and calls the Mersenne twister.
 */

///////////////////////////////////////////////////////////////////////////////////////
//



class randomGenerator {

    private :

    MTRand mtwist;                   ///< Mersenne twister object

    public :

    inline double getRand() {
        return mtwist.rand();
    }
    inline int getRandInt(unsigned n) {
        return mtwist.randInt(n);
    }


randomGenerator(const int aseed) : mtwist(proc.myid + aseed) {
    //randomGenerator(const int aseed) : mtwist(aseed) {
        if(aseed < 0)                // randomize the seed
            mtwist.seed();
    }
    randomGenerator()          : mtwist() {}
};



