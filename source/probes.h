/**
 * @file probes.h
 * Header file for class probes
 */

#pragma once

#include <vector>
#include <ostream>
#include <fstream>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing probes object
 *
 *  @author Marten Klein
 */

class probes {

    //////////////////// DATA MEMBERS //////////////////////

    public:

        domain                  *domn;          ///< pointer to domain object

        int                     nprb;           ///< number of probes

        int                     nof;            ///< output frequency in terms of ignored function calls

        vector<double>          pos;            ///< vector of probe position

        vector<ofstream*>       probeFiles;     ///< vector of file objects (ofstream)

    private:

        int                     ncall;          ///< output call counter for reduction of output frequency
        int                     istart;         ///< OPTIONAL: saved index to speed up ppos2cell index search (USE WITH CAUTION!)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        void openFiles();                       ///< open probe output files
        void closeFiles();                      ///< close probe output files
        void outputProbes(const double time, 
                          const bool LdoIt=false);   ///< dumb output of all probes
        void outputOneProbe(const double time, 
                            const int idx, 
                                  int &istart); ///< dumb output of selected probe; if in doubt use istart=0

    private:

        int findPos(const vector<double> &pos,
                    const double loc, 
                    const int istart=0);        ///< find the grid cell index left to the probe

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        probes(){};
        void init(domain *p_domn);
        ~probes(){};

};


////////////////////////////////////////////////////////////////////////////////


