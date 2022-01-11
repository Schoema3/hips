
#include "probes.h"
#include "domain.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** probes initialization function
 *  NOTE: call after io initialization!
 *
 * @param p_domn  \input set domain pointer with.
 */

void probes::init(domain *p_domn) {

    domn   = p_domn;

    nof    = domn->pram->modProbesCU;  ///< output freq. in terms of blind function calls

    ncall  = 0;         ///< init count of output calls

    istart = 0;         ///< init search index

    //----------------

    //MKlein 8/7/18 - TODO: HiPS support!
    // Problem is that HiPS uses the index array pLoc which complicates index searches, interpolations etc.
    if(domn->pram->LisHips) {
        *domn->io->ostrm << "\n\n*** WARNING: presently no probes are supported for HiPS" << endl << endl;
        nprb = 0;
        return;
    }

    //----------------

    pos.assign(domn->io->probePos.begin(), domn->io->probePos.end());   ///< fill the probe positions vector

    nprb = pos.size();  ///< number of probes

}

///////////////////////////////////////////////////////////////////////////////
/** Open probe files and write header corresponding to data. 
 *  Output file names: probe_<ID>.dat
 *  The Probe ID is simply the array index (0,1,...) and probes are ordered 
 *  as specified in the input file. No probes will be written when the probe 
 *  position vector is empty (zero length).
 */

void probes::openFiles() {

    if(probeFiles.size() != 0) {
        *domn->io->ostrm << "\n\n***************** ERROR: probes existing.  probeFiles.size() = " << probeFiles.size() << endl << endl;
        exit(0);
    }

    //-------- make room

    probeFiles.resize(nprb);

    //-------- setup vector of ofstream objects

    for(int id=0; id<nprb; id++) {

        /** For std=C++11 with at least GCC-5.x: probeFiles may be declared as vector<ofstream> probeFiles
        if(!*probeFiles.at(id)) {
            *domn->io->ostrm << "\n\n***************** ERROR: probe file already open.  id = " << id << endl << endl;
            exit(0);
        }
        **/

        //---- open file

        stringstream ss; ss.clear(); ss << setfill('0') << setw(2) << id;    // we won't have more than 100 probes 
        string fname = domn->io->dataDir + "probe_" + ss.str() + ".dat";

        //probeFiles.at(id).open(fname.c_str());    // now open     /** This may be used for std=C++11 with at least GCC-5.x **/
        probeFiles.at(id) = new ofstream(fname.c_str());  // now open 

        if(!*probeFiles.at(id)) {
            *domn->io->ostrm << "\n\n***************** ERROR: opening file " << fname << endl << endl;
            exit(0);
        }

        //---- write header

        *probeFiles.at(id) << "# Probe ID = " << id;
        *probeFiles.at(id) << "\n# Probe Position = " << pos.at(id);
        *probeFiles.at(id) << "\n# Left Boundary = " << domn->posf->d.at(0);
        if(!domn->pram->LisHips)
            *probeFiles.at(id) << "\n# Domain Size = " << domn->Ldomain();

        *probeFiles.at(id) << "\n#";
        *probeFiles.at(id) << setw(9) << 1 << "_time";
        for(int k=0,j=2; k<(int)domn->v.size(); k++) {
            if(domn->v.at(k)->L_output && domn->v.at(k)->var_name != "pos" && domn->v.at(k)->var_name != "posf") {
                *probeFiles.at(id) << setw(10) << j++ << "_" << domn->v.at(k)->var_name;
                *probeFiles.at(id) << setw( 6) << j++ << "_" << domn->v.at(k)->var_name << "Grad";
            }
        }
    }

}


///////////////////////////////////////////////////////////////////////////////
/** Close probe files.
 */

void probes::closeFiles() {

    //-------- close files

    for(int id=0; id<nprb; id++) {

        if(!*probeFiles.at(id)) {
            *domn->io->ostrm << "\n\n***************** ERROR: requested probe file is not open. id = " << id << endl << endl;
            //exit(0);
            continue;
        }

        probeFiles.at(id)->close();    // now close
    }

    //-------- clear the files vector

    probeFiles.clear();

}


///////////////////////////////////////////////////////////////////////////////
/** Dump all probes (in order)
 *
 *  @param time  \input time of the output
 */

void probes::outputProbes(const double time, const bool LdoIt) {

    ncall++;

    if(!LdoIt && (nprb == 0 || ncall < nof))   // do nothing for no probes and reduced output freq.
        return;

    ncall = 0;

    //-------- make sure the variables are up to date

    for(int i=0; i<(int)domn->v.size(); i++)
        domn->v.at(i)->setVar();

    //-------- dump each probe

    for(int id=0; id<nprb; id++) {
        istart = 0;        // enforce save probe search, but this is also less economic
        outputOneProbe(time, id, istart);   // USE WITH CAUTION! (requires monotonically increasing probe pos)
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Writes a data file of all properties for a single probe (in order)
 *  File name: probe_<id>.dat
 *  Format:    1_time 2_<var_name> ...
 *
 *  @param time     \input  time of the output
 *  @param idx      \input  probe index
 *  @param istart   \output start search from this index
 */

void probes::outputOneProbe(const double time, const int idx, int &istart) {

//    if(!probeFiles.size() > idx)
//        return;

//////////////////////////////
#ifndef INTERP_PROBES_CASE    /* allows over-ride during compile time */
#define INTERP_PROBES_CASE 2  /* select interpolation in this routine (1: linear, 2: quadratic) */
#endif
//////////////////////////////

    //-------- local variables

    double ppos = pos.at(idx);

    double val, grad;

    ofstream *p_ofile = probeFiles.at(idx);    // file handle
    *p_ofile << endl;   // cr
    *p_ofile << scientific;
    *p_ofile << setprecision(6);

    //-------- time

    *p_ofile << setw(15) << time;


//////////////////////////////
#if INTERP_PROBES_CASE == 1
//////////////////////////////
    // NOTE: boundary values may not coincide with BCs since these are NOT enforced here 
    //        (so those boundary values obtained here are meaningless)
    // NOTE: overestimates grad = 0 situations near the wall due to split-cell operations there
    // TODO: distinguish between inner and bdry probes

    //-------- indices 

    i = findPos(domn->pos->d, ppos, istart);  // cell index is left of probe and leaves room for i+1
    ip = i+1;

    istart = i>0 ? i-1 : 0;     // store start index for following calls; the save way is to comment this line

    //-------- interpolate and dump

    // lowest order
    for(int k=0; k<(int)domn->v.size(); k++) {
        if(domn->v.at(k)->L_output && domn->v.at(k)->var_name != "pos" && domn->v.at(k)->var_name != "posf") {

            //---- gradient by differencing interior cell-centered neighbours

            grad = (domn->v.at(k)->d.at(ip) - domn->v.at(k)->d.at(i)) / (domn->pos->d.at(ip) - domn->pos->d.at(i));

            //---- field value by linear interp/extrap of cell-centered variables

            //val = domn->v.at(k)->d.at(i) + (domn->v.at(k)->d.at(ip) - domn->v.at(k)->d.at(i)) * ( (ppos - domn->pos->d.at(i)) / (domn->pos->d.at(ip) - domn->pos->d.at(i)) );     // linear interp/extrap

            val = domn->v.at(k)->d.at(i) + grad * (ppos - domn->pos->d.at(i));   // interp/extrap with the gradient

            //---- dump (val, grad)

            *p_ofile << setw(15) << val;
            *p_ofile << setw(15) << grad;

        }
    }

//////////////////////////////
#elif INTERP_PROBES_CASE == 2
//////////////////////////////
    // NOTE: boundary values may not coincide with BCs since these are NOT enforced here 
    //        (so those boundary values obtained here are meaningless)
    // NOTE: quadratic extrap at walls can lead to large gradients when cells have been split

    //-------- indices

    int i = findPos(domn->pos->d, ppos, istart);  // cell index is left of probe and leaves room for i+1

    istart = i>0 ? i-1 : 0;     // store start index for following calls; the save way is to comment this line

    // higher-order exceptions
    if(i == 0)
        i++;        // left bdry: shift for higher-order
    else if(i == domn->ngrd-2)
        ;           // right bdry: index is save due to findPos implementation
    else
        if(ppos > domn->posf->d.at(i+1))
            i++;    // interior: shift if probe is actually behind the cell face

    int im = i-1;
    int ip = i+1;

    //-------- interpolate and dump

    double xim = domn->pos->d.at(im);
    double xi  = domn->pos->d.at(i);
    double xip = domn->pos->d.at(ip);
    double dxi  = xi - xim;
    double dxip = xip - xi;
    double qdxf = 1.0 / (dxi*dxip*(dxi+dxip));

    ppos -= xi;     // shift origin to reference point

    double vi, vip, vim;
    double a, b; //, c;

    for(int k=0; k<(int)domn->v.size(); k++) {
        if(domn->v.at(k)->L_output && domn->v.at(k)->var_name != "pos" && domn->v.at(k)->var_name != "posf") {

            // function values
            vim = domn->v.at(k)->d.at(im);
            vi  = domn->v.at(k)->d.at(i);
            vip = domn->v.at(k)->d.at(ip);

            // parabola coeffs
            a = ((vip-vi)*dxi     - (vi-vim)*dxip     ) * qdxf;
            b = ((vip-vi)*dxi*dxi + (vi-vim)*dxip*dxip) * qdxf;
            //c = vi;

            //---- gradient: second order approx on any grid (using cell-centered variables)

            //grad = (vip*dxi*dxi - vim*dxip*dxip + vi*(dxip*dxip - dxi*dxi)) / (dxi*dxip*(dxi+dxip));  // eval at stencil center
            grad = 2*a*ppos + b;    // eval at probe

            //---- field value: quadratic interp (using cell-centered variables)

            val = a*ppos*ppos + b*ppos + vi; //+ c;

            //---- dump (val, grad)

            *p_ofile << setw(15) << val;
            *p_ofile << setw(15) << grad;

        }
    }

//////////////////////////////
#else
//////////////////////////////

    *domn->io->ostrm << endl << "ERROR in probes::outputOneProbe: case '" << INTERP_PROBES_CASE << "' is invalid" << endl;
    exit(0);

//////////////////////////////
#endif
//////////////////////////////

    //--------

    p_ofile->flush();

}


////////////////////////////////////////////////////////////////////////////////
/** Given a probe position, find the grid position index to the left (use with pos)
 *  Similar to the function in the mesher.
 * 
 * @param pos    \input vector of cell positions
 * @param loc    \input location which is converted to index
 * @param istart \input OPTIONAL: start index for prode index search, default is 0
 * @return       index to the left of the probe position
 */
int probes::findPos(const vector<double> &pos, const double loc, const int istart) {

    if(loc <= pos.at(0))
        return 0;

    if(loc >= pos.at(pos.size()-1))
        return pos.size()-2;

    for(int i=istart; i<(int)pos.size(); i++)
        if(pos.at(i) >= loc)
            return i-1;

    // if here, the search has failed and all values must be considered
    for(int i=0; i<(int)pos.size(); i++)
        if(pos.at(i) >= loc)
            return i-1;

    // if here, something is just wrong
    *domn->io->ostrm << endl << "ERROR in probes::findPos: search failed" << endl;
    exit(0);
    return 0;

}
