/**
 * @file dv.cc
 * Source file for class \ref dv
 */


#include "dv.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/*! dv constructor function
 *
 * @param line     \input set domain pointer with.
 * @param s        \input set var_name with.
 * @param Lt       \input set L_transported with.
 * @param Lo       \input set L_output with.
 */

dv::dv(domain    *line,
       const      string s,
       const bool Lt,
       const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    LagSrc = false;

}

void dv::resize() {
    d.resize(domn->ngrd);
}




