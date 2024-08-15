//phaseblock.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector phaseblock(NumericVector pos, NumericVector bin, double max_phaseable_distance){

    // initialize variables    
    int n = pos.length();
    NumericVector block (n); // create a vector for each rows block value

    // get starting values before loop
    double currentbin = bin[0]; // starting bin
    double currentpos = pos[0]; // starting position
    double blockcount = 0;
    block[0] = 0;

    // loop over all positions for this bin
    for(int i=1; i<n; ++i){
        double nextpos = pos[i]; // new pos
        double nextbin = bin[i]; // new bin
        double posdiff = nextpos - currentpos; // distance from the previous pos

        if(posdiff > max_phaseable_distance || nextbin > currentbin) {
            blockcount = blockcount + 1; // increment the block count
        }

        // note down the block number for this position
        block[i] = blockcount; 

        // update the pos and bin variables
        currentpos = nextpos;
        currentbin = nextbin;
    }
    return(block);
}

/*
// [[Rcpp::export]]
NumericMatrix extrapolate_BAF(NumericMatrix b, NumericVector initial_vals, double na_val) {
    int nr = b.rows();
    int nc = b.cols();
    int diff;
    
    NumericMatrix b2(nr,nc);
    NumericVector current_vals(nc);
    NumericVector new_vals(nc);

    // initialize the current values
    for(int j=0; j<nc; j++) {
        current_vals[j] = initial_vals[j]; 
    }

    // initialize matrix for the result
    for(int i=0; i<nr; i++) {
        for(int j=0; j<nc; j++) {
            b2(i,j) = b(i,j);
        }
    }

    for(int i=0; i<nr; i++) {
        new_vals = b(i,_); // get new values from the BAF matrix          

        // check if any new values are different from current values
        diff = 0;
        for(int j=0; j<nc; j++) {
            if(new_vals[j] != current_vals[j]) {
                diff=diff+1;
            }                    
        }

        if(diff > 0) {
            // if anything changes, update the current values
            for(int j=0; j<nc; j++) {
                current_vals[j] = new_vals[j];
            }       

        } else {
             // otherwise, update the BAF matrix to use the last available valuesw
            for(int j=0; j<nc; j++) {
                b(i,j) = current_vals[j];
            }
        }
    }

    return(b2);
}

*/
