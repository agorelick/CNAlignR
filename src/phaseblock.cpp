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


