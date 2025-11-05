#include <Rcpp.h>
using namespace Rcpp;
using std::vector;
using std::size_t;

// [[Rcpp::export]]
IntegerVector pam_build_swap(NumericMatrix D, int K) {
  int TT = D.nrow();
  if (D.ncol() != TT) stop("D must be a square matrix.");
  
  // Step 1: BUILD phase — greedy selection of initial medoids
  vector<int> medoids;
  vector<bool> is_medoid(TT, false);
  
  for (int m = 0; m < K; ++m) {
    double best_cost = R_PosInf;
    int best_idx = -1;
    
    for (int i = 0; i < TT; ++i) {
      if (is_medoid[i]) continue;
      
      // Temporarily add i as medoid
      medoids.push_back(i);
      is_medoid[i] = true;
      
      // Compute total cost
      double cost = 0.0;
      for (int t = 0; t < TT; ++t) {
        double d_min = R_PosInf;
        for (int j : medoids) {
          if (D(t, j) < d_min)
            d_min = D(t, j);
        }
        cost += d_min;
      }
      
      // Update best if needed
      if (cost < best_cost) {
        best_cost = cost;
        best_idx = i;
      }
      
      // Remove i again
      medoids.pop_back();
      is_medoid[i] = false;
    }
    
    // Permanently add best medoid
    medoids.push_back(best_idx);
    is_medoid[best_idx] = true;
  }
  
  // Step 2: SWAP phase — try to improve medoids
  bool improved = true;
  while (improved) {
    improved = false;
    for (int i = 0; i < TT; ++i) {
      if (is_medoid[i]) continue;
      
      for (int m = 0; m < K; ++m) {
        int old_medoid = medoids[m];
        medoids[m] = i;
        
        // Recompute cost
        double cost = 0.0;
        for (int t = 0; t < TT; ++t) {
          double d_min = R_PosInf;
          for (int j : medoids) {
            if (D(t, j) < d_min)
              d_min = D(t, j);
          }
          cost += d_min;
        }
        
        // Check improvement
        if (cost < R_PosInf) {
          is_medoid[old_medoid] = false;
          is_medoid[i] = true;
          improved = true;
          break;
        } else {
          // Undo change
          medoids[m] = old_medoid;
        }
      }
      
      if (improved) break;
    }
  }
  
  return wrap(medoids);
}
