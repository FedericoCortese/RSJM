#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>
#include <cmath>
using namespace Rcpp;

// helper to compute median of a std::vector<double>
double median_vec(std::vector<double>& v) {
  std::sort(v.begin(), v.end());
  int n = v.size();
  if (n % 2 == 1) {
    return v[(n - 1) / 2];
  } else {
    return 0.5 * (v[n/2 - 1] + v[n/2]);
  }
}

// [[Rcpp::export]]
NumericMatrix WCD(const IntegerVector& s,
                  const NumericMatrix& Y,
                  int K) {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  // output: K x P
  NumericMatrix wcd(K, P);
  
  // compute sk[p] = IQR(Y[,p]) / 1.35
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    std::vector<double> col(TT);
    for (int i = 0; i < TT; ++i) col[i] = Y(i,p);
    std::sort(col.begin(), col.end());
    // type‐1 quantiles at 25% and 75%
    double q1 = col[(int)std::floor((TT - 1) * 0.25)];
    double q3 = col[(int)std::floor((TT - 1) * 0.75)];
    sk[p] = (q3 - q1) / 1.35;
    if (sk[p] == 0.0) sk[p] = 1.0;
  }
  
  // for each cluster i = 1..K
  for (int ci = 1; ci <= K; ++ci) {
    // collect row‐indices belonging to cluster ci
    std::vector<int> rows;
    for (int i = 0; i < TT; ++i) {
      if (s[i] == ci) rows.push_back(i);
    }
    int n = rows.size();
    if (n < 2) {
      continue;
    }
    
    // for each dimension p
    for (int p = 0; p < P; ++p) {
      // build the within‐cluster distance matrix for dimension p (Spk)
      NumericMatrix mat(n, n);
      for (int a = 0; a < n; ++a) {
        mat(a,a) = 0.0;
        for (int b = a+1; b < n; ++b) {
          double d = std::abs(Y(rows[a],p) - Y(rows[b],p)) / sk[p];
          mat(a,b) = mat(b,a) = d;
        }
      }
      // compute row‐medians (including the zero on the diagonal)
      double sum_meds = 0.0;
      std::vector<double> tmp(n);
      for (int a = 0; a < n; ++a) {
        for (int j = 0; j < n; ++j) tmp[j] = mat(a,j);
        sum_meds += median_vec(tmp);
      }
      wcd(ci-1, p) = sum_meds / n;
    }
  }
  
  return wcd;
}


// [[Rcpp::export]]
NumericMatrix weight_inv_exp_dist(
    const NumericMatrix& Y,
    const IntegerVector& s,
    const NumericMatrix& W,
    double zeta,
    Nullable<IntegerVector> medoids = R_NilValue
) {
  int T = Y.nrow();
  int P = Y.ncol();
  
  // 1) compute robust scales sk[p] from Y’s columns
  std::vector<double> sk(P);
  for (int p = 0; p < P; ++p) {
    std::vector<double> col(T);
    for (int i = 0; i < T; ++i) col[i] = Y(i, p);
    std::sort(col.begin(), col.end());
    double q1 = col[(int)std::floor((T - 1) * 0.25)];
    double q3 = col[(int)std::floor((T - 1) * 0.75)];
    sk[p] = (q3 - q1) / 1.35;
    if (sk[p] == 0.0) sk[p] = 1.0;  // avoid div by zero
  }
  
  // check s length
  if (s.size() != T)
    stop("Length of s must equal nrow(Y)");
  
  // branch on whether medoids is provided
  if (medoids.isNull()) {
    //  full Y vs Y: return T×T symmetric matrix 
    NumericMatrix D(T, T);
    for (int i = 0; i < T; ++i) {
      for (int j = i + 1; j < T; ++j) {
        double sum_exp = 0.0;
        int si = s[i] - 1;
        int sj = s[j] - 1;
        for (int p = 0; p < P; ++p) {
          double diff = std::abs(Y(i, p) - Y(j, p)) / sk[p];
          double wmax = std::max(W(si, p), W(sj, p));
          sum_exp += wmax * std::exp(-diff / zeta);
        }
        double val = -zeta * std::log(std::max(sum_exp, 1e-16));
        D(i, j) = D(j, i) = val;
      }
    }
    return D;
    
  } else {
    //  Y vs medoids: return T×K matrix 
    IntegerVector med(medoids);
    int K = med.size();
    // convert to 0-based and validate
    std::vector<int> midx(K);
    for (int k = 0; k < K; ++k) {
      if (med[k] < 1 || med[k] > T)
        stop("medoids must be 1..nrow(Y)");
      midx[k] = med[k] - 1;
    }
    
    NumericMatrix D(T, K);
    for (int i = 0; i < T; ++i) {
      int si = s[i] - 1;
      for (int k = 0; k < K; ++k) {
        int j = midx[k];
        int sj = s[j] - 1;
        double sum_exp = 0.0;
        for (int p = 0; p < P; ++p) {
          double diff = std::abs(Y(i, p) - Y(j, p)) / sk[p];
          double wmax = std::max(W(si, p), W(sj, p));
          sum_exp += wmax * std::exp(-diff / zeta);
        }
        D(i, k) = -zeta * std::log(std::max(sum_exp, 1e-16));
      }
    }
    return D;
  }
}

// [[Rcpp::export]]
NumericMatrix weight_inv_exp_dist_2(
    const NumericVector& darray,   // dim = c(T, T, P)
    const IntegerVector& s,         // length T, 1..K
    const NumericMatrix& W,         // K x P
    double zeta
) {
  // 0) dimension checks
  if (darray.attr("dim") == R_NilValue)
    stop("darray must have a dim attribute (T x T x P).");
  
  IntegerVector dd = darray.attr("dim");
  if (dd.size() != 3)
    stop("darray must be a 3D array with dim = c(T, T, P).");
  
  int T  = dd[0];
  int T2 = dd[1];
  int P  = dd[2];
  
  if (T2 != T)
    stop("darray must have dim[2] == dim[1] (square T x T slices).");
  
  if (s.size() != T)
    stop("Length of s must equal dim(darray)[1] (= T).");
  
  if (W.ncol() != P)
    stop("ncol(W) must equal dim(darray)[3] (= P).");
  
  // helper: access darray(i, j, p)
  // column-major layout: i + T*j + T*T*p
  auto darr = [&](int i, int j, int p) -> double {
    return darray[i + T * j + T * T * p];
  };
  
  // 1) full T x T symmetric distance matrix
  NumericMatrix D(T, T);
  
  for (int i = 0; i < T; ++i) {
    for (int j = i + 1; j < T; ++j) {
      
      int si = s[i] - 1;
      int sj = s[j] - 1;
      
      double sum_exp = 0.0;
      for (int p = 0; p < P; ++p) {
        double diff = darr(i, j, p);     // already scaled distance
        double wmax = std::max(W(si, p), W(sj, p));
        sum_exp += wmax * std::exp(-diff / zeta);
      }
      
      double val = -zeta * std::log(std::max(sum_exp, 1e-16));
      D(i, j) = D(j, i) = val;
    }
  }
  
  return D;
}




// [[Rcpp::export]]
NumericMatrix gower_dist(const NumericMatrix& Y,
                         const NumericMatrix& mu,
                         Nullable<IntegerVector> feat_type = R_NilValue,
                         std::string scale = "m" // default = "m" (max-min), "i"=IQR/1.35, "s"=std-dev
) {
  int n = Y.nrow();            // observations
  int p = Y.ncol();            // features
  int m = mu.nrow();           // prototypes
  
  // 1. Handle NULL feat_type -> all continuous
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if (ft.size() != p)
      stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(p, 0);
  }
  
  // 2. Compute s_p: scale for continuous features per 'scale' flag
  std::vector<double> s_p(p);
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 0) {
      // continuous feature
      std::vector<double> col(n);
      for (int i = 0; i < n; ++i) col[i] = Y(i, j);
      if (scale == "m") {
        // max-min
        double mn = col[0], mx = col[0];
        for (double v: col) {
          if (v < mn) mn = v;
          if (v > mx) mx = v;
        }
        s_p[j] = mx - mn;
      } else if (scale == "i") {
        // IQR / 1.35
        std::sort(col.begin(), col.end());
        double q1 = col[(int)std::floor((n - 1) * 0.25)];
        double q3 = col[(int)std::floor((n - 1) * 0.75)];
        s_p[j] = (q3 - q1) / 1.35;
      } else if (scale == "s") {
        // standard deviation
        double sum = 0.0;
        for (double v: col) sum += v;
        double mu_col = sum / n;
        double ss = 0.0;
        for (double v: col) ss += (v - mu_col) * (v - mu_col);
        s_p[j] = std::sqrt(ss / (n - 1));
      } else {
        stop("Invalid scale flag: must be 'm', 'i', or 's'");
      }
      if (s_p[j] == 0.0) s_p[j] = 1.0;
    } else {
      // categorical or ordinal
      s_p[j] = 1.0;
    }
  }
  
  // 3. Precompute ordinal levels and ranks
  std::vector< std::vector<int> > ord_rank_Y(p);
  std::vector< std::vector<int> > ord_rank_mu(p);
  std::vector<int> M(p, 1);
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 2) {
      std::vector<double> vals;
      vals.reserve(n + m);
      for (int i = 0; i < n; ++i) vals.push_back(Y(i, j));
      for (int u = 0; u < m; ++u) vals.push_back(mu(u, j));
      std::sort(vals.begin(), vals.end());
      vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
      int levels = vals.size();
      M[j] = levels > 1 ? levels : 1;
      ord_rank_Y[j].resize(n);
      ord_rank_mu[j].resize(m);
      for (int i = 0; i < n; ++i) {
        ord_rank_Y[j][i] = std::lower_bound(vals.begin(), vals.end(), Y(i, j)) - vals.begin();
      }
      for (int u = 0; u < m; ++u) {
        ord_rank_mu[j][u] = std::lower_bound(vals.begin(), vals.end(), mu(u, j)) - vals.begin();
      }
    }
  }
  
  // 4. Compute Gower distances
  NumericMatrix V(n, m);
  for (int i = 0; i < n; ++i) {
    for (int u = 0; u < m; ++u) {
      double acc = 0.0;
      for (int j = 0; j < p; ++j) {
        double diff;
        if (ft[j] == 0) {
          // continuous
          diff = std::abs(Y(i, j) - mu(u, j)) / s_p[j];
        } else if (ft[j] == 1) {
          // categorical
          diff = (Y(i, j) != mu(u, j)) ? 1.0 : 0.0;
        } else {
          // ordinal
          double denom = double(M[j] - 1);
          diff = denom > 0.0 ? std::abs(ord_rank_Y[j][i] - ord_rank_mu[j][u]) / denom : 0.0;
        }
        acc += diff;
      }
      V(i, u) = acc / p;  // mean over features
    }
  }
  
  return V;
}

// [[Rcpp::export]]
NumericVector gower_dist_array(const NumericMatrix& Y,
                               const NumericMatrix& mu,
                               Nullable<IntegerVector> feat_type = R_NilValue,
                               std::string scale = "m" // "m"=max-min, "i"=IQR/1.35, "s"=sd
) {
  int n = Y.nrow();   // observations
  int p = Y.ncol();   // features
  int m = mu.nrow();  // prototypes
  
  // 1) Handle NULL feat_type -> all continuous
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if (ft.size() != p) stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(p, 0);
  }
  
  // 2) Compute s_p: scale for continuous features per 'scale' flag
  std::vector<double> s_p(p);
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 0) {
      std::vector<double> col(n);
      for (int i = 0; i < n; ++i) col[i] = Y(i, j);
      
      if (scale == "m") {
        double mn = col[0], mx = col[0];
        for (double v : col) {
          if (v < mn) mn = v;
          if (v > mx) mx = v;
        }
        s_p[j] = mx - mn;
        
      } else if (scale == "i") {
        std::sort(col.begin(), col.end());
        double q1 = col[(int)std::floor((n - 1) * 0.25)];
        double q3 = col[(int)std::floor((n - 1) * 0.75)];
        s_p[j] = (q3 - q1) / 1.35;
        
      } else if (scale == "s") {
        double sum = 0.0;
        for (double v : col) sum += v;
        double mu_col = sum / n;
        double ss = 0.0;
        for (double v : col) ss += (v - mu_col) * (v - mu_col);
        s_p[j] = std::sqrt(ss / (n - 1));
        
      } else {
        stop("Invalid scale flag: must be 'm', 'i', or 's'");
      }
      
      if (s_p[j] == 0.0) s_p[j] = 1.0;
      
    } else {
      // categorical or ordinal
      s_p[j] = 1.0;
    }
  }
  
  // 3) Precompute ordinal levels and ranks
  std::vector< std::vector<int> > ord_rank_Y(p);
  std::vector< std::vector<int> > ord_rank_mu(p);
  std::vector<int> M(p, 1);
  
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 2) {
      std::vector<double> vals;
      vals.reserve(n + m);
      for (int i = 0; i < n; ++i) vals.push_back(Y(i, j));
      for (int u = 0; u < m; ++u) vals.push_back(mu(u, j));
      std::sort(vals.begin(), vals.end());
      vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
      
      int levels = (int)vals.size();
      M[j] = levels > 1 ? levels : 1;
      
      ord_rank_Y[j].resize(n);
      ord_rank_mu[j].resize(m);
      
      for (int i = 0; i < n; ++i) {
        ord_rank_Y[j][i] =
          (int)(std::lower_bound(vals.begin(), vals.end(), Y(i, j)) - vals.begin());
      }
      for (int u = 0; u < m; ++u) {
        ord_rank_mu[j][u] =
          (int)(std::lower_bound(vals.begin(), vals.end(), mu(u, j)) - vals.begin());
      }
    }
  }
  
  // 4) Compute per-feature Gower distances: array n x m x p
  NumericVector D(n * m * p);
  D.attr("dim") = IntegerVector::create(n, m, p);
  
  for (int j = 0; j < p; ++j) {
    for (int u = 0; u < m; ++u) {
      for (int i = 0; i < n; ++i) {
        
        double diff;
        if (ft[j] == 0) {
          diff = std::abs(Y(i, j) - mu(u, j)) / s_p[j];
        } else if (ft[j] == 1) {
          diff = (Y(i, j) != mu(u, j)) ? 1.0 : 0.0;
        } else {
          double denom = double(M[j] - 1);
          diff = (denom > 0.0) ? std::abs(ord_rank_Y[j][i] - ord_rank_mu[j][u]) / denom : 0.0;
        }
        
        // R is column-major: i + n*u + n*m*j corresponds to [i, u, j]
        D[i + n * u + n * m * j] = diff;
      }
    }
  }
  
  return D;
}


// [[Rcpp::export]]
IntegerVector initialize_states(const NumericMatrix& Y,
                                int K,
                                Nullable<IntegerVector> feat_type = R_NilValue,
                                int reps = 10,
                                std::string scale = "m") {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  // handle feat_type
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if ((int)ft.size() != P) stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(P, 0);
  }
  
  // Precompute full Gower distance matrix Y->Y
  NumericMatrix Dall = gower_dist(Y, Y, ft, scale);
  
  double best_sum = R_PosInf;
  IntegerVector best_assign(TT);
  
  // Repeat multiple times
  for (int rep = 0; rep < reps; ++rep) {
    // 1) Initialize centroids indices via kmeans++
    std::vector<int> centIdx;
    centIdx.reserve(K);
    // first centroid random
    int idx0 = std::floor(R::runif(0, TT));
    centIdx.push_back(idx0);
    
    // distances to nearest centroid
    std::vector<double> closestDist(TT);
    for (int j = 0; j < TT; ++j) closestDist[j] = Dall(idx0, j);
    
    // choose remaining centroids
    for (int k = 1; k < K; ++k) {
      // sample next centroid with prob proportional to closestDist
      double sumd = std::accumulate(closestDist.begin(), closestDist.end(), 0.0);
      if (sumd <= 0) {
        idx0 = std::floor(R::runif(0, TT));
      } else {
        double u = R::runif(0, sumd);
        double cum = 0;
        int idx = 0;
        for (; idx < TT; ++idx) {
          cum += closestDist[idx];
          if (cum >= u) break;
        }
        if (idx >= TT) idx = TT - 1;
        idx0 = idx;
      }
      centIdx.push_back(idx0);
      // update closestDist
      for (int j = 0; j < TT; ++j) {
        closestDist[j] = std::min(closestDist[j], Dall(centIdx[k], j));
      }
    }
    
    // Assign each point to nearest centroid
    // Compute distance Y->centroids via Dall
    double sum_intra = 0.0;
    IntegerVector assign(TT);
    for (int i = 0; i < TT; ++i) {
      int best_k = 0;
      double best_d = Dall(centIdx[0], i);
      for (int k = 1; k < K; ++k) {
        double d = Dall(centIdx[k], i);
        if (d < best_d) {
          best_d = d;
          best_k = k;
        }
      }
      assign[i] = best_k + 1;  // 1-based cluster
      sum_intra += best_d;
    }
    
    // 3) Keep the best initialization
    if (sum_intra < best_sum) {
      best_sum = sum_intra;
      best_assign = assign;
    }
  }
  
  return best_assign;
}

double mad_vec(const std::vector<double>& v, double med) {
  int n = v.size();
  if (n == 0) return NA_REAL;
  std::vector<double> absdev(n);
  for (int i = 0; i < n; ++i) absdev[i] = std::abs(v[i] - med);
  return median_vec(absdev);
}

// [[Rcpp::export]]
NumericVector v_1(const NumericMatrix& Y,
                  int knn = 10,
                  Nullable<NumericVector> c_param = R_NilValue,
                  Nullable<NumericVector> qt = R_NilValue,
                  Nullable<NumericVector> Mparam = R_NilValue,
                  Nullable<IntegerVector> feat_type = R_NilValue,
                  std::string scale = "m") {
  
  int T = Y.nrow();
  int P = Y.ncol();
  
  // handle feat_type
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if (ft.size() != P) stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(P, 0);
  }
  
  // handle Mparam
  bool useM = Mparam.isNotNull();
  double Mval = 0.0;
  if (useM) {
    NumericVector tmp = Mparam.get();
    if (tmp.size() != 1) stop("M must be length 1");
    Mval = tmp[0];
  }
  
  // handle c and qt
  bool use_c = c_param.isNotNull();
  double c_val = 0.0;
  double qt_val = 0.95;  // default quantile
  
  if (qt.isNotNull()) {
    NumericVector tmp_qt = qt.get();
    if (tmp_qt.size() != 1) stop("qt must be length 1");
    qt_val = tmp_qt[0];
  }
  
  // 1) Gower dissimilarity
  NumericMatrix D = gower_dist(Y, Y, feat_type, scale);
  
  // 2) k-distance and neighborhoods
  std::vector<double> d_knn(T);
  std::vector<std::vector<int>> N_knn(T);
  for (int i = 0; i < T; ++i) {
    std::vector<double> dists;
    dists.reserve(T - 1);
    for (int j = 0; j < T; ++j) {
      if (j != i) dists.push_back(D(i, j));
    }
    std::sort(dists.begin(), dists.end());
    d_knn[i] = dists[knn - 1];
    for (int j = 0; j < T; ++j) {
      if (j != i && D(i, j) <= d_knn[i]) N_knn[i].push_back(j);
    }
  }
  
  // 3) reachability distances
  std::vector<std::vector<double>> reach_dist(T);
  for (int i = 0; i < T; ++i) {
    for (int o : N_knn[i]) {
      reach_dist[i].push_back(std::max(d_knn[o], D(i, o)));
    }
  }
  
  // 4) local reachability density
  std::vector<double> lrd(T);
  for (int i = 0; i < T; ++i) {
    if (reach_dist[i].empty()) {
      lrd[i] = NA_REAL;
    } else {
      double sum = std::accumulate(reach_dist[i].begin(), reach_dist[i].end(), 0.0);
      lrd[i] = 1.0 / (sum / reach_dist[i].size());
    }
  }
  
  // 5) standard LOF
  std::vector<double> lof(T);
  for (int i = 0; i < T; ++i) {
    auto& neigh = N_knn[i];
    if (neigh.empty() || R_IsNA(lrd[i])) {
      lof[i] = NA_REAL;
    } else {
      double sum = 0.0;
      for (int o : neigh) sum += lrd[o] / lrd[i];
      lof[i] = sum / neigh.size();
    }
  }
  
  // 6) scaled LOF*
  std::vector<double> lof_star(T);
  for (int i = 0; i < T; ++i) {
    auto& neigh = N_knn[i];
    if (neigh.size() < 2 || R_IsNA(lof[i])) {
      lof_star[i] = NA_REAL;
    } else {
      std::vector<double> vals;
      vals.reserve(neigh.size());
      for (int o : neigh) vals.push_back(lof[o]);
      double mu = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
      double var = 0.0;
      for (double v : vals) var += (v - mu) * (v - mu);
      var = var / (vals.size() - 1);
      double sd = var > 0 ? std::sqrt(var) : 1.0;
      lof_star[i] = (lof[i] - mu) / sd;
    }
  }
  
  // ompute c from quantile if not provided 
  std::vector<double> lof_clean;
  lof_clean.reserve(T);
  for (double v : lof_star)
    if (!R_IsNA(v)) lof_clean.push_back(v);
    std::sort(lof_clean.begin(), lof_clean.end());
    
    if (use_c) {
      NumericVector tmp = c_param.get();
      if (tmp.size() != 1) stop("c must be length 1");
      c_val = tmp[0];
    } else {
      int pos = std::floor(qt_val * (lof_clean.size() - 1));
      c_val = lof_clean[pos];
    }
    
    // 7) modified score v
    NumericVector v(T);
    for (int i = 0; i < T; ++i) {
      double ls = lof_star[i];
      auto& neigh = N_knn[i];
      if (R_IsNA(ls) || neigh.empty()) {
        v[i] = NA_REAL;
        continue;
      }
      
      double M_i;
      if (useM) {
        M_i = Mval;
      } else {
        std::vector<double> nb;
        nb.reserve(neigh.size());
        for (int o : neigh) nb.push_back(lof_star[o]);
        double med = median_vec(nb);
        double mad = mad_vec(nb, med);
        M_i = med + mad;
      }
      
      if (ls <= M_i) {
        v[i] = 1.0;
      } else if (ls >= c_val) {
        v[i] = 0.0;
      } else {
        double t = (ls - M_i) / (c_val - M_i);
        v[i] = std::pow(1.0 - t * t, 2);
      }
    }
    
    return v;
}


// [[Rcpp::export]]
Rcpp::List E_step(const NumericMatrix& loss_by_state,
                  const NumericMatrix& Gamma) {
  int T = loss_by_state.nrow();
  int K = loss_by_state.ncol();
  
  if (Gamma.nrow() != K || Gamma.ncol() != K)
    stop("Gamma must be K x K with K = ncol(loss_by_state)");
  
  // 1) Forward pass: V[t,j] = loss[t,j] + min_i { V[t+1,i] + Gamma[i,j] }
  NumericMatrix V(T, K);
  
  // initialize with the immediate loss
  for (int t = 0; t < T; ++t)
    for (int j = 0; j < K; ++j)
      V(t, j) = loss_by_state(t, j);
  
  // fill rows T-2 ... 0
  for (int t = T - 2; t >= 0; --t) {
    for (int j = 0; j < K; ++j) {
      double m = V(t+1, 0) + Gamma(0, j);
      for (int i = 1; i < K; ++i) {
        double cand = V(t+1, i) + Gamma(i, j);
        if (cand < m) m = cand;
      }
      V(t, j) = loss_by_state(t, j) + m;
    }
  }
  
  // 2) Backtrack
  IntegerVector s(T);
  
  // first time‐point: pick argmin over j of V(0,j)
  {
    double m0 = V(0, 0);
    int idx = 0;
    for (int j = 1; j < K; ++j) {
      if (V(0, j) < m0) {
        m0  = V(0, j);
        idx = j;
      }
    }
    s[0] = idx + 1;  
  }
  
  // t = 1..T-1
  for (int t = 1; t < T; ++t) {
    int prev = s[t-1] - 1;  // zero-based
    double m = V(t, 0) + Gamma(prev, 0);
    int idx = 0;
    for (int j = 1; j < K; ++j) {
      double cand = V(t, j) + Gamma(prev, j);
      if (cand < m) {
        m   = cand;
        idx = j;
      }
    }
    s[t] = idx + 1;
  }
  
  // Return both
  return Rcpp::List::create(
    Rcpp::Named("s") = s,
    Rcpp::Named("V") = V
  );
}

// restituisce Q1 e Q3 usando median-halves
void median_halves_q1q3(const std::vector<double>& vals, double &q1, double &q3) {
  std::vector<double> v = vals;            // copia per ordinare
  std::sort(v.begin(), v.end());
  int n = v.size();
  if (n == 0) {
    q1 = R_NaReal;
    q3 = R_NaReal;
    return;
  }
  int half = n / 2;
  // costruisci lower e upper secondo la regola "median-halves"
  std::vector<double> lower, upper;
  if (n % 2 == 0) {
    // even: lower = v[0..half-1], upper = v[half..n-1]
    lower.assign(v.begin(), v.begin() + half);
    upper.assign(v.begin() + half, v.end());
  } else {
    // odd: exclude median at index half
    // lower = v[0..half-1], upper = v[half+1..n-1]
    lower.assign(v.begin(), v.begin() + half);
    upper.assign(v.begin() + half + 1, v.end());
  }
  
  if (lower.empty()) q1 = R_NaReal; else q1 = median_vec(lower);
  if (upper.empty()) q3 = R_NaReal; else q3 = median_vec(upper);
}

// tukey versione C++ (già fornita prima)
// [[Rcpp::export]]
Rcpp::NumericVector tukey_biw_vec_cpp(Rcpp::NumericVector u, double c = 4.685) {
  int n = u.size();
  Rcpp::NumericVector out(n);
  double c2_over6 = (c * c) / 6.0;
  
  for (int i = 0; i < n; ++i) {
    if (Rcpp::NumericVector::is_na(u[i])) {
      out[i] = NA_REAL;
    } else {
      double v = std::abs(u[i]);
      if (v <= c) {
        double w = v / c;
        double tmp = 1.0 - w * w;
        out[i] = c2_over6 * (1.0 - tmp * tmp * tmp);
      } else {
        out[i] = c2_over6;
      }
    }
  }
  return out;
}

// safe_scale_slice che usa median-halves per "i"
// [[Rcpp::export]]
Rcpp::NumericMatrix safe_scale_slice_median_cpp(Rcpp::NumericMatrix mat, std::string scale = "m") {
  int nr = mat.nrow();
  int nc = mat.ncol();
  Rcpp::NumericMatrix out = Rcpp::clone(mat);
  
  if (scale == "m") {
    double minv = R_PosInf;
    double maxv = R_NegInf;
    bool any_non_na = false;
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        double val = mat(i, j);
        if (!Rcpp::NumericVector::is_na(val)) {
          any_non_na = true;
          if (val < minv) minv = val;
          if (val > maxv) maxv = val;
        }
      }
    }
    if (!any_non_na) return out;
    double sc = maxv - minv;
    if (sc > 0 && R_finite(sc)) {
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          double val = out(i, j);
          if (!Rcpp::NumericVector::is_na(val)) out(i, j) = val / sc;
        }
      }
    }
    return out;
  }
  
  if (scale == "i") {
    std::vector<double> vals;
    vals.reserve(nr * nc);
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        double val = mat(i, j);
        if (!Rcpp::NumericVector::is_na(val)) vals.push_back(val);
      }
    }
    if (vals.empty()) return out;
    
    double q1 = NA_REAL, q3 = NA_REAL;
    median_halves_q1q3(vals, q1, q3);
    if (Rcpp::NumericVector::is_na(q1) || Rcpp::NumericVector::is_na(q3)) return out;
    
    double iq = q3 - q1;
    if (iq > 0 && R_finite(iq)) {
      double scale_factor = 1.35 / iq;
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          double val = out(i, j);
          if (!Rcpp::NumericVector::is_na(val)) out(i, j) = val * scale_factor;
        }
      }
    }
    return out;
  }
  
  return out;
}

// =======================
// MAD (constant = 1, na.rm=TRUE)
// =======================
double mad_cpp(Rcpp::NumericVector x) {
  std::vector<double> vals;
  int n = x.size();
  vals.reserve(n);
  
  for (int i = 0; i < n; ++i) {
    if (!Rcpp::NumericVector::is_na(x[i]))
      vals.push_back(x[i]);
  }
  
  if (vals.empty()) return NA_REAL;
  
  double med = median_vec(vals);
  
  std::vector<double> dev;
  dev.reserve(vals.size());
  for (double v : vals)
    dev.push_back(std::fabs(v - med));
  
  double mad = median_vec(dev);
  return mad;
}

// =======================
// Tukey scalar
// =======================
inline double tukey_scalar(double u, double c) {
  double v = std::fabs(u);
  double c2_over6 = (c * c) / 6.0;
  
  if (v <= c) {
    double w = v / c;
    double tmp = 1.0 - w * w;
    return c2_over6 * (1.0 - tmp * tmp * tmp);
  } else {
    return c2_over6;
  }
}

// =======================
// median-halves IQR
// =======================
double iqr_median_halves(std::vector<double> vals) {
  int n = vals.size();
  if (n < 2) return NA_REAL;
  
  std::sort(vals.begin(), vals.end());
  
  int half = n / 2;
  std::vector<double> lower, upper;
  
  if (n % 2 == 0) {
    lower.assign(vals.begin(), vals.begin() + half);
    upper.assign(vals.begin() + half, vals.end());
  } else {
    lower.assign(vals.begin(), vals.begin() + half);
    upper.assign(vals.begin() + half + 1, vals.end());
  }
  
  if (lower.empty() || upper.empty()) return NA_REAL;
  
  double q1 = median_vec(lower);
  double q3 = median_vec(upper);
  
  return q3 - q1;
}

// =======================
// MAIN FUNCTION
// =======================
// [[Rcpp::export]]
Rcpp::NumericVector compute_dttp_cpp(
    Rcpp::NumericMatrix Y,
    Rcpp::IntegerVector cont_feat,
    bool tukey = true,
    std::string scale = "m",
    double c_tukey = 4.685
) {
  
  int TT = Y.nrow();
  int P_cont = cont_feat.size();
  
  // output: TT x TT x P_cont
  Rcpp::NumericVector out(TT * TT * P_cont);
  out.attr("dim") = Rcpp::IntegerVector::create(TT, TT, P_cont);
  
  for (int p_idx = 0; p_idx < P_cont; ++p_idx) {
    
    int p = cont_feat[p_idx] - 1; // R -> C index
    Rcpp::NumericVector x = Y(_, p);
    
    double sc_mad = mad_cpp(x);
    if (!R_finite(sc_mad) || sc_mad <= 0)
      sc_mad = 1.0;
    
    // compute TT x TT slice
    for (int i = 0; i < TT; ++i) {
      for (int j = 0; j < TT; ++j) {
        
        double xi = x[i];
        double xj = x[j];
        
        double val;
        
        if (Rcpp::NumericVector::is_na(xi) ||
            Rcpp::NumericVector::is_na(xj)) {
          val = NA_REAL;
        } else {
          val = std::fabs(xi - xj) / sc_mad;
          if (tukey)
            val = tukey_scalar(val, c_tukey);
        }
        
        out[i + TT * j + TT * TT * p_idx] = val;
      }
    }
    
    // ===== scaling step per slice =====
    if (scale == "m") {
      
      double minv = R_PosInf;
      double maxv = R_NegInf;
      
      for (int i = 0; i < TT; ++i)
        for (int j = 0; j < TT; ++j) {
          double v = out[i + TT * j + TT * TT * p_idx];
          if (!Rcpp::NumericVector::is_na(v)) {
            if (v < minv) minv = v;
            if (v > maxv) maxv = v;
          }
        }
        
        double sc = maxv - minv;
      
      if (sc > 0 && R_finite(sc)) {
        for (int i = 0; i < TT; ++i)
          for (int j = 0; j < TT; ++j) {
            double &v = out[i + TT * j + TT * TT * p_idx];
            if (!Rcpp::NumericVector::is_na(v))
              v /= sc;
          }
      }
      
    } else if (scale == "i") {
      
      std::vector<double> vals;
      vals.reserve(TT * TT);
      
      for (int i = 0; i < TT; ++i)
        for (int j = 0; j < TT; ++j) {
          double v = out[i + TT * j + TT * TT * p_idx];
          if (!Rcpp::NumericVector::is_na(v))
            vals.push_back(v);
        }
        
        double iq = iqr_median_halves(vals);
      
      if (iq > 0 && R_finite(iq)) {
        double sf = 1.35 / iq;
        
        for (int i = 0; i < TT; ++i)
          for (int j = 0; j < TT; ++j) {
            double &v = out[i + TT * j + TT * TT * p_idx];
            if (!Rcpp::NumericVector::is_na(v))
              v *= sf;
          }
      }
    }
  }
  
  return out;
}





