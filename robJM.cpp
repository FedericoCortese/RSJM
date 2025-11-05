#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <set>
using namespace Rcpp;

// helper to compute median of a std::vector<double>
// includes the diagonal zero entry as in your R code
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
      // leave wcd(ci-1, *) as zeros
      continue;
    }
    
    // for each dimension p
    for (int p = 0; p < P; ++p) {
      // build the within‐cluster distance matrix for dimension p
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
      // average of the medians
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
  
  // 2) branch on whether medoids is provided
  if (medoids.isNull()) {
    // ——— full Y vs Y: return T×T symmetric matrix ———
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
    // ——— Y vs medoids: return T×K matrix ———
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
    
    // 2) Assign each point to nearest centroid
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
  
  // --- compute c from quantile if not provided ---
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
    s[0] = idx + 1;  // store 1-based
  }
  
  // subsequent t = 1..T-1
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
