#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List generate_ordPY_Cpp(int n, double theta, double alpha) {
  NumericVector x(n);
  NumericVector ms;
  NumericVector rs;
  NumericVector rs1;
  
  int k = 1;
  x[0] = 0;
  ms.push_back(1);
  rs.push_back(1);
  rs1.push_back(0);
  int m, temp_int;
  double lp1,lq1,temp_double,temp_double2;
  bool temp_bool;
  NumericVector lps, lqs;
  NumericVector temp_index,temp_index2;
  NumericVector temp_vec1,temp_vec2;
  LogicalVector temp_index_bool;
  NumericVector newms, newrs, newrs1;
  for(int i = 1; i < n; i++){
    m = i;
    lp1 = log(ms[0]-alpha)+log(ms[0]+1)-log(theta+m)-log(ms[0])+log(rs[0]/(rs[0]+1));
    if(k>1){
      temp_index = seq(2,k)-1;
      temp_vec1 = rs[temp_index];
      lp1 += sum(log(temp_vec1/(temp_vec1 + 1)));
      lps = log(temp_vec1/(temp_vec1+1));
      temp_vec1 = rs1[temp_index];
      temp_vec2 = ms[temp_index];
      lp1 += sum(log( alpha*(temp_vec1+1) + theta * temp_vec2 )) - sum(log( alpha*temp_vec1 + theta*temp_vec2 ));
      lps += log(temp_vec2-alpha) + log(alpha*temp_vec1 + theta*(temp_vec2+1)) - log(theta+m) - log(alpha*temp_vec1 + theta*temp_vec2);
      // lps length is k-1, so I need to access (j-1)-1 [one -1 because of 0-indexing, the other because it's long k-1]
      for(int j = 2; j<k; j++){
        temp_index = seq(j+1,k)-1;
        temp_vec1 = rs[temp_index];
        lps[j-2] += sum(log(temp_vec1/(temp_vec1 + 1)));
        temp_vec1 = rs1[temp_index];
        temp_vec2 = ms[temp_index];
        lps[j-2] += sum(log( alpha*(temp_vec1+1) + theta*temp_vec2 ) - log( alpha*temp_vec1 + theta*temp_vec2 ));
      }
      lps.push_front(lp1);
    } else {
      lps = lp1;
    }
    lq1 = log(alpha + theta*ms[0]) - log(theta+m) - log(ms[0]+1);
    // std::cout << "before lqs ";
    if(k>1){
      temp_index = seq(2,k)-1;
      temp_vec1 = rs[temp_index];
      lq1 += sum(log(temp_vec1/(temp_vec1 + 1)));
      
      temp_vec1 = rs1[temp_index];
      temp_vec2 = ms[temp_index];
      lq1 += sum(log( alpha*(temp_vec1+1) + theta * temp_vec2 )) - sum(log( alpha*temp_vec1 + theta*temp_vec2 ));
      lqs = log(theta + alpha * temp_vec1) - log(1+temp_vec1) - log(theta + m);
      // lqs length is k-1, so I need to access (j-1)-1 [one -1 because of 0-indexing, the other because it's long k-1]
      for(int j = 2; j<k+1; j++){
        temp_index = seq(j,k)-1;
        temp_vec1 = rs[temp_index];
        lqs[j-2] += sum(log(temp_vec1/(temp_vec1 + 1)));
        temp_vec1 = rs1[temp_index];
        temp_vec2 = ms[temp_index];
        lqs[j-2] += sum(log( alpha*(temp_vec1+1) + theta*temp_vec2 ) - log( alpha*temp_vec1 + theta*temp_vec2 ));
      }
      lqs.push_front(lq1);
    } else {
      lqs = lq1;
    }
    lqs.push_back( log(theta + alpha*rs[k-1])-log(1+rs[k-1])-log(theta+m) );
    temp_double = sum(exp(lps));
    temp_double2 = (double) Rcpp::runif(1)[0];
    temp_bool = (temp_double2 < temp_double);
    if(temp_bool){
      lps = exp(lps);
      // std::cout << i << " old ";
      temp_int = sample(k, 1, true, lps)[0]-1;
      // std::cout << temp_int << std::endl;
      x[i] = temp_int;
      ms[temp_int] += 1;
      temp_index = seq(temp_int+1,k)-1;
      temp_vec1 = rs[temp_index];
      temp_vec1 = temp_vec1+1;
      rs[temp_index] = temp_vec1;
      if(temp_int<k-1){
        temp_index = seq(temp_int+2,k)-1;
        temp_vec1 = rs1[temp_index];
        temp_vec1 = temp_vec1+1;
        rs1[temp_index] = temp_vec1;
      }
    } else {
      lqs = exp(lqs);
      // std::cout << i << " new ";
      temp_int = sample(k+1, 1, true, lqs)[0]-1;
      // std::cout << temp_int << std::endl;
      if(temp_int == k){
        x[i] = k;
        ms.push_back(1);
        rs.push_back(rs[k-1]+1);
        rs1.push_back(rs[k-1]);
      } else {
        temp_index_bool = (x >= temp_int);
        temp_vec1 = x[temp_index_bool];
        temp_vec1 = temp_vec1 + 1;
        x[temp_index_bool] = temp_vec1;
        x[i] = temp_int;
        if(temp_int == 0){
          ms.push_front(1);
          rs.push_front(0);
          temp_vec1 = rs + 1;
          rs = temp_vec1;
          if(k == 1){
            rs1 = NumericVector::create(0,1);
          } else {
            newrs1 = rs1[seq(1,temp_int+1)-1]; // 0 ... temp_int
            newrs1.push_back(rs1[temp_int]+1); // temp_int
            for(int h = temp_int+2; h <= k;h++){ // temp_int ... k-1
              newrs1.push_back(rs1[h-1]+1);
            }
            rs1 = newrs1;
          }
        } else {
          if(temp_int == k-1){
            newms = ms[seq(1,temp_int)-1];
            newms.push_back(1);
            for(int h=temp_int+1;h<=k;h++){
              newms.push_back(ms[h-1]);
            }
            ms = newms;
            
            newrs = rs[seq(1,temp_int)-1];
            newrs.push_back(rs[temp_int-1]+1);
            for(int h=temp_int+1;h<=k;h++){
              newrs.push_back(rs[h-1]+1);
            }
            rs = newrs;
            
            newrs1 = rs1[seq(1,temp_int+1)-1];
            newrs1.push_back(rs1[temp_int]+1);
            rs1 = newrs1;
          } else {
            newms = ms[seq(1,temp_int)-1];
            newms.push_back(1);
            for(int h=temp_int+1;h<=k;h++){
              newms.push_back(ms[h-1]);
            }
            ms = newms;
            
            newrs = rs[seq(1,temp_int)-1];
            newrs.push_back(rs[temp_int-1]+1);
            for(int h=temp_int+1;h<=k;h++){
              newrs.push_back(rs[h-1]+1);
            }
            rs = newrs;
            
            newrs1 = rs1[seq(1,temp_int+1)-1];
            newrs1.push_back(rs1[temp_int]+1);
            for(int h=temp_int+2;h<=k;h++){
              newrs1.push_back(rs1[h-1]+1);
            }
            rs1 = newrs1;
          }
        }
      }
      k +=1;
    }
  }
  
  List ret = List::create(Named("x") = x , _["ms"] = ms, 
                          _["rs"] = rs, _["rs1"] = rs1, _["k"] = k);
  return ret;
}

// [[Rcpp::export]]
List generate_ordPY_post_Cpp(List res, int M, double theta, double alpha) {
  NumericVector old_x = res["x"];
  int n = old_x.length();
  NumericVector x(n+M);
  for(int i = 0; i < n; i++){
    x[i] = old_x[i];
  }
  NumericVector oldms = res["ms"];
  NumericVector oldrs = res["rs"];
  NumericVector oldrs1 = res["rs1"];
  int k = oldms.length();
  NumericVector ms(k);
  NumericVector rs(k);
  NumericVector rs1(k);
  for(int h=0;h<k;h++){
    ms[h] = oldms[h];
    rs[h] = oldrs[h];
    rs1[h] = oldrs1[h];
  }
  
  int m, temp_int;
  double lp1,lq1,temp_double,temp_double2;
  bool temp_bool;
  NumericVector lps, lqs;
  NumericVector temp_index,temp_index2;
  NumericVector temp_vec1,temp_vec2;
  LogicalVector temp_index_bool;
  NumericVector newms, newrs, newrs1;
  for(int i = n; i < n+M; i++){
    m = i;
    
    lp1 = log(ms[0]-alpha)+log(ms[0]+1)-log(theta+m)-log(ms[0])+log(rs[0]/(rs[0]+1));
    if(k>1){
      temp_index = seq(2,k)-1;
      temp_vec1 = rs[temp_index];
      lp1 += sum(log(temp_vec1/(temp_vec1 + 1)));
      lps = log(temp_vec1/(temp_vec1+1));
      temp_vec1 = rs1[temp_index];
      temp_vec2 = ms[temp_index];
      lp1 += sum(log( alpha*(temp_vec1+1) + theta * temp_vec2 )) - sum(log( alpha*temp_vec1 + theta*temp_vec2 ));
      lps += log(temp_vec2-alpha) + log(alpha*temp_vec1 + theta*(temp_vec2+1)) - log(theta+m) - log(alpha*temp_vec1 + theta*temp_vec2);
      // lps length is k-1, so I need to access (j-1)-1 [one -1 because of 0-indexing, the other because it's long k-1]
      for(int j = 2; j<k; j++){
        temp_index = seq(j+1,k)-1;
        temp_vec1 = rs[temp_index];
        lps[j-2] += sum(log(temp_vec1/(temp_vec1 + 1)));
        temp_vec1 = rs1[temp_index];
        temp_vec2 = ms[temp_index];
        lps[j-2] += sum(log( alpha*(temp_vec1+1) + theta*temp_vec2 ) - log( alpha*temp_vec1 + theta*temp_vec2 ));
      }
      lps.push_front(lp1);
    } else {
      lps = lp1;
    }
    lq1 = log(alpha + theta*ms[0]) - log(theta+m) - log(ms[0]+1);
    if(k>1){
      temp_index = seq(2,k)-1;
      temp_vec1 = rs[temp_index];
      lq1 += sum(log(temp_vec1/(temp_vec1 + 1)));
      
      temp_vec1 = rs1[temp_index];
      temp_vec2 = ms[temp_index];
      lq1 += sum(log( alpha*(temp_vec1+1) + theta * temp_vec2 )) - sum(log( alpha*temp_vec1 + theta*temp_vec2 ));
      lqs = log(theta + alpha * temp_vec1) - log(1+temp_vec1) - log(theta + m);
      // lqs length is k-1, so I need to access (j-1)-1 [one -1 because of 0-indexing, the other because it's long k-1]
      for(int j = 2; j<k+1; j++){
        temp_index = seq(j,k)-1;
        temp_vec1 = rs[temp_index];
        lqs[j-2] += sum(log(temp_vec1/(temp_vec1 + 1)));
        temp_vec1 = rs1[temp_index];
        temp_vec2 = ms[temp_index];
        lqs[j-2] += sum(log( alpha*(temp_vec1+1) + theta*temp_vec2 ) - log( alpha*temp_vec1 + theta*temp_vec2 ));
      }
      lqs.push_front(lq1);
    } else {
      lqs = lq1;
    }
    lqs.push_back( log(theta + alpha*rs[k-1])-log(1+rs[k-1])-log(theta+m) );
    temp_double = sum(exp(lps));
    temp_double2 = (double) Rcpp::runif(1)[0];
    temp_bool = (temp_double2 < temp_double);
    if(temp_bool){
      lps = exp(lps);
      temp_int = sample(k, 1, true, lps)[0]-1;
      x[i] = temp_int;
      ms[temp_int] += 1;
      for(int h = temp_int; h<k;h++){ // will not start if temp_int == k-1
        rs[h] += 1;
      }
      for(int h = temp_int+1; h<k;h++){ // will not start if temp_int == k-1
        rs1[h] += 1;
      }
    } else {
      lqs = exp(lqs);
      temp_int = sample(k+1, 1, true, lqs)[0]-1;
      if(temp_int == k){
        x[i] = k;
        ms.push_back(1);
        rs.push_back(rs[k-1]+1);
        rs1.push_back(rs[k-1]);
      } else {
        temp_index_bool = (x >= temp_int);
        temp_vec1 = x[temp_index_bool];
        temp_vec1 = temp_vec1 + 1;
        x[temp_index_bool] = temp_vec1;
        x[i] = temp_int;
        
        newms = ms; // oldms
        ms.push_back(1);
        for(int h = 0; h < temp_int; h++){
          ms[h] = newms[h];
        }
        ms[temp_int] = 1;
        for(int h = temp_int; h < k; h++){
          ms[h+1] = newms[h];
        }
        
        newrs = rs; // oldms
        rs.push_back(1);
        for(int h = 0; h < temp_int; h++){
          rs[h] = newrs[h];
        }
        // rs[temp_int] = newrs[temp_int-1]+1;
        for(int h = temp_int-1; h < k; h++){
          rs[h+1] = newrs[h]+1;
        }
        
        newrs1 = rs1; // oldms
        rs1.push_back(1);
        for(int h = 0; h < temp_int+1; h++){
          rs1[h] = newrs1[h];
        }
        // rs1[temp_int+1] = newrs1[temp_int]+1;
        for(int h = temp_int; h < k; h++){
          rs1[h+1] = newrs1[h]+1;
        }
      }
      k +=1;
    }
  }
  
  List ret = List::create(Named("x") = x , _["ms"] = ms, 
                          _["rs"] = rs, _["rs1"] = rs1, _["k"] = k);
  return ret;
}
