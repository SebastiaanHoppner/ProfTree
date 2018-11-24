#include <Rcpp.h>
#include <queue>
#include <stack>
#include <iostream>
using namespace Rcpp;


double B(double x, double a, double b){
  return ((R::pbeta(x, a, b, 1, 0)) * (R::beta(a, b)));
} // end B


void FilterRch(NumericVector F1, NumericVector F0, int LIP, int RIP){
  for(int i = LIP+1; i < RIP; i++){
    if(F0[i] != 0){
      if(F0[RIP]+(F0[RIP]-F0[LIP])/(F1[RIP]-F1[LIP])*(F1[i]-F1[RIP]) > F0[i]){
        F0[i] = 0;
      }
    }
  }
} // end FilterRch


int FurthestPoint(NumericVector F1, NumericVector F0, int LIP, int RIP){
  double hdist = 0.0;
  int result = -1;
  for(int i = LIP+1; i < RIP; i++){
    if(F0[i] != 0){
      double xi = F1[i];
      double xl = F1[LIP];
      double xr = F1[RIP];
      double yi = F0[i];
      double yl = F0[LIP];
      double yr = F0[RIP];
      
      double distance = std::abs((yr-yl)*xi-(xr-xl)*yi+xr*yl-yr*xl) / sqrt(pow(yr-yl,2) + pow(xr-xl,2));
      if(distance > hdist){
        hdist = distance;
        result = i;
      }
    }
  }
  return result;
} // end FurthestPoint


IntegerVector order_asc(NumericVector v){
  typedef std::pair<double, int> Elt;
  std::priority_queue< Elt, std::vector<Elt>, std::greater<Elt> > pq;
  std::vector<int> result;
  
  for (int i = 0; i != v.size(); ++i){
    pq.push(Elt(v[i], i));
  }
  
  result.reserve(pq.size());
  while (!pq.empty()){
    result.push_back(pq.top().second + 1);
    pq.pop();
  }
  
  return wrap(result);
} // end order_asc


// [[Rcpp::export]]
List ROCconvexhull(NumericVector scores, IntegerVector classes){
  int posLabel = max(classes);
  int negLabel = min(classes);
  
  IntegerVector scoresOrderDes = order_asc(-scores);
  IntegerVector classesSorted = classes[scoresOrderDes-1];
  
  LogicalVector tmpPos1 = (classesSorted == posLabel);
  IntegerVector tmpPos2 = as<IntegerVector>(tmpPos1);
  LogicalVector tmpNeg1 = (classesSorted == negLabel);
  IntegerVector tmpNeg2 = as<IntegerVector>(tmpNeg1);
  
  IntegerVector tp = cumsum(tmpPos2);
  IntegerVector fp = cumsum(tmpNeg2);
  
  // remove fp & tp for duplicated predictions
  NumericVector scoresSortDes = scores[scoresOrderDes-1]; // scores are sorted in descending order
  LogicalVector dups = rev(duplicated(rev(scoresSortDes)));
  tp = tp[!dups];
  fp = fp[!dups];
  tp.push_front(0); // first element of tp must be 0
  fp.push_front(0); // first element of fp must be 0
  
  double n0 = sum(classes == posLabel);
  double n1 = classes.size()-n0;
  double pi0 = n0 / (n0 + n1);
  double pi1 = n1 / (n0 + n1);
  
  NumericVector F1 = as<NumericVector>(fp) / n1; // F1 = x-values of ROC curve
  NumericVector F0 = as<NumericVector>(tp) / n0; // F0 = y-values of ROC curve
  
  // keep only finite values
  LogicalVector finite_bool = is_finite(F1) & is_finite(F0);
  F1 = F1[finite_bool];
  F0 = F0[finite_bool];
  
  if(F1.size() < 2){
    List out;   // empty list
    std::cout << "Not enough distinct predictions to compute ROC convex hull.";
    return out;
  }
  
  // keep only points on the convex hull
  // add the point (0,0) and (1,1)
  F1.push_front(0.0);
  F0.push_front(0.0);
  F1.push_back(1.0);
  F0.push_back(1.0);
  
  int LP = 0; // index left point
  int RP = F1.size()-1; // index right point
  IntegerVector rchPoints;
  rchPoints.push_back(LP);
  rchPoints.push_back(RP);
  
  std::stack<IntegerVector> intervals;
  intervals.push(IntegerVector::create(LP,RP));
  
  while(!intervals.empty()){
    IntegerVector interval = intervals.top();
    intervals.pop();
    int LIP = interval[0];
    int RIP = interval[1];
    
    if(F1[LIP] != F1[RIP]){
      FilterRch(F1, F0, LIP, RIP);
      int NRchPoint = FurthestPoint(F1, F0, LIP, RIP);
      
      if(NRchPoint != -1){
        intervals.push(IntegerVector::create(LIP,NRchPoint));
        intervals.push(IntegerVector::create(NRchPoint,RIP));
        rchPoints.push_back(NRchPoint);
      }
    }
  }
  
  F1 = F1[rchPoints];
  F0 = F0[rchPoints];
  
  // keep only convex hull points above the diagonal, except (0,0) and (1,1)
  LogicalVector indUpperTriangle = (F1 < F0);
  F1 = F1[indUpperTriangle];
  F0 = F0[indUpperTriangle];
  F1.push_front(0.0);
  F0.push_front(0.0);
  F1.push_back(1.0);
  F0.push_back(1.0);
  // sort remaining points by ascending x(=F1) value
  IntegerVector F1OrderAsc = order_asc(F1);
  F1 = F1[F1OrderAsc-1];
  F0 = F0[F1OrderAsc-1];
  
  List out(4);
  out[0] = F1;
  out[1] = F0;
  out[2] = pi1;
  out[3] = pi0;
  
  return out;
} // end ROCconvexhull


// [[Rcpp::export]]
double empChurnCpp(NumericVector scores, IntegerVector classes, double alpha, double beta, double clv, double d, double f){
  List temp = ROCconvexhull(scores, classes);
  NumericVector F1 = temp[0];
  NumericVector F0 = temp[1];
  double pi1 = temp[2];
  double pi0 = temp[3];
  
  double delta = d/clv;
  double phi = f/clv;
  
  NumericVector diffF1 = diff(F1);
  NumericVector diffF0 = diff(F0);
  NumericVector gammaTemp = ( pi1*(delta+phi)*diffF1+pi0*phi*diffF0 )/( pi0*(1-delta)*diffF0 );
  
  NumericVector gamma;
  gamma.push_back(0);
  NumericVector contr0;
  NumericVector contr1;
  
  int index = 0;
  for(int i = 0; i < gammaTemp.size(); i++){
    if(gammaTemp[i] < 1){
      gamma.push_back(gammaTemp[i]);
      contr0.push_back(( clv*(1-delta)*pi0*F0[index])                          * (B(gamma[index+1],alpha+1,beta)-B(gamma[index],alpha+1,beta))/B(1,alpha,beta));
      contr1.push_back((-clv*(phi*pi0*F0[index] + (delta+phi)*pi1*F1[index]))  * (B(gamma[index+1],alpha  ,beta)-B(gamma[index],alpha,  beta))/B(1,alpha,beta));
      index += 1;
    }
  }
  gamma.push_back(1);
  contr0.push_back(( clv*(1-delta)*pi0*F0[index])                          * (B(gamma[index+1],alpha+1,beta)-B(gamma[index],alpha+1,beta))/B(1,alpha,beta));
  contr1.push_back((-clv*(phi*pi0*F0[index] + (delta+phi)*pi1*F1[index]))  * (B(gamma[index+1],alpha  ,beta)-B(gamma[index],alpha,  beta))/B(1,alpha,beta));
  
  return sum(contr0 + contr1); // = EMP
} // end empChurnCpp
