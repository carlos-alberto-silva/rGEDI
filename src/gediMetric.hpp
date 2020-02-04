#include <Rcpp.h>

Rcpp::DataFrame createMetricsDataFrame(control *);
void writeMetricsDataFrame(dataStruct *,control *,metStruct *,int ,float *,float *, Rcpp::DataFrame);
