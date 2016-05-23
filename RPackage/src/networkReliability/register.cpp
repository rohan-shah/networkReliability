#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
#include "crudeMC.h"
#include "approximateZeroVariance.h"
#include "exhaustiveSearch.h"
extern "C" const char* package_name = "networkReliability";
R_CallMethodDef callMethods[] = 
{
	{"crudeMC_igraph", (DL_FUNC)&networkReliability::crudeMC_igraph, 5},
	{"crudeMC_graphAM", (DL_FUNC)&networkReliability::crudeMC_graphAM, 5},
	{"crudeMC_graphNEL", (DL_FUNC)&networkReliability::crudeMC_graphNEL, 5},
	{"approximateZeroVariance_igraph", (DL_FUNC)&networkReliability::approximateZeroVariance_igraph, 5},
	{"approximateZeroVariance_graphAM", (DL_FUNC)&networkReliability::approximateZeroVariance_graphAM, 5},
	{"approximateZeroVariance_graphNEL", (DL_FUNC)&networkReliability::approximateZeroVariance_graphNEL, 5},
	{"exhaustiveSearch_igraph", (DL_FUNC)&networkReliability::exhaustiveSearch_igraph, 3},
	{"exhaustiveSearch_graphAM", (DL_FUNC)&networkReliability::exhaustiveSearch_graphAM, 3},
	{"exhaustiveSearch_graphNEL", (DL_FUNC)&networkReliability::exhaustiveSearch_graphNEL, 3},
	{NULL, NULL, 0}
};
RcppExport void R_init_networkReliability(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* packageCallMethods = callMethods;
	while(packageCallMethods->name != NULL) packageCallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, packageCallMethods);

	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
	init_Rcpp_cache();
}
