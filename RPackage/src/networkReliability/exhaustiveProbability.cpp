#include "exhaustiveProbability.h"
#include "includeMPFRNetworkReliability.h"
#include "exhaustiveProbabilityImpl.h"
namespace networkReliability
{
	SEXP exhaustiveProbability(SEXP searchObject_sexp, SEXP probability_sexp)
	{
	BEGIN_RCPP
		double probability = Rcpp::as<double>(probability_sexp);
		mpfr_class probability_mpfr = probability;

		Rcpp::Function asFunction("as");
		Rcpp::S4 searchObject = Rcpp::as<Rcpp::S4>(searchObject_sexp);
		Rcpp::RObject dataSlot = searchObject.slot("data");
		bool countDisconnected = Rcpp::as<bool>(searchObject.slot("countDisconnected"));
		Rcpp::CharacterVector dataAsCharacter = asFunction(dataSlot, "character");
		std::vector<mpfr_class> countData;
		for(Rcpp::CharacterVector::iterator i = dataAsCharacter.begin(); i != dataAsCharacter.end(); i++)
		{
			countData.push_back(mpfr_class(Rcpp::as<std::string>(*i)));
		}
		return Rcpp::wrap(exhaustiveProbability(countData, probability_mpfr, countDisconnected));
	END_RCPP
	}
}
