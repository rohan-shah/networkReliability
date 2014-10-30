#include <boost/program_options.hpp>
#include "Arguments.h"
#include "NetworkReliabilityObs.h"
#include "NetworkReliabilitySubObs.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include "includeMPFR.h"
#include "Turnip.h"
namespace networkReliability
{
	std::string toString(mpfr_class number)
	{
		mp_exp_t exponent;
		char* resultCStr = mpfr_get_str(NULL, &exponent, 10, 10, number.backend().data(), MPFR_RNDN);
		std::string resultStr  = resultCStr;
		free(resultCStr);
		return resultStr.substr(0, 1) + "." + resultStr.substr(1, resultStr.size() - 1) + "e" + boost::lexical_cast<std::string>(exponent - 1);
	}
	std::string toString(double number);
	int main(int argc, char **argv)
	{
		typedef NetworkReliabilityObs::conditioning_type calculation_type;
		mpfr_class::default_precision(20);

		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("opProbability", boost::program_options::value<double>(), "(float) The probability that an edge is operational")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("initialRadius", boost::program_options::value<int>(), "(int) The initial radius to use")
			("splittingFactor", boost::program_options::value<std::vector<float> >()->multitoken(), "(float) The splitting factor to use at every step")
			("distributionsCache", boost::program_options::value<std::vector<std::string> >()->multitoken(), "(path) The paths to files that contain truncated binomial distribution data")
			("n", boost::program_options::value<int>(), "(int) The number of graphs initially generated")
			("variance", "Output variance estimate")
			("relativeError", "Output relative error")
			("parts", "Output estimates for transition probabilities")
			("usePMC", "(Flag) Use PMC for the last step")
			("help", "Display this message");

		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << options << std::endl;
			return 0;
		}

		double opProbability;
		if(!readProbability(variableMap, opProbability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbability))
		{
			return 0;
		}

		bool estimateParts = variableMap.count("parts") > 0;
		bool variance = variableMap.count("variance") > 0;
		bool relativeError = variableMap.count("relativeError") > 0;
		bool usePMC = variableMap.count("usePMC") > 0;

		std::string message;
		int initialRadius;
		if(!readInitialRadius(variableMap, initialRadius, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		if(variableMap.count("splittingFactor") != 1)
		{
			std::cout << "Please enter a value for input splittingFactor" << std::endl;
			return 0;
		}

		std::vector<float> splittingFactors = variableMap["splittingFactor"].as<std::vector<float> >();
		if(splittingFactors.size() == 1)
		{
			splittingFactors.insert(splittingFactors.end(), initialRadius - 1, splittingFactors[0]);
		}
		if(splittingFactors.size() != initialRadius)
		{
			std::cout << "Wrong number of values entered for input splittingFactor" << std::endl;
			return 0;
		}
		bool validSplittingFactors = true;

		bool integerSplittingFactors = true;
		for(std::vector<float>::const_iterator i = splittingFactors.begin(); i != splittingFactors.end(); i++)
		{
			if(*i != std::floor(*i))
			{
				integerSplittingFactors = false;
			}
			if(*i < 1) validSplittingFactors = false;
		}
		if(!validSplittingFactors)
		{
			std::cout << "Input splittingFactors must contain numbers greater than 1" << std::endl;
			return 0;
		}
		//If we have non-integer splitting factors then we need to generate integer splitting factors for each run through. 
		std::vector<int> splittingFactorsPerRun(splittingFactors.size());
		int productSplittingFactorsPerRun = 1;
		if(integerSplittingFactors)
		{
			for(int i = 0; i < splittingFactors.size(); i++)
			{
				splittingFactorsPerRun[i] = (int)splittingFactors[i];
				productSplittingFactorsPerRun *= (int)splittingFactors[i];
			}
		}

		int n;
		if(!readN(variableMap, n))
		{
			return 0;
		}

		boost::mt19937 randomSource;
		if(variableMap.count("seed") > 0)
		{
			randomSource.seed(variableMap["seed"].as<int>());
		}

		std::vector<NetworkReliabilitySubObs> observations;
		std::vector<NetworkReliabilitySubObs> nextStepObservations;

		boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum, boost::accumulators::tag::count> > zeroInitialisedAccumulator(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		//mean of getting to the next level, once we've conditioned on having enough edges.
		std::vector<boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum, boost::accumulators::tag::count> > > probabilities(initialRadius + 1, zeroInitialisedAccumulator);
		//mean conditioningProbability
		std::vector<boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum, boost::accumulators::tag::count> > > conditioningProbabilities(initialRadius + 1, zeroInitialisedAccumulator);
		const std::vector<int>& interestVertices = context.getInterestVertices();

		//working data for algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		std::size_t totalSamples = 0;

		//additional working data for getRadius1ReducedGraph
		std::vector<int> edgeCounts;
		mpfr_class exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));
		std::vector<int> reducedGraphInterestVertices(interestVertices.size());

		boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > firstMomentSum(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > secondMomentSum(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		//If the usePMC flag is set, don't use splitting on the last step. Instead use PMC. 
		int finalSplittingStep;
		if (usePMC) finalSplittingStep = initialRadius - 1;
		else finalSplittingStep = initialRadius;
		TurnipInput turnipInput(randomSource, NULL, reducedGraphInterestVertices);

		for(int i = 0; i < n; i++)
		{

			//Used to work out the variance. This involves adding up the getConditioningProb() values from everything in observations, dividing by productSplittingFactorsPerRun, and then adding the value to the accumulator
			boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > perTopLevelSum(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);

			if(!integerSplittingFactors)
			{
				productSplittingFactorsPerRun = 1;
				for(int i = 0; i < splittingFactors.size(); i++)
				{
					splittingFactorsPerRun[i] = (int)splittingFactors[i];
					boost::random::bernoulli_distribution<float> additional(splittingFactors[i] - floor(splittingFactors[i]));
					if(additional(randomSource))
					{
						splittingFactorsPerRun[i]++;
					}
					productSplittingFactorsPerRun *= splittingFactorsPerRun[i];
				}
			}
			totalSamples += productSplittingFactorsPerRun;
			NetworkReliabilityObs currentObs = NetworkReliabilityObs::constructConditional(context, randomSource);
			conditioningProbabilities[0](currentObs.getConditioningProb());
			NetworkReliabilitySubObs subObs = currentObs.getSubObservation(initialRadius);
			
			if(isSingleComponent(context, subObs.getState(), components, stack, colorMap))
			{
				probabilities[0](0);
				firstMomentSum(0);
				secondMomentSum(0);
				continue;
			}
			probabilities[0](1);

			observations.clear();
			observations.push_back(std::move(subObs));

			for(int splittingLevel = 0; splittingLevel < finalSplittingStep; splittingLevel++)
			{
				nextStepObservations.clear();
				for(std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
				{
					for(int k = 0; k < splittingFactorsPerRun[splittingLevel]; k++)
					{
						NetworkReliabilityObs newObs = j->getObservation(randomSource);
						conditioningProbabilities[splittingLevel+1](newObs.getConditioningProb());
						NetworkReliabilitySubObs sub = newObs.getSubObservation(initialRadius - splittingLevel - 1);
						if(!isSingleComponent(context, sub.getState(), components, stack, colorMap))
						{
							probabilities[splittingLevel+1](newObs.getConditioningProb());
							if (splittingLevel == initialRadius - 1) perTopLevelSum(newObs.getConditioningProb());

							nextStepObservations.push_back(std::move(sub));
						}
					}
				}
				observations.swap(nextStepObservations);
			}
			if (usePMC)
			{
				for (std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
				{
					Context::internalGraph reducedGraph;
					j->getRadius1ReducedGraph(reducedGraph, edgeCounts, components, stack, colorMap);
					turnipInput.graph = &reducedGraph;
					const std::size_t nReducedVertices = boost::num_vertices(reducedGraph);
					const std::size_t nReducedEdges = boost::num_edges(reducedGraph);

					turnipInput.exponentialRates.clear();
					turnipInput.exponentialRates.resize(nReducedEdges);
					turnipInput.edges.clear();
					turnipInput.edges.resize(nReducedEdges);
					Context::internalGraph::edge_iterator current, end;
					boost::tie(current, end) = boost::edges(reducedGraph);
					for (; current != end; current++)
					{
						int edgeIndex = boost::get(boost::edge_index, reducedGraph, *current);
						turnipInput.exponentialRates[edgeIndex] = (edgeCounts[current->m_source + current->m_target * nReducedVertices] + edgeCounts[current->m_target + current->m_source * nReducedVertices]) * exponentialRate;
						turnipInput.edges[edgeIndex] = std::make_pair((int)current->m_source, (int)current->m_target);
					}
					turnipInput.n = *splittingFactorsPerRun.rbegin();
					for (int k = 0; k < interestVertices.size(); k++)
					{
						reducedGraphInterestVertices[k] = components[interestVertices[k]];
					}
					turnip(turnipInput);
					std::string part = j->getConditioningProb().str();
					part = turnipInput.estimateFirstMoment.str();
					perTopLevelSum(j->getConditioningProb() * turnipInput.n * turnipInput.estimateFirstMoment);
				}
			}
			calculation_type term = boost::accumulators::sum(perTopLevelSum) / calculation_type(productSplittingFactorsPerRun);
			firstMomentSum(term);
			secondMomentSum(term*term);
		}
		if (estimateParts)
		{
			calculation_type factor = boost::accumulators::sum(conditioningProbabilities[0]) / calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(conditioningProbabilities[0])));
			std::cout << "Step 0, conditioning on event " << toString(factor) << std::endl;
			factor = boost::accumulators::sum(probabilities[0]) / calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(conditioningProbabilities[0])));
			std::cout << "Step 0, event " << toString(factor) << std::endl;
			if (initialRadius > 0)
			{
				/*
				  The extra factor of
				  (boost::accumulators::sum(conditioningProbabilities[0]) / boost::accumulators::count(conditioningProbabilities[0]))
				  is there because the code multiplies the probabilities together as it goes, but in the formula we just want the average
				  of \P\(C_1 \lmid \lowerBoundExt = ...\).
				  Hence we have to divide by the probability of C_0 for this level. This doesn't come up again later on because
				  past this point it cancels in the numerator and denominator.
				  */
				factor = boost::accumulators::sum(conditioningProbabilities[1]) / (calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(conditioningProbabilities[1]))) * (boost::accumulators::sum(conditioningProbabilities[0]) / calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(conditioningProbabilities[0])))));
				std::cout << "Step 1, conditioning on event " << toString(factor) << std::endl;

				factor = boost::accumulators::sum(probabilities[1]) / (boost::accumulators::sum(conditioningProbabilities[1]));
				std::cout << "Step 1, event " << toString(factor) << std::endl;

				for (int i = 2; i < finalSplittingStep + 1; i++)
				{
					factor = (calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(probabilities[i - 1]))) / calculation_type(boost::lexical_cast<std::string>(boost::accumulators::count(conditioningProbabilities[i])))) * (boost::accumulators::sum(conditioningProbabilities[i]) / boost::accumulators::sum(probabilities[i - 1]));
					std::cout << "Step " << i << ", conditioning on event " << toString(factor) << std::endl;

					factor = boost::accumulators::sum(probabilities[i]) / (boost::accumulators::sum(conditioningProbabilities[i]));
					std::cout << "Step " << i << ", event " << toString(factor) << std::endl;
				}
			}
		}
		//Estimate as the product of the various parts
		calculation_type estimate = boost::accumulators::sum(probabilities[initialRadius]) / calculation_type(boost::lexical_cast<std::string>(totalSamples));
		//...but this is the one we're actually using
		estimate = boost::accumulators::sum(firstMomentSum) / n;
		//This is the variance of the distribution we're taking the empirical expectation of. Hence
		//the /n in everything below - Variance of the estimate decreases at rate 1/n
		calculation_type estimatedVariance = boost::accumulators::sum(secondMomentSum) / n - (estimate * estimate);
		if (variance)
		{
			std::cout << "Estimated variance is " << toString(estimatedVariance/n) << std::endl;
		}
		if (relativeError)
		{
			std::cout << "Estimated relative error is " << toString(boost::multiprecision::sqrt(estimatedVariance / n) / estimate) << std::endl;
		}
		std::cout << "Estimate is " << toString(boost::accumulators::sum(firstMomentSum) / n) << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}