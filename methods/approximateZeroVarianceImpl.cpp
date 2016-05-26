#include "approximateZeroVarianceImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
namespace networkReliability
{
	void approximateZeroVariance(approximateZeroVarianceArgs& args)
	{
		std::size_t n = args.n;
		if(n < 1)
		{
			throw std::runtime_error("Input n must be positive");
		}
		mpfr_class opProbability = args.contextObj.getOperationalProbability();
		mpfr_class inopProbability = 1 - opProbability;
		const std::size_t nEdges = args.contextObj.getNEdges();
		
		if(args.contextObj.getInterestVertices().size() > 2 && args.optimiseMinCut)
		{
			throw std::runtime_error("Can only specify option optimiseMinCut with 2-terminal reliability");
		}

		boost::random::uniform_01<float,float> uniformReal;
		//Vector used for mincut calculations
		std::vector<int> state(2*nEdges);
		//Sum over all the n simulations
		args.estimateFirstMoment = 0;
		args.estimateSecondMoment = 0;
		//likelihood ratio of current term
		mpfr_class currentLikelihoodRatio = 0;

		//Cache powers of the inopProbability 
		boost::scoped_array<mpfr_class> cachedInopPowers(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			cachedInopPowers[i] = boost::multiprecision::pow(inopProbability, i);
		}
		//Get out the vector that holds the flow
		const std::vector<int>& residualCapacityVector = args.contextObj.edgeResidualCapacityVector;
		//Cache the two mincut calculations that are the same every time - The first two. 
		mpfr_class minCutUpProb, minCutDownProb;
		std::fill(state.begin(), state.end(), 1);
		state[0] = state[1] = 0;
		int initialMinCutSizeDown = args.contextObj.getMinCut(state);
		minCutDownProb = cachedInopPowers[initialMinCutSizeDown];

		state[0] = state[1] = HIGH_CAPACITY;
		int initialMinCutSizeUp = args.contextObj.getMinCut(state);
		//This branch is only taken if there are only two vertices of interest, and the first edge directly connects them.
		if(initialMinCutSizeUp >= HIGH_CAPACITY)
		{
			minCutUpProb = 0;
		}
		else minCutUpProb = cachedInopPowers[initialMinCutSizeUp];

		mpfr_class qTilde;
		mpfr_class initialQTilde = inopProbability*minCutDownProb;
		initialQTilde = initialQTilde / (initialQTilde  + opProbability * minCutUpProb);
		for(std::size_t i = 0; i < n; i++)
		{
			//we might break out of this loop early because the value of the indicator function is already known.
			bool exitedEarly = false;
			int indicatorValue;

			currentLikelihoodRatio = 0;
			std::fill(state.begin(), state.end(), 1);

			float random = uniformReal(args.randomSource);
			std::size_t edgeCounter = 0;
			if(random < initialQTilde)
			{
				state[0] = state[1] = 0; 
				currentLikelihoodRatio = inopProbability / initialQTilde;
				if(initialMinCutSizeDown == 0)
				{
					indicatorValue = 1;
					exitedEarly = true;
				}
			}
			else 
			{
				state[0] = state[1] = HIGH_CAPACITY; 
				currentLikelihoodRatio = (1-inopProbability) / (1-initialQTilde);
				if(initialMinCutSizeUp >= HIGH_CAPACITY)
				{
					indicatorValue = 0;
					exitedEarly = true;
				}
			}
			if(!exitedEarly)
			{
				bool canReuseMinCut = false;
				int minCutSizeDown, minCutSizeUp;
				for(edgeCounter = 1; edgeCounter < nEdges; edgeCounter++)
				{
					state[2*edgeCounter] = state[2*edgeCounter+1] = 0;
					if(!canReuseMinCut || !args.optimiseMinCut) minCutSizeDown = args.contextObj.getMinCut(state);
					minCutDownProb = cachedInopPowers[minCutSizeDown];

					state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY;
					if(!canReuseMinCut || !args.optimiseMinCut) minCutSizeUp = args.contextObj.getMinCut(state);
					if(minCutSizeUp >= HIGH_CAPACITY)
					{
						minCutUpProb = 0;
					}
					else minCutUpProb = cachedInopPowers[minCutSizeUp];

					qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					float random = uniformReal(args.randomSource);
					if(random < qTilde)
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = 0; 
						currentLikelihoodRatio *= inopProbability / qTilde;
						if(minCutSizeDown == 0)
						{
							indicatorValue = 1;
							break;
						}
					}
					else 
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY; 
						currentLikelihoodRatio *= opProbability / (1-qTilde);
						if(minCutSizeUp >= HIGH_CAPACITY)
						{
							indicatorValue = 0;
							break;
						}
					}
					canReuseMinCut = (residualCapacityVector[2*edgeCounter] == 1) && (residualCapacityVector[2*edgeCounter+1] == 1);
				}
			}
			args.estimateFirstMoment += indicatorValue * currentLikelihoodRatio;
			args.estimateSecondMoment += indicatorValue * currentLikelihoodRatio * currentLikelihoodRatio;
		}
		args.estimateFirstMoment /= args.n;
		args.estimateSecondMoment /= args.n;
	}
}
