#include "commonOptions.h"
namespace networkReliability
{
	void addGridGraph(boost::program_options::options_description& options)
	{
		options.add_options()("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile. ");
	}
	void addCompensate(boost::program_options::options_description& options)
	{
		options.add_options()("compensateResampling", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Attempt to componsate for the imperfect splitting function by resampling?");
	}
}
