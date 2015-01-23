#include <boost/graph/graphml.hpp>
#include <boost/program_options.hpp>
#include <algorithm>
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridDimension", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. ")
			("extraIncrement", boost::program_options::value<int>()->default_value(0), "(int) An optional additional increment to apply to the edge_ids. Applied once for each grid column. ")
			("initialOrder", boost::program_options::value<int>()->default_value(0), "(int) The initial edge_id value. ");
		
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch (boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		int extraIncrement = variableMap["extraIncrement"].as<int>();
		int initialOrder = variableMap["initialOrder"].as<int>();
		if (variableMap.count("help") > 0)
		{
			std::cout <<
				"This program generates a grid graph of the specified dimension. \n\n"
				;
			std::cout << options << std::endl;
			return 0;
		}
		if (variableMap.count("gridDimension") == 0)
		{
			std::cout << "Please enter a value for input gridDimension" << std::endl;
			return 0;
		}
		int gridDimension = variableMap["gridDimension"].as<int>();
		if (gridDimension <= 1)
		{
			std::cout << "Please enter a value bigger than one for input gridDimension" << std::endl;
			return 0;
		}
		
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int> > graphType;
		graphType graph(gridDimension * gridDimension);
		boost::vector_property_map<float> xCoordinates(gridDimension * gridDimension), yCoordinates(gridDimension * gridDimension);
	
		int counter = initialOrder;
		for (int i = 0; i < gridDimension; i++)
		{
			for (int j = 0; j < gridDimension; j++)
			{
				xCoordinates[i + j * gridDimension] = (float)i * 100;
				yCoordinates[i + j * gridDimension] = (float)j * 100;
				if (i != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + 1 + j*gridDimension, counter++, graph);
				if (j != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + (j + 1)*gridDimension, counter++, graph);
			}
			counter += extraIncrement;
		}

		typedef boost::property_map<graphType, boost::edge_index_t>::type edgeIndexMapType;
		edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
		boost::dynamic_properties properties;
		properties.property("order", edgeIndexMap);
		properties.property("x", xCoordinates);
		properties.property("y", yCoordinates);

		boost::write_graphml(std::cout, graph, properties);
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
