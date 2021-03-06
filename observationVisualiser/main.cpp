#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include "arguments.h"
#include "argumentsMPFR.h"
#include "context.h"
#include "observationVisualiser.h"
#include <QApplication>
#include <fstream>
#if defined(_WIN32)
	#include <Windows.h>
	#if defined(_MSC_VER)
		#include "windowsConsoleOutput.h"
	#endif
#endif
namespace networkReliability
{
#if defined(_WIN32)
	void registerQTPluginDir()
	{
		static bool pluginDir = false;
		if(!pluginDir)
		{
			WCHAR pathArray[500];
			GetModuleFileNameW(NULL, pathArray, 500);
			int error = GetLastError();
			if(error != ERROR_SUCCESS) 
			{
				exit(-1);
			}
			std::wstring path(&(pathArray[0]));
			
			path.erase(std::find(path.rbegin(), path.rend(), '\\').base(), path.end());
			QApplication::addLibraryPath(QString::fromStdWString(path));
			pluginDir = true;
		}
	}
#endif
	int main(int argc, char** argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use.")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. ")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs")
			("pointSize", boost::program_options::value<float>(), "(float) The size of graph vertices. Defaults to 0.1")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("help", "Display this message");

#if defined(_WIN32) && defined(_MSC_VER)
		redirectConsoleOutput();
#endif
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			std::cerr << "Only one of gridGraph, graphFile and completeGraph can be specified" << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << options << std::endl;
			std::cout << "Only one of gridGraph, graphFile and completeGraph can be specified" << std::endl;
			return 0;
		}

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		mpfr_class probability;
		if(!readProbabilityString(variableMap, probability))
		{
			std::cout << "Please enter a single value for input `opProbability'" << std::endl;
			return 0;
		}

		context contextObj = context::emptyContext();
		if(!readContext(variableMap, contextObj, probability))
		{
			return 0;
		}
#if defined(_WIN32)
		registerQTPluginDir();
#endif
		
		float pointSize = 0.1f;
		if(variableMap.count("pointSize") >= 1)
		{
			pointSize = variableMap["pointSize"].as<float>();
		}
		
		QApplication app(argc, argv);
		observationVisualiser viewer(contextObj, randomSource, pointSize);
		viewer.show();
		app.exec();
		return 0;
	}
}
int main(int argc, char** argv)
{
	return networkReliability::main(argc, argv);
}
