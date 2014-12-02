#include "NetworkReliabilityObs.h"
#include "Context.h"
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <QApplication>
#include <windows.h>
#include "splittingVisualiser.h"
namespace networkReliability
{
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
	int main()
	{
		LPWSTR* argv;
		int argc;
		argv = CommandLineToArgvW(GetCommandLineW(), &argc);
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs")
			("pointSize", boost::program_options::value<float>(), "(float) The size of graph vertices. Defaults to 0.1")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("initialRadius", boost::program_options::value<int>(), "(int) The initial radius to use")
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
			LocalFree(argv);
			return -1;
		}
		LocalFree(argv);
		if(variableMap.count("help") > 0)
		{
			std::cout << options << std::endl;
			return 0;
		}

		int seed = 1;
		if(variableMap.count("seed") > 0)
		{
			seed = variableMap["seed"].as<int>();
		}

		mpfr_class probability;
		if(!readProbabilityString(variableMap, probability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, probability))
		{
			return 0;
		}

		std::string message;
		int initialRadius;
		if(!readInitialRadius(variableMap, initialRadius, message))
		{
			std::cout << message << std::endl;
			return 0;
		}

		registerQTPluginDir();
		argc = 0;
		
		float pointSize = 0.1f;
		if(variableMap.count("pointSize") >= 1)
		{
			pointSize = variableMap["pointSize"].as<float>();
		}
		
		QApplication app(argc, (char**)NULL);
		splittingVisualiser viewer(context, seed, pointSize, initialRadius);
		viewer.show();
		app.exec();
		return 0;
	}
}
int CALLBACK WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
	return networkReliability::main();
}
