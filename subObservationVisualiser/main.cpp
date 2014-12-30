#include "NetworkReliabilityObs.h"
#include "Context.h"
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <QApplication>
#include "subObservationVisualiser.h"
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
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
			("pointSize", boost::program_options::value<float>(), "(float) The size of graph vertices. Defaults to 0.1")
			("fromFile", boost::program_options::value<std::string>(), "(path) The file to take the first observation from")
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
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << options << std::endl;
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
		if(!variableMap.count("fromFile"))
		{
			std::cout << "Input `fromFile' is required" << std::endl;
			return 0;
		}
		std::ifstream inputStream(variableMap["fromFile"].as<std::string>().c_str(), std::ios_base::in);
		if(!inputStream)
		{
			std::cout << "Unable to open specified file" << std::endl;
			return 0;
		}
		boost::archive::text_iarchive archive(inputStream, boost::archive::no_codecvt);
		boost::shared_ptr<const Context> context;
		boost::shared_ptr<NetworkReliabilitySubObs> subObs;
		try
		{
			readNetworkReliabilitySubObs(archive, context, subObs);
		}
		catch(std::runtime_error& err)
		{
			std::cout << "Error loading SubObs from file: " << err.what() << std::endl;
			return 0;
		}
		subObservationVisualiser viewer(subObs, pointSize);
		viewer.show();
		app.exec();
		return 0;
	}
}
int main(int argc, char** argv)
{
	return networkReliability::main(argc, argv);
}
