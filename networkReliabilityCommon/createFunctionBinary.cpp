#include "createFunctionBinary.h"
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
namespace networkReliability
{
	bool createFunctionsBinary(std::string function, const std::vector<std::vector<int> >& edgeIDs, std::string& message, std::string& outputFile)
	{
		std::vector<std::string> functions;
		boost::split(functions, function, boost::is_any_of(","), boost::token_compress_on);
		if(functions.size() != edgeIDs.size())
		{
			message = "Inconsistent number of functions";
			return false;
		}
		char codeFile[] = "./codeFileXXXXXX";
		mkstemps(codeFile, 0);
		std::ofstream outputCodeStream(codeFile, std::ios::out);
		for(int i = 0; i  < functions.size(); i++)
		{
			outputCodeStream << "int compiledFunction" << i << "(long edgeMask)" << std::endl;
			outputCodeStream << "{" << std::endl;
			for(std::vector<int>::const_iterator relevantEdgeIterator = edgeIDs[i].begin(); relevantEdgeIterator != edgeIDs[i].end(); relevantEdgeIterator++)
			{
				outputCodeStream << "\tint e" << *relevantEdgeIterator << " = (edgeMask & (1 << " << *relevantEdgeIterator << ")) >> " << *relevantEdgeIterator << ";" << std::endl;
			}
			outputCodeStream << "\treturn " << functions[i] << ";" << std::endl;
			outputCodeStream << "}" << std::endl;
		}
		outputCodeStream.flush();
		outputCodeStream.close();
		char binaryFile[] = "./binaryFileXXXXXX.o";
		mkstemps(binaryFile, 2);

		std::string binaryFileString = binaryFile;
		//Compile .o file
		std::stringstream commandStream;
		commandStream << "gcc -fpic -x c -O3 -c " << codeFile << " -o " << binaryFileString;
		system(commandStream.str().c_str());
		commandStream.str("");
		commandStream.clear();
		//We don't need code file any more
		unlink(codeFile);

		char sharedObjectFile[] = "./binaryFileXXXXXX.so";
		mkstemps(sharedObjectFile, 3);
		std::string sharedObjectString = sharedObjectFile;
		//compile .so file
		commandStream << "gcc -shared -o " << sharedObjectString << " "  << binaryFileString;
		system(commandStream.str().c_str());
		//We don't need .o file any more
		unlink(binaryFileString.c_str());
		outputFile = sharedObjectString;
		return true;
	}
}
