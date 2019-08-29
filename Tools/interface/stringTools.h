/*
methods for std::string formatting 
*/

#ifndef stringTools_h
#define stringTools_h


//include c++ library classes
#include <string>
#include <utility>

namespace stringTools{

    //remove all trailing spaces and tabs 
    std::string removeBackSpaces( const std::string& s );

	//remove all leading spaces and tabs
    std::string removeFrontSpaces( const std::string&s );

	//remove all leading and trailing spaces and tabs
    std::string cleanSpaces( const std::string& s );

	//extract first word from string ( separated by spaces )
    std::string extractFirst( std::string& s );

    //add trailing / to directoryName if needed
    std::string formatDirectoryName( const std::string& );

    //check whether string contains substring
    bool stringContains( const std::string& s, const std::string& substring );

    //check whether string ends with substring 
    bool stringEndsWith( const std::string& s, const std::string& ending );

    //split file name and extentions
    std::pair< std::string, std::string > splitFileExtension( const std::string& );

    //return file name without extension
    std::string fileWithoutExtension( const std::string& );

    //split file base name and directory
    std::pair< std::string, std::string > splitDirectoryFileName( const std::string& );

    //return directory name 
    std::string directoryNameFromPath( const std::string& );

    //return file name 
    std::string fileNameFromPath( const std::string& );

    //convert double to string, with given precision
    std::string doubleToString(const double, const unsigned precision = 0);

    //remove occurences of substring 
    std::string removeOccurencesOf( const std::string& s, const std::string& substring );
}
#endif
