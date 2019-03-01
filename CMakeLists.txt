cmake_minimum_required(VERSION 3.7.2)

##################################################
set (VERSION_NUMBER "0.5.1")
project("SolARModuleCeres")
set (SOURCES 
    src/SolARModuleCeres.cpp
    src/SolARBundlerCeres.cpp
)

set (HEADERS
    interfaces/SolARCeresAPI.h
    interfaces/SolARModuleCeres_traits.h
    interfaces/SolARBundlerCeres.h
)

##################################################

# various macros
include("$ENV{BCOMDEVROOT}/bcomBuild/SolARFramework/solarmacros.cmake")
# config setup
setup()
# process packagedependencies.txt
processPackagedependencies()
# define targets (type library)
defineTargets("library" "")



