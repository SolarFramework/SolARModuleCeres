cmake_minimum_required(VERSION 3.7.2)

##################################################
set (VERSION_NUMBER "0.4.0")
project("SolARModuleNonFreeOpenCV")
set (SOURCES 
		src/SolARModuleCeres.cpp
		src/SolARBundlerCeres.cpp

	)

#install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/xpcf_SolARModuleNonFreeOpenCV_registry.xml" DESTINATION $ENV{BCOMDEVROOT}/.xpcf/SolAR/	)
#set (HEADERS
	#interfaces/SolARCeresAPI.h
	#interfaces/SolARBundlerCeres.h

#)

##################################################

# various macros
include("$ENV{BCOMDEVROOT}/bcomBuild/SolARFramework/solarmacros.cmake")
# config setup
setup()
# process packagedependencies.txt
processPackagedependencies()
# define targets (type library)
defineTargets("library" "")



