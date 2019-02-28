cmake_minimum_required(VERSION 3.7.2)

##################################################
project("SolARTestModuleCeres")
set (SOURCES main.cpp)
##################################################

# various macros
include("$ENV{BCOMDEVROOT}/bcomBuild/SolARFramework/solarmacros.cmake")
# config setup
setup()
# process packagedependencies.txt
processPackagedependencies()

# define the list of files to copy to build directory
if (UNIX)
        set (LIBPREFIX "lib")
        set (LIBGLUT "")
        set (LIBEXTENSION "so")
        set (OPENCV_VERSION "")
        set (FBOW_EXT "")
endif(UNIX)
if (WIN32)
        set (LIBPREFIX "")
        set (LIBGLUT "free")
        set (LIBEXTENSION "dll")
        set (OPENCV_VERSION "343")
        set (FBOW_EXT "001")
endif(WIN32)
set (BUILDCONFIG $<$<CONFIG:Debug>:debug>$<$<NOT:$<CONFIG:Debug>>:release>)

set(FILES_TO_COPY
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_filesystem.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_system.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_timer.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_log.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_chrono.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_thread.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_date_time.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_context.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_fiber.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/boost/1.68.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}boost_unit_test_framework.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/ceres/1.14.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}ceres.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/glog/0.3.5/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}glog.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/gflags/2.2.1/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}gflags.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/freeglut/3.0.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}${LIBGLUT}glut.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/thirdParties/xpcf/2.0.1/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}xpcf.${LIBEXTENSION}
    $ENV{BCOMDEVROOT}/bcomBuild/SolARFramework/0.4.0/lib/x86_64/shared/${BUILDCONFIG}/${LIBPREFIX}SolARFramework.${LIBEXTENSION}

        )
# define targets (library, install and uninstall)
defineTargets("executable" "${FILES_TO_COPY}")

