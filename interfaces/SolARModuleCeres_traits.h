#ifndef SOLARMODULECERES_TRAITS_H
#define SOLARMODULECERES_TRAITS_H



namespace SolAR {
namespace MODULES {
/**
 * @namespace SolAR::MODULES::CERES
 * @brief <B>Provides a bundle adjustment component based on Ceres library: http://ceres-solver.org/</B>
 * <TT>UUID: 09f4a367-c5bf-4a9f-9f3b-42424d52f717</TT>
 *
 */
namespace CERES {

class SolARBundlerCeres;
}
}
}

XPCF_DEFINE_COMPONENT_TRAITS(SolAR::MODULES::CERES::SolARBundlerCeres,
                             "4897fc13-682c-4e95-8aba-abd9f7a17193",
                             "SolARBundlerCeres",
                             "Applies a bundle adjustment to optimize a 3D map and keyframes.")


#endif // SOLARMODULECERES_TRAITS_H
