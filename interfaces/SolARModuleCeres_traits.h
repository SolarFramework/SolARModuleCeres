#ifndef SOLARMODULECERES_TRAITS_H
#define SOLARMODULECERES_TRAITS_H



namespace SolAR {
namespace MODULES {
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
