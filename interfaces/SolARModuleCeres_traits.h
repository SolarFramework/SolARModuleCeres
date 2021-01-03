/**
 * @copyright Copyright (c) 2017 B-com http://www.b-com.com/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SOLARMODULECERES_TRAITS_H
#define SOLARMODULECERES_TRAITS_H

#include "xpcf/api/IComponentManager.h"

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
