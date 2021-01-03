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

#ifndef FIXED3DNORMPARAMETRIZATION_H_
#define FIXED3DNORMPARAMETRIZATION_H_

#include <ceres/local_parameterization.h>
#include "api/solver/map/IBundler.h"

namespace xpcf = org::bcom::xpcf;

namespace SolAR {
namespace MODULES {
namespace CERES {
    class Fixed3DNormParametrization : public ceres::LocalParameterization
    {
    public:
        Fixed3DNormParametrization(double norm)
            : mFixedNorm(norm)
        {
        }
        virtual ~Fixed3DNormParametrization()
        {
        }

        virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;
        virtual bool ComputeJacobian(const double *x, double *jacobian) const;
        virtual int GlobalSize() const
        {
            return 3;
        }
        virtual int LocalSize() const
        {
            return 2;
        }

        /**
         * @brief Calculates two vectors that are orthogonal to X. It first picks a non-colinear point C then basis1=(X-C) x C and basis2=X x basis1
         */
        static void GetBasis(const double *x, double *basis1, double *basis2);

    protected:
        const double mFixedNorm;
    };
}
}
}

#endif


