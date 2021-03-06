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

#include "SolARCeresLocalParametrization.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
// Fixed3DNormParametrization
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Fixed3DNormParametrization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace SolAR {
using namespace datastructure;
namespace MODULES {
namespace CERES {
    void Fixed3DNormParametrization::GetBasis(const double *x, double *basis1, double *basis2)
    {
        const double kThreshold = 0.1;

        //Check that the point we use is not colinear with x
        if (x[1] > kThreshold || x[1] < -kThreshold || x[2] > kThreshold || x[2] < -kThreshold)
        {
            //Use C=[1,0,0]
            basis1[0] = 0;
            basis1[1] = x[2];
            basis1[2] = -x[1];

            basis2[0] = -(x[1] * x[1] + x[2] * x[2]);
            basis2[1] = x[0] * x[1];
            basis2[2] = x[0] * x[2];
        }
        else
        {
            //Use C=[0,1,0]
            basis1[0] = -x[2];
            basis1[1] = 0;
            basis1[2] = x[0];

            basis2[0] = x[0] * x[1];
            basis2[1] = -(x[0] * x[0] + x[2] * x[2]);
            basis2[2] = x[1] * x[2];
        }
        double norm;
        norm = sqrt(basis1[0] * basis1[0] + basis1[1] * basis1[1] + basis1[2] * basis1[2]);
        basis1[0] /= norm;
        basis1[1] /= norm;
        basis1[2] /= norm;

        norm = sqrt(basis2[0] * basis2[0] + basis2[1] * basis2[1] + basis2[2] * basis2[2]);
        basis2[0] /= norm;
        basis2[1] /= norm;
        basis2[2] /= norm;

        //	cv::Matx31f xmat(x[0],x[1],x[2]);
        //	cv::Matx33f umat;
        //	cv::Matx<float,1,1> wmat;
        //	cv::Matx<float,1,1> vtmat;
        //	cv::SVDecomp(xmat, wmat, umat, vtmat, cv::SVD::MODIFY_A | cv::SVD::FULL_UV);
        //
        //	basis1[0] = umat(0,1);
        //	basis1[1] = umat(1,1);
        //	basis1[2] = umat(2,1);
        //
        //	basis2[0] = umat(0,2);
        //	basis2[1] = umat(1,2);
        //	basis2[2] = umat(2,2);
    }

    bool Fixed3DNormParametrization::Plus(const double *x, const double *delta, double *x_plus_delta) const
    {
        double basis1[3];
        double basis2[3];

        //Translation is constrained
        GetBasis(x, basis1, basis2);

        x_plus_delta[0] = x[0] + delta[0] * basis1[0] + delta[1] * basis2[0];
        x_plus_delta[1] = x[1] + delta[0] * basis1[1] + delta[1] * basis2[1];
        x_plus_delta[2] = x[2] + delta[0] * basis1[2] + delta[1] * basis2[2];

        double norm = sqrt(
            x_plus_delta[0] * x_plus_delta[0] + x_plus_delta[1] * x_plus_delta[1] + x_plus_delta[2] * x_plus_delta[2]);
        double factor = mFixedNorm / norm;
        x_plus_delta[0] *= factor;
        x_plus_delta[1] *= factor;
        x_plus_delta[2] *= factor;

        return true;
    }

    bool Fixed3DNormParametrization::ComputeJacobian(const double *x, double *jacobian) const
    {
        Transform3Df jacobian_;
        double basis1[3];
        double basis2[3];

        //Translation is special
        GetBasis(x, basis1, basis2);

        jacobian_(0, 0) = basis1[0];
        jacobian_(1, 0) = basis1[1];
        jacobian_(2, 0) = basis1[2];

        jacobian_(0, 1) = basis2[0];
        jacobian_(1, 1) = basis2[1];
        jacobian_(2, 1) = basis2[2];
        return true;
    }
}
}
}
