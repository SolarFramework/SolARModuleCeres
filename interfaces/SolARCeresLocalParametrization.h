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


