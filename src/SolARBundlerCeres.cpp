#include "SolARBundlerCeres.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <utility>

#include "ceres/ceres.h"
#include "ceres/rotation.h"

using namespace std;
namespace xpcf  = org::bcom::xpcf;


XPCF_DEFINE_FACTORY_CREATE_INSTANCE(SolAR::MODULES::CERES::SolARBundlerCeres)

namespace SolAR {
    using namespace datastructure;
        namespace MODULES {
            namespace CERES {

            class BALProblem {
            public:
                ~BALProblem() {
                    delete[] point_index_;
                    delete[] camera_index_;
                    delete[] observations_;
                    delete[] parameters_;
                }
                int num_observations()       const { return num_observations_; }
                int get_parameters()       const { return num_parameters_;}
                int get_cameras()       const { return num_cameras_; }
                int get_points()       const { return num_points_; }
                const double* observations() const { return observations_; }
                double* mutable_cameras() { return parameters_; }
                double* mutable_points() { return parameters_ + 9 * num_cameras_; }
                double* mutable_camera_for_observation(int i) {
                    return mutable_cameras() + camera_index_[i] * 9;
                }
                double* mutable_point_for_observation(int i) {
                    return mutable_points() + point_index_[i] * 3;
                }
                bool LoadFile(const char* filename) {
                    FILE* fptr = fopen(filename, "r");
                    if (fptr == NULL) {
                        return false;
                    };
                    FscanfOrDie(fptr, "%d", &num_cameras_);
                    FscanfOrDie(fptr, "%d", &num_points_);
                    FscanfOrDie(fptr, "%d", &num_observations_);
                    point_index_ = new int[num_observations_];
                    camera_index_ = new int[num_observations_];
                    observations_ = new double[2 * num_observations_];
                    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
                    parameters_ = new double[num_parameters_];
                    for (int i = 0; i < num_observations_; ++i) {
                        FscanfOrDie(fptr, "%d", camera_index_ + i);
                        FscanfOrDie(fptr, "%d", point_index_ + i);
                        for (int j = 0; j < 2; ++j) {
                            FscanfOrDie(fptr, "%lf", observations_ + 2 * i + j);
                        }
                    }
                    for (int i = 0; i < num_parameters_; ++i) {
                        FscanfOrDie(fptr, "%lf", parameters_ + i);
                    }
                    return true;
                }
            private:
                template<typename T>
                void FscanfOrDie(FILE *fptr, const char *format, T *value) {
                    int num_scanned = fscanf(fptr, format, value);
                    if (num_scanned != 1) {
                        LOG(FATAL) << "Invalid UW data file.";
                    }
                }
                int num_cameras_;
                int num_points_;
                int num_observations_;
                int num_parameters_;
                int* point_index_;
                int* camera_index_;
                double* observations_;
                double* parameters_;
            };
            // Templated pinhole camera model for used with Ceres.  The camera is
            // parameterized using 9 parameters: 3 for rotation, 3 for translation, 1 for
            // focal length and 2 for radial distortion. The principal point is not modeled
            // (i.e. it is assumed be located at the image center).
            struct SnavelyReprojectionError {
                SnavelyReprojectionError(double observed_x, double observed_y)
                    : observed_x(observed_x), observed_y(observed_y) {}
                template <typename T>
                bool operator()(const T* const camera,
                    const T* const point,
                    T* residuals) const {
                    // camera[0,1,2] are the angle-axis rotation.
                    T p[3];
                    ceres::AngleAxisRotatePoint(camera, point, p);
                    // camera[3,4,5] are the translation.
                    p[0] += camera[3];
                    p[1] += camera[4];
                    p[2] += camera[5];
                    // Compute the center of distortion. The sign change comes from
                    // the camera model that Noah Snavely's Bundler assumes, whereby
                    // the camera coordinate system has a negative z axis.
                    T xp = -p[0] / p[2];
                    T yp = -p[1] / p[2];
                    // Apply second and fourth order radial distortion.
                    const T& l1 = camera[7];
                    const T& l2 = camera[8];
                    T r2 = xp*xp + yp*yp;
                    T distortion = 1.0 + r2  * (l1 + l2  * r2);
                    // Compute final projected point position.
                    const T& focal = camera[6];
                    T predicted_x = focal * distortion * xp;
                    T predicted_y = focal * distortion * yp;
                    // The error is the difference between the predicted and observed position.
                    residuals[0] = predicted_x - observed_x;
                    residuals[1] = predicted_y - observed_y;
                    return true;
                }
                // Factory to hide the construction of the CostFunction object from
                // the client code.
                static ceres::CostFunction* Create(const double observed_x,
                    const double observed_y) {
                    return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 2, 9, 3>(
                        new SnavelyReprojectionError(observed_x, observed_y)));
                }
                double observed_x;
                double observed_y;
            };


                SolARBundlerCeres::SolARBundlerCeres():ComponentBase(xpcf::toUUID<SolARBundlerCeres>())
                {
                     addInterface<IBundler>(this);
                #ifdef DEBUG
                    std::cout << " SolARBundlerCeres constructor" << std::endl;
                #endif
                }

                bool SolARBundlerCeres::adjustBundle(const std::string&path_bundle,
                                                     std::vector<SRef<CloudPoint>>& cloud_before,
                                                     std::vector<SRef<CloudPoint>>& cloud_after){
                    std::string filename = path_bundle;
                    BALProblem bal_problem;
                    if (!bal_problem.LoadFile(filename.c_str())) {
                        std::cerr << "ERROR: unable to open file " << filename << "\n";
                        return 1;
                    }
                    const double* observations = bal_problem.observations();
                    // Create residuals for each observation in the bundle adjustment problem. The
                    // parameters for cameras and points are added automatically.
                    ceres::Problem problem;
                    for (int i = 0; i < bal_problem.num_observations(); ++i) {
                        // Each Residual block takes a point and a camera as input and outputs a 2
                        // dimensional residual. Internally, the cost function stores the observed
                        // image location and compares the reprojection against the observation.
                        ceres::CostFunction* cost_function =
                            SnavelyReprojectionError::Create(observations[2 * i + 0],
                                observations[2 * i + 1]);
                        problem.AddResidualBlock(cost_function,
                            NULL /* squared loss */,
                            bal_problem.mutable_camera_for_observation(i),
                            bal_problem.mutable_point_for_observation(i));
                    }
                    // Make Ceres automatically detect the bundle structure. Note that the
                    // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
                    // for standard bundle adjustment problems.
                    ceres::Solver::Options options;

                    //options.linear_solver_type = ceres::DENSE_SCHUR;
                 //   options.minimizer_progress_to_stdout = true;

                    double * parameters = bal_problem.mutable_cameras();
                    int cameras_no = bal_problem.get_cameras();
                    int param_no = bal_problem.get_parameters();
                    int points_no = bal_problem.get_points();
                    cloud_after.resize(points_no);
                    cloud_before.resize(points_no);


                    int idx0 = 0;
                    for (int j = 0; j < cameras_no; ++j) {
                        for (int jj = 0; jj < 9; ++jj) {
                            ++idx0;
                        }
                    }

                    int idx1 = idx0;
                    for (int j = 0; j < points_no; ++j) {
                        double x = parameters[idx1];
                        ++idx1;
                       double y = parameters[idx1];
                        ++idx1;
                       double z =parameters[idx1];
                        ++idx1;
                        double reprj_err = 0;  std::vector<int>visibility = std::vector<int>(50, -1);
                        cloud_before[j] = xpcf::utils::make_shared<CloudPoint>(x, y, z,0.0,0.0,0.0,reprj_err,visibility);
                    }


                    options.linear_solver_type = ceres::ITERATIVE_SCHUR;
                    options.preconditioner_type = ceres::SCHUR_JACOBI;
                    options.dense_linear_algebra_library_type = ceres::LAPACK;
                    options.max_num_iterations = 20;
                    options.num_threads = 8;
                    options.num_linear_solver_threads = 8;
                    options.logging_type = ceres::SILENT;


                    ceres::Solver::Summary summary;
                    ceres::Solve(options, &problem, &summary);
                    std::cout << summary.FullReport() << "\n";



                    std::cout << " number of parameters: " << param_no << std::endl;
                    std::cout << " number of cameras: " << cameras_no << std::endl;
                    std::cout << " number of points: " << points_no << std::endl;

                    double * new_3dpoints = bal_problem.mutable_cameras();
                    int idx2 = 0;
                    for (int j = 0; j < cameras_no; ++j) {
                        for (int jj = 0; jj < 9; ++jj) {
                            ++idx2;
                        }
                    }

                    int idx3 = idx2;
                    for (int j = 0; j < points_no; ++j) {
                        double x = new_3dpoints[idx3];
                        ++idx3;
                       double y = new_3dpoints[idx3];
                        ++idx3;
                       double z = new_3dpoints[idx3];
                        ++idx3;
                        double reprj_err = 0;  std::vector<int>visibility = std::vector<int>(50, -1);
                        cloud_after[j] = xpcf::utils::make_shared<CloudPoint>(x, y, z,0.0,0.0,0.0,reprj_err,visibility);
                    }
                    return true;
                }

            }
       }
}
