
#include "SolARBundlerCeres.h"
#include <core/Log.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>





using namespace std;
namespace xpcf  = org::bcom::xpcf;


XPCF_DEFINE_FACTORY_CREATE_INSTANCE(SolAR::MODULES::CERES::SolARBundlerCeres)

namespace SolAR {
    using namespace datastructure;
        namespace MODULES {
            namespace CERES {

            template <typename T>
            inline void SolARRadialDistorsion(const T &focal_length_x,
                                              const T &focal_length_y,
                                              const T &principal_point_x,
                                              const T &principal_point_y,
                                              const T &k1,
                                              const T &k2,
                                              const T &k3,
                                              const T &p1,
                                              const T &p2,
                                              const T &normalized_x,
                                              const T &normalized_y,
                                              T *image_x,
                                              T *image_y)
                                        {
                T x = normalized_x;
                T y = normalized_y;

                T r2 = x * x + y * y;
                T r4 = r2 * r2;
                T r6 = r4 * r2;

                T r_coeff = 1.0 + k1 * r2 + k2 * r4 + k3 * r6;

                T xd = x * r_coeff + 2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x);
                T yd = y * r_coeff + 2.0 * p2 * x * y + p1 * (r2 + 2.0 * y * y);

                *image_x = focal_length_x * xd + principal_point_x;
                *image_y = focal_length_y * yd + principal_point_y;
            }

            struct SolARReprojectionError {
                SolARReprojectionError(double observed_x, double observed_y): observed_x(observed_x), observed_y(observed_y) {}
                template <typename T>
                bool operator()(const T* const cameraIntr,
                                const T* const cameraExtr,
                                const T* const point,
                                T* residuals) const {
                    T p[3];
                    ceres::AngleAxisRotatePoint(cameraExtr, point, p);

                    p[0] += cameraExtr[3];
                    p[1] += cameraExtr[4];
                    p[2] += cameraExtr[5];

                    T xp = +p[0] / p[2];
                    T yp = +p[1] / p[2];


                    const T& fx = cameraIntr[0];
                    const T& fy = cameraIntr[1];

                    const T& cx = cameraIntr[2];
                    const T& cy = cameraIntr[3];

                    const T& k1 = cameraIntr[4];
                    const T& k2 = cameraIntr[5];

                    const T& p1 = cameraIntr[6];
                    const T& p2 = cameraIntr[7];

                    const T& k3 = cameraIntr[8];


                    T predicted_x, predicted_y;
                    SolARRadialDistorsion(fx,
                                          fy,
                                          cx,
                                          cy,
                                          k1, k2, k3,
                                          p1, p2,
                                          xp, yp,
                                          &predicted_x,
                                          &predicted_y);

                    residuals[0] = predicted_x - observed_x;
                    residuals[1] = predicted_y - observed_y;
                    return true;
                }
                static ceres::CostFunction* create(const double observed_x,
                    const double observed_y) {
                    return (new ceres::AutoDiffCostFunction<SolARReprojectionError, 2, 9, 6, 3>(
                        new SolARReprojectionError(observed_x, observed_y)));
                }
                double observed_x;
                double observed_y;
            };
            struct ceresObserv{
                    int cIdx;
                    int pIdx;
                    Point2Df oPt;
                    void show(){
                        std::cout<<" obervation: "<<std::endl;
                        std::cout<<"    # cam idx: "<<cIdx<<" #3d: "<<pIdx<<" #2d: "<<oPt.getX()<<" "<<oPt.getY()<<std::endl;
                    }
                    ceresObserv(){

                    };
                };
                SolARBundlerCeres::SolARBundlerCeres():ConfigurableBase(xpcf::toUUID<SolARBundlerCeres>())
                {
                     addInterface<IBundler>(this);
                     declareProperty("iterationsCount", m_iterationsNo);
                     declareProperty("fixedMap", m_fixedMap);
                     declareProperty("fixedExtrinsics", m_fixedExtrinsics);
                     declareProperty("fixedIntrinsics", m_fixedIntrinsics);
                     declareProperty("fixedFirstPose", m_holdFirstPose);
                     m_parameters=NULL;
                     LOG_DEBUG(" SolARBundlerCeres constructor");
                }


                double SolARBundlerCeres::solve(std::vector<SRef<Keyframe>> & framesToAdjust,
                                                std::vector<CloudPoint> & mapToAdjust,
                                                CamCalibration &K,
                                                CamDistortion &D,
                                                const std::vector<int>&selectKeyframes){
                    LOG_DEBUG("0. INIT CERES PROBLEM");
                    double reproj_error = 0.f;
                    initCeresProblem();
                    LOG_DEBUG("ITERATIONS NO: {}", m_iterationsNo);
                    LOG_DEBUG("MAP FIXED ? {}", m_fixedMap);
                    LOG_DEBUG("EXTRINSICS FIXED ? {}", m_fixedExtrinsics);
                    LOG_DEBUG("INTRINSICS FIXED ? {}", m_fixedIntrinsics);
                    LOG_DEBUG("HOLD FIRST POSE ? {}", m_holdFirstPose);
                    LOG_DEBUG("1. FILL CERES PROBLEM");
                    fillCeresProblem(framesToAdjust,
                                            mapToAdjust,
                                            K,
                                            D,
                                            selectKeyframes);
                   LOG_DEBUG("2. SOLVE CERES PROBLEM");
                   reproj_error = solveCeresProblem();
                   LOG_DEBUG("3. UPDATE  CERES PROBLEM");
                   updateCeresProblem(framesToAdjust,
                                      mapToAdjust,
                                      K,
                                      D);
                    return reproj_error;
                }

                void SolARBundlerCeres::initCeresProblem(){
                    m_options.use_nonmonotonic_steps = true;
                    m_options.preconditioner_type = ceres::SCHUR_JACOBI;      
                    m_options.linear_solver_type = ceres::ITERATIVE_SCHUR;                                        
                    m_options.use_inner_iterations = true;
                    m_options.max_num_iterations = m_iterationsNo;
                    m_options.minimizer_progress_to_stdout = false;
                }
                void SolARBundlerCeres::fillCeresProblem(std::vector<SRef<Keyframe>> & framesToAdjust,
                                                         std::vector<CloudPoint> & mapToAdjust,
                                                         CamCalibration &K,
                                                         CamDistortion &D,
                                                         const std::vector<int>&selectedKeyframes){



                    std::vector<ceresObserv>observations_temp;
                    int mapToBundleSize = 0;
                    if(selectedKeyframes.size()> 0){
                        LOG_DEBUG("#### LOCAL BUNDLER");
                        int minVisibleViews  = 2;
                        for(int i = 0; i < mapToAdjust.size(); ++i){
                            std::map<unsigned int, unsigned int> visibility = mapToAdjust[i].getVisibility();
                            map<unsigned int, unsigned int>::iterator it;
                            int minViz = 0;
                                for (unsigned int c = 0; c < selectedKeyframes.size(); ++c){
                                    it =  visibility.find(selectedKeyframes[c]);
                                    if(it!=visibility.end()) ++ minViz;
                               }
                            if(minViz<minVisibleViews)continue;
                            ++ mapToBundleSize;
                            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it){
                                int idxCam0 = it->first;
                                for (unsigned int c = 0; c < selectedKeyframes.size(); ++c){ // seen by at least two cameras...!
                                    int idxCam1 = selectedKeyframes[c];
                                    if(idxCam0 == idxCam1){
                                        int idxPoint = i;
                                        int idxLoc = it->second;
                                        ceresObserv v;
                                        v.oPt  = Point2Df(framesToAdjust[idxCam0]->getKeypoints()[idxLoc].getX(),
                                                          framesToAdjust[idxCam0]->getKeypoints()[idxLoc].getY());
                                        v.cIdx = idxCam0;
                                        v.pIdx = idxPoint;
                                        observations_temp.push_back(v);

                                    }
                                 }
                            }
                        }
                    }
                    else{
                        for(int i = 0; i < mapToAdjust.size(); ++i){
                            std::map<unsigned int, unsigned int> visibility = mapToAdjust[i].getVisibility();
                            ++mapToBundleSize;
                            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it){
                                if(it->second  != -1){
                                    ceresObserv v;
                                    int idxCam = it->first;
                                    int idxLoc = it->second;
                                    int idxPoint = i;
                                    v.oPt  = Point2Df(framesToAdjust[idxCam]->getKeypoints()[idxLoc].getX(),
                                                      framesToAdjust[idxCam]->getKeypoints()[idxLoc].getY());
                                    v.cIdx = idxCam;
                                    v.pIdx = idxPoint;
                                    observations_temp.push_back(v);
                                }
                             }
                        }
                    }

                    LOG_DEBUG("number of points to bundle:{} ",mapToBundleSize);
                    m_observationsNo = observations_temp.size();
                    m_camerasNo      = framesToAdjust.size();
                    m_pointsNo       = mapToAdjust.size();

                    m_observations = new double[OBSERV_DIM * m_observationsNo];

                    m_pointIndex     = new int[m_observationsNo];
                    m_extrinsicIndex = new int[m_observationsNo];
                    m_intrinsicIndex = new int[m_observationsNo];


                    m_parametersNo = (EXT_DIM + INT_DIM) * m_camerasNo + POINT_DIM * m_pointsNo;
                    if(!m_parameters)
                        delete[] m_parameters;
                    m_parameters = new double[m_parametersNo];

                    for (int i = 0; i < m_observationsNo; ++i) {
                        m_extrinsicIndex[i] = observations_temp[i].cIdx;
                        m_intrinsicIndex[i] = observations_temp[i].cIdx;
                        m_pointIndex[i]     = observations_temp[i].pIdx;

                        m_observations[OBSERV_DIM*i + 0] = observations_temp[i].oPt.getX();
                        m_observations[OBSERV_DIM*i + 1] = observations_temp[i].oPt.getY();

                    }
                    for(int  i = 0; i < framesToAdjust.size(); ++i){
                            Vector3f r, t;
                            Transform3Df pose_temp = framesToAdjust[i]->getPose();

                            toRodrigues(pose_temp, r);

                            pose_temp = pose_temp.inverse();

                            t[0] = pose_temp(0, 3);
                            t[1] = pose_temp(1, 3);
                            t[2] = pose_temp(2, 3);

                            float fc = -1.0;
                            m_parameters[EXT_DIM*i + 0] = r[0] * fc;
                            m_parameters[EXT_DIM*i + 1] = r[1] * fc;
                            m_parameters[EXT_DIM*i + 2] = r[2] * fc;
                            m_parameters[EXT_DIM*i + 3] = t[0];
                            m_parameters[EXT_DIM*i + 4] = t[1];
                            m_parameters[EXT_DIM*i + 5] = t[2];


                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 0] = K(0, 0);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 1] = K(1, 1);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 2] = K(0, 2);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 3] = K(1, 2);                            
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 4] = D(0);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 5] = D(1);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 6] = D(2);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 7] = D(3);
                            m_parameters[INT_DIM*i + m_camerasNo * EXT_DIM + 8] = D(4);

                    }
                    for(int i = 0; i < mapToAdjust.size(); ++i){
                        m_parameters[POINT_DIM*i + m_camerasNo * (EXT_DIM  + INT_DIM) + 0] = mapToAdjust[i].getX();
                        m_parameters[POINT_DIM*i + m_camerasNo * (EXT_DIM  + INT_DIM) + 1] = mapToAdjust[i].getY();
                        m_parameters[POINT_DIM*i + m_camerasNo * (EXT_DIM  + INT_DIM) + 2] = mapToAdjust[i].getZ();
                    }

                }
                double SolARBundlerCeres::solveCeresProblem(){
                    ceres::Problem m_problem;
                    for (int i = 0; i < num_observations(); ++i) {
                            ceres::CostFunction* cost_function = SolARReprojectionError::create(m_observations[OBSERV_DIM * i + 0],
                                                                                                m_observations[OBSERV_DIM * i + 1]);

                            m_problem.AddResidualBlock(cost_function, NULL, mutable_intrinsic_for_observation(i),
                                                                            mutable_extrinsic_for_observation(i),
                                                                            mutable_point_for_observation(i));
                   }

                    if(m_fixedExtrinsics){
                       for (int i = 0; i < m_camerasNo; ++i)
                          m_problem.SetParameterBlockConstant(mutable_extrinsic_for_observation(i));
                    }else if(m_holdFirstPose){
                          m_problem.SetParameterBlockConstant(mutable_extrinsic_for_observation(0));
                    }
                    if(m_fixedIntrinsics){
                        for (int i = 0; i < m_camerasNo; ++i)
                                m_problem.SetParameterBlockConstant(mutable_intrinsic_for_observation(i));
                    }
                    if(m_fixedMap){
                        for (int i = 0; i <num_observations(); ++i)
                                m_problem.SetParameterBlockConstant(mutable_point_for_observation(i));
                    }

                    ceres::Solve(m_options, &m_problem, &m_summary);
                    std::cout << m_summary.FullReport() << "\n";
                    return ((double)m_summary.final_cost/(double)m_pointsNo);
                }


                void SolARBundlerCeres::updateMap(std::vector<CloudPoint> & mapToAdjust){
                    for (int j = 0; j < get_points(); ++j) {
                        double x = m_parameters[(j * 3 + 0) + ((EXT_DIM + INT_DIM) * get_cameras())];
                        double y = m_parameters[(j * 3 + 1) + ((EXT_DIM + INT_DIM) * get_cameras())];
                        double z = m_parameters[(j * 3 + 2) + ((EXT_DIM + INT_DIM) * get_cameras())];
                        double reprj_err = mapToAdjust[j].getReprojError();
                        std::map<unsigned int, unsigned int>visibility = mapToAdjust[j].getVisibility();
                        mapToAdjust[j] = CloudPoint(x, y, z,0.0,0.0,0.0,reprj_err,visibility);
                     }
                }
                void SolARBundlerCeres::updateExtrinsic(std::vector<SRef<Keyframe>>&framesToAdjust){
                    for (int j = 0; j < m_camerasNo; ++j) {
                        Vector3d r,t, f;
                        Transform3Df pose_temp;
                        r[0] = m_parameters[EXT_DIM * j + 0];
                        r[1] = m_parameters[EXT_DIM * j + 1];
                        r[2] = m_parameters[EXT_DIM * j + 2];

                        iRodrigues(r,pose_temp);

                        pose_temp(0,3) = m_parameters[EXT_DIM * j + 3];
                        pose_temp(1,3) = m_parameters[EXT_DIM * j + 4];
                        pose_temp(2,3) = m_parameters[EXT_DIM * j + 5];

                        pose_temp(3,0) = 0.0;
                        pose_temp(3,1) = 0.0;
                        pose_temp(3,2) = 0.0;
                        pose_temp(3,3) = 1.0;


                       pose_temp = pose_temp.inverse();
                       framesToAdjust[j]->setPose(pose_temp);
                    }
                }

                void SolARBundlerCeres::updateIntrinsic(CamCalibration &Knew,CamDistortion &Dnew){

                    int idx = 0;
                    Knew(0, 0)  = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 0];
                    Knew(1, 1)  = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 1];
                    Knew(0, 2)  = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 2];
                    Knew(1, 2)  = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 3];

                    Dnew(0)     = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 4];
                    Dnew(1)     = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 5];
                    Dnew(2)     = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 6];
                    Dnew(3)     = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 7];
                    Dnew(4)     = m_parameters[INT_DIM*idx + m_camerasNo * EXT_DIM + 8];
                }




                void SolARBundlerCeres::updateCeresProblem(std::vector<SRef<Keyframe>> & framesToAdjust,
                                                           std::vector<CloudPoint> & mapToAdjust,
                                                           CamCalibration &K,
                                                           CamDistortion &D){
                    if(!m_fixedMap)
                       updateMap(mapToAdjust);
                    if(!m_fixedExtrinsics)
                        updateExtrinsic(framesToAdjust);
                    if(!m_fixedIntrinsics)
                        updateIntrinsic(K,D);
                }


            }
       }
}
