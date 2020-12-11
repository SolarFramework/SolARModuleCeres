
#include "SolARBundlerCeres.h"
#include "SolARCeresLocalParametrization.h"
#include <core/Log.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>

using namespace std;
namespace xpcf = org::bcom::xpcf;


XPCF_DEFINE_FACTORY_CREATE_INSTANCE(SolAR::MODULES::CERES::SolARBundlerCeres)

namespace SolAR {
using namespace datastructure;
using namespace api::storage;
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
		SolARReprojectionError(double observed_x, double observed_y) : observed_x(observed_x), observed_y(observed_y) {}
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

		double squaredError() {	
			double res[2];
			operator()(intr, ext, pt, res);
			return res[0] * res[0] + res[1] * res[1];
		}

        double observed_x;
        double observed_y;
		double * intr;
		double * ext;
		double * pt;
    };
    struct ceresObserv {
        int cIdx;
        int pIdx;
        Point2Df oPt;
        void show() {
            std::cout << " obervation: " << std::endl;
            std::cout << "    # cam idx: " << cIdx << " #3d: " << pIdx << " #2d: " << oPt.getX() << " " << oPt.getY() << std::endl;
        }
        ceresObserv() {

        };
    };
    SolARBundlerCeres::SolARBundlerCeres() :ConfigurableBase(xpcf::toUUID<SolARBundlerCeres>())
    {
        addInterface<IBundler>(this);
        declareInjectable<IPointCloudManager>(m_pointCloudManager);
        declareInjectable<IKeyframesManager>(m_keyframesManager);
		declareInjectable<ICovisibilityGraph>(m_covisibilityGraph);
        declareProperty("iterationsCount", m_iterationsNo);
        declareProperty("fixedMap", m_fixedMap);
        declareProperty("fixedKeyframes", m_fixedKeyframes);
        declareProperty("fixedIntrinsics", m_fixedIntrinsics);
        declareProperty("fixedFirstPose", m_fixedFirstPose);
        declareProperty("fixedNeighbourKeyframes", m_fixedNeighbourKeyframes);
        declareProperty("nbMaxFixedKeyframes", m_nbMaxFixedKeyframes);
		declareProperty("useSpanningTree", m_useSpanningTree);
        LOG_DEBUG(" SolARBundlerCeres constructor");
    }

	SolARBundlerCeres::~SolARBundlerCeres()
	{
		LOG_DEBUG(" SolARBundlerCeres destructor")
	}

	FrameworkReturnCode SolARBundlerCeres::setMapper(const SRef<api::solver::map::IMapper>& map)
	{
		map->getPointCloudManager(m_pointCloudManager);
		map->getKeyframesManager(m_keyframesManager);
		map->getCovisibilityGraph(m_covisibilityGraph);
		return FrameworkReturnCode::_SUCCESS;
	}

	double SolARBundlerCeres::bundleAdjustment(CamCalibration & K, CamDistortion & D, const std::vector<uint32_t> & selectedKeyframes) {

        // Init Ceres Problem
        // ------------------

        LOG_DEBUG("0. INIT CERES PROBLEM");
        ceres::Solver::Options options;
        options.use_nonmonotonic_steps = true;
        options.preconditioner_type = ceres::SCHUR_JACOBI;
        options.linear_solver_type = ceres::ITERATIVE_SCHUR;
        options.use_inner_iterations = true;
        options.max_num_iterations = m_iterationsNo;
        options.minimizer_progress_to_stdout = false;

        LOG_DEBUG("ITERATIONS NO: {}", m_iterationsNo);
        LOG_DEBUG("MAP FIXED ? {}", (bool)m_fixedMap);
        LOG_DEBUG("EXTRINSICS FIXED ? {}", (bool)m_fixedKeyframes);
        LOG_DEBUG("INTRINSICS FIXED ? {}", (bool)m_fixedIntrinsics);
        LOG_DEBUG("FIRST POSE FIXED ? {}", (bool)m_fixedFirstPose);
        LOG_DEBUG("NEIGHBOUR KEYFRAMES FIXED ? {}", (bool)m_fixedNeighbourKeyframes);


        // Fill Ceres Problem
        // ------------------
        LOG_DEBUG("1. FILL CERES PROBLEM");
        vector<ceresObserv> observationsTemp;

        // Local KeyFrames
        std::vector< SRef<Keyframe>> localKeyframes;
        std::map<uint32_t, uint32_t> solar2CeresKeyframeID;

        // Local KeyPoints
        std::vector<SRef<CloudPoint>> localCloudPoints;
        std::map<uint32_t, uint32_t> solar2CeresCloudPointID;
        std::map<uint32_t, uint32_t> ceres2SolarCloudPointID;

		/// Local bundle adjustment
        if (selectedKeyframes.size() > 0)
        {
            // Fill the local keyframeid set
            for (const auto& skf_id : selectedKeyframes) {
                SRef<Keyframe> localKeyframe;
                m_keyframesManager->getKeyframe(skf_id, localKeyframe);
                int ceresKeyframeID = (int)localKeyframes.size();
                solar2CeresKeyframeID.insert(std::pair<uint32_t, uint32_t>(skf_id, ceresKeyframeID));
                localKeyframes.push_back(localKeyframe);
                const std::map<uint32_t, uint32_t>& mapPointVisibility = localKeyframe->getVisibility();

                for (auto const &it_kpId2cpIdVisibility : mapPointVisibility) {
                    SRef<CloudPoint> localCloudPoint;
                    m_pointCloudManager->getPoint(it_kpId2cpIdVisibility.second, localCloudPoint);
                    // Check if the solar point has already been added the set of cloud point
                    std::map<uint32_t, uint32_t>::iterator itCeresCPIdCorresp = solar2CeresCloudPointID.find(it_kpId2cpIdVisibility.second);
                    uint32_t ceresCloudPointID;
                    if (itCeresCPIdCorresp==solar2CeresCloudPointID.end())
                    {
                        //insert its id in the SolAR/Ceres cloud point id correspondences map
                        ceresCloudPointID = (uint32_t)localCloudPoints.size();
                        solar2CeresCloudPointID.insert(std::pair<uint32_t, uint32_t>(it_kpId2cpIdVisibility.second, ceresCloudPointID));
						ceres2SolarCloudPointID.insert(std::pair<uint32_t, uint32_t>(ceresCloudPointID, it_kpId2cpIdVisibility.second));

                        //add the point cloud to the vector of PointCloud
                        localCloudPoints.push_back(localCloudPoint);
                    }
                    else
                        ceresCloudPointID = itCeresCPIdCorresp->second;

                    //create and add the observation
                    ceresObserv observTemp;
                    observTemp.cIdx = ceresKeyframeID;
                    observTemp.pIdx = ceresCloudPointID;
                    Keypoint kp = localKeyframe->getKeypoint(it_kpId2cpIdVisibility.first);
                    observTemp.oPt = Point2Df(kp.getX(), kp.getY());
                    observationsTemp.push_back(observTemp);
                }
            }
        }
		/// Global bundle adjustment based on spanning tree
		else if (m_useSpanningTree) {
			// get all keyframes
			std::vector<SRef<Keyframe>> keyframes;
			m_keyframesManager->getAllKeyframes(keyframes);
			for (auto localKeyframe : keyframes) {
				int ceresKeyframeID = (int)localKeyframes.size();
				solar2CeresKeyframeID.insert(std::pair<uint32_t, uint32_t>(localKeyframe->getId(), ceresKeyframeID));
				localKeyframes.push_back(localKeyframe);
			}

			// get the maximal spanning tree
			std::vector<std::tuple<uint32_t, uint32_t, float>> edgesSpanningTree;
			float totalWeights;
			m_covisibilityGraph->maximalSpanningTree(edgesSpanningTree, totalWeights);

			// get cloud points belong to maximal spanning tree	
			std::set<uint32_t> idxCloudPoints;
			for (const auto &edge : edgesSpanningTree) {
				SRef<Keyframe> kf1, kf2;
				m_keyframesManager->getKeyframe(std::get<0>(edge), kf1);
				m_keyframesManager->getKeyframe(std::get<1>(edge), kf2);
				std::map<uint32_t, uint32_t> kf1_visibilites = kf1->getVisibility();
				std::map<uint32_t, uint32_t> kf2_visibilites = kf2->getVisibility();
				std::map<uint32_t, int> countNbSeenCP;
				// get common cloud points of two keyframes
				for (const auto &it : kf1_visibilites)
					countNbSeenCP[it.second]++;
				for (const auto &it : kf2_visibilites)
					countNbSeenCP[it.second]++;
				for (const auto &it : countNbSeenCP)
					if (it.second >= 2)
						idxCloudPoints.insert(it.first);
			}
			for (const auto &it : idxCloudPoints) {
				SRef<CloudPoint> localCloudPoint;
				m_pointCloudManager->getPoint(it, localCloudPoint);
				uint32_t ceresCloudPointID = (uint32_t)localCloudPoints.size();
				solar2CeresCloudPointID.insert(std::pair<uint32_t, uint32_t>(it, ceresCloudPointID));
				ceres2SolarCloudPointID.insert(std::pair<uint32_t, uint32_t>(ceresCloudPointID, it));
				localCloudPoints.push_back(localCloudPoint);

				std::map<uint32_t, uint32_t> cpVisibilities = localCloudPoint->getVisibility();
				for (const auto &vis : cpVisibilities) {
					//create and add the observation
					ceresObserv observTemp;
					observTemp.cIdx = solar2CeresKeyframeID[vis.first];
					observTemp.pIdx = ceresCloudPointID;
					Keypoint kp = localKeyframes[solar2CeresKeyframeID[vis.first]]->getKeypoint(vis.second);
					observTemp.oPt = Point2Df(kp.getX(), kp.getY());
					observationsTemp.push_back(observTemp);
				}
			}
		}
		/// Global bundle adjustment on all keyframe and all point cloud
		else
        {
            // Fill the local keyframeid set
            vector<SRef<Keyframe>> keyframes;
            m_keyframesManager->getAllKeyframes(keyframes);
            for (auto localKeyframe : keyframes) {
                int ceresKeyframeID = (int)localKeyframes.size();
                solar2CeresKeyframeID.insert(std::pair<uint32_t, uint32_t>(localKeyframe->getId(), ceresKeyframeID));
                localKeyframes.push_back(localKeyframe);
                const std::map<uint32_t, uint32_t>& mapPointVisibility = localKeyframe->getVisibility();

                for (auto const &it_kpId2cpIdVisibility : mapPointVisibility) {
                    SRef<CloudPoint> localCloudPoint;
                    m_pointCloudManager->getPoint(it_kpId2cpIdVisibility.second, localCloudPoint);
                    // Check if the solar point has already been added the set of cloud point
                    std::map<uint32_t, uint32_t>::iterator itCeresCPIdCorresp = solar2CeresCloudPointID.find(it_kpId2cpIdVisibility.second);
                    uint32_t ceresCloudPointID;
                    if (itCeresCPIdCorresp==solar2CeresCloudPointID.end())
                    {
                        //insert its id in the SolAR/Ceres cloud point id correspondences map
                        ceresCloudPointID = (uint32_t)localCloudPoints.size();
                        solar2CeresCloudPointID.insert(std::pair<uint32_t, uint32_t>(it_kpId2cpIdVisibility.second, ceresCloudPointID));
						ceres2SolarCloudPointID.insert(std::pair<uint32_t, uint32_t>(ceresCloudPointID, it_kpId2cpIdVisibility.second));
                        //add the point cloud to the vector of PointCloud
                        localCloudPoints.push_back(localCloudPoint);
                    }
                    else
                        ceresCloudPointID = itCeresCPIdCorresp->second;

                    //create and add the observation
                    ceresObserv observTemp;
                    observTemp.cIdx = ceresKeyframeID;
                    observTemp.pIdx = ceresCloudPointID;
                    Keypoint kp = localKeyframe->getKeypoint(it_kpId2cpIdVisibility.first);
                    observTemp.oPt = Point2Df(kp.getX(), kp.getY());
                    observationsTemp.push_back(observTemp);
                }
            }
        }
        uint32_t nbLocalKeyframes = (uint32_t) localKeyframes.size();

        // Fixed Keyframes. Keyframes in the neighbourhood of the Local keyframes that see Local MapPoints but that are not Local Keyframes
        std::map<uint32_t, uint32_t> solar2CeresNeighbourKeyframeID;
        vector<SRef<Keyframe>> neighbourKeyframes;
        vector<ceresObserv> neighbourObservationsTemp;

        if (m_fixedNeighbourKeyframes && (selectedKeyframes.size() > 0))
        {
            // For al cloudPoint to optimize, search for keyframes viewing them, but not is the selected keyframes
            for (auto const cp : localCloudPoints)
            {
                const std::map<uint32_t, uint32_t> &kpVisibility = cp->getVisibility();
                for (auto const &it_kfIDkpIDVisibility : kpVisibility)
                {
                    uint32_t ceresNeighbourKFid;
                    if (solar2CeresKeyframeID.find(it_kfIDkpIDVisibility.first) == solar2CeresKeyframeID.end())
                     {
                        ceresNeighbourKFid = nbLocalKeyframes + (uint32_t)solar2CeresNeighbourKeyframeID.size();
                        solar2CeresNeighbourKeyframeID.insert(std::pair<uint32_t, uint32_t>(it_kfIDkpIDVisibility.first, ceresNeighbourKFid));

                        SRef<Keyframe> neighbourKF;
                        m_keyframesManager->getKeyframe(it_kfIDkpIDVisibility.first, neighbourKF);
                        //create and add the observation
                        std::map<uint32_t, uint32_t>::iterator it_solar2CeresCPid = solar2CeresCloudPointID.find(cp->getId());
                        if (it_solar2CeresCPid != solar2CeresCloudPointID.end())
                        {
                            ceresObserv neighbourObservTemp;
                            neighbourObservTemp.cIdx = ceresNeighbourKFid;
                            neighbourObservTemp.pIdx = it_solar2CeresCPid->second;
                            Keypoint kp = neighbourKF->getKeypoint(it_kfIDkpIDVisibility.second);
                            neighbourObservTemp.oPt = Point2Df(kp.getX(), kp.getY());

                            neighbourObservationsTemp.push_back(neighbourObservTemp);
                        }
                    }
                }
            }
        }

        // Fill in Local Keyframe observations
        uint32_t nbObservations = (uint32_t)observationsTemp.size() + (uint32_t) neighbourObservationsTemp.size();
        double* observations = new double[OBSERV_DIM * nbObservations];
        int*    pointIndex = new int[nbObservations];
        int*    extrinsicIndex = new int[nbObservations];
        int*    intrinsicIndex = new int[nbObservations];

        // Add Local observations to Ceres
        for (uint32_t i = 0; i < (uint32_t) observationsTemp.size(); i++)
        {
            extrinsicIndex[i] = observationsTemp[i].cIdx;
            intrinsicIndex[i] = observationsTemp[i].cIdx;
            pointIndex[i] = observationsTemp[i].pIdx;
            observations[OBSERV_DIM*i + 0] = observationsTemp[i].oPt.getX();
            observations[OBSERV_DIM*i + 1] = observationsTemp[i].oPt.getY();
        }

        uint32_t nbLocalObservations = (uint32_t) observationsTemp.size();
        uint32_t nbNeighbourObservations = (uint32_t) neighbourObservationsTemp.size();
        // Add also neighbourhoud observations to Ceres if active
        if (m_fixedNeighbourKeyframes && (selectedKeyframes.size() > 0))
        {
            for (uint32_t i = 0; i < (uint32_t) neighbourObservationsTemp.size(); i++)
            {
                extrinsicIndex[nbLocalObservations+i] = neighbourObservationsTemp[i].cIdx;
                intrinsicIndex[nbLocalObservations+i] = neighbourObservationsTemp[i].cIdx;
                pointIndex[nbLocalObservations+i] = neighbourObservationsTemp[i].pIdx;
                observations[OBSERV_DIM*(nbLocalObservations+i) + 0] = neighbourObservationsTemp[i].oPt.getX();
                observations[OBSERV_DIM*(nbLocalObservations+i) + 1] = neighbourObservationsTemp[i].oPt.getY();
            }
        }


        // Fill the table of parameters to optimize by Ceres (Keyframe poses and intrinsics + cloud points)
        uint32_t nbNeighbourKeyframes = (uint32_t) neighbourKeyframes.size();
        uint32_t nbKeyframes = nbLocalKeyframes + nbNeighbourKeyframes;

        uint32_t nbCloudPoint = (uint32_t) localCloudPoints.size();
        uint32_t parametersSize = (EXT_DIM + INT_DIM) * nbKeyframes + POINT_DIM * nbCloudPoint;
        double* parameters =  new double[parametersSize];

        // Fill local Keyframes
        for (uint32_t i = 0; i < nbLocalKeyframes; i++)
        {
            Transform3Df kfpose = localKeyframes[i]->getPose();
            Vector3f r,t;

            toRodrigues(kfpose,r);
            kfpose = kfpose.inverse();

            t[0] = kfpose(0, 3);
            t[1] = kfpose(1, 3);
            t[2] = kfpose(2, 3);

            float fc = -1.0;
            parameters[EXT_DIM*i + 0] = r[0] * fc;
            parameters[EXT_DIM*i + 1] = r[1] * fc;
            parameters[EXT_DIM*i + 2] = r[2] * fc;
            parameters[EXT_DIM*i + 3] = t[0];
            parameters[EXT_DIM*i + 4] = t[1];
            parameters[EXT_DIM*i + 5] = t[2];


            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 0] = K(0, 0);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 1] = K(1, 1);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 2] = K(0, 2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 3] = K(1, 2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 4] = D(0);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 5] = D(1);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 6] = D(2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 7] = D(3);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*i + 8] = D(4);
         }

        // Fill neightbour Keyframes
        for (uint32_t i = 0; i < nbNeighbourKeyframes; i++)
        {
            Transform3Df kfpose = neighbourKeyframes[i]->getPose();
            Vector3f r,t;

            toRodrigues(kfpose,r);
            kfpose = kfpose.inverse();

            t[0] = kfpose(0, 3);
            t[1] = kfpose(1, 3);
            t[2] = kfpose(2, 3);

            float fc = -1.0;
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 0] = r[0] * fc;
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 1] = r[1] * fc;
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 2] = r[2] * fc;
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 3] = t[0];
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 4] = t[1];
            parameters[EXT_DIM*(i+nbLocalKeyframes) + 5] = t[2];


            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 0] = K(0, 0);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 1] = K(1, 1);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 2] = K(0, 2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 3] = K(1, 2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 4] = D(0);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 5] = D(1);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 6] = D(2);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 7] = D(3);
            parameters[nbKeyframes * EXT_DIM + INT_DIM*(i+nbLocalKeyframes) + 8] = D(4);
         }

        // Fill cloud point
        for (uint32_t i = 0; i < nbCloudPoint; ++i)
        {
            parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 0] = localCloudPoints[i]->getX();
            parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 1] = localCloudPoints[i]->getY();
            parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 2] = localCloudPoints[i]->getZ();
        }


        // Solve Ceres Problem
        //--------------------
        LOG_DEBUG("2. SOLVE CERES PROBLEM");
        ceres::Problem problem;
		std::vector<SolARReprojectionError*> vReprojError;
        for (uint32_t i = 0; i < nbObservations; ++i) {
			SolARReprojectionError* reprojError = new SolARReprojectionError(observations[OBSERV_DIM * i + 0],
				observations[OBSERV_DIM * i + 1]);

			ceres::CostFunction*cost_function = new ceres::AutoDiffCostFunction<SolARReprojectionError, 2, 9, 6, 3>(reprojError);		
            problem.AddResidualBlock(cost_function,
                                     NULL,
                                     parameters + nbKeyframes * EXT_DIM + intrinsicIndex[i] * INT_DIM, // mutable_intrinsic_for_observation(i),
                                     parameters + extrinsicIndex[i] * EXT_DIM, //mutable_extrinsic_for_observation(i),
                                     parameters + (EXT_DIM + INT_DIM) * nbKeyframes + pointIndex[i] * POINT_DIM); //mutable_point_for_observation(i));
			reprojError->intr = parameters + nbKeyframes * EXT_DIM + intrinsicIndex[i] * INT_DIM;
			reprojError->ext = parameters + extrinsicIndex[i] * EXT_DIM;
			reprojError->pt = parameters + (EXT_DIM + INT_DIM) * nbKeyframes + pointIndex[i] * POINT_DIM;
			vReprojError.push_back(reprojError);
        }

		LOG_DEBUG("Nb keyframes: {}", localKeyframes.size());
		LOG_DEBUG("Nb neighbor keyframes: {}", neighbourKeyframes.size());
		LOG_DEBUG("Nb cloud points: {}", localCloudPoints.size());
		LOG_DEBUG("Nb observations: {}", nbObservations);

        if (m_fixedKeyframes) {
            for (uint32_t i = 0; i < nbLocalKeyframes; ++i)
                problem.SetParameterBlockConstant(parameters + extrinsicIndex[i] * EXT_DIM);
        }
        else if (m_fixedFirstPose) {
            problem.SetParameterBlockConstant(parameters + extrinsicIndex[solar2CeresKeyframeID[0]] * EXT_DIM);
        }
        if (m_fixedNeighbourKeyframes) {
            for (uint32_t i = 0; i < nbNeighbourObservations; ++i)
                problem.SetParameterBlockConstant(parameters + extrinsicIndex[i + nbLocalObservations] * EXT_DIM);
        }
        if (m_fixedIntrinsics) {
            for (uint32_t i = 0; i < nbKeyframes; ++i)
                problem.SetParameterBlockConstant(parameters + nbKeyframes * EXT_DIM + intrinsicIndex[i] * INT_DIM);
        }
        if (m_fixedMap) {
            for (uint32_t i = 0; i < nbObservations; ++i)
                problem.SetParameterBlockConstant(parameters + (EXT_DIM + INT_DIM) * nbKeyframes + pointIndex[i] * POINT_DIM);
        }

        ceres::Solver::Summary summary;

        ceres::Solve(options, &problem, &summary);
        LOG_DEBUG("Ceres Report: {}", summary.FullReport());		
        // Update Ceres Problem
        //--------------------
        LOG_DEBUG("3. UPDATE CERES PROBLEM");
		
        // Update cloud points
        if (!m_fixedMap)
        {
            for (uint32_t i = 0; i < (uint32_t)localCloudPoints.size(); i++) {
                SRef<CloudPoint> cp = localCloudPoints[i];
                cp->setX(parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 0]);
                cp->setY(parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 1]);
                cp->setZ(parameters[nbKeyframes * (EXT_DIM + INT_DIM) + POINT_DIM*i + 2]);
            }
        }

        if (!m_fixedKeyframes)
        {
            // Update keyframes
            for (uint32_t i = 0; i < (uint32_t)localKeyframes.size(); i++)
            {
                SRef<Keyframe> kf = localKeyframes[i];
                Vector3d r, t, f;
                Transform3Df kf_pose;
                r[0] = parameters[EXT_DIM*i + 0];
                r[1] = parameters[EXT_DIM*i + 1];
                r[2] = parameters[EXT_DIM*i + 2];

                iRodrigues(r,kf_pose);

                kf_pose(0,3) = parameters[EXT_DIM*i + 3];
                kf_pose(1,3) = parameters[EXT_DIM*i + 4];
                kf_pose(2,3) = parameters[EXT_DIM*i + 5];

                kf_pose(3,0) = 0.0f;
                kf_pose(3,1) = 0.0f;
                kf_pose(3,2) = 0.0f;
                kf_pose(3,3) = 1.0f;

                kf_pose = kf_pose.inverse();

                kf->setPose(kf_pose);
            }
        }

        if (!m_fixedIntrinsics)
        {
            K(0,0) = parameters[nbKeyframes * EXT_DIM + 0];
            K(1,1) = parameters[nbKeyframes * EXT_DIM + 1];
            K(0,2) = parameters[nbKeyframes * EXT_DIM + 2];
            K(1,2) = parameters[nbKeyframes * EXT_DIM + 3];

            D(0) = parameters[nbKeyframes * EXT_DIM + 4];
            D(1) = parameters[nbKeyframes * EXT_DIM + 5];
            D(2) = parameters[nbKeyframes * EXT_DIM + 6];
            D(3) = parameters[nbKeyframes * EXT_DIM + 7];
            D(4) = parameters[nbKeyframes * EXT_DIM + 8];
        }

		// get the final mean squared reprojection error
		double errorReproj(0.0);
		for (const auto &it : vReprojError) {
			errorReproj += it->squaredError();
		}

		//Update re-projection error of point cloud
		std::map<uint32_t, std::vector<double>> projErrors;
		for (int i = 0; i < vReprojError.size(); i++) {
			uint32_t cloudPoint_id = ceres2SolarCloudPointID[pointIndex[i]];
			projErrors[cloudPoint_id].push_back(std::sqrt(vReprojError[i]->squaredError()));
		}
		for (const auto &it : projErrors) {
			SRef<CloudPoint> mapPoint;
			m_pointCloudManager->getPoint(it.first, mapPoint);
			mapPoint->setReprojError(std::accumulate(it.second.begin(), it.second.end(), 0.0) / it.second.size());
		}
        return errorReproj / nbObservations;
    }
}
}
}
