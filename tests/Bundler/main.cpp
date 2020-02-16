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


#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <boost/log/core.hpp>

#include "xpcf/xpcf.h"
#include "api/display/I3DPointsViewer.h"
#include "api/solver/map/IBundler.h"
#include "core/Log.h"

using namespace SolAR;
using namespace SolAR::datastructure;
using namespace SolAR::api;
namespace xpcf = org::bcom::xpcf;


///@brief: Bundle problem loader struct:
/// Loads:
///   a) load 2D points.
///   b) load 3D points and visibility (to estbalih 2D/3D correspondances).
///   c) load camera extrinsics.
///   d) load camera intrinsics.
///Shows : all the loaded parameters.
///
struct SolARBALoader{
	std::vector<SRef<Keyframe>>m_keyframes;
	std::vector<std::vector<Keypoint>>m_points2d;
	std::vector<Transform3Df>m_poses;
    std::vector<CloudPoint>m_points3d;
    std::vector<SRef<Image>> m_views;
    std::vector<SRef<DescriptorBuffer>> m_descriptors;
	CamCalibration  m_intrinsic;
	CamDistortion   m_distorsion;

	const std::vector<CloudPoint> get3DPoints() {
		return m_points3d;
	}

	const std::vector<Keypoint> get2DPoints(int i) {
		return m_points2d[i];
	}


	const std::vector<SRef<Keyframe>> getKeyframes() {
		return m_keyframes;
	}

	const SRef<Keyframe> &getKeyframe(int i) {
		return m_keyframes[i];
	}

	const CamCalibration getCamCalibration() {
		return m_intrinsic;
	}
	const CamDistortion getCamDistorsion() {
		return m_distorsion;
	}

	const std::vector<Transform3Df> getPoses() {
		return m_poses;
	}
	void fillKeyframes() {
		m_keyframes.resize(m_poses.size());
		m_descriptors.resize(m_poses.size());
		m_views.resize(m_poses.size());
		for (unsigned int i = 0; i < m_keyframes.size(); ++i) {
			m_keyframes[i] = xpcf::utils::make_shared<Keyframe>(get2DPoints(i),
				m_descriptors[i],
				m_views[i],
				m_poses[i]);

		}
	}

    bool load2DPoints(const std::string & path_measures) {
        int N;
        std::ifstream ox(path_measures);
        if (!ox.is_open()) {
            std::cerr << " can't read measurements file from: " << path_measures << std::endl;
            return false;
        }
        else {
            std::cout<<" LOADING 2D POINTS: ";
            ox >> N;
            m_points2d.resize(N);
            for (int i = 0; i < N; ++i) {
                std::cout<<i<<" ";
                std::string path_measure;
                ox >> path_measure;
                std::filesystem::path path2DPoints(path_measures);
                std::ifstream ox(path2DPoints.parent_path().string() + path_measure);
                if (!ox.is_open()) {
                    std::cerr << " can't find observation file from: " << path_measure << std::endl;
                    return false;
                }
                else {
                    int kp_no;
                    ox >> kp_no;
                    m_points2d[i].resize(kp_no);
                    for (int j = 0; j < kp_no; ++j) {
                        float x,y;
                        ox >>x;
                        ox >>y;
                        m_points2d[i][j] = Keypoint(x,y,0.0,0.0,0.0,0.0,0);
                    }
                }   
            }
            std::cout<<" done"<<std::endl;
            return true;
        }
    }
    bool load3DPoints(const std::string & path_obs) {
        std::ifstream ox(path_obs);
        if (!ox.is_open()) {
            std::cerr << "can't find cloud from: " << path_obs << std::endl;
            return false;
        }
        else{
            std::cout<<"LOADING 3D POINTS: ";
            int obs_no;
            ox >> obs_no;
            m_points3d.resize(obs_no);
            for (int i = 0; i < obs_no; ++i) {
                double x,y,z;
                ox >> x;
                ox >> y;
                ox >> z;

                std::map<unsigned int, unsigned int> visibility_temp;

               m_points3d[i] = CloudPoint(x, y, z,0.0,0.0,0.0,0.0,visibility_temp);
               int viz_no; ox >> viz_no;
               for(int j = 0; j < viz_no; ++j) {
                   int idxView,idxLoc;
                    ox >>idxView;
                    ox >>idxLoc;
                    m_points3d[i].visibilityAddKeypoint(idxView, idxLoc);
                }
            }
            std::cout<<" done"<<std::endl;
        }
        return true;
        
    }
    
    bool loadIntrinsic(const std::string&path_calib){
        std::cout<<"loading intrinsics: ";
        std::ifstream ox(path_calib);
        if(ox.is_open()){
           for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    ox>>m_intrinsic(i,j);
                }
            }
           std::cout<<"done"<<std::endl;
        }
        else{
            LOG_INFO("can't read calirabtion file from", path_calib);
           return false;
        }
        return true;
    }
    bool loadDistorsions(const std::string&path_dist) {
        std::ifstream ox(path_dist);
        if (!ox.is_open()) {
            LOG_INFO("can't read distorsion file from", path_dist);
            return false;
        }
        else {
            std::cout<<"loading intrinsic: ";
                for (int i = 0; i < 5; ++i) {
                    ox >> m_distorsion[i];
            }
                std::cout<<" done"<<std::endl;
        }
        return true;
    }

    bool loadExtrinsics(const std::string & path_poses) {
        std::ifstream ox(path_poses);
        if (!ox.is_open()) {
            std::cerr << "can't find poses file from: " << path_poses << std::endl;
            return false;
        }
        else{
            std::cout<<" loading poses: ";
            int N;
            ox >> N;
            m_poses.resize(N);
            for (unsigned int i = 0; i < N; ++i) {
                for (int ii = 0; ii < 3; ++ii) {
                    for (int jj = 0; jj < 4; ++jj) {
                        ox >> m_poses[i](ii, jj);
                    }
                }
                m_poses[i](3,0) = 0.0;m_poses[i](3,1) = 0.0;m_poses[i](3,2) = 0.0;m_poses[i](3,3) = 1.0;
            }
        }
        std::cout<<" done"<<std::endl;
        return true;
    }
    
    void showExtrinsics()const {
        int idx = 0;
        for (const auto &p : m_poses) {
            std::cout << " EXTRINSIC: " << idx << std::endl;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 4; ++j) {
                    std::cout << p(i, j) << " ";
                }
                std::cout << std::endl;
            }
            ++idx;
            std::cout << std::endl;
        }
    }
    
    void show3DPoints()const {
        int idx = 0;
        std::cout << "<3D POINTS>: " << std::endl;
        std::cout<<"    ->size: "<<m_points3d.size()<<std::endl;
        for (unsigned int i = 0; i < m_points3d.size(); ++i){
            std::cout << "p: " << m_points3d[i].getX() << " " <<  m_points3d[i].getY() << " " <<  m_points3d[i].getZ() << "  ";

            std::map<unsigned int, unsigned int> visibility = m_points3d[i].getVisibility();
            int idxFrame = 0;
            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it){
                std::cout<<it->first<<" "<<it->second<<" ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    void show2Dpoints()const {
        std::cout << "<2D POINTS>: " << std::endl;
        for (int i = 0; i < m_points2d.size(); ++i) {
            std::cout << "	<2D POINTS from view: " << i << ">:" << std::endl;
            for (int j = 0; j < m_points2d[i].size(); ++j) {
                std::cout << m_points2d[i][j].getX() << " " <<  m_points2d[i][j].getY() << std::endl;
            }
        }
    }
    
    void  showIntrinsics()const {
        std::cout << "<INTRINSIC>: " << std::endl;
            for (int ii = 0; ii < 3; ++ii) {
                for (int jj = 0; jj < 3; ++jj) {
                    std::cout << m_intrinsic(ii, jj) << " ";
                }
                std::cout << std::endl;
            }
    }
    
    void  showDistorsions()const {
        std::cout << "<DISTORSION>: " << std::endl;
            for (int ii = 0; ii < 5; ++ii) {
                    std::cout << m_distorsion[ii] << " ";
            }
            std::cout << std::endl;
    }
    
};


int run_bundle(std::string & scene){
    LOG_ADD_LOG_TO_CONSOLE();

#if NDEBUG
    boost::log::core::get()->set_logging_enabled(false);
#endif

    try {

        SolARBALoader *ba = new SolARBALoader();
        const std::string path_poses        = "../../" + scene + "Bundle/" + scene + "Poses.txt";
        const std::string path_points3d     = "../../" + scene + "Bundle/" + scene + "Pts3D.txt";;
        const std::string path_points2d     = "../../" + scene + "Bundle/" + scene + "Pts2D.txt";
        const std::string path_calibration  = "../../" + scene + "Bundle/" + scene + "Calibration.txt";
        const std::string path_distorison   = "../../" + scene + "Bundle/" + scene + "Distorsion.txt";

        LOG_INFO("-<SolAR BA PROBLEM LOADING>-");
        ba->load3DPoints(path_points3d);
        ba->load2DPoints(path_points2d);
        ba->loadExtrinsics(path_poses);
        ba->loadIntrinsic(path_calibration);
        ba->loadDistorsions(path_distorison);
        ba->fillKeyframes();

        const std::string path_config = "SolARCeresBundler_conf.xml";
        SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
        if(xpcfComponentManager->load(path_config.c_str())!=org::bcom::xpcf::_SUCCESS)
        {
            LOG_ERROR("Failed to load the configuration {}",path_config)
            return -1;
        }
        LOG_INFO("-<SolARBundlerCeres: >-");
        SRef < solver::map::IBundler> bundler = xpcfComponentManager->resolve<api::solver::map::IBundler>();
        LOG_INFO("-<SolAR3DPointsViewerOpengl: >-");
        SRef<display::I3DPointsViewer> viewer3DPoints = xpcfComponentManager->resolve<display::I3DPointsViewer>();

        std::vector<int>selectedKeyframes; // = { 0,1,2,3 };
        std::vector<CloudPoint>correctedCloud;
        std::vector<Transform3Df>correctedPoses;


        correctedPoses.resize(ba->getKeyframes().size());
        CamCalibration correctedCalib;
        CamDistortion correctedDist;

        double reproj_errorFinal = 0.f;
        reproj_errorFinal = bundler->solve(ba->getKeyframes(),
                                            ba->get3DPoints(),
                                            ba->getCamCalibration(),
                                            ba->getCamDistorsion(),
                                            selectedKeyframes,
                                            correctedPoses,
                                            correctedCloud,
                                            correctedCalib,
                                            correctedDist);

       LOG_INFO("reprojection error final: {}",reproj_errorFinal);

        std::vector<Transform3Df>framePoses;
        Transform3Df pp = Transform3Df::Identity();
        while (true) {
            if (viewer3DPoints->display(ba->get3DPoints(),
                pp,
                ba->getPoses(),
                framePoses,
                correctedCloud,
                correctedPoses) == FrameworkReturnCode::_STOP) {
                return 0;
            }

        }
        delete ba;
    }
    catch (xpcf::Exception e)
    {
        LOG_ERROR ("The following exception has been catched: {}", e.what());
        return -1;
    }
    return 0;
}
int main(int argc, char ** argv){
    std::string scene_name = "room15";
    run_bundle(scene_name);
    return 0;
}



