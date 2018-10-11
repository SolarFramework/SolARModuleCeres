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

//#define USE_FREE


#include <iostream>
#include <string>
#include <vector>

#include <boost/log/core.hpp>
// ADD COMPONENTS HEADERS HERE

#include "SolARModuleOpencv_traits.h"
#include "SolARModuleOpengl_traits.h"
#include "SolARModuleNonFreeOpencv_traits.h"
#include "SolARModuleCeres_traits.h"
#include "SolARModuleTools_traits.h"

#include "xpcf/xpcf.h"

#include "api/input/devices/ICamera.h"
#include "api/features/IKeypointDetector.h"
#include "api/features/IDescriptorsExtractor.h"
#include "api/features/IDescriptorMatcher.h"
#include "api/solver/pose/I3DTransformFinderFrom2D2D.h"
#include "api/solver/map/ITriangulator.h"
#include "api/solver/map/IMapper.h"
#include "api/solver/map/IMapFilter.h"
#include "api/solver/pose/I2D3DCorrespondencesFinder.h"
#include "api/solver/pose/I3DTransformFinderFrom2D3D.h"
#include "api/features/IMatchesFilter.h"
#include "api/display/ISideBySideOverlay.h"
#include "api/display/I2DOverlay.h"
#include "api/display/I3DOverlay.h"
#include "api/display/IImageViewer.h"
#include "api/display/I3DPointsViewer.h"

#include "api/image/IImageLoader.h"
#include "api/solver/map/IBundler.h"

#include "SolAROpenCVHelper.h"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"


using namespace SolAR;
using namespace SolAR::datastructure;
using namespace SolAR::api;
using namespace SolAR::MODULES::OPENCV;
using namespace SolAR::MODULES::NONFREEOPENCV;
using namespace SolAR::MODULES::OPENGL;
using namespace SolAR::MODULES::CERES;
using namespace SolAR::MODULES::TOOLS;

namespace xpcf = org::bcom::xpcf;

struct streamConfig{

  std::string m_dir;
  int m_viewsNo;

  void show(){
      std::cout<<"  ->stream configuration: "<<std::endl;
      std::cout<<"      # views directory: "<<m_dir<<std::endl;
            std::cout<<"      # views number: "<<m_viewsNo<<std::endl;
  }

  bool load(std::string&path){
      std::ifstream ox(path);
      if(!ox.is_open()){
          std::cerr<<" can't read stream config file from: "<<path<<std::endl;
          return false;
      }
      else{
          std::string v[2];
          for(int i = 0; i < 2; ++i)
              ox>>v[i];
          m_dir = v[0];
          m_viewsNo = std::stoi(v[1]);
          show();

          return true;
      }
  }
  streamConfig(){
      m_dir = ""; m_viewsNo = 0;
  }
  streamConfig(std::string a, int b): m_dir(a), m_viewsNo(b){}
};

void basic_mapFiltering( std::vector<SRef<CloudPoint>>&in, double thresh, double&reproj_err){
     std::vector<SRef<CloudPoint>> out;
     out.reserve(in.size());
    reproj_err = 0.0;
    for(const auto &p: in){
        if(p->getReprojError()<thresh){
            out.push_back(p);
            reproj_err+= p->getReprojError();
        }
    }
    in = out;
    reproj_err/=(double)out.size();
}
int run_bundle(){

    // stream config need data folder (ask for it).
    std::string path_stream = "stream_config.txt";
    streamConfig* mStream = new streamConfig();
    mStream->load(path_stream);
    LOG_ADD_LOG_TO_CONSOLE();

    /* instantiate component manager*/
    /* this is needed in dynamic mode */
    SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
    if(xpcfComponentManager->load("bundle_config.xml")!=org::bcom::xpcf::_SUCCESS)
    {
        LOG_ERROR("Failed to load the configuration file bundle_config.xml")
        return -1;
    }
    LOG_INFO("Start creating components");
    auto camera =xpcfComponentManager->create<SolARCameraOpencv>()->bindTo<input::devices::ICamera>();
#ifdef USE_FREE
    auto keypointsDetector =xpcfComponentManager->create<SolARKeypointDetectorOpencv>()->bindTo<features::IKeypointDetector>();
    auto descriptorExtractor =xpcfComponentManager->create<SolARDescriptorsExtractorAKAZE2Opencv>()->bindTo<features::IDescriptorsExtractor>();
#else
   auto  keypointsDetector = xpcfComponentManager->create<SolARKeypointDetectorNonFreeOpencv>()->bindTo<features::IKeypointDetector>();
   auto descriptorExtractor = xpcfComponentManager->create<SolARDescriptorsExtractorSURF64Opencv>()->bindTo<features::IDescriptorsExtractor>();
#endif

    auto matcher =xpcfComponentManager->create<SolARDescriptorMatcherKNNOpencv>()->bindTo<features::IDescriptorMatcher>();
    auto poseFinderFrom2D2D =xpcfComponentManager->create<SolARPoseFinderFrom2D2DOpencv>()->bindTo<solver::pose::I3DTransformFinderFrom2D2D>();
    auto mapper =xpcfComponentManager->create<SolARSVDTriangulationOpencv>()->bindTo<solver::map::ITriangulator>();
    auto poseGraph =xpcfComponentManager->create<SolARMapper>()->bindTo<solver::map::IMapper>();
    auto viewer3DPoints =xpcfComponentManager->create<SolAR3DPointsViewerOpengl>()->bindTo<display::I3DPointsViewer>();
    auto bundler =xpcfComponentManager->create<SolARBundlerCeres>()->bindTo<api::solver::map::IBundler>();

    auto  mapFilter =xpcfComponentManager->create<SolARMapFilter>()->bindTo<solver::map::IMapFilter>();



    std::vector<SRef<Image>>                            views;
    SRef<Keyframe>                                      keyframe[2];
    std::vector<SRef<Keypoint>>                         keypoints[2];
    SRef<DescriptorBuffer>                              descriptors[2];
    std::vector<DescriptorMatch>                        matches;
    std::vector<SRef<CloudPoint>>                       cloud, filtredCloud;
    std::vector<Transform3Df>                           keyframePoses;

    SRef<Map> map;
    views.resize(mStream->m_viewsNo);
    cv::Mat cvView;
    LOG_INFO("loading views: ");
    int c = 0;
    for(unsigned int i = 0; i < mStream->m_viewsNo; ++i){
        std::string path_temp = mStream->m_dir + std::to_string(i) + ".png";
        cvView = cv::imread(path_temp);
        if(!cvView.empty()) ++c;
        SolAROpenCVHelper::convertToSolar(cvView, views[i]);
    }
    std::cout<<c <<"/"<<mStream->m_viewsNo<<" loaded correctly"<<std::endl;


    keyframePoses.resize(2);
    poseFinderFrom2D2D->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());
    mapper->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());

    LOG_DEBUG("Intrincic parameters : \n {}", camera->getIntrinsicsParameters());

    if (camera->start() != FrameworkReturnCode::_SUCCESS) // videoFile
    {
        LOG_ERROR("Camera cannot start");
        return -1;
    }

    for(unsigned int v = 0; v < 2; ++v){
        keypointsDetector->detect(views[v], keypoints[v]);
        descriptorExtractor->extract(views[v], keypoints[v], descriptors[v]);
    }
    matcher->match(descriptors[0], descriptors[1], matches);

    int nbOriginalMatches = matches.size();
    keyframePoses[0] = Transform3Df::Identity();

    // Estimate the pose of of the second frame (the first frame being the reference of our coordinate system)
    poseFinderFrom2D2D->estimate(keypoints[0], keypoints[1], keyframePoses[0], keyframePoses[1], matches);

    LOG_INFO("Nb matches for triangulation: {}\\{}", matches.size(), nbOriginalMatches);
    LOG_INFO("Estimate pose of the camera for the frame 2: \n {}", keyframePoses[1].matrix());

    // Triangulate
    double reproj_error = mapper->triangulate(keypoints[0],
                                              keypoints[1],
                                              matches,
                                              std::make_pair(0, 1),
                                              keyframePoses[0],
                                              keyframePoses[1],
                                              cloud);


    for(unsigned int v = 0; v < 2; ++v){
        keyframe[v] = xpcf::utils::make_shared<Keyframe>(keypoints[v],
                                                         descriptors[v],
                                                         views[v],
                                                         keyframePoses[v]);
    }

    mapFilter->filter(keyframePoses[0], keyframePoses[1], cloud, filtredCloud);


    poseGraph->update(map, keyframe[0]);
    poseGraph->update(map,
                       keyframe[1],
                       filtredCloud,
                       matches);


    std::vector<int>selectedKeyframes = {0,1};

    std::vector<SRef<CloudPoint>>cloud_before_ba;
    std::vector<SRef<CloudPoint>>cloud_after_ba;

    cloud_before_ba = *poseGraph->getMap()->getPointCloud();

    std::vector<SRef<Keyframe>> kf=poseGraph->getKeyframes();
    bundler->adjustBundle(kf,
                          *poseGraph->getMap()->getPointCloud(),
                          camera->getIntrinsicsParameters(),
                          camera->getDistorsionParameters(),
                          selectedKeyframes);



    cloud_after_ba = *poseGraph->getMap()->getPointCloud();

    std::vector<Transform3Df>KeyframePoses_after;
    Transform3Df p0,p1;

    p0 = (poseGraph->getKeyframes()[0]->getPose());
    p1 = (poseGraph->getKeyframes()[1]->getPose());

    KeyframePoses_after.push_back(p0);
    KeyframePoses_after.push_back(p1);

    std::cout<<" after: "<<cloud_after_ba.size()<<std::endl;

    std::vector<float>color_noba = {1.0,0.0,0.0}; // color for cloud before
    std::vector<float>color_withba = {0.0,1.0,0.0}; // color for cloud after

    // here the drawing cameras!

    while(true){
        if (viewer3DPoints->displayCloudsAndPoses(cloud_before_ba,
                                                  cloud_after_ba,
                                                  keyframePoses,
                                                  KeyframePoses_after,
                                                  color_noba,
                                                  color_withba) == FrameworkReturnCode::_STOP){
            return 0;
        }
    }

    return 0;
}
int main(){
  run_bundle();
  return 0;
}



