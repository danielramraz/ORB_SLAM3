/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#include "GeometricTools.h"

#include "KeyFrame.h"
#include "PolyBase.h"
#include "PolyAbs.h"
#include "../../opencv/include/opencv2/opencv.hpp"

namespace ORB_SLAM3
{

Eigen::Matrix3f GeometricTools::ComputeF12(KeyFrame* &pKF1, KeyFrame* &pKF2)
{
    Sophus::SE3<float> Tc1w = pKF1->GetPose();
    Sophus::Matrix3<float> Rc1w = Tc1w.rotationMatrix();
    Sophus::SE3<float>::TranslationMember tc1w = Tc1w.translation();

    Sophus::SE3<float> Tc2w = pKF2->GetPose();
    Sophus::Matrix3<float> Rc2w = Tc2w.rotationMatrix();
    Sophus::SE3<float>::TranslationMember tc2w = Tc2w.translation();

    Sophus::Matrix3<float> Rc1c2 = Rc1w * Rc2w.transpose();
    Eigen::Vector3f tc1c2 = -Rc1c2 * tc2w + tc1w;

    Eigen::Matrix3f tc1c2x = Sophus::SO3f::hat(tc1c2);

    const Eigen::Matrix3f K1 = pKF1->mpCamera->toK_();
    const Eigen::Matrix3f K2 = pKF2->mpCamera->toK_();

    return K1.transpose().inverse() * tc1c2x * Rc1c2 * K2.inverse();
}

bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, 
                                Eigen::Vector3f &x_c2,
                                Eigen::Matrix<float,3,4> &Tc1w ,
                                Eigen::Matrix<float,3,4> &Tc2w , 
                                Eigen::Vector3f &x3D)
{
    //---------------- linear triangulation - orb slam 3 origin ----------------------
    /*
    Eigen::Matrix4f A;
    A.block<1,4>(0,0) = x_c1(0) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(0,0);
    A.block<1,4>(1,0) = x_c1(1) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(1,0);
    A.block<1,4>(2,0) = x_c2(0) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(0,0);
    A.block<1,4>(3,0) = x_c2(1) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(1,0);

    Eigen::JacobiSVD<Eigen::Matrix4f> svd(A, Eigen::ComputeFullV);

    Eigen::Vector4f x3Dh = svd.matrixV().col(3);

    if(x3Dh(3)==0)
        return false;

    // Euclidean coordinates
    x3D = x3Dh.head(3)/x3Dh(3);

    return true;
    */


    //---------------- poly abs triangulation - orb slam 3 michal and daniel ----------------------
     /*   
    std::vector<double> params;
    params = {x_c1(0), x_c1(1), x_c1(2), x_c2(0), x_c2(1), x_c2(2)}
    result = Triangulation::PolyAbs::PreparePolyCoeffs(params);
    
    std::cout << result << std::endl;

    if (result == 0 )
    {
        return false;
    }
    else
    {
        x3D = result;
        return true;
    }
    */

    //---------------- poly abs triangulation - orb slam 3 ----------------------
    cv::Mat P0 = cv::Mat::eye(3,4, CV_64F);
    cv::Mat P1 = cv::Mat::eye(3,4, CV_64F);

    P0 = [ Tc1w(0,0), Tc1w(0,1), Tc1w(0,2), Tc1w(0,3),
           Tc1w(1,0), Tc1w(1,1), Tc1w(1,2), Tc1w(1,3),
           Tc1w(2,0), Tc1w(2,1), Tc1w(2,2), Tc1w(2,3)];
    
    P1 = [ Tc2w(0,0), Tc2w(0,1), Tc2w(0,2), Tc2w(0,3),
           Tc2w(1,0), Tc2w(1,1), Tc2w(1,2), Tc2w(1,3),
           Tc2w(2,0), Tc2w(2,1), Tc2w(2,2), Tc2w(2,3)];

    //cv::eigen2cv(Tc1w, &P0);
    //cv::eigen2cv(Tc2w, &P1);

    cv::Point2d point1 = cv::Point2d(x_c1(0), x_c1(1));
    cv::Point2d point2 = cv::Point2d(x_c2(0), x_c2(1));

    cv::Point3d result = Triangulation::PolyAbs(P0, P1).triangulate(point1, point2);
    Eigen::Vector3f result_eigen = Eigen::Vector3f(static_cast<float>(result.x),
                                                   static_cast<float>(result.y),
                                                   static_cast<float>(result.z));

    x3D = result_eigen;
    //std::cout << result << std::endl;
    return true;

}

} //namespace ORB_SLAM
