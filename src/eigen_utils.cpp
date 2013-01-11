/*
 * Copyright (c) 2013, Willow Garage, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Willow Garage, Inc. nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Mario Prats
 */

#include <eigen_utils/eigen_utils.h>
#include <Eigen/SVD>

namespace eigen_utils {

void pseudoinverse(const Eigen::MatrixXd &M, Eigen::MatrixXd &Minv, double tolerance)
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfM(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::MatrixXd U = svdOfM.matrixU();
  const Eigen::MatrixXd V = svdOfM.matrixV();
  const Eigen::VectorXd S = svdOfM.singularValues();

  Eigen::VectorXd Sinv = S;
  double maxsv = 0 ;
  for (std::size_t i = 0; i < S.rows(); ++i)
    if (fabs(S(i)) > maxsv) maxsv = fabs(S(i));
  for (std::size_t i = 0; i < S.rows(); ++i)
  {
    //Those singular values smaller than a percentage of the maximum singular value are removed
    if ( fabs(S(i)) > maxsv * tolerance )
      Sinv(i) = 1.0 / S(i);
    else Sinv(i) = 0;
  }

  Minv = V * Sinv.asDiagonal() * U.transpose();
}

Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &M, double tolerance)
{
  Eigen::MatrixXd Minv;
  pseudoinverse(M, Minv, tolerance);
  return Minv;
}

void transformToTwist(const Eigen::Affine3d &M, Eigen::VectorXd &pose)
{
  pose.resize(6);

  //fill translation
  pose.matrix().block(0,0,3,1) = M.matrix().block(0,3,3,1);

  //Compute and fill rotation (theta-u convention)
  //Most of this code is adapted from ViSP vpThetaUVector::build_from(vpRotationMatrix)
  const Eigen::Matrix3d R = M.matrix().block(0, 0, 3, 3);

  double s,c,theta,sinc;
  s = (R(1,0) - R(0,1)) * (R(1,0) - R(0,1))
    + (R(2,0) - R(0,2)) * (R(2,0) - R(0,2))
    + (R(2,1) - R(1,2)) * (R(2,1) - R(1,2));
  s = sqrt(s) / 2.0;
  c = (R(0,0) + R(1,1) + R(2,2) - 1.0) / 2.0;
  theta = atan2(s,c);  /* theta in [0, PI] since s > 0 */

  // General case when theta != pi. If theta=pi, c=-1
  static const double minimum = 0.0001;
  if ( (1 + c) > minimum) // Since -1 <= c <= 1, no fabs(1+c) is required
  {
    static const double threshold = 1.0e-8;
    if (fabs(theta) < threshold) sinc = 1.0 ;
    else  sinc = (s / theta) ;

    pose(3) = (R(2,1) - R(1,2)) / (2*sinc);
    pose(4) = (R(0,2) - R(2,0)) / (2*sinc);
    pose(5) = (R(1,0) - R(0,1)) / (2*sinc);
  }
  else /* theta near PI */
  {
    if ( (R(0,0) - c) < std::numeric_limits<double>::epsilon() )
      pose(3) = 0.;
    else
      pose(3) = theta * (sqrt((R(0,0) - c) / (1 - c)));
    if ((R(2,1) - R(1,2)) < 0) pose(3) = -pose(3);

    if ( (R(1,1) - c) < std::numeric_limits<double>::epsilon() )
      pose(4) = 0.;
    else
      pose(4) = theta * (sqrt((R(1,1) - c) / (1 - c)));

    if ((R(0,2) - R(2,0)) < 0) pose(4) = -pose(4);

    if ( (R(2,2) - c) < std::numeric_limits<double>::epsilon() )
      pose(5) = 0.;
    else
      pose(5) = theta * (sqrt((R(2,2) - c) / (1 - c)));

    if ((R(1,0) - R(0,1)) < 0) pose(5) = -pose(5);
  }
}

Eigen::VectorXd transformToTwist(const Eigen::Affine3d &M) 
{
  Eigen::VectorXd twist;
  transformToTwist(M, twist);
  return twist;
}

} // namespace
