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

#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <tf/transform_listener.h>

namespace eigen_utils {

/** Computes the pseudoinverse of a matrix via SVD 
 * @param M the input matrix
 * @param Minv[out] the output pseudoinverse matrix
 * @param tolerance Those singular values smaller than tolerance * max_singular_value are removed
 */
void pseudoinverse(const Eigen::MatrixXd &M, Eigen::MatrixXd &Minv, double tolerance = 1.e-06);

/** Computes the pseudoinverse of a matrix via SVD 
 * @param M the input matrix
 * @param tolerance Those singular values smaller than tolerance * max_singular_value are removed
 * @returns the output pseudoinverse matrix
 */
Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &M, double tolerance = 1.e-06);

/** Computes a pose vector in (t utheta) notation from an homogeneous transformation matrix
 * @param M an homogeneous transformation matrix
 * @param pose[out] an output 6x1 vector containing the translation and rotation (utheta notation) described by matrix M
 */
void transformToPoseVector(const Eigen::Affine3d &M, Eigen::VectorXd &pose);

/** Computes a pose vector in (t utheta) notation from an homogeneous transformation matrix
 * @param M an homogeneous transformation matrix
 * @returns a 6x1 vector containing the translation and rotation (utheta notation) described by matrix M
 */
Eigen::VectorXd transformToPoseVector(const Eigen::Affine3d &M);

/** Computes a rotation matrix corresponding to a utheta rotation
 * @param u a 3x1 vector containing a rotation in utheta notation 
 * @returns an homogeneous matrix containing a rotation corresponding to u
 */
Eigen::Affine3d UThetaToAffine3d(const Eigen::Vector3d &u);

/** Computes the direct exponential map for a twist acting during a time
 * @param v a 6x1 twist vector in (t utheta) notation
 * @param delta_t the time during which the twist v is applied
 * @returns an homogeneous matrix that contains the displacement after applying the velocity v during the time interval delta_t
 */
Eigen::Affine3d direct_exponential_map(const Eigen::VectorXd &v, double delta_t);

/** Queries tf for a transform and returns it as an eigen affine3d object 
 * @param listener the tf listener object
 * @param target the target frame name (t)
 * @param source the source frame name (s)
 * @param[out] tMs an homogeneous matrix describing the source frame in target coordinates
 * @param timestamp the desired timestamp
 * @param timeout the maximum time to wait for the transform
 * @returns true if the transform exists, false otherwise
 */
bool getTransform(const tf::TransformListener &listener, const std::string &target, const std::string source, Eigen::Affine3d &tMs,
                  const ros::Time &timestamp = ros::Time::now(), const ros::Duration &timeout = ros::Duration(5.0));


} // namespace

#endif
