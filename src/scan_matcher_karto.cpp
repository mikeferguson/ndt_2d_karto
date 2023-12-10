/*
 * Copyright 2023 Michael Ferguson
 * All Rights Reserved
 */

#include <open_karto/Karto.h>
#include <ndt_2d_karto/scan_matcher_karto.hpp>

namespace ndt_2d
{

void ScanMatcherKarto::initialize(const std::string & name,
                                  rclcpp::Node * node, double range_max)
{
  *params_.m_pCorrelationSearchSpaceDimension =
    node->declare_parameter<double>(name + ".correlation_search_space_dimension", 0.3);
  *params_.m_pCorrelationSearchSpaceResolution =
    node->declare_parameter<double>(name + ".correlation_search_space_resolution", 0.01);
  *params_.m_pCorrelationSearchSpaceSmearDeviation =
    node->declare_parameter<double>(name + ".correlation_search_space_smear_deviation", 0.03);
  *params_.m_pDistanceVariancePenalty =
    node->declare_parameter<double>(name + ".distance_variance_penalty", 0.09);
  *params_.m_pAngleVariancePenalty =
    node->declare_parameter<double>(name + ".angle_variance_penalty", 0.1218);
  *params_.m_pFineSearchAngleOffset =
    node->declare_parameter<double>(name + ".fine_search_angle_offset", 0.003491);
  *params_.m_pCoarseSearchAngleOffset =
    node->declare_parameter<double>(name + ".course_search_angle_offset", 0.3491);
  *params_.m_pCoarseAngleResolution =
    node->declare_parameter<double>(name + ".course_angle_resolution", 0.03491);
  *params_.m_pMinimumAnglePenalty =
    node->declare_parameter<double>(name + ".minimum_angle_penalty", 0.9);
  *params_.m_pMinimumDistancePenalty =
    node->declare_parameter<double>(name + ".minimum_distance_penalty", 0.5);
  *params_.m_pUseResponseExpansion =
    node->declare_parameter<bool>(name + ".use_response_expansion", false);

  range_max_ = range_max;
  karto_matcher_ = nullptr;
}

void ScanMatcherKarto::addScans(const std::vector<ScanPtr>::const_iterator & begin,
                                const std::vector<ScanPtr>::const_iterator & end)
{
  // Scan Matcher is now invalid
  if (karto_matcher_) delete karto_matcher_;
  karto_matcher_ = nullptr;

  // Add scans to our candidates
  for (auto scan = begin; scan != end; ++scan)
  {
    karto::LocalizedRangeScan * range_scan = makeKartoScan(*scan);
    candidates_.push_back(range_scan);
  }
}

double ScanMatcherKarto::matchScan(const ScanPtr & scan, Pose2d & pose,
                                   Eigen::Matrix3d & covariance,
                                   size_t /*scan_points_to_use*/) const
{
  karto::LocalizedRangeScan * karto_scan = makeKartoScan(scan);

  karto::Matrix3 karto_covariance;
  karto_covariance.SetToIdentity();

  karto::Pose2 best_pose;
  karto::ScanMatcher * matcher =
    karto::ScanMatcher::Create(&params_,
                               params_.m_pCorrelationSearchSpaceDimension->GetValue(),
                               params_.m_pCorrelationSearchSpaceResolution->GetValue(),
                               params_.m_pCorrelationSearchSpaceSmearDeviation->GetValue(),
                               range_max_);

  double score = matcher->MatchScan(karto_scan, candidates_, best_pose, karto_covariance);
  delete matcher;
  delete karto_scan;

  // Set correction pose
  pose.x = best_pose.GetX() - scan->pose.x;
  pose.y = best_pose.GetY() - scan->pose.y;
  pose.theta = best_pose.GetHeading() - scan->pose.theta;

  // Set covariance
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      covariance(i, j) = karto_covariance(i, j);
    }
  }

  return -score;
}

double ScanMatcherKarto::scoreScan(const ScanPtr & scan) const
{
  Pose2d pose;
  return scoreScan(scan, pose);
}

double ScanMatcherKarto::scoreScan(const ScanPtr & scan, const Pose2d & pose) const
{
  karto::LocalizedRangeScan * karto_scan = makeKartoScan(scan);
  karto::Pose2 karto_pose(scan->pose.x + pose.x, scan->pose.y + pose.y,
                          scan->pose.theta + pose.theta);
  karto_scan->SetOdometricPose(karto_pose);
  karto_scan->SetCorrectedPose(karto_pose);

  karto::ScanMatcher * matcher =
    karto::ScanMatcher::Create(&params_,
                               params_.m_pCorrelationSearchSpaceDimension->GetValue(),
                               params_.m_pCorrelationSearchSpaceResolution->GetValue(),
                               params_.m_pCorrelationSearchSpaceSmearDeviation->GetValue(),
                               range_max_);

  double score = matcher->ComputeScore(karto_scan, candidates_);

  delete matcher;
  delete karto_scan;

  return -score;
}

double ScanMatcherKarto::scorePoints(const std::vector<Point> & points, const Pose2d & pose) const
{
  karto::LocalizedRangeScan * karto_scan = makeKartoScan(points);
  karto::Pose2 karto_pose(pose.x, pose.y, pose.theta);
  karto_scan->SetOdometricPose(karto_pose);
  karto_scan->SetCorrectedPose(karto_pose);

  karto::ScanMatcher * matcher =
    karto::ScanMatcher::Create(&params_,
                               params_.m_pCorrelationSearchSpaceDimension->GetValue(),
                               params_.m_pCorrelationSearchSpaceResolution->GetValue(),
                               params_.m_pCorrelationSearchSpaceSmearDeviation->GetValue(),
                               range_max_);

  double score = matcher->ComputeScore(karto_scan, candidates_);

  delete matcher;
  delete karto_scan;

  return -score;
}

void ScanMatcherKarto::reset()
{
  // Clean up matcher
  if (karto_matcher_) delete karto_matcher_;
  karto_matcher_ = nullptr;

  // Clean up candidates_
  for (auto candidate : candidates_)
  {
    delete candidate;
  }
  candidates_.clear();
}

karto::LocalizedRangeScan * ScanMatcherKarto::makeKartoScan(const ScanPtr & scan) const
{
  karto::LocalizedRangeScan * range_scan = new karto::LocalizedRangeScan(karto::Name("karto"));
  
  range_scan->raw_points.resize(scan->points.size());
  for (size_t i = 0; i < scan->points.size(); ++i)
  {
    auto & p = scan->points[i];
    range_scan->raw_points[i] = karto::Vector2<kt_double>(p.x, p.y);
  }

  karto::Pose2 karto_pose(scan->pose.x, scan->pose.y, scan->pose.theta);
  range_scan->SetOdometricPose(karto_pose);
  range_scan->SetCorrectedPose(karto_pose);

  return range_scan;
}

karto::LocalizedRangeScan * ScanMatcherKarto::makeKartoScan(const std::vector<Point> & points) const
{
  karto::LocalizedRangeScan * range_scan = new karto::LocalizedRangeScan(karto::Name("karto"));

  range_scan->raw_points.resize(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    auto & p = points[i];
    range_scan->raw_points[i] = karto::Vector2<kt_double>(p.x, p.y);
  }

  karto::Pose2 karto_pose(0.0, 0.0, 0.0);
  range_scan->SetOdometricPose(karto_pose);
  range_scan->SetCorrectedPose(karto_pose);

  return range_scan;
}

}  // namespace ndt_2d

#include <pluginlib/class_list_macros.hpp>
PLUGINLIB_EXPORT_CLASS(ndt_2d::ScanMatcherKarto, ndt_2d::ScanMatcher)
