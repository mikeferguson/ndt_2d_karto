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
  correlation_search_space_dimension_ =
    node->declare_parameter<double>(name + ".correlation_search_space_dimension", 0.3);
  correlation_search_space_resolution_ =
    node->declare_parameter<double>(name + ".correlation_search_space_resolution", 0.01);
  correlation_search_space_smear_deviation_ =
    node->declare_parameter<double>(name + ".correlation_search_space_smear_deviation", 0.03);

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
                                   size_t scan_points_to_use) const
{
  karto::LocalizedRangeScan * karto_scan = makeKartoScan(scan);

  karto::Matrix3 karto_covariance;
  karto_covariance.SetToIdentity();

  karto::Pose2 best_pose;
  karto::ScanMatcher * matcher =
    karto::ScanMatcher::Create(&params_, correlation_search_space_dimension_,
                               correlation_search_space_resolution_,
                               correlation_search_space_smear_deviation_,
                               range_max_);

  double score = matcher->MatchScan(karto_scan, candidates_, best_pose, karto_covariance);
  std::cout << score << std::endl;

  // Set correction pose
  pose.x = best_pose.GetX() - scan->pose.x;
  pose.y = best_pose.GetY() - scan->pose.y;
  pose.theta = best_pose.GetHeading() - scan->pose.theta;
  std::cout << pose.x << ", " << pose.y << ". " << pose.theta << std::endl;

  // Set covariance
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      covariance(i, j) = karto_covariance(i, j);
    }
    std::cout << covariance(i, 0) << " " << covariance(i, 1) << " " << covariance(i, 2) << std::endl;
  }

  return score;
}

double ScanMatcherKarto::scoreScan(const ScanPtr & scan) const
{


  return 0.0;
}

double ScanMatcherKarto::scoreScan(const ScanPtr & scan, const Pose2d & pose) const
{
  return 0.0;
}

double ScanMatcherKarto::scorePoints(const std::vector<Point> & points, const Pose2d & pose) const
{
  return 0.0;
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

}  // namespace ndt_2d

#include <pluginlib/class_list_macros.hpp>
PLUGINLIB_EXPORT_CLASS(ndt_2d::ScanMatcherKarto, ndt_2d::ScanMatcher)