/*
 * Copyright 2023 Michael Ferguson
 * All Rights Reserved
 */

#ifndef NDT_2D_KARTO__SCAN_MATCHER_KARTO_HPP_
#define NDT_2D_KARTO__SCAN_MATCHER_KARTO_HPP_

#include <memory>
#include <string>
#include <vector>
#include <open_karto/Mapper.h>  // ScanMatcher, ScanMatcherParams
#include <ndt_2d/scan_matcher.hpp>
#include <rclcpp/rclcpp.hpp>

namespace ndt_2d
{

class ScanMatcherKarto : public ScanMatcher
{
public:
  virtual ~ScanMatcherKarto() = default;

  /**
   * @brief Initialize an NDT scan matcher instance.
   * @param name Name for ths scan matcher instance.
   * @param node Node instance to use for getting parameters.
   * @param range_max Maximum range of laser scanner.
   */
  virtual void initialize(const std::string & name,
                          rclcpp::Node * node, double range_max);

  /**
   * @brief Add scans to the internal NDT map.
   * @param begin Starting iterator of scans for NDT map building.
   * @param end Ending iterator of scans for NDT map building.
   */
  virtual void addScans(const std::vector<ScanPtr>::const_iterator & begin,
                        const std::vector<ScanPtr>::const_iterator & end);

  /**
   * @brief Match a scan against the internal NDT map.
   * @param scan Scan to match against internal NDT map.
   * @param pose The corrected pose that best matches scan to NDT map.
   * @param covariance Covariance matrix for the match.
   * @param scan_points_to_use Number of points to match from the scan.
   * @returns The likelihood score when scan is at corrected pose.
   */
  virtual double matchScan(const ScanPtr & scan, Pose2d & pose,
                           Eigen::Matrix3d & covariance,
                           size_t scan_points_to_use) const;

  /**
   * @brief Score a scan against the internal NDT map.
   * @param scan Scan to score against internal NDT map.
   */
  virtual double scoreScan(const ScanPtr & scan) const;

  /**
   * @brief Score a scan against the internal NDT map.
   * @param scan Scan to score against internal NDT map.
   * @param pose The pose of the scan within the internal NDT map.
   */
  virtual double scoreScan(const ScanPtr & scan, const Pose2d & pose) const;

  /**
   * @brief Score a set of points against the internal NDT map.
   * @param points Points to score against internal NDT map.
   * @param pose The pose of the points within the internal NDT map.
   */
  virtual double scorePoints(const std::vector<Point> & points, const Pose2d & pose) const;

  /**
   * @brief Reset the internal NDT map, removing all scans.
   */
  virtual void reset();

private:
  karto::LocalizedRangeScan * makeKartoScan(const ScanPtr & scan) const;

  mutable karto::ScanMatcherParams params_;
  karto::ScanMatcher * karto_matcher_;
  karto::LocalizedRangeScanVector candidates_;

  double correlation_search_space_dimension_;
  double correlation_search_space_resolution_;
  double correlation_search_space_smear_deviation_;
  double range_max_;
};

}  // namespace ndt_2d

#endif  // NDT_2D_KARTO__SCAN_MATCHER_KARTO_HPP_
