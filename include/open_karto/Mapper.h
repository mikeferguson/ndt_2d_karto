/*
 * Copyright 2010 SRI International
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPEN_KARTO_MAPPER_H
#define OPEN_KARTO_MAPPER_H

#include <map>
#include <vector>
#include <queue>

#include <open_karto/Karto.h>

namespace karto
{
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  class ScanMatcher;
  class ScanMatcherParams;

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Implementation of a correlation grid used for scan matching
   */
  class CorrelationGrid : public Grid<kt_int8u>
  {
  public:
    /**
     * Destructor
     */
    virtual ~CorrelationGrid()
    {
      delete [] m_pKernel;
    }

  public:
    /**
     * Create a correlation grid of given size and parameters
     * @param width
     * @param height
     * @param resolution
     * @param smearDeviation
     * @return correlation grid
     */
    static CorrelationGrid* CreateGrid(kt_int32s width,
                                       kt_int32s height,
                                       kt_double resolution,
                                       kt_double smearDeviation)
    {
      assert(resolution != 0.0);

      // +1 in case of roundoff
      kt_int32u borderSize = GetHalfKernelSize(smearDeviation, resolution) + 1;

      CorrelationGrid* pGrid = new CorrelationGrid(width, height, borderSize, resolution, smearDeviation);

      return pGrid;
    }

    /**
     * Gets the index into the data pointer of the given grid coordinate
     * @param rGrid
     * @param boundaryCheck
     * @return grid index
     */
    virtual kt_int32s GridIndex(const Vector2<kt_int32s>& rGrid, kt_bool boundaryCheck = true) const
    {
      kt_int32s x = rGrid.GetX() + m_Roi.GetX();
      kt_int32s y = rGrid.GetY() + m_Roi.GetY();

      return Grid<kt_int8u>::GridIndex(Vector2<kt_int32s>(x, y), boundaryCheck);
    }

    /**
     * Get the Region Of Interest (ROI)
     * @return region of interest
     */
    inline const Rectangle2<kt_int32s>& GetROI() const
    {
      return m_Roi;
    }

    /**
     * Sets the Region Of Interest (ROI)
     * @param roi
     */
    inline void SetROI(const Rectangle2<kt_int32s>& roi)
    {
      m_Roi = roi;
    }

    /**
     * Smear cell if the cell at the given point is marked as "occupied"
     * @param rGridPoint
     */
    inline void SmearPoint(const Vector2<kt_int32s>& rGridPoint)
    {
      assert(m_pKernel != NULL);

      int gridIndex = GridIndex(rGridPoint);
      if (GetDataPointer()[gridIndex] != GridStates_Occupied)
      {
        return;
      }

      kt_int32s halfKernel = m_KernelSize / 2;

      // apply kernel
      for (kt_int32s j = -halfKernel; j <= halfKernel; j++)
      {
        kt_int8u* pGridAdr = GetDataPointer(Vector2<kt_int32s>(rGridPoint.GetX(), rGridPoint.GetY() + j));

        kt_int32s kernelConstant = (halfKernel) + m_KernelSize * (j + halfKernel);

        // if a point is on the edge of the grid, there is no problem
        // with running over the edge of allowable memory, because
        // the grid has margins to compensate for the kernel size
        for (kt_int32s i = -halfKernel; i <= halfKernel; i++)
        {
          kt_int32s kernelArrayIndex = i + kernelConstant;

          kt_int8u kernelValue = m_pKernel[kernelArrayIndex];
          if (kernelValue > pGridAdr[i])
          {
            // kernel value is greater, so set it to kernel value
            pGridAdr[i] = kernelValue;
          }
        }
      }
    }

  protected:
    /**
     * Constructs a correlation grid of given size and parameters
     * @param width
     * @param height
     * @param borderSize
     * @param resolution
     * @param smearDeviation
     */
    CorrelationGrid(kt_int32u width, kt_int32u height, kt_int32u borderSize,
                    kt_double resolution, kt_double smearDeviation)
      : Grid<kt_int8u>(width + borderSize * 2, height + borderSize * 2)
      , m_SmearDeviation(smearDeviation)
      , m_pKernel(NULL)
    {
      GetCoordinateConverter()->SetScale(1.0 / resolution);

      // setup region of interest
      m_Roi = Rectangle2<kt_int32s>(borderSize, borderSize, width, height);

      // calculate kernel
      CalculateKernel();
    }

    /**
     * Sets up the kernel for grid smearing.
     */
    virtual void CalculateKernel()
    {
      kt_double resolution = GetResolution();

      assert(resolution != 0.0);
      assert(m_SmearDeviation != 0.0);

      // min and max distance deviation for smearing;
      // will smear for two standard deviations, so deviation must be at least 1/2 of the resolution
      const kt_double MIN_SMEAR_DISTANCE_DEVIATION = 0.5 * resolution;
      const kt_double MAX_SMEAR_DISTANCE_DEVIATION = 10 * resolution;

      // check if given value too small or too big
      if (!math::InRange(m_SmearDeviation, MIN_SMEAR_DISTANCE_DEVIATION, MAX_SMEAR_DISTANCE_DEVIATION))
      {
        std::stringstream error;
        error << "Mapper Error:  Smear deviation too small:  Must be between "
              << MIN_SMEAR_DISTANCE_DEVIATION
              << " and "
              << MAX_SMEAR_DISTANCE_DEVIATION;
        throw std::runtime_error(error.str());
      }

      // NOTE:  Currently assumes a two-dimensional kernel

      // +1 for center
      m_KernelSize = 2 * GetHalfKernelSize(m_SmearDeviation, resolution) + 1;

      // allocate kernel
      m_pKernel = new kt_int8u[m_KernelSize * m_KernelSize];
      if (m_pKernel == NULL)
      {
        throw std::runtime_error("Unable to allocate memory for kernel!");
      }

      // calculate kernel
      kt_int32s halfKernel = m_KernelSize / 2;
      for (kt_int32s i = -halfKernel; i <= halfKernel; i++)
      {
        for (kt_int32s j = -halfKernel; j <= halfKernel; j++)
        {
#ifdef WIN32
          kt_double distanceFromMean = _hypot(i * resolution, j * resolution);
#else
          kt_double distanceFromMean = hypot(i * resolution, j * resolution);
#endif
          kt_double z = exp(-0.5 * pow(distanceFromMean / m_SmearDeviation, 2));

          kt_int32u kernelValue = static_cast<kt_int32u>(math::Round(z * GridStates_Occupied));
          assert(math::IsUpTo(kernelValue, static_cast<kt_int32u>(255)));

          int kernelArrayIndex = (i + halfKernel) + m_KernelSize * (j + halfKernel);
          m_pKernel[kernelArrayIndex] = static_cast<kt_int8u>(kernelValue);
        }
      }
    }

    /**
     * Computes the kernel half-size based on the smear distance and the grid resolution.
     * Computes to two standard deviations to get 95% region and to reduce aliasing.
     * @param smearDeviation
     * @param resolution
     * @return kernel half-size based on the parameters
     */
    static kt_int32s GetHalfKernelSize(kt_double smearDeviation, kt_double resolution)
    {
      assert(resolution != 0.0);

      return static_cast<kt_int32s>(math::Round(2.0 * smearDeviation / resolution));
    }

  private:
    /**
     * The point readings are smeared by this value in X and Y to create a smoother response.
     * Default value is 0.03 meters.
     */
    kt_double m_SmearDeviation;

    // Size of one side of the kernel
    kt_int32s m_KernelSize;

    // Cached kernel for smearing
    kt_int8u* m_pKernel;

    // region of interest
    Rectangle2<kt_int32s> m_Roi;
  };  // CorrelationGrid

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Scan matcher
   */
  class KARTO_EXPORT ScanMatcher
  {
  public:
    /**
     * Destructor
     */
    virtual ~ScanMatcher();

  public:
    /**
     * Create a scan matcher with the given parameters
     */
    static ScanMatcher* Create(ScanMatcherParams* pParams,
                               kt_double searchSize,
                               kt_double resolution,
                               kt_double smearDeviation,
                               kt_double rangeThreshold);

    /**
     * Match given scan against set of scans
     * @param pScan scan being scan-matched
     * @param rBaseScans set of scans whose points will mark cells in grid as being occupied
     * @param rMean output parameter of mean (best pose) of match
     * @param rCovariance output parameter of covariance of match
     * @param doPenalize whether to penalize matches further from the search center
     * @param doRefineMatch whether to do finer-grained matching if coarse match is good (default is true)
     * @return strength of response
     */
    kt_double MatchScan(LocalizedRangeScan* pScan,
                        const LocalizedRangeScanVector& rBaseScans,
                        Pose2& rMean, Matrix3& rCovariance,
                        kt_bool doPenalize = true,
                        kt_bool doRefineMatch = true);

    /**
     * Finds the best pose for the scan centering the search in the correlation grid
     * at the given pose and search in the space by the vector and angular offsets
     * in increments of the given resolutions
     * @param pScan scan to match against correlation grid
     * @param rSearchCenter the center of the search space
     * @param rSearchSpaceOffset searches poses in the area offset by this vector around search center
     * @param rSearchSpaceResolution how fine a granularity to search in the search space
     * @param searchAngleOffset searches poses in the angles offset by this angle around search center
     * @param searchAngleResolution how fine a granularity to search in the angular search space
     * @param doPenalize whether to penalize matches further from the search center
     * @param rMean output parameter of mean (best pose) of match
     * @param rCovariance output parameter of covariance of match
     * @param doingFineMatch whether to do a finer search after coarse search
     * @return strength of response
     */
    kt_double CorrelateScan(LocalizedRangeScan* pScan,
                            const Pose2& rSearchCenter,
                            const Vector2<kt_double>& rSearchSpaceOffset,
                            const Vector2<kt_double>& rSearchSpaceResolution,
                            kt_double searchAngleOffset,
                            kt_double searchAngleResolution,
                            kt_bool doPenalize,
                            Pose2& rMean,
                            Matrix3& rCovariance,
                            kt_bool doingFineMatch);

    /**
     * Compute the score
     */
    double ComputeScore(LocalizedRangeScan * scan,
                        const LocalizedRangeScanVector& rBaseScans);

    /**
     * Computes the positional covariance of the best pose
     * @param rBestPose
     * @param bestResponse
     * @param rSearchCenter
     * @param rSearchSpaceOffset
     * @param rSearchSpaceResolution
     * @param searchAngleResolution
     * @param rCovariance
     */
    void ComputePositionalCovariance(const Pose2& rBestPose,
                                     kt_double bestResponse,
                                     const Pose2& rSearchCenter,
                                     const Vector2<kt_double>& rSearchSpaceOffset,
                                     const Vector2<kt_double>& rSearchSpaceResolution,
                                     kt_double searchAngleResolution,
                                     Matrix3& rCovariance);

    /**
     * Computes the angular covariance of the best pose
     * @param rBestPose
     * @param bestResponse
     * @param rSearchCenter
     * @param searchAngleOffset
     * @param searchAngleResolution
     * @param rCovariance
     */
    void ComputeAngularCovariance(const Pose2& rBestPose,
                                  kt_double bestResponse,
                                  const Pose2& rSearchCenter,
                                  kt_double searchAngleOffset,
                                  kt_double searchAngleResolution,
                                  Matrix3& rCovariance);

    /**
     * Gets the correlation grid data (for debugging)
     * @return correlation grid
     */
    inline CorrelationGrid* GetCorrelationGrid() const
    {
      return m_pCorrelationGrid;
    }

  private:
    /**
     * Marks cells where scans' points hit as being occupied
     * @param rScans scans whose points will mark cells in grid as being occupied
     * @param viewPoint do not add points that belong to scans "opposite" the view point
     */
    void AddScans(const LocalizedRangeScanVector& rScans, Vector2<kt_double> viewPoint);

    /**
     * Marks cells where scans' points hit as being occupied.  Can smear points as they are added.
     * @param pScan scan whose points will mark cells in grid as being occupied
     * @param viewPoint do not add points that belong to scans "opposite" the view point
     * @param doSmear whether the points will be smeared
     */
    void AddScan(LocalizedRangeScan* pScan, const Vector2<kt_double>& rViewPoint, kt_bool doSmear = true);

    /**
     * Compute which points in a scan are on the same side as the given viewpoint
     * @param pScan
     * @param rViewPoint
     * @return points on the same side
     */
    PointVectorDouble FindValidPoints(LocalizedRangeScan* pScan, const Vector2<kt_double>& rViewPoint) const;

    /**
     * Get response at given position for given rotation (only look up valid points)
     * @param angleIndex
     * @param gridPositionIndex
     * @return response
     */
    kt_double GetResponse(kt_int32u angleIndex, kt_int32s gridPositionIndex) const;

  protected:
    /**
     * Default constructor
     */
    ScanMatcher(ScanMatcherParams* pParams)
      : m_pParams(pParams)
      , m_pCorrelationGrid(NULL)
      , m_pSearchSpaceProbs(NULL)
      , m_pGridLookup(NULL)
    {
    }

  private:
    ScanMatcherParams* m_pParams;

    CorrelationGrid* m_pCorrelationGrid;
    Grid<kt_double>* m_pSearchSpaceProbs;

    GridIndexLookup<kt_int8u>* m_pGridLookup;
  };  // ScanMatcher

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  class KARTO_EXPORT ScanMatcherParams : public Module
  {
    friend class ScanMatcher;

  public:
    ScanMatcherParams();
    ~ScanMatcherParams();
    void Reset();

    //////////////////////////////////////////////////////////////////////////////
    //    CorrelationParameters correlationParameters;

    /**
     * The size of the search grid used by the matcher.
     * Default value is 0.3 meters which tells the matcher to use a 30cm x 30cm grid.
     */
    Parameter<kt_double>* m_pCorrelationSearchSpaceDimension;

    /**
     * The resolution (size of a grid cell) of the correlation grid.
     * Default value is 0.01 meters.
     */
    Parameter<kt_double>* m_pCorrelationSearchSpaceResolution;

    /**
     * The point readings are smeared by this value in X and Y to create a smoother response.
     * Default value is 0.03 meters.
     */
    Parameter<kt_double>* m_pCorrelationSearchSpaceSmearDeviation;

    //////////////////////////////////////////////////////////////////////////////
    // ScanMatcherParameters;

    // Variance of penalty for deviating from odometry when scan-matching.
    // The penalty is a multiplier (less than 1.0) is a function of the
    // delta of the scan position being tested and the odometric pose
    Parameter<kt_double>* m_pDistanceVariancePenalty;
    Parameter<kt_double>* m_pAngleVariancePenalty;

    // The range of angles to search during a coarse search and a finer search
    Parameter<kt_double>* m_pFineSearchAngleOffset;
    Parameter<kt_double>* m_pCoarseSearchAngleOffset;

    // Resolution of angles to search during a coarse search
    Parameter<kt_double>* m_pCoarseAngleResolution;

    // Minimum value of the penalty multiplier so scores do not
    // become too small
    Parameter<kt_double>* m_pMinimumAnglePenalty;
    Parameter<kt_double>* m_pMinimumDistancePenalty;

    // whether to increase the search space if no good matches are initially found
    Parameter<kt_bool>* m_pUseResponseExpansion;
  };
}  // namespace karto

#endif  // OPEN_KARTO_MAPPER_H
