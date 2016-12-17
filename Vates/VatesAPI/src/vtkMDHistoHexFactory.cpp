#include "MantidAPI/IMDWorkspace.h"
#include "MantidKernel/CPUTimer.h"
#include "MantidDataObjects/MDHistoWorkspace.h"
#include "MantidVatesAPI/vtkMDHWSignalArray.h"
#include "MantidVatesAPI/Common.h"
#include "MantidVatesAPI/Normalization.h"
#include "MantidVatesAPI/ProgressAction.h"
#include "MantidVatesAPI/vtkNullStructuredGrid.h"
#include "MantidVatesAPI/vtkMDHistoHexFactory.h"
#include "MantidAPI/NullCoordTransform.h"
#include "MantidKernel/ReadLock.h"

#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkSMPTools.h"
#include "vtkStructuredGrid.h"
#include "vtkUnsignedCharArray.h"

#include <cmath>

using Mantid::API::IMDWorkspace;
using Mantid::API::IMDHistoWorkspace;
using Mantid::Kernel::CPUTimer;
using namespace Mantid::DataObjects;
using Mantid::Kernel::ReadLock;

namespace Mantid {
namespace VATES {

vtkMDHistoHexFactory::vtkMDHistoHexFactory(
    ThresholdRange_scptr thresholdRange,
    const VisualNormalization normalizationOption)
    : m_normalizationOption(normalizationOption),
      m_thresholdRange(thresholdRange) {}

/**
Assigment operator
@param other : vtkMDHistoHexFactory to assign to this instance from.
@return ref to assigned current instance.
*/
vtkMDHistoHexFactory &vtkMDHistoHexFactory::
operator=(const vtkMDHistoHexFactory &other) {
  if (this != &other) {
    this->m_normalizationOption = other.m_normalizationOption;
    this->m_thresholdRange = other.m_thresholdRange;
    this->m_workspace = other.m_workspace;
  }
  return *this;
}

/**
Copy Constructor
@param other : instance to copy from.
*/
vtkMDHistoHexFactory::vtkMDHistoHexFactory(const vtkMDHistoHexFactory &other) {
  this->m_normalizationOption = other.m_normalizationOption;
  this->m_thresholdRange = other.m_thresholdRange;
  this->m_workspace = other.m_workspace;
}

void vtkMDHistoHexFactory::initialize(Mantid::API::Workspace_sptr workspace) {
  m_workspace = doInitialize<MDHistoWorkspace, 3>(workspace);

  // Setup range values according to whatever strategy object has been injected.
  m_thresholdRange->setWorkspace(workspace);
  m_thresholdRange->calculate();
}

void vtkMDHistoHexFactory::validateWsNotNull() const {

  if (NULL == m_workspace.get()) {
    throw std::runtime_error("IMDWorkspace is null");
  }
}

void vtkMDHistoHexFactory::validate() const { validateWsNotNull(); }

namespace {
struct Worker {
  vtkMDHWSignalArray<double> *m_signal;
  vtkUnsignedCharArray *m_cga;
  Worker(vtkMDHWSignalArray<double> *signal, vtkUnsignedCharArray *cga)
      : m_signal(signal), m_cga(cga) {}
  void operator()(vtkIdType begin, vtkIdType end) {
    for (vtkIdType index = begin; index < end; ++index) {
      if (!std::isfinite(m_signal->GetValue(index))) {
        m_cga->SetValue(index, m_cga->GetValue(index) |
                                   vtkDataSetAttributes::HIDDENCELL);
      }
    }
  }
};

struct Worker2 {
  vtkPoints *m_pts;
  coord_t incrementX, incrementY, incrementZ;
  coord_t minX, minY, minZ;
  vtkIdType nPointsX, nPointsY;
  Worker2(Mantid::DataObjects::MDHistoWorkspace &ws, vtkPoints *pts)
      : m_pts(pts) {
    int nBinsX = static_cast<int>(ws.getXDimension()->getNBins());
    int nBinsY = static_cast<int>(ws.getYDimension()->getNBins());
    int nBinsZ = static_cast<int>(ws.getZDimension()->getNBins());

    minX = ws.getXDimension()->getMinimum();
    minY = ws.getYDimension()->getMinimum();
    minZ = ws.getZDimension()->getMinimum();
    coord_t maxX = ws.getXDimension()->getMaximum();
    coord_t maxY = ws.getYDimension()->getMaximum();
    coord_t maxZ = ws.getZDimension()->getMaximum();

    incrementX = (maxX - minX) / static_cast<coord_t>(nBinsX);
    incrementY = (maxY - minY) / static_cast<coord_t>(nBinsY);
    incrementZ = (maxZ - minZ) / static_cast<coord_t>(nBinsZ);

    nPointsX = nBinsX + 1;
    nPointsY = nBinsY + 1;
  }
  void operator()(vtkIdType begin, vtkIdType end) {
    float in[3];
    vtkIdType pos = begin * nPointsX * nPointsY;
    for (int z = static_cast<int>(begin); z < static_cast<int>(end); ++z) {
      in[2] = minZ + static_cast<coord_t>(z) * incrementZ;
      for (int y = 0; y < nPointsY; ++y) {
        in[1] = minY + static_cast<coord_t>(y) * incrementY;
        for (int x = 0; x < nPointsX; ++x) {
          in[0] = minX + static_cast<coord_t>(x) * incrementX;
          m_pts->SetPoint(pos, in);
          ++pos;
        }
      }
    }
  }
};
} // end anon namespace

/** Method for creating a 3D or 4D data set
 *
 * @param timestep :: index of the time step (4th dimension) in the workspace.
 *        Set to 0 for a 3D workspace.
 * @param progressUpdate: Progress updating. passes progress information up the
 *stack.
 * @return the vtkDataSet created
 */
vtkSmartPointer<vtkDataSet>
vtkMDHistoHexFactory::create3Dor4D(size_t timestep,
                                   ProgressAction &progressUpdate) const {
  // Acquire a scoped read-only lock to the workspace (prevent segfault from
  // algos modifying ws)
  ReadLock lock(*m_workspace);

  const size_t nDims = m_workspace->getNonIntegratedDimensions().size();

  std::vector<size_t> indexMultiplier(nDims, 0);

  // For quick indexing, accumulate these values
  // First multiplier
  indexMultiplier[0] = m_workspace->getDimension(0)->getNBins();
  for (size_t d = 1; d < nDims; d++) {
    indexMultiplier[d] =
        indexMultiplier[d - 1] * m_workspace->getDimension(d)->getNBins();
  }

  const int nBinsX = static_cast<int>(m_workspace->getXDimension()->getNBins());
  const int nBinsY = static_cast<int>(m_workspace->getYDimension()->getNBins());
  const int nBinsZ = static_cast<int>(m_workspace->getZDimension()->getNBins());

  const int imageSize = (nBinsX) * (nBinsY) * (nBinsZ);

  auto visualDataSet = vtkSmartPointer<vtkStructuredGrid>::New();
  visualDataSet->SetDimensions(nBinsX + 1, nBinsY + 1, nBinsZ + 1);

  // Array with true where the voxel should be shown
  // double progressFactor = 0.5 / static_cast<double>(imageSize);

  std::size_t offset = 0;
  if (nDims == 4) {
    offset = timestep * indexMultiplier[2];
  }

  vtkNew<vtkMDHWSignalArray<double>> signal;

  signal->SetName(vtkDataSetFactory::ScalarName.c_str());
  signal->InitializeArray(m_workspace.get(), m_normalizationOption, offset);
  visualDataSet->GetCellData()->SetScalars(signal.GetPointer());
  auto cga = visualDataSet->AllocateCellGhostArray();
  auto start = std::chrono::high_resolution_clock::now();

  Worker func(signal.GetPointer(), cga);
  vtkSMPTools::For(0, imageSize, func);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time1: " << elapsed_seconds.count() << "s\n";
  vtkNew<vtkPoints> points;

  const vtkIdType nPointsX = nBinsX + 1;
  const vtkIdType nPointsY = nBinsY + 1;
  const vtkIdType nPointsZ = nBinsZ + 1;

  points->SetNumberOfPoints(nPointsX * nPointsY * nPointsZ);

  start = std::chrono::high_resolution_clock::now();
  Worker2 func2(*m_workspace, points.GetPointer());
  vtkSMPTools::For(0, nPointsZ, func2);
  end = std::chrono::high_resolution_clock::now();
  elapsed_seconds = end - start;
  std::cout << "elapsed time2: " << elapsed_seconds.count() << "s\n";

  visualDataSet->SetPoints(points.GetPointer());
  visualDataSet->Register(NULL);
  visualDataSet->Squeeze();

  // Hedge against empty data sets
  if (visualDataSet->GetNumberOfPoints() <= 0) {
    vtkNullStructuredGrid nullGrid;
    visualDataSet = nullGrid.createNullData();
  }

  vtkSmartPointer<vtkDataSet> dataset = visualDataSet;
  return dataset;
}

/**
Create the vtkStructuredGrid from the provided workspace
@param progressUpdating: Reporting object to pass progress information up the
stack.
@return fully constructed vtkDataSet.
*/
vtkSmartPointer<vtkDataSet>
vtkMDHistoHexFactory::create(ProgressAction &progressUpdating) const {
  auto product =
      tryDelegatingCreation<MDHistoWorkspace, 3>(m_workspace, progressUpdating);
  if (product != NULL) {
    return product;
  } else {
    // Create in 3D mode
    return this->create3Dor4D(0, progressUpdating);
  }
}

vtkMDHistoHexFactory::~vtkMDHistoHexFactory() {}
}
}
