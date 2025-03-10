// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "ALCInterface.h"

#include "ALCBaselineModellingView.h"
#include "ALCDataLoadingView.h"

#include "ALCPeakFittingView.h"

#include "ALCDataLoadingPresenter.h"
#include "ALCPeakFittingPresenter.h"

#include "ALCBaselineModellingModel.h"
#include "ALCDataLoadingModel.h"
#include "ALCPeakFittingModel.h"

#include "QInputDialog"
#include <QCloseEvent>

#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidKernel/Logger.h"
#include "MantidQtWidgets/Common/QtJobRunner.h"

#include <algorithm>

using namespace Mantid::API;

namespace {
Mantid::Kernel::Logger logger("ALC Interface");

auto &ads = AnalysisDataService::Instance();

MatrixWorkspace_sptr getWorkspace(const std::string &workspaceName) {
  if (ads.doesExist(workspaceName)) {
    return ads.retrieveWS<MatrixWorkspace>(workspaceName);
  } else {
    logger.warning("Workspace " + workspaceName + " was not found");
    return nullptr;
  }
}

QHash<QString, QVariant> createPointKwargs() {
  QHash<QString, QVariant> kwargs{{"marker", "."}, {"linestyle", "None"}};
  return kwargs;
}

QHash<QString, QVariant> createLineKwargs() {
  QHash<QString, QVariant> kwargs{{"marker", "None"}};
  return kwargs;
}

std::vector<std::optional<QHash<QString, QVariant>>> createPointAndLineKwargs() {
  std::vector<std::optional<QHash<QString, QVariant>>> kwargs;
  kwargs.emplace_back(createPointKwargs());
  kwargs.emplace_back(createLineKwargs());
  return kwargs;
}

} // namespace

namespace MantidQt::CustomInterfaces {
DECLARE_SUBWINDOW(ALCInterface)

const QStringList ALCInterface::STEP_NAMES = {"Data loading", "Baseline modelling", "Peak fitting"};

// %1 - current step no., %2 - total no. of steps, %3 - current step label
const QString ALCInterface::LABEL_FORMAT = "Step %1/%2 - %3";

// Custom close event - only allows the window to close if loading is not taking
// place.
void ALCInterface::closeEvent(QCloseEvent *event) {
  // check if currently loading data
  if (m_dataLoading) {
    if (m_dataLoading->isLoading()) {
      m_dataLoading->cancelLoading();
      event->ignore();
    } else {
      event->accept();
    }
  } else {
    event->accept();
  }
}
// namespace CustomInterfaces
ALCInterface::ALCInterface(QWidget *parent)
    : UserSubWindow(parent), m_ui(), m_peakFittingView(nullptr), m_dataLoading(nullptr), m_baselineModelling(nullptr),
      m_peakFitting(nullptr), m_externalPlotter(std::make_unique<Widgets::MplCpp::ExternalPlotter>()) {}

void ALCInterface::initLayout() {
  m_ui.setupUi(this);

  connect(m_ui.nextStep, &QPushButton::clicked, this, &ALCInterface::nextStep);
  connect(m_ui.previousStep, &QPushButton::clicked, this, &ALCInterface::previousStep);
  connect(m_ui.exportResults, &QPushButton::clicked, this, &ALCInterface::exportResults);
  connect(m_ui.importResults, &QPushButton::clicked, this, &ALCInterface::importResults);
  connect(m_ui.externalPlotButton, &QPushButton::clicked, this, &ALCInterface::externalPlotRequested);

  auto dataLoadingModel = std::make_unique<ALCDataLoadingModel>();
  auto dataLoadingView = new ALCDataLoadingView(m_ui.dataLoadingView);

  m_dataLoading = std::make_unique<ALCDataLoadingPresenter>(dataLoadingView, std::move(dataLoadingModel));
  m_dataLoading->initialize();

  auto baselineModellingModel = std::make_unique<ALCBaselineModellingModel>();
  auto *baselineModellingView = new ALCBaselineModellingView(m_ui.baselineModellingView);

  m_baselineModelling =
      std::make_unique<ALCBaselineModellingPresenter>(baselineModellingView, std::move(baselineModellingModel));
  m_baselineModelling->initialize();

  m_peakFittingView = new ALCPeakFittingView(m_ui.peakFittingView);
  auto jobRunner = std::make_unique<MantidQt::API::QtJobRunner>(true);
  auto algorithmRunner = std::make_unique<MantidQt::API::AlgorithmRunner>(std::move(jobRunner));
  m_peakFittingModel = std::make_shared<ALCPeakFittingModel>(std::move(algorithmRunner));
  m_peakFitting = std::make_unique<ALCPeakFittingPresenter>(m_peakFittingView, m_peakFittingModel.get());
  m_peakFitting->initialize();

  m_dataLoading->setSubscriber(this);
  m_baselineModelling->setSubscriber(this);

  assert(m_ui.stepView->count() == STEP_NAMES.count()); // Should have names for all steps

  switchStep(0); // We always start from the first step
}

void ALCInterface::loadedDataChanged() {

  // Make sure we do have some data
  if (m_dataLoading->loadedData()) {

    // Send the data to BaselineModelling
    m_baselineModelling->setData(m_dataLoading->loadedData());

    // If we have a fitting function and a fitting range
    // we can update the baseline model
    if ((!m_baselineModelling->function().empty()) && (m_baselineModelling->noOfSectionRows() > 0)) {

      // Fit the data
      m_baselineModelling->fit();
    }
  }
}

void ALCInterface::correctedDataChanged() {

  // Make sure we do have some data
  if (m_baselineModelling->correctedData()) {

    // Send the data to PeakFitting
    m_peakFittingModel->setData(m_baselineModelling->correctedData());

    // If we have a fitting function
    if (m_peakFittingView->function("")) {

      // Fit the data
      m_peakFittingView->fitRequested();
    }
  }
}

void ALCInterface::nextStep() {
  int next = m_ui.stepView->currentIndex() + 1;

  switchStep(next);
}

void ALCInterface::previousStep() {
  int previous = m_ui.stepView->currentIndex() - 1;

  switchStep(previous);
}

void ALCInterface::switchStep(int newStepIndex) {
  // Should be disallowed by disabling buttons
  assert(newStepIndex >= 0);
  assert(newStepIndex < m_ui.stepView->count());

  m_ui.label->setText(LABEL_FORMAT.arg(newStepIndex + 1).arg(STEP_NAMES.count()).arg(STEP_NAMES[newStepIndex]));

  int nextStepIndex = newStepIndex + 1;
  int prevStepIndex = newStepIndex - 1;

  bool hasNextStep = (nextStepIndex < m_ui.stepView->count());
  bool hasPrevStep = (prevStepIndex >= 0);

  m_ui.previousStep->setVisible(hasPrevStep);

  // On last step - hide next step button, but show "Export results..."
  m_ui.nextStep->setVisible(hasNextStep);

  if (hasPrevStep) {
    m_ui.previousStep->setText("< " + STEP_NAMES[prevStepIndex]);
  }

  if (hasNextStep) {
    m_ui.nextStep->setText(STEP_NAMES[nextStepIndex] + " >");
  }

  m_ui.stepView->setCurrentIndex(newStepIndex);
}

void ALCInterface::exportResults() {

  bool ok;
  QString label = QInputDialog::getText(this, "Results label", "Label to assign to the results: ", QLineEdit::Normal,
                                        "ALCResults", &ok);

  if (!ok) // Cancelled
  {
    return;
  }

  std::string groupName = label.toStdString();

  using namespace Mantid::API;

  std::map<std::string, Workspace_sptr> results;

  if (const auto loadedData = m_dataLoading->exportWorkspace())
    results["Loaded_Data"] = loadedData->clone();
  if (const auto baseline = m_baselineModelling->exportWorkspace())
    results["Baseline_Workspace"] = baseline->clone();
  if (const auto baselineSections = m_baselineModelling->exportSections())
    results["Baseline_Sections"] = baselineSections->clone();
  if (const auto baselineModel = m_baselineModelling->exportModel())
    results["Baseline_Model"] = baselineModel->clone();
  if (const auto peaksWorkspace = m_peakFittingModel->exportWorkspace())
    results["Peaks_Workspace"] = peaksWorkspace->clone();
  if (const auto peaksResults = m_peakFittingModel->exportFittedPeaks())
    results["Peaks_FitResults"] = peaksResults->clone();

  // Check if any of the above is not empty
  bool nothingToExport =
      std::none_of(results.cbegin(), results.cend(), [](const auto &result) { return result.second; });

  // There is something to export
  if (!nothingToExport) {

    // Add output group to the ADS
    AnalysisDataService::Instance().addOrReplace(groupName, std::make_shared<WorkspaceGroup>());

    for (auto &result : results) {
      if (result.second) {
        std::string wsName = groupName + "_" + result.first;
        AnalysisDataService::Instance().addOrReplace(wsName, result.second);
        AnalysisDataService::Instance().addToGroup(groupName, wsName);
      }
    }
  } else {
    // Nothing to export, show error message
    QMessageBox::critical(this, "Error", "Nothing to export");
  }
}

void ALCInterface::importResults() {
  bool okClicked;
  const auto workspaceName =
      QInputDialog::getText(this, "Results label", "Label to assign to the results: ", QLineEdit::Normal, "ALCResults",
                            &okClicked)
          .toStdString();

  if (!okClicked) {
    return;
  }
  if (!ads.doesExist(workspaceName)) {
    QMessageBox::critical(this, "Error", "Workspace " + QString::fromStdString(workspaceName) + " could not be found.");
    return;
  }

  if (ads.retrieveWS<WorkspaceGroup>(workspaceName)) {
    importLoadedData(workspaceName + "_Loaded_Data");
    importBaselineData(workspaceName + "_Baseline_Workspace");
    importPeakData(workspaceName + "_Peaks_Workspace");
  } else if (ads.retrieveWS<MatrixWorkspace>(workspaceName)) {
    importLoadedData(workspaceName);
    importBaselineData(workspaceName);
  }
}

void ALCInterface::importLoadedData(const std::string &workspaceName) {
  if (const auto dataWS = getWorkspace(workspaceName)) {
    m_dataLoading->setData(dataWS);
  }
}

void ALCInterface::importBaselineData(const std::string &workspaceName) {
  if (const auto baselineWS = getWorkspace(workspaceName)) {
    m_baselineModelling->setData(baselineWS);
    m_baselineModelling->setCorrectedData(baselineWS);
  }
}

void ALCInterface::importPeakData(const std::string &workspaceName) {
  if (const auto peaksWS = getWorkspace(workspaceName)) {
    m_peakFittingModel->setData(peaksWS);
  }
}

/**
 * Handles when External Plot is pressed on the ALC interface
 */
void ALCInterface::externalPlotRequested() {
  // Get current step to determine what data to externally plot
  switch (m_ui.stepView->currentIndex()) {
  case DataLoading:
    externalPlotDataLoading();
    break;
  case BaselineModel:
    externalPlotBaselineModel();
    break;
  case PeakFitting:
    externalPlotPeakFitting();
    break;
  }
}

/**
 * Plots in workbench the single workspace from the data given
 * @param data The workspace to add to the ADS before plotting
 * @param workspaceName The name of workspace to plot
 * @param workspaceIndices String list of indices to plot (e.g.
 * '0-2,5,7-10')
 * @param errorBars Boolean to add/remove error bars to plot
 * @param kwargs The kwargs used when plotting the workspace
 */
void ALCInterface::externallyPlotWorkspace(MatrixWorkspace_sptr &data, std::string const &workspaceName,
                                           std::string const &workspaceIndices, bool errorBars,
                                           std::optional<QHash<QString, QVariant>> const &kwargs) {
  AnalysisDataService::Instance().addOrReplace(workspaceName, data);
  m_externalPlotter->plotSpectra(workspaceName, workspaceIndices, errorBars, kwargs);
}

/**
 * Plots in workbench all the provided workspaces from the data given
 * @param data The workspace to add to the ADS before plotting
 * @param workspaceNames List of names of workspaces to plot
 * @param workspaceIndices List of indices to plot
 * @param errorBars List of booleans to add/remove error bars to each line individually
 * @param kwargs The kwargs used when plotting each of the workspaces
 */
void ALCInterface::externallyPlotWorkspaces(MatrixWorkspace_sptr &data, std::vector<std::string> const &workspaceNames,
                                            std::vector<int> const &workspaceIndices,
                                            std::vector<bool> const &errorBars,
                                            std::vector<std::optional<QHash<QString, QVariant>>> const &kwargs) {
  AnalysisDataService::Instance().addOrReplace(workspaceNames[0], data);
  m_externalPlotter->plotCorrespondingSpectra(workspaceNames, workspaceIndices, errorBars, kwargs);
}

/**
 * Handle Data Loading external plot requested. Will plot the loaded data if available
 */
void ALCInterface::externalPlotDataLoading() {
  if (auto data = m_dataLoading->exportWorkspace()) {
    externallyPlotWorkspace(data, "ALC_External_Plot_Loaded_Data", "0", true, createPointKwargs());
  } else
    logger.warning("Load some data before externally plotting");
}

/**
 * Handle Baseline Model external plot requested. Will plot the baseline model data if available otherwise will plot the
 * loaded data if available
 */
void ALCInterface::externalPlotBaselineModel() {
  if (auto data = m_baselineModelling->exportWorkspace()) {
    externallyPlotWorkspaces(data, std::vector<std::string>{2, "ALC_External_Plot_Baseline_Workspace"},
                             std::vector<int>{0, 1}, std::vector<bool>{true, false}, createPointAndLineKwargs());
  } else {
    // If we don't have a baseline model workspace, try to plot the raw data from the data loading tab
    externalPlotDataLoading();
  }
}

/**
 * Handle Baseline Model external plot requested. Will plot the peak fitting data if available otherwise will plot the
 * corrected data from the baseline model if available
 */
void ALCInterface::externalPlotPeakFitting() {
  if (auto data = m_peakFittingModel->exportWorkspace()) {
    externallyPlotWorkspaces(data, std::vector<std::string>{2, "ALC_External_Plot_Peaks_Workspace"},
                             std::vector<int>{0, 1}, std::vector<bool>{true, false}, createPointAndLineKwargs());
  } else {
    // If we don't have a peaks fit workspace, try to plot the raw peak data from the baseline model workspace (diff
    // spec (2))
    if (auto data = m_baselineModelling->exportWorkspace()) {
      // Plot the diff spec from the baseline model workspace
      externallyPlotWorkspace(data, "ALC_External_Plot_Baseline_Workspace", "2", true, createPointKwargs());
    } else
      logger.warning("Perform a baseline fit before externally plotting");
  }
}

} // namespace MantidQt::CustomInterfaces
