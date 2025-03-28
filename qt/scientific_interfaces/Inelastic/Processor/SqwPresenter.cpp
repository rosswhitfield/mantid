// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "SqwPresenter.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidQtWidgets/Common/UserInputValidator.h"
#include "MantidQtWidgets/Common/WorkspaceUtils.h"
#include "MantidQtWidgets/Plotting/AxisID.h"
#include "MantidQtWidgets/Spectroscopy/DataValidationHelper.h"
#include "MantidQtWidgets/Spectroscopy/InterfaceUtils.h"

#include "MantidGeometry/IComponent.h"
#include "MantidGeometry/Instrument.h"

#include <QFileInfo>

using namespace DataValidationHelper;
using namespace Mantid::API;
using MantidQt::API::BatchAlgorithmRunner;

using namespace MantidQt::MantidWidgets::WorkspaceUtils;
using namespace MantidQt::CustomInterfaces::InterfaceUtils;

namespace {
Mantid::Kernel::Logger g_log("S(Q,w)");

} // namespace

namespace MantidQt::CustomInterfaces {
//----------------------------------------------------------------------------------------------
/** Constructor
 */
SqwPresenter::SqwPresenter(QWidget *parent, std::unique_ptr<MantidQt::API::IAlgorithmRunner> algorithmRunner,
                           ISqwView *view, std::unique_ptr<ISqwModel> model)
    : DataProcessor(parent, std::move(algorithmRunner)), m_view(view), m_model(std::move(model)) {
  m_view->subscribePresenter(this);
  setRunWidgetPresenter(std::make_unique<RunPresenter>(this, m_view->getRunView()));
  setOutputPlotOptionsPresenter(m_view->getPlotOptions(), PlotWidget::SpectraSliceSurface);
}

/**
 * Handles the event of data being loaded. Validates the loaded data.
 *
 */
void SqwPresenter::handleDataReady(std::string const &dataName) {
  if (!m_view->validate()) {
    return;
  }
  m_model->setInputWorkspace(dataName);
  if (auto const eFixed = getEFixed(m_model->inputWorkspace())) {
    m_model->setEFixed(*eFixed);
  } else {
    m_view->showMessageBox("An 'Efixed' value could not be found in the provided workspace.");
    return;
  }

  plotRqwContour();
  m_view->setDefaultQAndEnergy();
}

void SqwPresenter::handleValidation(IUserInputValidator *validator) const {
  m_model->validate(validator, m_view->getQRangeFromPlot(), m_view->getERangeFromPlot());
}

/**
 * Handles plotting the S(Q, w) workspace when the algorithm chain is finished.
 *
 * @param error If the algorithm chain failed
 */
void SqwPresenter::runComplete(Mantid::API::IAlgorithm_sptr const algorithm, bool const error) {
  (void)algorithm;
  if (!error) {
    setOutputPlotOptionsWorkspaces({m_model->getOutputWorkspace()});
  }
  m_view->setEnableOutputOptions(!error);
}

/**
 * Plots the data as a contour plot
 *
 * Creates a colour 2D plot of the data
 */
void SqwPresenter::plotRqwContour() {
  try {
    auto const rqwWorkspace = m_model->getRqwWorkspace();
    if (rqwWorkspace)
      m_view->plotRqwContour(rqwWorkspace);
  } catch (std::exception const &ex) {
    g_log.warning(ex.what());
    m_view->showMessageBox("Invalid file. Please load a valid reduced workspace.");
  }
}

void SqwPresenter::setFileExtensionsByName(bool filter) {
  QStringList const noSuffixes{""};
  auto const tabName("Sqw");
  m_view->setFBSuffixes(filter ? getSampleFBSuffixes(tabName) : getExtensions(tabName));
  m_view->setWSSuffixes(filter ? getSampleWSSuffixes(tabName) : noSuffixes);
}

void SqwPresenter::setLoadHistory(bool doLoadHistory) { m_view->setLoadHistory(doLoadHistory); }

void SqwPresenter::handleRun() {
  clearOutputPlotOptionsWorkspaces();
  std::deque<API::IConfiguredAlgorithm_sptr> algoQueue = {};
  if (m_model->isRebinInEnergy())
    algoQueue.emplace_back(m_model->setupRebinAlgorithm());
  algoQueue.emplace_back(m_model->setupSofQWAlgorithm());
  algoQueue.emplace_back(m_model->setupAddSampleLogAlgorithm());
  m_view->setEnableOutputOptions(false);
  m_algorithmRunner->execute(std::move(algoQueue));
}

void SqwPresenter::handleSaveClicked() {
  if (checkADSForPlotSaveWorkspace(m_model->getOutputWorkspace(), false))
    m_algorithmRunner->execute(setupSaveAlgorithm(m_model->getOutputWorkspace()));
}

void SqwPresenter::handleQLowChanged(double const value) { m_model->setQMin(value); }

void SqwPresenter::handleQWidthChanged(double const value) { m_model->setQWidth(value); }

void SqwPresenter::handleQHighChanged(double const value) { m_model->setQMax(value); }

void SqwPresenter::handleELowChanged(double const value) { m_model->setEMin(value); }

void SqwPresenter::handleEWidthChanged(double const value) { m_model->setEWidth(value); }

void SqwPresenter::handleEHighChanged(double const value) { m_model->setEMax(value); }

void SqwPresenter::handleRebinEChanged(int const value) { m_model->setRebinInEnergy(value != 0); }
} // namespace MantidQt::CustomInterfaces
