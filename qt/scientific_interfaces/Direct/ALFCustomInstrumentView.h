// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllConfig.h"
#include "MantidQtWidgets/Common/FileFinderWidget.h"
#include "MantidQtWidgets/Common/ObserverPattern.h"
#include "MantidQtWidgets/InstrumentView/BaseCustomInstrumentView.h"
#include "MantidQtWidgets/InstrumentView/PlotFitAnalysisPaneView.h"

#include <QObject>
#include <QSplitter>
#include <QString>
#include <string>

namespace MantidQt {
namespace CustomInterfaces {

class MANTIDQT_DIRECT_DLL IALFCustomInstrumentView : public virtual MantidWidgets::IBaseCustomInstrumentView {
public:
  IALFCustomInstrumentView(){};
  virtual ~IALFCustomInstrumentView(){};
  virtual void addSpectrum(const std::string &wsName) = 0;
};

class MANTIDQT_DIRECT_DLL ALFCustomInstrumentView : public MantidWidgets::BaseCustomInstrumentView,
                                                    public IALFCustomInstrumentView {
  Q_OBJECT

public:
  explicit ALFCustomInstrumentView(const std::string &instrument, QWidget *parent = nullptr);
  void setUpInstrument(const std::string &fileName,
                       std::vector<std::function<bool(std::map<std::string, bool>)>> &binders) override final;

  void addSpectrum(const std::string &wsName) override final;

public slots:
  void selectWholeTube();
  void extractSingleTube();
  void averageTube();

private:
  QAction *m_extractAction;
  QAction *m_averageAction;
  MantidWidgets::IPlotFitAnalysisPaneView *m_analysisPane;
};
} // namespace CustomInterfaces
} // namespace MantidQt
