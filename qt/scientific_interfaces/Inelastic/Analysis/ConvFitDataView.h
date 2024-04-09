// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "FitDataView.h"

#include "DllConfig.h"

#include <QTabWidget>

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

/**
Presenter for a table of convolution fitting data.
*/
class MANTIDQT_INELASTIC_DLL ConvFitDataView : public FitDataView {
  Q_OBJECT
public:
  ConvFitDataView(QWidget *parent = nullptr);
  void addTableEntry(size_t row, FitDataRow newRow) override;

protected:
  ConvFitDataView(const QStringList &headers, QWidget *parent = nullptr);

protected slots:
  void showAddWorkspaceDialog() override;
};

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt