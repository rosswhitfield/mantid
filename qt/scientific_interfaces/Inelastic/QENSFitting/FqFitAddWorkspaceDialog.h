// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllConfig.h"
#include "MantidQtWidgets/Common/IAddWorkspaceDialog.h"
#include "ui_FqFitAddWorkspaceDialog.h"

#include <QDialog>

namespace MantidQt {
namespace CustomInterfaces {
namespace Inelastic {

class MANTIDQT_INELASTIC_DLL FqFitAddWorkspaceDialog : public QDialog, public MantidWidgets::IAddWorkspaceDialog {
  Q_OBJECT
public:
  explicit FqFitAddWorkspaceDialog(QWidget *parent);

  std::string workspaceName() const override;
  std::string parameterType() const;
  int parameterNameIndex() const;

  void setParameterTypes(const std::vector<std::string> &types);
  void setParameterNames(const std::vector<std::string> &names);
  void setWSSuffices(const QStringList &suffices) override;
  void setFBSuffices(const QStringList &suffices) override;

  void enableParameterSelection();
  void disableParameterSelection();

  void updateSelectedSpectra() override{};

signals:
  void addData(MantidWidgets::IAddWorkspaceDialog *dialog);
  void workspaceChanged(FqFitAddWorkspaceDialog *dialog, const std::string &workspace);
  void parameterTypeChanged(FqFitAddWorkspaceDialog *dialog, const std::string &type);

private slots:
  void emitWorkspaceChanged(const QString &name);
  void emitParameterTypeChanged(const QString &index);
  void emitAddData();
  void handleAutoLoaded();

private:
  Ui::FqFitAddWorkspaceDialog m_uiForm;
};

} // namespace Inelastic
} // namespace CustomInterfaces
} // namespace MantidQt
