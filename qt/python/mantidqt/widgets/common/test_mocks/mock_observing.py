# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench.
from __future__ import (absolute_import, division, print_function)

from mock import Mock

from mantidqt.widgets.common.observing_presenter import ObservingPresenter
from mantidqt.widgets.common.observing_view import ObservingView
from mantidqt.widgets.common.test_mocks.mock_qt import MockQtSignal


class MockObservingView(ObservingView):
    def __init__(self, _):
        ObservingView.__init__(self, _)
        self.close_signal = MockQtSignal()
        self.rename_signal = MockQtSignal()
        self.presenter = Mock(spec={'clear_observer'})
        self.setWindowTitle = Mock()


class MockObservingPresenter(ObservingPresenter):
    def __init__(self, workspaces_are_equal=True):
        self.ads_observer = Mock()
        self.view = MockObservingView(None)
        self.model = Mock()
        self.model.workspace_equals = Mock(return_value=workspaces_are_equal)
