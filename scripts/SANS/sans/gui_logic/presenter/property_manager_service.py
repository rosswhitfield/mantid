""" The property manager service.

The property manager service serializes a SANS state object into a PropertyManager object and places it on the
 PropertyManagerDataService. It is also used to retrieve the state from this service.
"""

from __future__ import (absolute_import, division, print_function)

from mantid.api import Algorithm
from mantid.kernel import (PropertyManagerProperty, PropertyManagerDataService)

from sans.state.state_base import create_deserialized_sans_state_from_property_manager


class PropertyManagerService(object):
    sans_property_manager_prefix = "SANS_PROPERTY_MANAGER_THIS_NEEDS_TO_BE_UNIQUE_"

    class ConverterAlgorithm(Algorithm):
        def PyInit(self):
            self.declareProperty(PropertyManagerProperty("Args"))

        def PyExec(self):
            pass

    def __init__(self):
        super(PropertyManagerService, self).__init__()
        self._convert_alg = PropertyManagerService.ConverterAlgorithm()
        self._convert_alg.initialize()

    def add_states_to_pmds(self, states):
        # 1. Remove all property managers which belong to the sans property manager type
        self.remove_sans_property_managers()

        # 2. Convert states to property managers
        property_managers = self._convert_states_to_property_managers(states)

        # 2. Add all property managers
        self._add_property_managers_to_pmds(property_managers)

    def get_single_state_from_pmds(self, index_to_retrieve):
        # 1. Find all sans state names
        sans_property_managers = {}
        for name in PropertyManagerDataService.getObjectNames():
            if name.startswith(self.sans_property_manager_prefix):
                property_manager = PropertyManagerDataService.retrieve(name)
                index = self._get_index_from_name(name)
                sans_property_managers.update({index: property_manager})

        # 2. Convert property managers to states
        if index_to_retrieve not in list(sans_property_managers.keys()):
            return []
        sans_property_manager = sans_property_managers[index_to_retrieve]
        states_map = self._convert_property_manager_to_state({index_to_retrieve: sans_property_manager})

        # 3. Create a sequence container
        return self._get_states_list(states_map)

    def get_states_from_pmds(self):
        # 1. Find all sans state names
        sans_property_managers = {}
        for name in PropertyManagerDataService.getObjectNames():
            if name.startswith(self.sans_property_manager_prefix):
                property_manager = PropertyManagerDataService.retrieve(name)
                index = self._get_index_from_name(name)
                sans_property_managers.update({index: property_manager})

        # 2. Convert property managers to states
        states_map = self._convert_property_manager_to_state(sans_property_managers)

        # 3. Create a sequence container
        return self._get_states_list(states_map)

    def remove_sans_property_managers(self):
        property_manager_names_to_delete = []
        for name in PropertyManagerDataService.getObjectNames():
            if name.startswith(self.sans_property_manager_prefix):
                property_manager_names_to_delete.append(name)

        for element in property_manager_names_to_delete:
            PropertyManagerDataService.remove(element)

    def _convert_states_to_property_managers(self, states):
        property_managers = {}
        for index, state in states.items():
            property_manager = self._convert_state_to_property_manager(state)
            property_managers.update({index: property_manager})
        return property_managers

    def _convert_state_to_property_manager(self, state):
        self._convert_alg.setProperty("Args", state.property_manager)
        return self._convert_alg.getProperty("Args").value

    def _add_property_managers_to_pmds(self, property_managers):
        for index, property_manager in property_managers.items():
            name = self.sans_property_manager_prefix + str(index)
            PropertyManagerDataService.addOrReplace(name, property_manager)

    def _get_index_from_name(self, name):
        return int(name.replace(self.sans_property_manager_prefix, ""))

    @staticmethod
    def _convert_property_manager_to_state(property_managers):
        states = {}
        for key, property_manager in property_managers.items():
            state = create_deserialized_sans_state_from_property_manager(property_manager)
            states.update({key: state})
        return states

    @staticmethod
    def _get_states_list(states_map):
        states = []
        indices = list(states_map.keys())
        indices.sort()
        for index in indices:
            states.append(states_map[index])
        return states
