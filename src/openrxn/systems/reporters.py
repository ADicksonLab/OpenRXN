"""Reporters are attached to system objects and control
what information is retained and attached to system.results
during system dynamics."""

import numpy as np

class Reporter(object):
    """Base class for Reporters"""

    def __init__(self, freq=None):
        """The following attributes are necessary for Reporter
        objects.

        freq:  int
        the frequency with which results are reported
        as measured in number of timesteps

        """
        self.freq = freq
        self._reports = []

    def report(self,current_time,current_state_vec):
        """Function for returning results that are attached
        to system.results.  current_state_vec is the state vector 
        of the system."""

        return NotImplementedError

    def reports(self):
        """Return the list of reports that have been collected."""

        return self._reports

class AllReporter(Reporter):
    """Reports the entire state_vec at some specified frequency."""

    def __init__(self, freq=1):
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        self._reports.append({'t': current_time,
                              'report': current_state_vec})

class SelectionReporter(Reporter):
    """Reports a subsection of the full state_vec at some specified
    frequency.
    
    selection_idxs : list of int
    A list of indexes to report.  This needs to be compatible 
    with the elements of the system.state.q_vec vector.
    """

    def __init__(self, selection_idxs, freq=1):
        self.selection_idxs = selection_idxs
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        self._reports.append({'t': current_time,
                              'report': current_state_vec[self.selection_idxs]})

class SumReporter(Reporter):
    """Reports the sum of values for a subselection of the 
    full state_vec at some specified frequency.
    
    selection_idxs : list of int
    A list of indexes to report.  This needs to be compatible 
    with the elements of the system.state.q_vec vector.
    """

    def __init__(self, selection_idxs, freq=1):
        self.selection_idxs = selection_idxs
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        self._reports.append({'t': current_time,
                              'report': np.sum(current_state_vec[self.selection_idxs])})

class AvgReporter(Reporter):
    """Reports the average of values for a subselection of the 
    full state_vec at some specified frequency.
    
    selection_idxs : list of int
    A list of indexes to report.  This needs to be compatible 
    with the elements of the system.state.q_vec vector.
    """

    def __init__(self, selection_idxs, freq=1):
        self.selection_idxs = selection_idxs
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        self._reports.append({'t': current_time,
                              'report': np.sum(current_state_vec[self.selection_idxs])/len(self.selection_idxs)})

class MaxReporter(Reporter):
    """Reports the maxmimum value over a subselection of the 
    full state_vec variables at some specified frequency.
    
    selection_idxs : list of int
    A list of indexes to report.  This needs to be compatible 
    with the elements of the system.state.q_vec vector.

    returns:

    (max, argmax)
    """

    def __init__(self, selection_idxs, freq=1):
        self.selection_idxs = selection_idxs
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        tmp = current_state_vec[self.selection_idxs]
        self._reports.append({'t': current_time,
                              'report': (np.max(tmp),np.argmax(tmp))})

class MinReporter(Reporter):
    """Reports the minimum value over a subselection of the 
    full state_vec at some specified frequency.
    
    selection_idxs : list of int
    A list of indexes to report.  This needs to be compatible 
    with the elements of the system.state.q_vec vector.

    returns:

    (min, argmin)
    """

    def __init__(self, selection_idxs, freq=1):
        self.selection_idxs = selection_idxs
        super().__init__(freq=freq)

    def report(self,current_time,current_state_vec):
        tmp = current_state_vec[self.selection_idxs]
        self._reports.append({'t': current_time,
                              'report': (np.max(tmp),np.argmax(tmp))})

