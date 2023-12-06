from typing import Sequence
import numpy as np

class Utilities:

    @staticmethod
    def convert_time_to_scan(time_value: float|Sequence[float], scan_rate: float) -> list[int]:

        '''
        Takes in single time or sequence of two times in minutes, returns the closest scan index or a sequence of two
        indices.
        '''

        if not type(time_value) in [np.ndarray, list, tuple]: time_value = [time_value]
        indices = []
        for value in time_value:
            value = int(np.round(value / scan_rate))
            indices.append(value)
        if len(indices) == 1: indices = indices[0]
        return indices

    @staticmethod
    def check_interval(value:float, midpoint:float, tolerance:float) -> bool:

        '''
        Takes in a value, midpoint, and tolerance and returns True if the value is within the interval defined by the
        midpoint and tolerance.
        '''

        if (midpoint-tolerance) < value < (midpoint+tolerance):
            return True
        else:
            return False

