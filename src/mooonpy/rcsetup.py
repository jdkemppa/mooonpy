# -*- coding: utf-8 -*-


# mooonpy/_config.py
from collections.abc import MutableMapping


# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.weight"] = "bold"
# plt.rcParams['axes.titleweight'] = "bold"
# plt.rcParams['figure.titleweight'] = "bold"
# plt.rcParams['axes.labelweight'] = "bold"
# plt.rcParams['figure.dpi'] = 163
# plt.rcParams['figure.figsize'] = (15,8.43)


# Matplotlib GitHub:
#    https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/rcsetup.py
#    https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/__init__.py
# Names:
#  rc = runtime configure
#  rcPhase
#  rcParams
#  rcDefaults
class RCParams(MutableMapping):
    def __init__(self, defaults):
        self._params = dict(defaults)

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, value):
        self._params[key] = value

    def __delitem__(self, key):
        del self._params[key]

    def __iter__(self):
        return iter(self._params)

    def __len__(self):
        return len(self._params)

    def update(self, new_params):
        self._params.update(new_params)

    def reset(self):
        self._params = dict(self._defaults)

    def __str__(self):
        return str(self._params)

# Default parameters
_defaults = {'molspace.read.dsect': ['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'],
             'molspace.read.steps': None,
             'molspace.write.data.astyle': 'full',
             'molspace.astyles': ['all'],
             
             'molspace.C.radii.ff.ReaxFF': 1.7,
             
             'thermospace.read': 'all',
             
             'xrdspace.read': 'all',
             
             'guis.size': 'large'
             }


rcParams = RCParams(_defaults)


if __name__ == "__main__": 
    print(rcParams.get('help'))
    print(rcParams.get('molspace.write.data.astyle'))