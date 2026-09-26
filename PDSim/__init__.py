__version__ = "2.15"
__git_revision__ = 'bdb3f251695a1cf43bac32cfc6e5df7ed5326650'
__git_branch__ = 'master'

from PDSim._abi_check import check_coolprop_abi as _check_coolprop_abi
_check_coolprop_abi()
del _check_coolprop_abi
