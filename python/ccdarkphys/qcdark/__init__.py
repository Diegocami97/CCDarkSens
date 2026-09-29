# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  __init__.py -- QCDark subpackage; it re-exports compute_dRdE and
#  repo_qcdark_data_dir.
# ============================================================================

from ccdarkphys.qcdark.entry import compute_dRdE, repo_qcdark_data_dir

__all__ = ["compute_dRdE", "repo_qcdark_data_dir"]
