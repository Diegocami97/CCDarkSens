"""
Load crystal |F|^2 tables produced by the standalone QCDark workflow (HDF5 layout
compatible with ``dark_matter_rates.form_factor`` in the reference QCDark code).

CCDarkSens does not ship these large files; set ``CCDARK_SENS_QCDARK_FORM_FACTOR``
or pass ``form_factor_h5`` to ``compute_dRdE``.
"""

from __future__ import annotations

import os

from ccdarkphys.common import constants as QEC


class CrystalFormFactor:
    """
    Crystal form factor on a (q, E) grid from an HDF5 file.

    Attributes mirror the reference QCDark reader: ``dq``, ``dE``, ``mCell``,
    ``ff`` (|F|^2), ``band_gap``.
    """

    def __init__(self, filename: str) -> None:
        try:
            import h5py
        except ImportError as exc:
            raise ImportError(
                "QCDark crystal tables require h5py (pip install h5py)."
            ) from exc

        path = os.path.abspath(os.path.expanduser(filename))
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        with h5py.File(path, "r") as data:
            self._source_path = path
            self.lattice = data["run_settings/a"][...].copy()
            self.atom = data["run_settings"].attrs["atom"]
            self.basis = data["run_settings"].attrs["basis"]
            self.ecp = data["run_settings"].attrs["ecp"]
            self.dft_rcut = float(data["run_settings/rcut"][...])
            self.dft_precision = float(data["run_settings/precision"][...])
            try:
                self.dark_Rvec = data["run_settings/Rvec"][...].copy()
            except KeyError:
                self.dark_Rvec = None
            self.num_con = data["run_settings"].attrs["numcon"]
            self.num_val = data["run_settings"].attrs["numval"]
            self.dft_xc = data["run_settings"].attrs["xc"]
            self.dft_density_fitting = data["run_settings"].attrs["df"]
            self.kpts = data["run_settings/kpts"][...].copy()
            self.VCell = data["results"].attrs["VCell"]
            self.mCell = float(data["results"].attrs["mCell"])
            bohr_inv_to_ev = QEC.alpha * QEC.me_eV
            self.dq = float(data["results"].attrs["dq"]) * bohr_inv_to_ev
            self.dE = float(data["results"].attrs["dE"])
            self.ff = data["results/f2"][...].copy()
            self.band_gap = float(data["results"].attrs["bandgap"])
            try:
                self.scissor_corrected = bool(data["results"].attrs["scissor"])
            except KeyError:
                self.scissor_corrected = False
        if self.ecp == "None":
            self.ecp = None
