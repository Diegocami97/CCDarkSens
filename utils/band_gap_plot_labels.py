# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: band_gap_plot_labels.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  band_gap_plot_labels.py -- Axis-label and colour helpers for band-gap
#  pheno plots
# ============================================================================
"""
Shared band-gap pheno plot labels (E_gap, epsilon_h).

Internal scenario ids (D-equal, B-thresh) stay in configs/paths; plot text uses
the same convention as utils/plot_band_gap_2d_heatmap.py.
"""

from __future__ import annotations

EH_B_THRESH_EV = 3.8


# ----------------------------------------------------------------------------
# eh_for_scenario
#   epsilon_h of a scenario: "D-equal" uses eps_h = E_gap, "B-thresh" the fixed B-thresh value; anything else raises ValueError.
# ----------------------------------------------------------------------------
def eh_for_scenario(gap_eV: float, scenario: str) -> float:
    if scenario == "D-equal":
        return gap_eV
    if scenario in ("B-thresh", "B_thresh"):
        return EH_B_THRESH_EV
    raise ValueError(f"unknown scenario: {scenario!r}")


def pheno_param_label_mpl(gap_eV: float, eh_eV: float) -> str:
    """Matplotlib / LaTeX: E_gap and epsilon_h."""
    return rf"$E_{{\mathrm{{gap}}}} = {gap_eV:g}$ eV, $\varepsilon_h = {eh_eV:g}$ eV"


def pheno_param_label_root(gap_eV: float, eh_eV: float) -> str:
    """ROOT TLatex for ccdarksens_plot_dmelectron_limit legends (plotter adds ' (q-map)')."""
    # Use E_{gap} not E_{#mathrm{gap}} — nested #mathrm inside _{} does not render in TLegend.
    return f"E_{{gap}} = {gap_eV:g} eV, #varepsilon_{{h}} = {eh_eV:g} eV"


# ----------------------------------------------------------------------------
# limit_curve_label
#   Legend label of a limit curve for the given gap and scenario.
# ----------------------------------------------------------------------------
def limit_curve_label(gap_eV: float, scenario: str) -> str:
    return pheno_param_label_root(gap_eV, eh_for_scenario(gap_eV, scenario))


# ----------------------------------------------------------------------------
# limit_sweep_title
#   Title of a Phase C limit-sweep plot for the heavy or light mediator, with the eps_h condition of the tier appended (ROOT or matplotlib syntax).
# ----------------------------------------------------------------------------
def limit_sweep_title(
    mediator: str = "heavy",
    *,
    tier: str | None = None,
    root: bool = False,
) -> str:
    med = "light" if mediator == "light" else "heavy"
    title = f"Band-gap pheno Phase C ({med} mediator)"
    if tier == "B-thresh":
        if root:
            title += f", #varepsilon_{{h}} = {EH_B_THRESH_EV:g} eV"
        else:
            title += rf", $\varepsilon_h = {EH_B_THRESH_EV:g}$ eV"
    elif tier == "D-equal":
        if root:
            title += ", #varepsilon_{h} = E_{gap}"
        else:
            title += r", $\varepsilon_h = E_{\mathrm{gap}}$"
    return title


# ----------------------------------------------------------------------------
# gap_title_root
#   ROOT-syntax title fragment "E_gap = <gap> eV".
# ----------------------------------------------------------------------------
def gap_title_root(gap_eV: float) -> str:
    return f"E_{{gap}} = {gap_eV:g} eV"


def si_reference_label_root() -> str:
    """Legend label for Si reference overlay (E_gap=1.2 eV, epsilon_h=3.8 eV)."""
    return pheno_param_label_root(1.2, 3.8) + " [Si ref]"
