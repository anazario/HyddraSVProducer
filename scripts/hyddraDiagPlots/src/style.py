"""
style.py — CMS label, canvas, and grid-line utilities for hyddraDiagPlots.
"""
import math
import ROOT


def draw_cms_label(right_label="Leptonic HYDDRA"):
    """Draw CMS / Simulation Preliminary label with right-side descriptor."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(0.057)
    latex.DrawLatex(0.16, 0.91, "CMS")
    latex.SetTextFont(52)
    latex.SetTextSize(0.044)
    latex.DrawLatex(0.25, 0.91, "Simulation Preliminary")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.90, 0.91, right_label)


def make_canvas(name, title, logy=True, width=800, height=600):
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.16)
    c.SetRightMargin(0.10)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.10)
    if logy:
        c.SetLogy(1)
    c.SetGridx(1)
    c.SetGridy(1)
    c.SetTicks(1, 1)
    return c


def draw_grid_lines(canvas, n_stages, min_y, max_y, logy=True):
    """Draw dotted vertical stage separators and horizontal y-grid lines."""
    canvas._vlines, canvas._hlines = [], []
    for i in range(1, n_stages):
        ln = ROOT.TLine(i + 0.5, min_y, i + 0.5, max_y)
        ln.SetLineStyle(3)
        ln.SetLineColor(ROOT.kGray + 2)
        ln.Draw()
        canvas._vlines.append(ln)
    if logy:
        lo = int(math.floor(math.log10(max(min_y, 1e-3))))
        hi = int(math.ceil(math.log10(max(max_y, 1))))
        h_vals = [10**exp for exp in range(lo, hi)]
    else:
        step = 10 ** math.floor(math.log10(max(max_y - min_y, 1e-9)) - 1)
        h_vals = [
            round(min_y + k * step, 10)
            for k in range(int((max_y - min_y) / step) + 2)
            if min_y < round(min_y + k * step, 10) < max_y
        ]
    for y in h_vals:
        ln = ROOT.TLine(0.5, y, n_stages + 0.5, y)
        ln.SetLineStyle(3)
        ln.SetLineColor(ROOT.kGray + 2)
        ln.Draw()
        canvas._hlines.append(ln)
