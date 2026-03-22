"""
style.py — CMS label, canvas, and grid-line utilities for hyddraDiagPlots.
"""
import math
import ROOT


def draw_cms_label(right_label="Leptonic HYDDRA", x_pos_cms=0.16, x_pos_label=0.9, y_pos=0.91):
    """Draw CMS / Simulation Preliminary label with right-side descriptor."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(0.057)
    latex.DrawLatex(x_pos_cms, 0.91, "CMS")
    latex.SetTextFont(52)
    latex.SetTextSize(0.044)
    latex.DrawLatex(x_pos_cms+0.09, y_pos, "Simulation Preliminary")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.DrawLatex(x_pos_label, y_pos, right_label)


def make_canvas(name, title, logy=True, width=800, height=600):
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.16)
    c.SetRightMargin(0.10)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.10)
    if logy:
        c.SetLogy(1)
    c.SetTicks(1, 1)
    return c


def _nice_ticks(lo, hi):
    """Return a list of nice round tick positions strictly inside [lo, hi]."""
    span = hi - lo
    if span <= 0:
        return []
    raw = 10 ** math.floor(math.log10(span) - 0.5)
    for f in [1, 2, 2.5, 5, 10]:
        step = raw * f
        if span / step <= 8:
            break
    first = math.ceil(lo / step + 1e-9) * step
    ticks, v = [], first
    while v < hi - 1e-9:
        ticks.append(round(v, 12))
        v += step
    return ticks


def draw_colz_grid(h2):
    """
    Draw dotted grid lines on top of a 2D COLZ histogram.
    Must be called after h2.Draw("COLZ") and canvas.Update().
    Returns the list of TLine objects (keep them alive on the canvas).
    """
    x_lo = h2.GetXaxis().GetXmin()
    x_hi = h2.GetXaxis().GetXmax()
    y_lo = h2.GetYaxis().GetXmin()
    y_hi = h2.GetYaxis().GetXmax()
    col  = ROOT.kGray + 1
    lines = []
    for x in _nice_ticks(x_lo, x_hi):
        ln = ROOT.TLine(x, y_lo, x, y_hi)
        ln.SetLineStyle(3); ln.SetLineColor(col); ln.Draw()
        lines.append(ln)
    for y in _nice_ticks(y_lo, y_hi):
        ln = ROOT.TLine(x_lo, y, x_hi, y)
        ln.SetLineStyle(3); ln.SetLineColor(col); ln.Draw()
        lines.append(ln)
    return lines


def draw_axis_grid(h_ax, logy=False):
    """
    Draw dotted grid lines for a plot with a continuous axis.
    Must be called AFTER all Draw("HIST SAME") calls so lines appear on top.
    Returns the list of TLine objects (keep them alive on the canvas).
    """
    x_lo = h_ax.GetXaxis().GetXmin()
    x_hi = h_ax.GetXaxis().GetXmax()
    y_lo = h_ax.GetMinimum()
    y_hi = h_ax.GetMaximum()
    col  = ROOT.kGray + 1
    lines = []

    # Horizontal (y) grid lines
    if logy:
        lo_e = int(math.floor(math.log10(max(y_lo, 1e-12))))
        hi_e = int(math.ceil(math.log10(max(y_hi, 1e-12))))
        y_ticks = [10**e for e in range(lo_e, hi_e) if y_lo < 10**e < y_hi]
    else:
        y_ticks = _nice_ticks(y_lo, y_hi)
    for y in y_ticks:
        ln = ROOT.TLine(x_lo, y, x_hi, y)
        ln.SetLineStyle(3); ln.SetLineColor(col); ln.Draw()
        lines.append(ln)

    # Vertical (x) grid lines
    for x in _nice_ticks(x_lo, x_hi):
        ln = ROOT.TLine(x, y_lo, x, y_hi)
        ln.SetLineStyle(3); ln.SetLineColor(col); ln.Draw()
        lines.append(ln)

    return lines


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
