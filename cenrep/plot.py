"""S8: centromere overview: one panel per genome (A, B, ...; two columns), one track per chromosome; a borrowed
(contig) genome, or one with more than MAX_ALL sequences, shows the MAX_TRACKS called sequences with the most
centromeric satellite. Every track has the same height (capped so the figure stays <= ~13 in tall).
Track = satellite bp per 100-kb bin (primary families coloured, cen-associated grey), scaled gene density line,
and a box for the centromere call (dashed for fallback calls). No titles; panel letters upper-left.
Palette: Okabe-Ito without yellow (fails the lightness band on white; dataviz validator), fixed order, never
cycled: a 7th+ centromeric family folds into one labelled 'other centromeric' grey."""
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from .util import union_bp  # noqa: E402

OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]
OTHER_CEN = "#555555"
ASSOC = "#AAAAAA"
GENE = "#999999"        # context: grey, drawn under the satellite bars
BIN = 100_000
H = 0.8
MAX_TRACKS = 20
MAX_ALL = 50            # chromosome-level genomes with up to this many sequences show them all
TRACK_IN, TRACKS_MAX_IN = 0.32, 12.0


def bin_bp(arrays, L: int, bin_size: int = BIN) -> np.ndarray:
    out = np.zeros(max(1, int(np.ceil(L / bin_size))))
    for s, e in zip(arrays.start, arrays.end):
        for b in range(int(s) // bin_size, (int(e) - 1) // bin_size + 1):
            out[b] += min(e, (b + 1) * bin_size) - max(s, b * bin_size)
    return out


def stack_heights(layers: list, bin_size: int = BIN) -> list:
    """Scale stacked per-bin bp so the total never exceeds one bin (overlapping arrays of several families would
    otherwise spill into the next chromosome track)."""
    if not layers:
        return layers
    tot = np.sum(layers, axis=0)
    f = np.where(tot > bin_size, bin_size / np.maximum(tot, 1), 1.0)
    return [x * f for x in layers]


def family_colors(fams: list) -> dict:
    return {f: OKABE_ITO[i] if i < len(OKABE_ITO) else OTHER_CEN for i, f in enumerate(fams)}


def grid_shape(n: int) -> tuple:
    return (1, 1) if n == 1 else (math.ceil(n / 2), 2)


def track_height(heights: list) -> float:
    """Inches per track: TRACK_IN, shrunk so all tracks together take at most TRACKS_MAX_IN."""
    return min(TRACK_IN, TRACKS_MAX_IN / max(1, sum(heights)))


def panel_tracks(pn: dict, max_tracks: int = MAX_TRACKS) -> tuple:
    """(tracks [(chrom, length)], y label). Chromosome-level genomes show every sequence; a borrowed (contig)
    genome, or one with more than MAX_ALL sequences, shows the max_tracks called ones with the most
    centromeric-satellite bp, longest first."""
    chroms = list(pn["chroms"].items())
    if len(chroms) <= max_tracks or (not pn.get("borrowed") and len(chroms) <= MAX_ALL):
        return chroms, pn["genome"]
    called = [c for c, _ in chroms if pn["regions"].get(c, {}).get("evidence", "none") != "none"]
    a = pn["arrays"][pn["arrays"].cls == "centromeric"]
    bp = {c: union_bp(zip(d.start, d.end)) for c, d in a.groupby("chrom")}
    top = sorted(called, key=lambda c: (-bp.get(c, 0), c))[:max_tracks]
    tracks = sorted(((c, pn["chroms"][c]) for c in top), key=lambda x: (-x[1], x[0]))
    return tracks, f"{pn['genome']}: {len(tracks)} of {len(called)}\nsequences with a call\n(most centromeric bp)"


def plot_overview(path, panels: list) -> int:
    """Writes the figure; returns the number of panels drawn."""
    plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "ps.fonttype": 42})
    fams = sorted({f for pn in panels for f in pn["arrays"].loc[pn["arrays"].cls == "centromeric", "family"]})
    color = family_colors(fams)
    shown = [panel_tracks(pn) for pn in panels]
    nrows, ncols = grid_shape(len(panels))
    heights = [max([1] + [len(t) for t, _ in shown[r * ncols:(r + 1) * ncols]]) for r in range(nrows)]
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 1.2 + track_height(heights) * sum(heights)),
                             squeeze=False, gridspec_kw={"height_ratios": heights}, constrained_layout=True)
    row_n = [heights[i // ncols] for i in range(len(panels))]
    letters = [chr(65 + i) if i < 26 else chr(64 + i // 26) + chr(65 + i % 26) for i in range(len(panels))]
    for ax, pn, (chroms, label), letter, rn in zip(axes.flat, panels, shown, letters, row_n):
        n = len(chroms)
        for row, (c, L) in enumerate(chroms):
            y = n - 1 - row
            ax.add_patch(Rectangle((0, y), L / 1e6, H, fill=False, ec="#BBBBBB", lw=0.4, zorder=1))
            a = pn["arrays"][pn["arrays"].chrom == c]
            x = (np.arange(int(np.ceil(L / BIN))) + 0.5) * BIN / 1e6
            bottom = np.zeros(len(x))
            layers = [(a[(a.family == f) & (a.cls == "centromeric")], color[f]) for f in fams]
            layers.append((a[a.cls == "cen_associated"], ASSOC))
            layers = [(bin_bp(sel, L)[:len(x)], col) for sel, col in layers if not sel.empty]
            for h, (_, col) in zip(stack_heights([bp for bp, _ in layers]), layers):
                h = h / BIN * H
                ax.bar(x, h, width=BIN / 1e6, bottom=y + bottom, color=col, linewidth=0, zorder=2)
                bottom += h
            genes = pn["genes"]
            if genes is not None and len(genes.get(c, [])):
                cnt = np.bincount(genes[c] // BIN, minlength=len(x))[:len(x)]
                ax.plot(x, y + cnt / max(cnt.max(), 1) * H, color=GENE, lw=0.3, zorder=1.5)
            r = pn["regions"].get(c)
            if r is not None and r["evidence"] != "none":
                ax.add_patch(Rectangle((r["cen_start"] / 1e6, y - 0.06), (r["cen_end"] - r["cen_start"]) / 1e6,
                                       H + 0.12, fill=False, ec="black", lw=0.7,
                                       ls="-" if r["evidence"] == "satellite" else "--", zorder=4))
        ax.set_yticks([n - 1 - i + H / 2 for i in range(n)])
        ax.set_yticklabels([c for c, _ in chroms])
        ax.set_ylim(n - rn - 0.3, max(n, 1) + 0.1)             # tracks top-aligned, same height in a row
        if n:
            ax.set_xlim(0, max(L for _, L in chroms) / 1e6 * 1.01)
        else:
            ax.text(0.5, 0.5, "no centromere calls", transform=ax.transAxes, ha="center", va="center",
                    color="#555555")
        ax.set_xlabel("Position (Mb)")
        ax.set_ylabel(label)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.text(-0.2, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    handles = [Patch(color=color[f], label=f) for f in fams[:len(OKABE_ITO)]]
    if len(fams) > len(OKABE_ITO):
        handles.append(Patch(color=OTHER_CEN, label="other centromeric families"))
    handles += [Patch(color=ASSOC, label="cen-associated satellite"),
                Line2D([], [], color=GENE, lw=0.6, label="gene density (scaled)"),
                Patch(fill=False, ec="black", label="centromere call (dashed: fallback)")]
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return len(panels)
