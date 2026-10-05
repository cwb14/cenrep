import numpy as np
import pandas as pd

from cenrep.annotate import ARR_COLS
from cenrep.plot import OKABE_ITO, OTHER_CEN, bin_bp, family_colors, plot_overview, stack_heights


def _panel(name, with_fam=True):
    chroms = {f"{name}01": 3_000_000, f"{name}02": 2_000_000}
    rows = [(0, f"{name}01", 1_000_000, 1_250_000, "+", "Cx_sat01_178", 1400, 0.9, 178)] if with_fam else []
    a = pd.DataFrame(rows, columns=ARR_COLS).assign(cls="centromeric")
    reg = {f"{name}01": dict(evidence="satellite", cen_start=1_000_000, cen_end=1_250_000),
           f"{name}02": dict(evidence="none")} if with_fam else {c: dict(evidence="none") for c in chroms}
    genes = {f"{name}01": np.arange(0, 3_000_000, 20_000, dtype=np.int64)}
    return dict(genome=name, chroms=chroms, arrays=a, regions=reg, genes=genes)


def test_bin_bp_splits_across_bins():
    a = pd.DataFrame({"start": [50_000], "end": [250_000]})
    assert bin_bp(a, 300_000, 100_000).tolist() == [50_000, 100_000, 50_000]


def test_family_colors_fixed_order_never_cycled():
    fams = [f"F{i}" for i in range(8)]
    col = family_colors(fams)
    assert [col[f] for f in fams[:6]] == OKABE_ITO and col["F6"] == col["F7"] == OTHER_CEN


def test_plot_writes_pdf(tmp_path):
    out = tmp_path / "o.pdf"
    plot_overview(out, [_panel(n) for n in ("Aa_A", "Ab_B", "Ac_A", "Ad_B")])
    assert out.stat().st_size > 5_000


def test_plot_without_families(tmp_path):
    out = tmp_path / "o.pdf"
    plot_overview(out, [_panel(n, with_fam=False) for n in ("Aa_A", "Ab_B", "Ac_A", "Ad_B")])
    assert out.exists()


def test_stacked_bars_never_exceed_the_bin():
    """Review M8 / real figures: overlapping arrays of several families summed past 100% of a bin and spilled
    into the next chromosome track. Stacked heights are scaled so their total is at most one bin."""
    layers = [np.array([80_000.0, 10_000.0]), np.array([70_000.0, 5_000.0])]
    out = stack_heights(layers, 100_000)
    assert np.allclose(sum(out), [100_000.0, 15_000.0]) and np.allclose(out[0][1], 10_000.0)


def test_grid_shape_for_any_number_of_panels():
    from cenrep.plot import grid_shape
    assert [grid_shape(n) for n in (1, 2, 4, 5, 7)] == [(1, 1), (1, 2), (2, 2), (3, 2), (4, 2)]


def _contig_panel(name, n=40, called=35):
    chroms = {f"{name}_{i:02d}": 1_000_000 + i * 10_000 for i in range(n)}
    rows = [(i, c, 0, 10_000 * (i + 1), "+", "Cx_sat01_178", 10, 0.9, 178)
            for i, c in enumerate(list(chroms)[:called])]
    a = pd.DataFrame(rows, columns=ARR_COLS).assign(cls="centromeric")
    reg = {c: (dict(evidence="satellite", cen_start=0, cen_end=10_000 * (i + 1)) if i < called
               else dict(evidence="none")) for i, c in enumerate(chroms)}
    return dict(genome=name, chroms=chroms, arrays=a, regions=reg, genes=None, borrowed=True)


def test_many_sequences_show_the_called_ones_with_most_satellite():
    from cenrep.plot import panel_tracks
    pn = _contig_panel("Ctg")
    tracks, label = panel_tracks(pn, max_tracks=20)
    assert len(tracks) == 20 and {c for c, _ in tracks} == {f"Ctg_{i:02d}" for i in range(15, 35)}
    assert [L for _, L in tracks] == sorted((L for _, L in tracks), reverse=True)
    assert "20 of 35" in label and "Ctg" in label
    few = _panel("Aa_A")
    assert panel_tracks(few, max_tracks=20) == (list(few["chroms"].items()), "Aa_A")
    chrom_level = dict(pn, borrowed=False)                      # 40 chromosomes: all shown (review M2)
    assert len(panel_tracks(chrom_level, max_tracks=20)[0]) == 40


def test_plot_seven_panels_including_contig_genomes(tmp_path):
    out = tmp_path / "o.pdf"
    panels = [_panel(n) for n in ("Aa_A", "Ab_B", "Ac_A", "Ad_B", "Ae_A")]
    panels += [_contig_panel("Ctg"), _contig_panel("Nocall", called=0)]
    assert plot_overview(out, panels) == 7
    assert out.stat().st_size > 5_000


def test_track_height_is_capped_for_tall_figures():
    from cenrep.plot import track_height
    assert track_height([7, 7]) == 0.32 and track_height([7, 7, 20, 20]) * 54 <= 12.0 + 1e-9
