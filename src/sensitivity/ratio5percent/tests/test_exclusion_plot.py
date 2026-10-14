"""Tests for the additional pure-pion exclusion sensitivity plot."""

from pathlib import Path
import tempfile
import unittest

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import uproot

from plot_exclusion_sensitivity import build_figure, plot_exclusion, OUTPUT_NAME


class ExclusionPlotTests(unittest.TestCase):
    def test_style_and_sorted_data(self):
        fig = build_figure([0.20, 0.01, 0.10], [4e4, 1e7, 2e5])
        try:
            ax = fig.axes[0]
            np.testing.assert_allclose(ax.lines[0].get_xdata(), [0.01, 0.10, 0.20])
            self.assertEqual(ax.get_yscale(), "log")
            self.assertEqual(ax.get_ylim(), (3e4, 2e7))
            self.assertEqual(ax.get_xlim(), (0.01, 0.20))
            np.testing.assert_allclose(ax.collections[0].get_facecolor()[0],
                                       colors.to_rgba("#f6f6ffff"))
            self.assertEqual(ax.texts[0].get_text(), "Kaon not significant")
            self.assertEqual(ax.texts[0].get_color(), "#8a8af5ff")
            self.assertTrue(np.isfinite(ax.texts[0].get_rotation()))
        finally:
            plt.close(fig)

    def test_invalid_input(self):
        for ratios, counts in [([0.01, 0.2], [0, 1]), ([0.1, 0.1], [1, 2]),
                               ([0.01, 0.2], [np.nan, 1]), ([0.01, 0.05], [1, 2])]:
            with self.assertRaises(ValueError):
                build_figure(ratios, counts)

    def test_root_to_images(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            source = folder / "sensitivity.root"
            with uproot.recreate(source) as output:
                output.mktree("sensitivityResult", {
                    "ratio": np.array([0.01, 0.1, 0.2]),
                    "targetStatistics": np.array([1e7, 2e5, 4e4]),
                })
            result = plot_exclusion(source, folder)
            self.assertEqual(result["points"], 3)
            for extension in ("pdf", "png"):
                self.assertGreater((folder / f"{OUTPUT_NAME}.{extension}").stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
