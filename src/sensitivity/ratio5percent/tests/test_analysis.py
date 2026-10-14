"""Unit tests plus optional regression against the local, untracked dev reference."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
import uproot

from analysis import generate_distributions, mix_templates, weighted_chi2, normalize, scan_distributions
from plot_ratio import find_intersections, plot_scan

REPO_ROOT = Path(__file__).resolve().parents[4]


class AnalysisTests(unittest.TestCase):
    def test_mixture_and_statistics_scaling(self):
        pion, kaon = [1, 3], [3, 1]
        shapes = mix_templates(pion, kaon, [0, 1])
        np.testing.assert_allclose(shapes, [[0.25, 0.75], [0.5, 0.5]])
        chi2 = weighted_chi2(shapes[0], shapes, 100)
        self.assertEqual(chi2[0], 0)
        np.testing.assert_allclose(weighted_chi2(shapes[0], shapes, 200), chi2 * 2)
        self.assertAlmostEqual(chi2[1], 40 / 3)

    def test_invalid_histograms_and_events(self):
        for counts in ([0, 0], [-1, 2], [np.nan, 1]):
            with self.assertRaises(ValueError):
                normalize(counts)
        with self.assertRaises(ValueError):
            weighted_chi2([1, 2], [[1, 2]], 0)
        with self.assertRaises(ValueError):
            mix_templates([1, 2], [1, 2], [-0.1])
        np.testing.assert_allclose(weighted_chi2([0, 1], [[1, 1]], 100), [200 / 3])
        np.testing.assert_allclose(weighted_chi2([0, 1], [[0, 1]], 100), [0])

    def test_intersections(self):
        x = np.linspace(0.01, 0.1, 1001)
        y = ((x - 0.05) / 0.003) ** 2
        _, target, _, roots = find_intersections(x, y)
        np.testing.assert_allclose(roots, [0.05 - 0.003 * np.sqrt(target),
                                           0.05 + 0.003 * np.sqrt(target)], atol=1e-10)
        with self.assertRaises(ValueError):
            find_intersections(x, y, 1000)

    def test_root_roundtrip(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            templates = folder / "templates.root"
            distributions, scan = folder / "distributions.root", folder / "scan.root"
            with uproot.recreate(templates) as output:
                output["piMichelDis"] = (np.array([1., 3.]), np.array([-1., 0., 1.]))
                output["kaMichelDis"] = (np.array([3., 1.]), np.array([-1., 0., 1.]))
            physical = generate_distributions(templates, distributions)
            x, y = scan_distributions(distributions, scan)
            np.testing.assert_allclose(x[[0, 444, -1]], [0.00991, 0.04987, 0.09982])
            self.assertLess(y[444], 1e-20)
            summary = plot_scan(scan, folder)
            self.assertEqual(len(summary["refined_intersections"]), 2)
            self.assertTrue((folder / "measure_ratio.pdf").is_file())
            self.assertTrue((folder / "measure_ratio.png").is_file())
            corrected, _ = scan_distributions(distributions, scan, ratio_axis="physical")
            np.testing.assert_allclose(corrected, physical)

    @unittest.skipUnless((REPO_ROOT / "dev/megaMichelDis.root").exists()
                         and (REPO_ROOT / "dev/2601reviseArticle/measure_ratio_chi2test.root").exists(),
                         "Local dev reference data are not tracked in git")
    def test_notebook_regression(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            generate_distributions(REPO_ROOT / "dev/megaMichelDis.root", folder / "michel.root")
            x, y = scan_distributions(folder / "michel.root", folder / "scan.root")
            with uproot.open(REPO_ROOT / "dev/2601reviseArticle/measure_ratio_chi2test.root") as source:
                reference = source["tree"].arrays(library="np")
            np.testing.assert_allclose(x, reference["ratio"], rtol=0, atol=1e-16)
            np.testing.assert_allclose(y, reference["result"], rtol=1e-10, atol=1e-20)
            _, _, markers, _ = find_intersections(x, y)
            np.testing.assert_allclose(markers, [0.045415, 0.054415], atol=1e-12)


if __name__ == "__main__":
    unittest.main()
