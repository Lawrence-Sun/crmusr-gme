"""Synthetic tests for exclusion calculation and the upstream/downstream interface."""

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot

from analysis import mix_templates, required_events, weighted_chi2
from run_end_to_end import HERE, UPSTREAM, call_root


class ExclusionTests(unittest.TestCase):
    def test_exact_threshold_and_zero_discrimination(self):
        pion, kaon = [1, 3], [3, 1]
        ratios = [0.01, 0.05, 0.1, 0.2]
        events = required_events(pion, kaon, ratios)
        self.assertTrue(np.all(np.diff(events) < 0))
        for ratio, count in zip(ratios, events):
            np.testing.assert_allclose(weighted_chi2(pion, mix_templates(pion, kaon, [ratio]), count), [9])
        with self.assertRaises(ValueError):
            required_events(pion, pion, ratios)
        with self.assertRaises(ValueError):
            required_events(pion, kaon, [0, 0.1])

    @unittest.skipUnless((HERE.parents[2] / "dev/2601reviseArticle/sensitivity_revise/eigen_sensitiviy_output.root").exists()
                         and (HERE.parents[2] / "dev/2601reviseArticle/sensitivity_revise/2601_michelElectronWithWeight.root").exists(),
                         "Local historical exclusion data are not tracked")
    def test_historical_exclusion_regression(self):
        from analysis import load_templates
        folder = HERE.parents[2] / "dev/2601reviseArticle/sensitivity_revise"
        pion, kaon, _ = load_templates(folder / "2601_michelElectronWithWeight.root", "piAntiMuon", "KAntiMuon")
        with uproot.open(folder / "eigen_sensitiviy_output.root") as source:
            reference = source["sensitivityResult"].arrays(library="np")
        events = required_events(pion, kaon, reference["ratio"])
        # Legacy ROOT Brent solve and TNtuple float32 storage limit precision.
        np.testing.assert_allclose(events, reference["targetStatistics"], rtol=1e-6)

    def test_templates_to_both_plots_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            templates = folder / "templates.root"
            with uproot.recreate(templates) as output:
                for name, values in {"piMuon": [1., 3.], "KMuon": [3., 1.],
                                     "piAntiMuon": [3., 1.], "KAntiMuon": [1., 3.]}.items():
                    output[name] = (np.array(values), np.array([-1., 0., 1.]))
            result = folder / "result"
            command = [sys.executable, str(HERE / "run_end_to_end.py"),
                       "--templates", str(templates), "--output-dir", str(result)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            manifest = json.loads((result / "manifest.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["parameters"]["ratio_axis"], "physical")
            self.assertTrue((result / "measurement/measure_ratio.pdf").is_file())
            self.assertTrue((result / "exclusion/241013kpiRatioNSensitivity.png").is_file())
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            with uproot.open(result / "exclusion/eigen_sensitiviy_output.root") as source:
                self.assertEqual(source["sensitivityResult"].num_entries, 9)

    @unittest.skipUnless(shutil.which("root"), "ROOT is required for raw-data integration")
    def test_raw_to_both_plots(self):
        # Reuse the upstream raw schema fixture rather than defining a second schema.
        spec = importlib.util.spec_from_file_location("upstream_test", UPSTREAM / "test_pipeline.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        fixture = module.FIXTURE.replace(
            "helicity = source ? 0.5f : -0.5f;",
            "helicity = (species == 211 ? -0.8f : 0.8f);", 1
        )
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            macro = folder / "make_fixture.C"
            macro.write_text(fixture)
            raw = folder / "raw"
            call_root(macro, [str(raw)], folder / "fixture.log")
            output = folder / "output"
            subprocess.run([sys.executable, str(HERE / "run_end_to_end.py"),
                            "--proton", str(raw / "proton.root"),
                            "--helium", str(raw / "helium.root"),
                            "--output-dir", str(output), "--template-events", "10000"],
                           check=True, capture_output=True, text=True)
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(len(manifest["root_commands"]), 2)
            self.assertIn("include_antimuon=true", manifest["template_provenance"])
            with uproot.open(output / "upstream/michelData.root") as source:
                for name in ("piMuon", "KMuon", "KLongMuon", "piAntiMuon", "KAntiMuon", "KLongAntiMuon"):
                    self.assertEqual(source[name].values().sum(), 10000)
            # Keep the first three species and random stream backward compatible.
            legacy = folder / "muon_only"
            call_root(UPSTREAM / "generate_michel.C",
                      [str(output / "upstream/polarizationData.root"), str(legacy), 10000, 20241025],
                      folder / "muon_only.log")
            with uproot.open(legacy / "michelData.root") as old, uproot.open(output / "upstream/michelData.root") as new:
                for name in ("piMuon", "KMuon", "KLongMuon"):
                    np.testing.assert_array_equal(old[name].values(), new[name].values())


if __name__ == "__main__":
    unittest.main()
