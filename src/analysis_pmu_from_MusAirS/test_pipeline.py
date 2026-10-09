#!/usr/bin/env python3
"""Small synthetic ROOT integration tests; never touches the real raw data."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
FIXTURE = r'''
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <string>
void make_fixture(const char *directory) {
    for (int source = 0; source < 2; ++source) {
        std::string path = std::string(directory) + (source ? "/helium.root" : "/proton.root");
        TFile file(path.c_str(), "RECREATE");
        auto folder = file.mkdir("G4Run0");
        folder->cd();
        for (int id : {13, -13}) {
            auto name = "ReactionChain(" + std::to_string(id) + ")";
            TTree tree(name.c_str(), name.c_str());
            int pdg = id, parent = 211, parentTrack = 1, event = 0;
            float energy = 200, zenith = 0, helicity = 0;
            std::string kill = "<0|";
            tree.Branch("PDGID", &pdg, "PDGID/I");
            tree.Branch("ParPDGID", &parent, "ParPDGID/I");
            tree.Branch("ParTrkID", &parentTrack, "ParTrkID/I");
            tree.Branch("EvtID", &event, "EvtID/I");
            tree.Branch("Ek", &energy, "Ek/F");
            tree.Branch("Zenith", &zenith, "Zenith/F");
            tree.Branch("Helicity", &helicity, "Helicity/F");
            tree.Branch("KillProc", &kill);
            for (int species : {211, 321, 130}) {
                ++event;
                parent = id == 13 ? -species : species;
                helicity = -0.1f;
                tree.Fill(); // Same event/parent: overwritten in legacy mode.
                helicity = source ? 0.5f : -0.5f;
                tree.Fill();
                kill = "other"; helicity = 0.8f; tree.Fill(); kill = "<0|";
                zenith = 0.5f; tree.Fill(); zenith = 0;
            }
            // A merged tree may repeat an earlier event non-contiguously.
            event = 1; parent = id == 13 ? -211 : 211;
            helicity = source ? 0.5f : -0.5f;
            tree.Fill();
            tree.Write();
        }
        file.Close();
    }
}
'''


@unittest.skipUnless(shutil.which("root"), "ROOT is required")
class PipelineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix="pmu-test-")
        cls.directory = Path(cls.workspace.name)
        fixture = cls.directory / "make_fixture.C"
        fixture.write_text(FIXTURE)
        cls.root_call(fixture, [str(cls.directory)])
        cls.root_call(HERE / "analyze_pmu.C", [str(cls.directory / "proton.root"),
                     str(cls.directory / "helium.root"), str(cls.directory / "legacy"), 4])

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    @staticmethod
    def root_call(macro, arguments, success=True):
        expression = f"{macro}({','.join(json.dumps(arg) for arg in arguments)})"
        result = subprocess.run(["root", "-l", "-b", "-q", expression],
                                capture_output=True, text=True, timeout=120)
        if success and result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        if not success and result.returncode == 0:
            raise AssertionError("Expected failure: " + result.stdout + result.stderr)
        return result

    def test_weighting_errors_and_legacy_deduplication(self):
        with (self.directory / "legacy/polarizationData.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 6 * 4)
        weight = 9.704844263299673
        for name in {row["histogram"] for row in rows}:
            histogram = [row for row in rows if row["histogram"] == name]
            self.assertAlmostEqual(sum(float(row["y"]) for row in histogram), 1)
            for row in histogram:
                index = int(row["bin"])
                expected = weight / (weight + 100) if index == 2 else (
                    100 / (weight + 100) if index == 4 else 0)
                self.assertAlmostEqual(float(row["y"]), expected)
                self.assertAlmostEqual(float(row["yerr"]), expected)
        self.assertTrue((self.directory / "legacy/polarizationData.root").exists())
        self.assertIn("legacy_deduplicate=1",
                      (self.directory / "legacy/provenance.txt").read_text())

    def test_all_muons_mode(self):
        self.root_call(HERE / "analyze_pmu.C", [str(self.directory / "proton.root"),
                       str(self.directory / "helium.root"), str(self.directory / "all"), 4, False])
        with (self.directory / "all/polarizationData.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        kaon = [row for row in rows if row["histogram"] == "KMuonPmu"]
        self.assertAlmostEqual(float(kaon[1]["y"]),
                               (100 + 2 * 9.704844263299673) / (2 * (100 + 9.704844263299673)))

    def test_refuses_overwrite_and_missing_input(self):
        self.root_call(HERE / "analyze_pmu.C", [str(self.directory / "proton.root"),
                       str(self.directory / "helium.root"), str(self.directory / "legacy"), 4], False)
        self.root_call(HERE / "analyze_pmu.C", ["/nonexistent/pmu.root",
                       str(self.directory / "helium.root"), str(self.directory / "missing")], False)
        self.assertFalse((self.directory / "missing").exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is required")
    def test_plot_and_overwrite_protection(self):
        import sys
        output = self.directory / "figure.pdf"
        command = [sys.executable, str(HERE / "plot_pmu.py"),
                   str(self.directory / "legacy/polarizationData.csv"), "-o", str(output)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(output.read_bytes().startswith(b"%PDF"))
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Refusing to overwrite", result.stderr)

    def test_csv_validation(self):
        spec = importlib.util.spec_from_file_location("plot_pmu", HERE / "plot_pmu.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.read_data(self.directory / "legacy/polarizationData.csv")
        invalid = self.directory / "invalid.csv"
        invalid.write_text("histogram,bin,x,y,yerr\npiMuonPmu,1,0,-1,0\n")
        with self.assertRaises(ValueError):
            module.read_data(invalid)


if __name__ == "__main__":
    unittest.main()
