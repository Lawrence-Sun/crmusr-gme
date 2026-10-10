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

    def test_michel_generation_is_reproducible(self):
        arguments = [str(self.directory / "legacy/polarizationData.root")]
        for name in ("michel_a", "michel_b"):
            self.root_call(HERE / "generate_michel.C",
                           arguments + [str(self.directory / name), 10000, 12345])
        first = self.directory / "michel_a/michelData.csv"
        self.assertEqual(first.read_bytes(), (self.directory / "michel_b/michelData.csv").read_bytes())
        with first.open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 300)
        self.assertEqual({row["histogram"] for row in rows}, {"piMuon", "KMuon", "KLongMuon"})
        for name in ("piMuon", "KMuon", "KLongMuon"):
            histogram = [row for row in rows if row["histogram"] == name]
            total = sum(float(row["y"]) for row in histogram)
            self.assertGreater(total, 2 / 3 * 10000)
            self.assertLess(total, 4 / 3 * 10000)
            self.assertEqual(total, 10000)
            self.assertTrue(all(float(row["yerr"]) > 0 for row in histogram))
        metadata = (self.directory / "michel_a/michel_provenance.txt").read_text()
        self.assertIn("seed=12345", metadata)
        self.assertIn("events_per_component=10000", metadata)
        self.assertIn("mode=independent_components", metadata)
        self.assertNotIn("100:15:5", metadata)
        self.assertIn("legacy_deduplicate=1", metadata)
        self.root_call(HERE / "generate_michel.C",
                       arguments + [str(self.directory / "michel_a"), 10000, 12345], False)
        self.root_call(HERE / "generate_michel.C",
                       arguments + [str(self.directory / "invalid_seed"), 100, 0], False)

    def test_michel_component_isolation(self):
        # Changing only K polarization must not change the pi or KLong curves.
        macro = self.directory / "change_kaon.C"
        macro.write_text(r'''
#include <TFile.h>
#include <TH1D.h>
void change_kaon(const char *original, const char *modified) {
    TFile source(original), target(modified,"RECREATE");
    for(auto name:{"piMuonPmu","KMuonPmu","KLongMuonPmu"}) {
        auto h=(TH1D*)source.Get(name);
        if(std::string(name)=="KMuonPmu") {h->Reset();h->Fill(0.9);}
        target.cd();h->Write();
    }
}
''')
        original = str(self.directory / "legacy/polarizationData.root")
        modified = str(self.directory / "modified_kaon.root")
        self.root_call(macro, [original, modified])
        data = []
        for source, name in ((original, "isolation_original"), (modified, "isolation_modified")):
            output = self.directory / name
            self.root_call(HERE / "generate_michel.C", [source, str(output), 10000, 42])
            with (output / "michelData.csv").open() as stream:
                data.append(list(csv.DictReader(stream)))
        for name in ("piMuon", "KLongMuon"):
            self.assertEqual([row for row in data[0] if row["histogram"] == name],
                             [row for row in data[1] if row["histogram"] == name])
        self.assertNotEqual([row for row in data[0] if row["histogram"] == "KMuon"],
                            [row for row in data[1] if row["histogram"] == "KMuon"])

    def test_direct_michel_sampling(self):
        inverse_test = self.directory / "check_inverse.C"
        inverse_test.write_text('#include ' + json.dumps(str(HERE / "generate_michel.C")) + r'''
void check_inverse() {
    for(double k:{-1.0/3,0.0,1.0/3}) {
        double previous=-2;
        for(int i=0;i<=1000;i++) {
            double u=i/1000.0, x=michel_direct::inverseCDF(u,k);
            double cdf=(x+1)/2+k*(x*x-1)/4;
            if(!std::isfinite(x)||x < -1-1e-14||x > 1+1e-14||x<previous||fabs(cdf-u)>1e-14)
                gSystem->Exit(1);
            previous=x;
        }
    }
}
''')
        self.root_call(inverse_test, [])
        polarization = self.directory / "legacy/polarizationData.csv"
        with polarization.open() as stream:
            polar_rows = list(csv.DictReader(stream))
        output = self.directory / "direct_a"
        arguments = [str(self.directory / "legacy/polarizationData.root")]
        for name, macro in (("direct_a", "generate_michel.C"),
                            ("direct_b", "generate_michel_direct.C")):
            self.root_call(HERE / macro, arguments + [str(self.directory / name), 100000, 42])
        self.assertEqual((output / "michelData.csv").read_bytes(),
                         (self.directory / "direct_b/michelData.csv").read_bytes())
        with (output / "michelData.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        for name in ("piMuon", "KMuon", "KLongMuon"):
            histogram = [row for row in rows if row["histogram"] == name]
            source = [row for row in polar_rows if row["histogram"] == name + "Pmu"]
            mean = sum(float(row["x"]) * float(row["y"]) for row in source) / sum(float(row["y"]) for row in source)
            self.assertEqual(sum(float(row["y"]) for row in histogram), 100000)
            chi_square = 0
            for row in histogram:
                count, error = float(row["y"]), float(row["yerr"])
                self.assertTrue(count.is_integer())
                self.assertAlmostEqual(error * error, count, places=8)
                expected = 1000 * (1 + mean * float(row["x"]) / 3)
                chi_square += (count - expected) ** 2 / expected
            self.assertLess(chi_square, 180)
        metadata = (output / "michel_provenance.txt").read_text()
        self.assertIn("direct unweighted inverse-CDF", metadata)
        self.root_call(HERE / "generate_michel_direct.C", arguments + [str(output), 100000, 42], False)
        for count, seed in ((0, 42), (10000000001, 42), (100, 0)):
            self.root_call(HERE / "generate_michel_direct.C",
                           arguments + [str(self.directory / "invalid_direct"), count, seed], False)
        self.root_call(HERE / "generate_michel_direct.C", ["/nonexistent/root", str(output)], False)
        self.root_call(HERE / "generate_michel_direct.C", ["/nonexistent/root", str(self.directory / "missing_direct")], False)

    def test_direct_and_weighted_shapes_agree(self):
        datasets = []
        for method in ("generate_michel.C", "generate_michel_weighted.C"):
            output = self.directory / ("shape_" + method)
            self.root_call(HERE / method,
                           [str(self.directory / "legacy/polarizationData.root"), str(output), 100000, 123])
            with (output / "michelData.csv").open() as stream:
                datasets.append(list(csv.DictReader(stream)))
        for name in ("piMuon", "KMuon", "KLongMuon"):
            first = [row for row in datasets[0] if row["histogram"] == name]
            second = [row for row in datasets[1] if row["histogram"] == name]
            chi_square = sum((float(a["y"]) - float(b["y"])) ** 2 /
                             (float(a["yerr"]) ** 2 + float(b["yerr"]) ** 2)
                             for a, b in zip(first, second))
            self.assertLess(chi_square, 200)

    def test_michel_weight_matches_original_formula(self):
        macro = self.directory / "check_weight.C"
        macro.write_text('#include ' + json.dumps(str(HERE / "generate_michel_weighted.C")) + r'''
void check_weight() {
    TRandom3 random(42);
    for(int i=0; i<10000; ++i) {
        double alpha=random.Uniform(0, 6.283185307179586);
        double p=random.Uniform(-1, 1), c=random.Uniform(-1, 1);
        double phi=random.Uniform(0, 6.283185307179586);
        double ca=cos(alpha), sa=sin(alpha), s=sqrt(1-c*c), spin=sqrt(1-p*p);
        // The unsimplified expression used by GenMichel::getWeight.
        double numerator=p*c+pow(abs(ca),2)/ca*spin*s*(sin(phi)-cos(phi)*tan(alpha));
        double a=ca*cos(alpha-phi)*s+pow(abs(ca),2)/ca*(c*spin+p*s*sin(alpha-phi))*tan(alpha);
        double b2=pow(abs(ca),4)*pow(-c/ca*spin+pow(abs(1/ca),2)*cos(alpha-phi)*sa*s
                                    +p*s*(sin(phi)-cos(phi)*tan(alpha)),2);
        double expected=1+numerator/(3*sqrt(a*a+b2+numerator*numerator));
        double actual=michel_weighted::weight(alpha,p,c,phi);
        // Before marginalization, the vector norm is identically one.
        double simplified=1+(p*c-spin*s*sin(alpha-phi))/3;
        if(abs(actual-simplified)>1e-12) gSystem->Exit(1);
        if(abs(actual-expected)>1e-11 || actual<2.0/3 || actual>4.0/3) gSystem->Exit(1);
    }
    // The simplified expression also remains finite at cos(alpha)=0.
    if(!std::isfinite(michel_weighted::weight(acos(-1.0)/2,0.5,0.2,0.7))) gSystem->Exit(1);
}
''')
        self.root_call(macro, [])

    @unittest.skipUnless(importlib.util.find_spec("matplotlib") and importlib.util.find_spec("scipy"),
                         "matplotlib and scipy are required")
    def test_full_one_command_pipeline(self):
        import os
        import sys
        output = self.directory / "full pipeline"
        command = ["bash", str(HERE / "run_pipeline.sh"), str(self.directory / "proton.root"),
                   str(self.directory / "helium.root"), str(output), "10000", "42"]
        environment = dict(os.environ, PYTHON=sys.executable)
        result = subprocess.run(command, capture_output=True, text=True, env=environment, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for filename in ("combinedPmu_revise.pdf", "michelComponents.pdf"):
            self.assertTrue((output / filename).read_bytes().startswith(b"%PDF"))
        original = (output / "michelData.csv").read_bytes()
        result = subprocess.run(command, capture_output=True, text=True, env=environment, timeout=120)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Refusing to overwrite", result.stderr)
        self.assertEqual(original, (output / "michelData.csv").read_bytes())
        # An existing final figure is detected before any intermediate is written.
        blocked = self.directory / "blocked"
        blocked.mkdir()
        (blocked / "michelComponents.pdf").write_bytes(b"keep")
        command[4] = str(blocked)
        result = subprocess.run(command, capture_output=True, text=True, env=environment, timeout=120)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((blocked / "polarizationData.root").exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib") and importlib.util.find_spec("scipy"),
                         "matplotlib and scipy are required")
    def test_michel_periodic_interpolation(self):
        import numpy as np
        import sys
        sys.path.insert(0, str(HERE))
        try:
            import plot_michel
            rows = [(i + 1, -0.99 + i * 0.02, 1 + 0.1 * (-0.99 + i * 0.02), 0.1)
                    for i in range(100)]
            theta, values = plot_michel.polar_curve(rows, 501)
            self.assertAlmostEqual(values[0], values[-1])
            np.testing.assert_allclose(values, values[::-1], atol=1e-12)
            self.assertAlmostEqual(theta[-1], 2 * np.pi)
            self.assertTrue(np.all(values > 0))
        finally:
            sys.path.pop(0)

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
