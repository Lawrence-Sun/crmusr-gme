// Run with ROOT; see README.md. No absolute paths or dependency on dev/.
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace pmu_analysis {
constexpr double protonWeight = 9.704844263299673;
constexpr double heliumWeight = 100.0;
const std::array<std::string, 6> names = {
    "piMuonPmu", "piAntiMuonPmu", "KMuonPmu", "KAntiMuonPmu",
    "KLongMuonPmu", "KLongAntiMuonPmu"};
using Histograms = std::array<std::unique_ptr<TH1D>, 6>;

Histograms makeHistograms(int bins, const std::string &prefix) {
    Histograms histograms;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto name = prefix + names[i];
        histograms[i] = std::make_unique<TH1D>(name.c_str(), name.c_str(), bins, -1, 1);
        histograms[i]->SetDirectory(nullptr);
        histograms[i]->Sumw2();
    }
    return histograms;
}

struct Candidate {
    int pdg, parent;
    float energy, polarization;
    double weight;
};

void fill(Histograms &histograms, const Candidate &muon) {
    int species = -1;
    const auto parent = std::abs(muon.parent);
    if (parent == 211) species = 0;
    if (parent == 321) species = 1;
    if (parent == 130) species = 2;
    if (species < 0) return;
    const bool inEnergyRange = species == 2
        ? muon.energy > 0 && muon.energy < 1e9
        : muon.energy > 100 && muon.energy < 500;
    if (inEnergyRange) {
        histograms[2 * species + (muon.pdg == -13)]->Fill(muon.polarization, muon.weight);
    }
}

void readRaw(const std::string &path, Histograms &histograms,
             bool legacyDeduplicate, bool useWeights) {
    std::unique_ptr<TFile> input(TFile::Open(path.c_str(), "READ"));
    if (!input || input->IsZombie()) throw std::runtime_error("Cannot open input: " + path);
    const double zenithLimit = 1 - std::sqrt(3.0) / 2;
    for (const int muonId : {13, -13}) {
        const std::string treeName = "G4Run0/ReactionChain(" + std::to_string(muonId) + ")";
        if (!input->Get(treeName.c_str())) {
            throw std::runtime_error("Missing tree " + treeName + " in " + path);
        }
        std::cout << "Reading " << path << ": " << treeName << std::endl;
        ROOT::RDataFrame frame(treeName, path);
        const auto columns = frame.GetColumnNames();
        for (const std::string column : {"PDGID", "KillProc", "ParPDGID", "Ek",
                                        "ParTrkID", "EvtID", "Zenith", "Helicity"}) {
            if (std::find(columns.begin(), columns.end(), column) == columns.end()) {
                throw std::runtime_error("Missing column " + column + " in " + treeName);
            }
        }
        if (useWeights && std::find(columns.begin(), columns.end(), "Weight") == columns.end()) {
            throw std::runtime_error("Weighted mode requires column Weight in " + treeName + " in " + path);
        }
        // Cast numeric Float_t/Double_t (and other scalar numeric types) to one
        // callback type. Unit-weight mode deliberately ignores any raw Weight.
        auto events = frame.Define("__pmu_event_weight", useWeights ? "static_cast<double>(Weight)" : "1.0");
        // Historical polarizationEnergyRelation() keeps the LAST selected muon
        // for each (EvtID, ParTrkID) across the entire tree. Merged raw files
        // can repeat event IDs non-contiguously, so event-local flushing is wrong.
        // Store only plotting fields, not parent energies or whole raw rows.
        std::unordered_map<std::uint64_t, Candidate> selectedMuons;
        events.Foreach([&](int pdg, const std::string &kill, int parent, float energy,
                           int parentTrack, int event, float zenith, float polarization, double weight) {
            if (pdg != muonId || kill != "<0|" || !(zenith < zenithLimit)) return;
            if (!std::isfinite(energy) || !std::isfinite(polarization)) {
                throw std::runtime_error("Non-finite energy/polarization in " + treeName);
            }
            if (!std::isfinite(weight) || weight < 0) {
                throw std::runtime_error("Raw Weight must be finite and nonnegative in " + treeName + " in " + path);
            }
            // Keep the final candidate's weight together with its polarization.
            const Candidate candidate{pdg, parent, energy, polarization, weight};
            if (legacyDeduplicate) {
                const auto key = (std::uint64_t(std::uint32_t(event)) << 32)
                               | std::uint32_t(parentTrack);
                selectedMuons.insert_or_assign(key, candidate);
            } else {
                fill(histograms, candidate);
            }
        }, {"PDGID", "KillProc", "ParPDGID", "Ek", "ParTrkID", "EvtID", "Zenith", "Helicity", "__pmu_event_weight"});
        for (const auto &entry : selectedMuons) fill(histograms, entry.second);
    }
}

void run(const std::string &proton, const std::string &helium,
         const std::string &outputDirectory, int bins, bool legacyDeduplicate, bool overwrite,
         bool useWeights, double protonScale, double heliumScale) {
    if (proton.empty() || helium.empty()) {
        throw std::runtime_error("Provide proton and helium ROOT paths; see README.md.");
    }
    if (bins <= 0 || bins > 100000) throw std::runtime_error("bins must be in [1, 100000]");
    if (outputDirectory.empty()) throw std::runtime_error("Output directory must not be empty");
    if (!std::isfinite(protonScale) || !std::isfinite(heliumScale) ||
        protonScale < 0 || heliumScale < 0 || !(protonScale > 0 || heliumScale > 0)) {
        throw std::runtime_error("Source scales must be finite, nonnegative and not both zero");
    }
    namespace fs = std::filesystem;
    const fs::path directory(outputDirectory);
    const auto rootPath = directory / "polarizationData.root";
    const auto csvPath = directory / "polarizationData.csv";
    const auto metadataPath = directory / "provenance.txt";
    for (const auto &output : {rootPath, csvPath, metadataPath}) {
        if (!overwrite && fs::exists(output)) {
            throw std::runtime_error("Refusing to overwrite " + output.string());
        }
        for (const auto &source : {proton, helium}) {
            if (fs::weakly_canonical(output) == fs::weakly_canonical(source)) {
                throw std::runtime_error("Output must not replace an input file");
            }
        }
    }
    // Callbacks mutate histograms/event state, so processing must be sequential.
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
    auto protonHists = makeHistograms(bins, "P_");
    auto heliumHists = makeHistograms(bins, "He_");
    auto totalHists = makeHistograms(bins, "");
    std::cout << "Raw event weights: " << (useWeights ? "Weight branch" : "unit weights") << std::endl;
    readRaw(proton, protonHists, legacyDeduplicate, useWeights);
    readRaw(helium, heliumHists, legacyDeduplicate, useWeights);
    for (std::size_t i = 0; i < names.size(); ++i) {
        totalHists[i]->Add(protonHists[i].get(), protonScale);
        totalHists[i]->Add(heliumHists[i].get(), heliumScale);
        for (int bin = 0; bin <= bins + 1; ++bin) {
            if (!std::isfinite(totalHists[i]->GetBinContent(bin)) ||
                !std::isfinite(totalHists[i]->GetBinError(bin))) {
                throw std::runtime_error("Non-finite weighted histogram " + names[i] + "; check raw weights/source scales");
            }
        }
        const double integral = totalHists[i]->Integral();
        if (!std::isfinite(integral)) throw std::runtime_error("Non-finite histogram integral: " + names[i]);
        if (integral > 0) totalHists[i]->Scale(1.0 / integral);
        else std::cerr << "Warning: empty in-range histogram " << names[i] << '\n';
        if (totalHists[i]->GetBinContent(0) || totalHists[i]->GetBinContent(bins + 1)) {
            std::cerr << "Warning: " << names[i]
                      << " has underflow/overflow (Pmu=+1 is overflow, as in the original).\n";
        }
    }
    fs::create_directories(directory);
    std::ostringstream metadata;
    metadata << std::setprecision(17)
             << "proton=" << fs::absolute(proton).string() << '\n'
             << "helium=" << fs::absolute(helium).string() << '\n'
             << "proton_bytes=" << fs::file_size(proton) << '\n'
             << "helium_bytes=" << fs::file_size(helium) << '\n'
             << "ROOT=" << gROOT->GetVersion() << '\n'
             << "bins=" << bins << "\nrange=[-1,1)\n"
             << "legacy_deduplicate=" << legacyDeduplicate << '\n'
             << "selection=PDGID==+/-13; KillProc==<0|; Zenith<1-sqrt(3)/2\n"
             << "pi_K_energy=100<Ek<500; KLong_energy=0<Ek<1e9\n"
             << "use_raw_weights=" << useWeights << '\n'
             << "raw_weight_column=" << (useWeights ? "Weight" : "none (unit weights)") << '\n'
             << "proton_weight=" << protonScale << "\nhelium_weight=" << heliumScale << '\n'
             << "effective_weight=source scale * raw Weight (or 1 in unit-weight mode)\n"
             << "normalization=sum of regular bins; empty histograms remain zero\n";
    TFile output(rootPath.string().c_str(), overwrite ? "RECREATE" : "CREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create " + rootPath.string());
    for (const auto &histogram : totalHists) {
        if (histogram->Write() <= 0) throw std::runtime_error("Failed to write histogram");
    }
    TNamed provenance("provenance", metadata.str().c_str());
    provenance.Write();
    output.Close();
    std::ofstream csv(csvPath);
    csv << std::setprecision(17) << "histogram,bin,x,y,yerr\n";
    for (const auto &histogram : totalHists) {
        for (int bin = 1; bin <= bins; ++bin) {
            csv << histogram->GetName() << ',' << bin << ',' << histogram->GetBinCenter(bin)
                << ',' << histogram->GetBinContent(bin) << ',' << histogram->GetBinError(bin) << '\n';
        }
    }
    csv.close();
    std::ofstream info(metadataPath);
    info << metadata.str();
    info.close();
    if (!csv || !info) throw std::runtime_error("Failed to write CSV/provenance");
    std::cout << "Saved " << rootPath << "\nSaved " << csvPath
              << "\nSaved " << metadataPath << std::endl;
}
} // namespace pmu_analysis

void analyze_pmu(const char *proton = "", const char *helium = "",
                 const char *outputDirectory = "pmu_output", int bins = 50,
                 bool legacyDeduplicate = true, bool overwrite = false,
                 bool useWeights = false, double protonScale = pmu_analysis::protonWeight,
                 double heliumScale = pmu_analysis::heliumWeight) {
    try {
        pmu_analysis::run(proton, helium, outputDirectory, bins, legacyDeduplicate, overwrite,
                          useWeights, protonScale, heliumScale);
    } catch (const std::exception &error) {
        std::cerr << "PMU analysis failed: " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
