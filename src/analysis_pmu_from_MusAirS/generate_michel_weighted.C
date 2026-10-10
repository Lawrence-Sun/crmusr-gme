// Optional reference generator; run_pipeline.sh uses direct sampling instead.
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace michel_weighted {
// Algebraically identical to GenMichel::getWeight in src/getMichel.cc.
// Cancel cos(alpha) factors analytically to avoid artificial tan/division poles.
inline double weight(double alpha, double polarization, double cosine, double phi) {
    const double s = std::sqrt(std::max(0.0, 1 - cosine * cosine));
    const double spin = std::sqrt(std::max(0.0, 1 - polarization * polarization));
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    const double difference = alpha - phi;
    const double a = ca * std::cos(difference) * s
                   + sa * (cosine * spin + polarization * s * std::sin(difference));
    const double b = -ca * cosine * spin + sa * std::cos(difference) * s
                   - polarization * s * ca * std::sin(difference);
    const double c = polarization * cosine - spin * s * std::sin(difference);
    const double norm = std::sqrt(a * a + b * b + c * c);
    if (!(norm > 0)) throw std::runtime_error("Undefined Michel angular weight");
    return 1 + c / (3 * norm);
}

void run(const std::string &inputPath, const std::string &outputDirectory,
         long long eventsPerComponent, unsigned int seed, bool overwrite) {
    if (eventsPerComponent < 1 || eventsPerComponent > 10000000000LL) {
        throw std::runtime_error("eventsPerComponent must be in [1, 10000000000]");
    }
    if (seed == 0) throw std::runtime_error("Seed must be nonzero for reproducibility");
    if (outputDirectory.empty()) throw std::runtime_error("Output directory must not be empty");
    namespace fs = std::filesystem;
    const fs::path directory(outputDirectory);
    const auto rootPath = directory / "michelData.root";
    const auto csvPath = directory / "michelData.csv";
    const auto metadataPath = directory / "michel_provenance.txt";
    for (const auto &path : {rootPath, csvPath, metadataPath}) {
        if (!overwrite && fs::exists(path)) throw std::runtime_error("Refusing to overwrite " + path.string());
        if (fs::weakly_canonical(path) == fs::weakly_canonical(inputPath)) {
            throw std::runtime_error("Output must not replace input");
        }
    }
    std::unique_ptr<TFile> input(TFile::Open(inputPath.c_str(), "READ"));
    if (!input || input->IsZombie()) throw std::runtime_error("Cannot open " + inputPath);
    std::array<std::unique_ptr<TH1D>, 3> sources;
    const std::array<std::string, 3> sourceNames = {"piMuonPmu", "KMuonPmu", "KLongMuonPmu"};
    for (std::size_t i = 0; i < sources.size(); ++i) {
        auto histogram = dynamic_cast<TH1D *>(input->Get(sourceNames[i].c_str()));
        if (!histogram || histogram->GetXaxis()->GetXmin() != -1 ||
            histogram->GetXaxis()->GetXmax() != 1) {
            throw std::runtime_error("Missing or invalid [-1,1] TH1D: " + sourceNames[i]);
        }
        if (i && histogram->GetNbinsX() != sources[0]->GetNbinsX()) {
            throw std::runtime_error("Polarization histogram binning must match");
        }
        for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            if (!std::isfinite(histogram->GetBinContent(bin)) || histogram->GetBinContent(bin) < 0) {
                throw std::runtime_error("Invalid polarization bin content");
            }
        }
        const double integral = histogram->Integral();
        if (!(integral > 0) || !std::isfinite(integral)) {
            throw std::runtime_error("Empty polarization histogram: " + sourceNames[i]);
        }
        sources[i].reset(static_cast<TH1D *>(histogram->Clone(("source_" + sourceNames[i]).c_str())));
        sources[i]->SetDirectory(nullptr);
        sources[i]->Scale(1.0 / integral);
    }
    std::string upstream;
    if (auto provenance = dynamic_cast<TNamed *>(input->Get("provenance"))) upstream = provenance->GetTitle();
    input->Close();
    // Independent parent components, matching the michelElectron.root plot:
    // piMuonPmu -> piMuon; KMuonPmu -> KMuon; KLongMuonPmu -> KLongMuon.
    const std::array<std::string, 3> names = {"piMuon", "KMuon", "KLongMuon"};
    std::array<std::unique_ptr<TH1D>, 3> electrons;
    TRandom3 random(seed);
    const double twoPi = 2 * std::acos(-1.0);
    for (std::size_t i = 0; i < names.size(); ++i) {
        electrons[i] = std::make_unique<TH1D>(names[i].c_str(), names[i].c_str(), 100, -1, 1);
        electrons[i]->SetDirectory(nullptr);
        electrons[i]->Sumw2();
        std::cout << "Generating " << names[i] << " from " << sourceNames[i]
                  << ": " << eventsPerComponent << " samples" << std::endl;
        for (long long event = 0; event < eventsPerComponent; ++event) {
            const double polarization = sources[i]->GetRandom(&random);
            const double alpha = random.Uniform(0, twoPi);
            const double cosine = random.Uniform(-1, 1);
            const double phi = random.Uniform(0, twoPi);
            const double value = weight(alpha, polarization, cosine, phi);
            if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid angular weight");
            electrons[i]->Fill(cosine, value);
        }
    }
    std::ostringstream metadata;
    metadata << "input=" << fs::absolute(inputPath).string() << '\n'
             << "ROOT=" << gROOT->GetVersion() << "\nseed=" << seed
             << "\nevents_per_component=" << eventsPerComponent
             << "\nmode=independent_components; no parent mixture scaling\n"
             << "sources=piMuon:piMuonPmu; KMuon:KMuonPmu; KLongMuon:KLongMuonPmu\n"
             << "bins=100; range=[-1,1); variable=cos(theta)\n"
             << "model=GenMichel::getWeight algebraic equivalent; weighted angular distribution only\n"
             << "normalization=none; CSV y/yerr are weighted counts\n"
             << "upstream_provenance_begin\n" << upstream << "upstream_provenance_end\n";
    fs::create_directories(directory);
    TFile output(rootPath.string().c_str(), overwrite ? "RECREATE" : "CREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create " + rootPath.string());
    for (std::size_t i = 0; i < names.size(); ++i) {
        // Preserve each normalized input shape for explicit component provenance.
        sources[i]->SetName(sourceNames[i].c_str());
        if (sources[i]->Write() <= 0 || electrons[i]->Write() <= 0) {
            throw std::runtime_error("Failed to write Michel histogram");
        }
    }
    TNamed provenance("provenance", metadata.str().c_str());
    provenance.Write();
    output.Close();
    std::ofstream csv(csvPath);
    csv << std::setprecision(17) << "histogram,bin,x,y,yerr\n";
    for (const auto &histogram : electrons) {
        for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            csv << histogram->GetName() << ',' << bin << ',' << histogram->GetBinCenter(bin)
                << ',' << histogram->GetBinContent(bin) << ',' << histogram->GetBinError(bin) << '\n';
        }
    }
    csv.close();
    std::ofstream info(metadataPath);
    info << metadata.str();
    info.close();
    if (!csv || !info) throw std::runtime_error("Failed to write Michel CSV/provenance");
    std::cout << "Saved " << rootPath << "\nSaved " << csvPath << std::endl;
}
} // namespace michel_weighted

void generate_michel_weighted(const char *input = "pmu_output/polarizationData.root",
                     const char *outputDirectory = "pmu_output", long long eventsPerComponent = 100000000LL,
                     unsigned int seed = 20241025, bool overwrite = false) {
    try {
        michel_weighted::run(input, outputDirectory, eventsPerComponent, seed, overwrite);
    } catch (const std::exception &error) {
        std::cerr << "Michel generation failed: " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
