// Direct unweighted sampling of the existing generator's marginal cos(theta) law.
// Not a replacement for its joint auxiliary-variable distribution or weighted Sumw2.
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace michel_direct {
// For p(x)=(1+k*x)/2: F(x)=(x+1)/2+k*(x*x-1)/4.
// Rationalized inverse avoids cancellation and works at k=0 without a branch.
inline double inverseCDF(double uniform, double k) {
    const double y = 2 * uniform - 1;
    return (k + 2 * y) / (1 + std::sqrt(1 + k * k + 2 * k * y));
}

double meanPolarization(const TH1D &histogram) {
    if (histogram.GetXaxis()->GetXmin() != -1 || histogram.GetXaxis()->GetXmax() != 1)
        throw std::runtime_error("Polarization histogram must span [-1,1]");
    double integral = 0, moment = 0;
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        const double value = histogram.GetBinContent(bin);
        if (!std::isfinite(value) || value < 0)
            throw std::runtime_error("Invalid polarization bin content");
        integral += value;
        // GetRandom samples uniformly within each bin, so use the bin center,
        // NOT TH1::GetMean(), which may retain original unbinned fill statistics.
        moment += value * histogram.GetBinCenter(bin);
    }
    if (!(integral > 0) || !std::isfinite(integral) || !std::isfinite(moment))
        throw std::runtime_error("Invalid polarization integral or moment");
    return moment / integral;
}

void run(const std::string &inputPath, const std::string &outputDirectory,
         long long events, unsigned int seed, bool overwrite) {
    namespace fs = std::filesystem;
    if (events < 1 || events > 10000000000LL || seed == 0 || outputDirectory.empty())
        throw std::runtime_error("Require 1..10000000000 events, nonzero seed, nonempty output directory");
    const fs::path directory(outputDirectory);
    const auto rootPath = directory / "michelData.root";
    const auto csvPath = directory / "michelData.csv";
    const auto provenancePath = directory / "michel_provenance.txt";
    for (const auto &path : {rootPath, csvPath, provenancePath}) {
        if (!overwrite && fs::exists(path)) throw std::runtime_error("Refusing to overwrite " + path.string());
        if (fs::weakly_canonical(path) == fs::weakly_canonical(inputPath))
            throw std::runtime_error("Output must not replace input");
    }
    std::unique_ptr<TFile> input(TFile::Open(inputPath.c_str(), "READ"));
    if (!input || input->IsZombie()) throw std::runtime_error("Cannot open " + inputPath);
    const std::array<std::string, 3> names = {"piMuon", "KMuon", "KLongMuon"};
    std::array<double, 3> means;
    for (std::size_t i = 0; i < names.size(); ++i) {
        auto source = dynamic_cast<TH1D *>(input->Get((names[i] + "Pmu").c_str()));
        if (!source) throw std::runtime_error("Missing polarization histogram: " + names[i] + "Pmu");
        means[i] = meanPolarization(*source);
    }
    std::string upstream;
    if (auto source = dynamic_cast<TNamed *>(input->Get("provenance"))) upstream = source->GetTitle();
    input->Close();
    std::array<std::unique_ptr<TH1D>, 3> histograms;
    std::ostringstream provenance;
    provenance << std::setprecision(17) << "input=" << fs::absolute(inputPath).string()
               << "\nROOT=" << gROOT->GetVersion() << "\nseed=" << seed
               << "\nevents_per_component=" << events
               << "\nmode=independent_components; direct unweighted inverse-CDF sampling\n"
               << "model=p(x)=(1+mean(Pmu)*x/3)/2; x=cos(theta)\n"
               << "normalization=none; CSV y/yerr are counts/sqrt(counts), not legacy weighted Sumw2\n"
               << "bins=100; range=[-1,1); ROOT moments use bin centers\n";
    TRandom3 random(seed);
    for (std::size_t i = 0; i < names.size(); ++i) {
        const double k = means[i] / 3;
        std::cout << "Direct sampling " << names[i] << ": " << events
                  << " events; mean(Pmu)=" << std::setprecision(10) << means[i] << std::endl;
        std::array<std::uint64_t, 100> counts{};
        for (long long event = 0; event < events; ++event) {
            const double x = inverseCDF(random.Rndm(), k);
            const int bin = std::clamp(static_cast<int>((x + 1) * 50), 0, 99);
            ++counts[bin];
        }
        auto &histogram = histograms[i];
        histogram = std::make_unique<TH1D>(names[i].c_str(), names[i].c_str(), 100, -1, 1);
        histogram->SetDirectory(nullptr);
        histogram->Sumw2();
        double stats[4] = {double(events), double(events), 0, 0};
        for (int bin = 1; bin <= 100; ++bin) {
            const double count = counts[bin - 1], x = histogram->GetBinCenter(bin);
            histogram->SetBinContent(bin, count);
            histogram->SetBinError(bin, std::sqrt(count));
            stats[2] += count * x;
            stats[3] += count * x * x;
        }
        histogram->PutStats(stats);
        histogram->SetEntries(events);
        provenance << names[i] << ": source=" << names[i] << "Pmu; mean_polarization=" << means[i] << '\n';
    }
    provenance << "upstream_provenance_begin\n" << upstream << "upstream_provenance_end\n";
    fs::create_directories(directory);
    TFile output(rootPath.string().c_str(), overwrite ? "RECREATE" : "CREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create " + rootPath.string());
    std::ofstream csv(csvPath), metadata(provenancePath);
    if (!csv || !metadata) throw std::runtime_error("Cannot create CSV or provenance");
    csv << std::setprecision(17) << "histogram,bin,x,y,yerr\n";
    for (const auto &histogram : histograms) {
        if (histogram->Write() <= 0) throw std::runtime_error("Failed to write histogram");
        for (int bin = 1; bin <= 100; ++bin)
            csv << histogram->GetName() << ',' << bin << ',' << histogram->GetBinCenter(bin)
                << ',' << histogram->GetBinContent(bin) << ',' << histogram->GetBinError(bin) << '\n';
    }
    TNamed record("provenance", provenance.str().c_str());
    if (record.Write() <= 0) throw std::runtime_error("Failed to write provenance");
    output.Close();
    metadata << provenance.str();
    csv.close(); metadata.close();
    if (!csv || !metadata) throw std::runtime_error("Failed to write CSV or provenance");
    std::cout << "Saved direct-sampling ROOT, CSV and provenance in " << directory << std::endl;
}
} // namespace michel_direct

void generate_michel(const char *input = "pmu_output/polarizationData.root",
                     const char *outputDirectory = "pmu_output",
                     long long eventsPerComponent = 100000000LL,
                            unsigned int seed = 20241025, bool overwrite = false) {
    try {
        michel_direct::run(input, outputDirectory, eventsPerComponent, seed, overwrite);
    } catch (const std::exception &error) {
        std::cerr << "generate_michel: " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
