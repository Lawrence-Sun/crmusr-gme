//  ██████╗██████╗ ███╗   ███╗██╗   ██╗███████╗██████╗
// ██╔════╝██╔══██╗████╗ ████║██║   ██║██╔════╝██╔══██╗
// ██║     ██████╔╝██╔████╔██║██║   ██║███████╗██████╔╝
// ██║     ██╔══██╗██║╚██╔╝██║██║   ██║╚════██║██╔══██╗
// ╚██████╗██║  ██║██║ ╚═╝ ██║╚██████╔╝███████║██║  ██║
//  ╚═════╝╚═╝  ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝

#include <iostream>
#include <vector>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TNtuple.h"

auto evaluate(TH1D *exp_hist, TH1D *fit_hist)
{
    double likelihood{0.0};
    double chi2{0.0};

    for (int i = 1; i <= exp_hist->GetNbinsX(); ++i)
    {
        double observed = exp_hist->GetBinContent(i);
        double expected = fit_hist->GetBinContent(i);
        double error = exp_hist->GetBinError(i);

        // 计算likelihood
        if (expected > 0)
        {
            likelihood += observed * std::log(expected) - expected;
        }

        // 计算chi^2
        if (error > 0)
        {
            chi2 += std::pow((observed - expected) / error, 2);
        }
    }

    double value = chi2;
    return value;
}

auto measure_ratio()
{
    int number_of_events{3000000}; // 事例数 30000 per m^2 per year

    // 读取文件并且获取实验测量分布
    TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/measure_ratio/michelDis.root", "READ");
    // TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/measure_ratio/michelDis_bad.root", "READ");

    auto experiment_michel_dis = (TH1D *)eigen_michel_dis.Get("simulated_data");
    experiment_michel_dis->Scale(double(number_of_events) / experiment_michel_dis->Integral());
    for (int i = 1; i <= experiment_michel_dis->GetNbinsX(); ++i)
    {
        double bin_content = experiment_michel_dis->GetBinContent(i);
        double bin_error = std::sqrt(bin_content); // 假设泊松误差
        experiment_michel_dis->SetBinError(i, bin_error);
    }

    std::vector<std::pair<double, double>> ratios{};

    // 遍历ROOT文件中的所有直方图
    TIter next(eigen_michel_dis.GetListOfKeys());
    TKey *key;
    bool is_first_histogram = true;
    while ((key = (TKey *)next()))
    {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TH1"))
        {
            if (is_first_histogram)
            {
                is_first_histogram = false;
                continue; // 跳过第一个直方图
            }

            TH1D *hist = (TH1D *)obj;
            std::cout << "Histogram name: " << hist->GetName() << std::endl;

            // 提取直方图名字的最后一位或两位整数
            std::string hist_name = hist->GetName();
            int extracted_number = -1;
            if (!hist_name.empty())
            {
                size_t pos = hist_name.find_last_not_of("0123456789");
                if (pos != std::string::npos && pos + 1 < hist_name.size())
                {
                    extracted_number = std::stoi(hist_name.substr(pos + 1));
                }
            }
            std::cout << "Extracted number: " << extracted_number << std::endl;

            hist->Scale(double(number_of_events) / hist->Integral());
            for (int i = 1; i <= hist->GetNbinsX(); ++i)
            {
                double bin_content = hist->GetBinContent(i);
                double bin_error = std::sqrt(bin_content); // 假设泊松误差
                hist->SetBinError(i, bin_error);
            }

            double ratio = 0.15 + (extracted_number - 1) * 0.0001; // 假设ratio从0.2开始，每个直方图增加0.001
            double result = evaluate(experiment_michel_dis, hist);

            ratios.push_back(std::make_pair(ratio, result));
            std::cout << "Evaluation result: " << result << std::endl;
        }
    }

    TFile *outputFile = TFile::Open("measure_ratio_mega.root", "RECREATE");
    // TFile *outputFile = TFile::Open("measure_ratio_bad.root", "RECREATE");
    outputFile->cd();

    TTree *tree = new TTree("tree", "Tree with ratios");
    double r1, r2;
    tree->Branch("ratio", &r1, "ratio/D");
    tree->Branch("result", &r2, "result/D");

    for (const auto &pair : ratios)
    {
        r1 = pair.first;
        r2 = pair.second;
        tree->Fill();
        std::cout << "Ratio: " << r1 << ", Result: " << r2 << std::endl;
    }
    tree->Write();
    outputFile->Close();
    eigen_michel_dis.Close();
}