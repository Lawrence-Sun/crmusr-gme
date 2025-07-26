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

auto measure_ratio_bad()
{
    long number_of_events{30000000}; // 事例数 100 个每天 3W  （一台一年3W)

    // 读取文件并且获取实验测量分布
    // TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/measure_ratio/michelDis.root", "READ");
    TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/measure_ratio/michelDis_bad.root", "READ");

    auto experiment_michel_dis = (TH1D *)eigen_michel_dis.Get("simulated_data");
    experiment_michel_dis->Scale(double(number_of_events) / experiment_michel_dis->Integral());
    for (int i = 1; i <= experiment_michel_dis->GetNbinsX(); ++i)
    {
        double bin_content = experiment_michel_dis->GetBinContent(i);
        double bin_error = std::sqrt(bin_content); // 假设泊松误差
        experiment_michel_dis->SetBinError(i, bin_error);
    }

    std::vector<std::tuple<int, int, double>> ratios{};

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

            std::string hist_name = hist->GetName();
            int extracted_number1 = -1, extracted_number2 = -1;
            size_t last_underscore = hist_name.find_last_of('_');
            if (last_underscore != std::string::npos)
            {
                size_t second_last_underscore = hist_name.find_last_of('_', last_underscore - 1);
                if (second_last_underscore != std::string::npos)
                {
                    std::string number1_str = hist_name.substr(second_last_underscore + 1, last_underscore - second_last_underscore - 1);
                    std::string number2_str = hist_name.substr(last_underscore + 1);
                    try
                    {
                        extracted_number1 = std::stoi(number1_str);
                        extracted_number2 = std::stoi(number2_str);
                    }
                    catch (const std::invalid_argument &)
                    {
                        std::cerr << "Error: Unable to extract numbers from histogram name." << std::endl;
                    }
                }
            }
            // std::cout << "Extracted numbers: " << extracted_number1 << ", " << extracted_number2 << std::endl;

            hist->Scale(double(number_of_events) / hist->Integral());
            for (int i = 1; i <= hist->GetNbinsX(); ++i)
            {
                double bin_content = hist->GetBinContent(i);
                double bin_error = std::sqrt(bin_content); // 假设泊松误差
                hist->SetBinError(i, bin_error);
            }

            double ratio_N = 0.15 + (extracted_number1 - 1) * 0.0001; // 假设ratio从0.2开始，每个直方图增加0.001
            double ratio_P = 0.15 + (extracted_number2 - 1) * 0.0001; // 假设ratio从0.2开始，每个直方图增加0.001
            double result = evaluate(experiment_michel_dis, hist);

            ratios.push_back(std::make_tuple(extracted_number1, extracted_number2, result));
            std::cout << "Evaluation result: " << result << std::endl;
        }
    }

    TH2D *hist2D = new TH2D("hist2D", "Ratio vs Result", 100, 0.15, 0.25, 100, 0.15, 0.25);
    for (const auto &pair : ratios)
    {
        int number1 = std::get<0>(pair);
        int number2 = std::get<1>(pair);
        double result = std::get<2>(pair);

        double ratio_N = 0.15 + (number1 - 1) * 0.0001;
        double ratio_P = 0.15 + (number2 - 1) * 0.0001;

        hist2D->Fill(ratio_N, ratio_P, result);
    }

    // // TFile *outputFile = TFile::Open("measure_ratio.root", "RECREATE");
    TFile *outputFile = TFile::Open("measure_ratio_bad.root", "RECREATE");
    outputFile->cd();
    hist2D->Write();
    outputFile->Close();

    // TTree *tree = new TTree("tree", "Tree with ratios");
    // int r1, r2;
    // double r3;
    // tree->Branch("ratioN", &r1, "ratioN/I");
    // tree->Branch("ratioP", &r2, "ratioP/I");
    // tree->Branch("result", &r3, "result/D");

    // for (const auto &tuple : ratios)
    // {
    //     r1 = std::get<0>(tuple);
    //     r2 = std::get<1>(tuple);
    //     r3 = std::get<2>(tuple);
    //     tree->Fill();
    // }
    // tree->Write();
    // outputFile->Close();
    // eigen_michel_dis.Close();
}