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

auto generate_michel()
{
    // 读取文件获取四个分布
    TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/megaMichelDis.root", "READ");
    if (!eigen_michel_dis.IsOpen())
    {
        std::cerr << "Error: Could not open file megaMichelDis.root!" << std::endl;
        return;
    }

    auto pi_N_Michel = (TH1D *)eigen_michel_dis.Get("piMichelDis");
    auto pi_P_Michel = (TH1D *)eigen_michel_dis.Get("piMichelDisAnti");
    auto ka_N_Michel = (TH1D *)eigen_michel_dis.Get("kaMichelDis");
    auto ka_P_Michel = (TH1D *)eigen_michel_dis.Get("kaMichelDisAnti");

    pi_N_Michel->Scale(1 / pi_N_Michel->Integral());
    pi_P_Michel->Scale(1 / pi_P_Michel->Integral());
    ka_N_Michel->Scale(1 / ka_N_Michel->Integral());
    ka_P_Michel->Scale(1 / ka_P_Michel->Integral());

    TFile outputFile("michelDis.root", "RECREATE");
    outputFile.cd();

    // 生成假设ratio为0.2的分布
    TH1D simulated_data("simulated_data", "simulated_data", 100, -1, 1);
    simulated_data.Add(ka_N_Michel, pi_N_Michel, 0.2, 1);
    simulated_data.Write();

    // 从0.2到0.3均匀分布100个点
    int num_points = 1000;
    double start_ratio = 0.15;
    double end_ratio = 0.25;
    double step = (end_ratio - start_ratio) / (num_points - 1);

    for (int i = 0; i < num_points; ++i)
    {
        double ratio = start_ratio + i * step;
        TH1D simulated_data_point(Form("simulated_data_%d", i), "simulated_data", 100, -1, 1);
        simulated_data_point.Add(ka_N_Michel, pi_N_Michel, ratio, 1);
        simulated_data_point.Write();
    }
}