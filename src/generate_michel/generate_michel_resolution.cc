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
#include "TF1.h"
#include "TVector3.h"

using namespace std;

// 现在我们首先需要确定的是最终的数据格式，应为这个数据格式会决定我们怎样处理数据并且后续的分析过程要如何进行
// 首先我们现在的两个变量分别是宇生缪子的动量分辨以及Michel电子的空间分辨
// 我们最终要对比的是这两个分辨率的劣化会如何影响到最终的结果

auto add_resolution_to_alpha(TRandom &random_generator, double alpha_michel_muon, double LGA_resolution, double PDR_resolution)
{
    // TRandom random_generator;
    // random_generator.SetSeed(time(0));

    double MDD1_h{200.0}; // 第一层距离靶中心的高度
    double MDD2_h{400.0}; // 第二层距离靶中心的高度
    double PDR_R{150.0};  // PDR的半径

    // 首先我们是要按照天顶角分布的形式选定一个天顶角
    TF1 zenith_angle_distribution("zenith_angle_distribution", "x*x", 0, 1);
    double zenith_angle = zenith_angle_distribution.GetRandom();
    // 然后我们需要得到的这个粒子击中靶的坐标
    double hit_x = random_generator.Uniform(-100, 100);
    double hit_y = random_generator.Uniform(-100, 100);
    double muon_azimuth = random_generator.Uniform(0, 2 * M_PI);

    // 接下来要求解的就是缪子击中两层的坐标
    double MDD1_hit_x = hit_x + MDD1_h / zenith_angle * std::sqrt(1 - zenith_angle * zenith_angle) * std::cos(muon_azimuth);
    double MDD1_hit_y = hit_y + MDD1_h / zenith_angle * std::sqrt(1 - zenith_angle * zenith_angle) * std::sin(muon_azimuth);
    double MDD2_hit_x = hit_x + MDD2_h / zenith_angle * std::sqrt(1 - zenith_angle * zenith_angle) * std::cos(muon_azimuth);
    double MDD2_hit_y = hit_y + MDD2_h / zenith_angle * std::sqrt(1 - zenith_angle * zenith_angle) * std::sin(muon_azimuth);

    // 然后我们要生成Michel电子击中PDR的坐标
    double gamma = random_generator.Uniform(0, 2 * M_PI); // 定义圆锥曲面相对方位角
    // 我们需要选取和缪子动量方向的正交基
    TVector3 muon_momentum_direction(sqrt(1 - pow(zenith_angle, 2)) * cos(muon_azimuth),
                                     sqrt(1 - pow(zenith_angle, 2)) * sin(muon_azimuth),
                                     zenith_angle);
    TVector3 u, v; // 正交基向量
    if (zenith_angle != 0)
    {
        u = TVector3(-sqrt(1 - pow(zenith_angle, 2)) * sin(muon_azimuth),
                     sqrt(1 - pow(zenith_angle, 2)) * cos(muon_azimuth),
                     0);
        u = u.Unit();                         // 归一化
        v = muon_momentum_direction.Cross(u); // 计算正交基
        v = v.Unit();                         // 归一化
    }
    else
    {
        u = TVector3(1, 0, 0);
        v = muon_momentum_direction.Cross(u);
        v = v.Unit();
    }

    double b_x = muon_momentum_direction.X() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.X() * cos(gamma) + v.X() * sin(gamma));
    double b_y = muon_momentum_direction.Y() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.Y() * cos(gamma) + v.Y() * sin(gamma));
    double b_z = muon_momentum_direction.Z() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.Z() * cos(gamma) + v.Z() * sin(gamma));

    double t_gamma = (sqrt(pow(hit_x * b_x + hit_y * b_y, 2) - (pow(b_x, 2) + pow(b_y, 2)) * (pow(hit_x, 2) + pow(hit_y, 2) - pow(PDR_R, 2))) - (hit_x * b_x + hit_y * b_y)) / (b_x * b_x + b_y * b_y);

    double PDR_hit_x = hit_x + t_gamma * (muon_momentum_direction.X() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.X() * cos(gamma) + v.X() * sin(gamma)));
    double PDR_hit_y = hit_y + t_gamma * (muon_momentum_direction.Y() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.Y() * cos(gamma) + v.Y() * sin(gamma)));
    double PDR_hit_z = t_gamma * (muon_momentum_direction.Z() * alpha_michel_muon + sqrt(1 - pow(zenith_angle, 2)) * (u.Z() * cos(gamma) + v.Z() * sin(gamma)));

    // 到现在为止我们获得了缪子击中两层LGA的坐标以及Michel电子击中PDR的坐标
    // 接下来就是要给LGA的击中点和PDR的击中点加上一个分辨率
    MDD1_hit_x += random_generator.Gaus(0, LGA_resolution);
    MDD1_hit_y += random_generator.Gaus(0, LGA_resolution);
    MDD2_hit_x += random_generator.Gaus(0, LGA_resolution);
    MDD2_hit_y += random_generator.Gaus(0, LGA_resolution);
    PDR_hit_x += random_generator.Gaus(0, PDR_resolution);
    PDR_hit_y += random_generator.Gaus(0, PDR_resolution);
    PDR_hit_z += random_generator.Gaus(0, PDR_resolution);
    // evaluate hit point with resolution
    auto t = (MDD1_h) / (MDD1_h - MDD2_h);
    auto reconstructed_target_hit_x = MDD1_hit_x + t * (MDD2_hit_x - MDD1_hit_x);
    auto reconstructed_target_hit_y = MDD1_hit_y + t * (MDD2_hit_y - MDD1_hit_y);
    // 定义缪子动量方向的向量和Michel电子动量方向的向量
    TVector3 reconstruct_muon_momentum(MDD1_hit_x - MDD2_hit_x, MDD1_hit_y - MDD2_hit_y, MDD1_h - MDD2_h);
    TVector3 reconstruct_michel_momentum(PDR_hit_x - reconstructed_target_hit_x, PDR_hit_y - reconstructed_target_hit_y, PDR_hit_z);

    double reconstructed_cos_alpha{};
    if (abs(PDR_hit_z) < 150)
    {
        reconstructed_cos_alpha = reconstruct_michel_momentum.Dot(reconstruct_muon_momentum) /
                                  (reconstruct_michel_momentum.Mag() * reconstruct_muon_momentum.Mag());
    }
    else
    {
        reconstructed_cos_alpha = NAN;
    }

    return reconstructed_cos_alpha;
}

// 我们可以得到一个带有分辨率的Michel电子空间分布
auto get_single_resolution_michel(double LGA_resolution, double PDR_resolution)
{
    TFile eigen_michel_dis("/Users/sunmingchen/codeDIR/GME/dev/megaMichelDis.root", "READ");
    auto pi_N_Michel = (TH1D *)eigen_michel_dis.Get("piMichelDis");
    auto pi_P_Michel = (TH1D *)eigen_michel_dis.Get("piMichelDisAnti");
    auto ka_N_Michel = (TH1D *)eigen_michel_dis.Get("kaMichelDis");
    auto ka_P_Michel = (TH1D *)eigen_michel_dis.Get("kaMichelDisAnti");

    // int random_number{1000000000}; // 1e9
    int random_number{10000}; // 1e3

    TH1D *pi_N_Michel_resolution = new TH1D("pi_N_Michel_resolution", "pi_N_Michel_resolution", 100, -1, 1);
    TH1D *pi_P_Michel_resolution = new TH1D("pi_P_Michel_resolution", "pi_P_Michel_resolution", 100, -1, 1);
    TH1D *ka_N_Michel_resolution = new TH1D("ka_N_Michel_resolution", "ka_N_Michel_resolution", 100, -1, 1);
    TH1D *ka_P_Michel_resolution = new TH1D("ka_P_Michel_resolution", "ka_P_Michel_resolution", 100, -1, 1);

    for (int i = 0; i < random_number; i++)
    {
        TRandom random_generate;
        random_generate.SetSeed(time(0) + i);
        double pi_N_cos_alpha = pi_N_Michel->GetRandom();
        double pi_P_cos_alpha = pi_P_Michel->GetRandom();
        double ka_N_cos_alpha = ka_N_Michel->GetRandom();
        double ka_P_cos_alpha = ka_P_Michel->GetRandom();

        double pi_N_reconstruct_cos_alpha = add_resolution_to_alpha(random_generate, pi_N_cos_alpha, LGA_resolution, PDR_resolution);
        double pi_P_reconstruct_cos_alpha = add_resolution_to_alpha(random_generate, pi_P_cos_alpha, LGA_resolution, PDR_resolution);
        double ka_N_reconstruct_cos_alpha = add_resolution_to_alpha(random_generate, ka_N_cos_alpha, LGA_resolution, PDR_resolution);
        double ka_P_reconstruct_cos_alpha = add_resolution_to_alpha(random_generate, ka_P_cos_alpha, LGA_resolution, PDR_resolution);

        if (pi_N_reconstruct_cos_alpha != NAN)
        {
            pi_N_Michel_resolution->Fill(pi_N_reconstruct_cos_alpha);
        }
        if (pi_P_reconstruct_cos_alpha != NAN)
        {
            pi_P_Michel_resolution->Fill(pi_P_reconstruct_cos_alpha);
        }
        if (ka_N_reconstruct_cos_alpha != NAN)
        {
            ka_N_Michel_resolution->Fill(ka_N_reconstruct_cos_alpha);
        }
        if (ka_P_reconstruct_cos_alpha != NAN)
        {
            ka_P_Michel_resolution->Fill(ka_P_reconstruct_cos_alpha);
        }

        // if (i % (random_number / 100) == 0)
        // {
        //     std::cout << "\rProgress: " << (i * 100 / random_number) << "% " << std::flush;
        // }
        // if (i == random_number - 1)
        // {
        //     std::cout << "\rProgress: 100% " << std::endl;
        // }
    }

    TFile output_file("single_resolution_michel.root", "RECREATE");
    pi_N_Michel_resolution->Write();
    pi_P_Michel_resolution->Write();
    ka_N_Michel_resolution->Write();
    ka_P_Michel_resolution->Write();
    output_file.Close();
}

// 计算两个直方图之间chi2的函数
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

// 对于带有分辨率的Michel电子空间分布，接下来就是要计算对应的chi2值并且得到对应的1sigma区间的宽度
auto generate_single_resolution_chi2(double LGA_resolution, double PDR_resolution)
{
    get_single_resolution_michel(LGA_resolution, PDR_resolution);
    // 指定事例数
    int number_of_events{300000};

    TFile eigen_michel_dis("single_resolution_michel.root", "READ");
    auto pi_P_Michel_resolution = (TH1D *)eigen_michel_dis.Get("pi_P_Michel_resolution");
    auto ka_P_Michel_resolution = (TH1D *)eigen_michel_dis.Get("ka_P_Michel_resolution");

    pi_P_Michel_resolution->Scale(1 / pi_P_Michel_resolution->Integral());
    ka_P_Michel_resolution->Scale(1 / ka_P_Michel_resolution->Integral());

    // 产生假想实验下的分布
    TH1D *experimental_michel_distribution = new TH1D("experimental_michel_distribution", "experimental_michel_distribution", 100, -1, 1);
    experimental_michel_distribution->Add(pi_P_Michel_resolution, 1.0);
    experimental_michel_distribution->Add(ka_P_Michel_resolution, 0.2);
    experimental_michel_distribution->Scale(number_of_events / experimental_michel_distribution->Integral());
    for (int i = 1; i <= experimental_michel_distribution->GetNbinsX(); i++)
    {
        double bin_content = experimental_michel_distribution->GetBinContent(i);
        double bin_error = sqrt(bin_content); // 这里假设误差为平方根
        experimental_michel_distribution->SetBinError(i, bin_error);
    }

    vector<pair<double, double>> chi2_values; // 存储所有的 ratio 对应的 chi2 值

    // 计算chi2值
    for (int i = 0; i < 1000; i++)
    {
        double test_ratio = 0.15 + i * 0.001; // 测试所用的 k-pi ratio
        TH1D *test_distribution = new TH1D(Form("test_distribution_%d", i), Form("test_distribution_%d", i), 100, -1, 1);
        test_distribution->Add(pi_P_Michel_resolution, 1);
        test_distribution->Add(ka_P_Michel_resolution, test_ratio);
        test_distribution->Scale(number_of_events / test_distribution->Integral());
        for (int i = 1; i <= test_distribution->GetNbinsX(); i++)
        {
            double bin_content = test_distribution->GetBinContent(i);
            double bin_error = sqrt(bin_content); // 这里假设误差为平方根
            test_distribution->SetBinError(i, bin_error);
        }

        double chi2_value = evaluate(experimental_michel_distribution, test_distribution);
        chi2_values.emplace_back(test_ratio, chi2_value);
        delete test_distribution; // 释放内存
    }

    auto min_chi2 = min_element(chi2_values.begin(), chi2_values.end(),
                                [](const pair<double, double> &a, const pair<double, double> &b)
                                {
                                    return a.second < b.second;
                                });

    double min_ratio{};
    double max_ratio{};
    auto item_with_chi2 = chi2_values.begin();
    // 找到 min_chi2 附近的 1-sigma 区间
    bool found_min = false, found_max = false;
    for (size_t i = 1; i < chi2_values.size(); ++i)
    {
        double prev_chi2 = chi2_values[i - 1].second;
        double curr_chi2 = chi2_values[i].second;
        double threshold = min_chi2->second + 1.0;

        // min_ratio: 从左往右，第一个穿过 threshold 的点
        if (!found_min && prev_chi2 > threshold && curr_chi2 <= threshold)
        {
            // 线性插值
            double prev_ratio = chi2_values[i - 1].first;
            double curr_ratio = chi2_values[i].first;
            min_ratio = prev_ratio + (threshold - prev_chi2) * (curr_ratio - prev_ratio) / (curr_chi2 - prev_chi2);
            found_min = true;
        }
        // max_ratio: 从右往左，第一个穿过 threshold 的点
        if (!found_max && prev_chi2 < threshold && curr_chi2 >= threshold)
        {
            double prev_ratio = chi2_values[i - 1].first;
            double curr_ratio = chi2_values[i].first;
            max_ratio = prev_ratio + (threshold - prev_chi2) * (curr_ratio - prev_ratio) / (curr_chi2 - prev_chi2);
            found_max = true;
        }
    }

    eigen_michel_dis.Close();
    double one_sigma_width = max_ratio - min_ratio;
    // cout << "LGA resolution: " << LGA_resolution << ", PDR resolution: " << PDR_resolution
    //      << ", One sigma width: " << one_sigma_width << endl;
    return one_sigma_width;
}

// auto generate_michel_resolution()
int main()
{
    // 现在要生成的对应分辨率，LGA的位置分辨应该是从3mm到30mm，PDR的分辨应该是从10mm到50mm
    int grid_size{10};                     // 分辨率的网格大小
    auto LGA_setp = (30 - 3) / grid_size;  // 从3mm到30mm，分成50个点
    auto PDR_setp = (50 - 10) / grid_size; // 从10mm到50mm，分成50个点
    // 生成放置结果的二维直方图
    TH2D *resolution = new TH2D("resolution", "resolution", grid_size, 3, 30, grid_size, 10, 50);
    for (int i = 0; i < grid_size; i++)
    {
        double LGA_resolution = 3 + i * LGA_setp; // 从3mm到30mm
        for (int j = 0; j < grid_size; j++)
        {
            double PDR_resolution = 10 + j * PDR_setp; // 从10mm到50mm
            auto one_sigma_width = generate_single_resolution_chi2(LGA_resolution, PDR_resolution);
            resolution->SetBinContent(i + 1, j + 1, one_sigma_width);

            // 打印进度条
            int total = grid_size * grid_size;
            int current = i * grid_size + j + 1;
            int barWidth = grid_size;
            double progress = (double)current / total;
            std::cout << "\r[";
            int pos = barWidth * progress;
            for (int k = 0; k < barWidth; ++k)
                std::cout << (k < pos ? "=" : " ");
            std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")" << std::flush;
        }
    }

    TFile output_file("michel_resolution.root", "RECREATE");
    resolution->Write();
    output_file.Close();

    generate_single_resolution_chi2(10, 30);

    return 0;
}