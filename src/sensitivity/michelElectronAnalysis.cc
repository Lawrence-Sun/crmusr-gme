//  ██████╗██████╗ ███╗   ███╗██╗   ██╗███████╗██████╗                ██████╗ ███╗   ███╗███████╗
// ██╔════╝██╔══██╗████╗ ████║██║   ██║██╔════╝██╔══██╗              ██╔════╝ ████╗ ████║██╔════╝
// ██║     ██████╔╝██╔████╔██║██║   ██║███████╗██████╔╝    █████╗    ██║  ███╗██╔████╔██║█████╗
// ██║     ██╔══██╗██║╚██╔╝██║██║   ██║╚════██║██╔══██╗    ╚════╝    ██║   ██║██║╚██╔╝██║██╔══╝
// ╚██████╗██║  ██║██║ ╚═╝ ██║╚██████╔╝███████║██║  ██║              ╚██████╔╝██║ ╚═╝ ██║███████╗
//  ╚═════╝╚═╝  ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝               ╚═════╝ ╚═╝     ╚═╝╚══════╝
// this code is used to analysis michel electron distribution with different expectation
// ==============================================================================================================
// version: 0.0.1
// date: 2026-01-15
// description: initial version

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/BrentRootFinder.h"
#include "Math/WrappedFunction.h"

class michelElectronAnalysis
{
private:
public:
    // 这部分分析应该是从一个假设开始，首先我们应该已经有了一个足以支撑分析的asimov dataset其中包含了击中本征的Michel电子分布
    // 有两种比较常用的工作模式：
    // 1. 生成假设分布并且讨论在一定统计量下排斥假设的能力（类似于sensitivity分析）
    // 2. 生成假设粉笔并且讨论在一定统计量下测量参数的精度（类似于fit analysis）

    // 生成asimov dataset的函数
    void genPurePiAsimovDataset(std::string eigenFileName, TFile &outputFileName);
    void genIsoAsomovDataset(int binNumber, double cosAlphaMin, double cosAlphaMax, int numEvents, TFile &outputFileName);
    void genCustomAsimovDataset(TFile &outputFile, double kpiRatio, std::string eigenFileName, std::string outputHistName);

    // 分析asimov dataset的函数
    double chi2Calculate(TH1D &hypothesisHist, TH1D &testHist, double targetStatistics = -1.0);
    double excludeHypothesis(TFile *outputFile, std::string hypothesisHistName, std::string testHistName);
};

void michelElectronAnalysis::genPurePiAsimovDataset(std::string eigenFileName, TFile &outputFileName)
{
    TFile *eigenFile = TFile::Open(eigenFileName.c_str(), "READ");
    TH1D *pionAntiMichelDis = (TH1D *)eigenFile->Get("piAntiMuon");
    pionAntiMichelDis->SetDirectory(&outputFileName);
    outputFileName.cd();
    pionAntiMichelDis->Write();
}

void michelElectronAnalysis::genIsoAsomovDataset(int binNumber, double cosAlphaMin, double cosAlphaMax, int numEvents, TFile &outputFileName)
{
    TH1D isoMichelDis("isoMichelDis", "isoMichelDis", binNumber, cosAlphaMin, cosAlphaMax);
    TRandom3 random(0);
    for (int i = 0; i < numEvents; ++i)
    {
        isoMichelDis.Fill(random.Uniform(cosAlphaMin, cosAlphaMax));
    }
    isoMichelDis.SetDirectory(&outputFileName);
    outputFileName.cd();
    isoMichelDis.Write();
}

void michelElectronAnalysis::genCustomAsimovDataset(TFile &outputFile, double kpiRatio, std::string eigenFileName, std::string outputHistName)
{
    TFile *eigenFile = TFile::Open(eigenFileName.c_str(), "READ");
    TH1D *pionAntiMichelDis = (TH1D *)eigenFile->Get("piAntiMuon");
    TH1D *kaonAntiMichelDis = (TH1D *)eigenFile->Get("KAntiMuon");

    TH1D *outputHist = (TH1D *)pionAntiMichelDis->Clone(outputHistName.c_str());
    outputHist->Add(kaonAntiMichelDis, kpiRatio);
    outputHist->SetDirectory(&outputFile);
    outputFile.cd();
    outputHist->Write();
}

double michelElectronAnalysis::chi2Calculate(TH1D &hypothesisHist, TH1D &testHist, double targetStatistics)
{
    // Scale to target statistics
    hypothesisHist.Scale(targetStatistics / hypothesisHist.Integral());
    testHist.Scale(targetStatistics / testHist.Integral());

    // Reset errors to Poisson errors: error = sqrt(N)
    // This is crucial for Chi2Test to work correctly with the scaled statistics
    for (int i = 1; i <= hypothesisHist.GetNbinsX(); ++i)
    {
        double c = hypothesisHist.GetBinContent(i);
        hypothesisHist.SetBinError(i, std::sqrt(c > 0 ? c : 0));
    }
    for (int i = 1; i <= testHist.GetNbinsX(); ++i)
    {
        double c = testHist.GetBinContent(i);
        testHist.SetBinError(i, std::sqrt(c > 0 ? c : 0));
    }

    // Calculate Chi2 and p-value using ROOT functions
    // "WW": Weighted-Weighted comparison. Uses error bars from both histograms.
    // "CHI2": Returns Chi2 value.
    // "P": Prints the result (optional).
    Double_t chi2Val = hypothesisHist.Chi2Test(&testHist, "WW CHI2");
    return chi2Val;
}

double michelElectronAnalysis::excludeHypothesis(TFile *outputFile, std::string hypothesisHistName, std::string testHistName)
{
    // read hypothesis histogram
    TH1D *hypothesisHist = (TH1D *)outputFile->Get(hypothesisHistName.c_str());
    TH1D *testHist = (TH1D *)outputFile->Get(testHistName.c_str());
    const ROOT::Math::WrappedFunction chi2function = [&](double n)
    {
        constexpr double targetChi2 = 9;
        return chi2Calculate(*hypothesisHist, *testHist, n) - targetChi2;
    };

    ROOT::Math::BrentRootFinder rootFinder;
    rootFinder.SetFunction(chi2function, 1e3, 1e9);
    rootFinder.SetLogScan(true);
    rootFinder.Solve();
    double targetStatistics = rootFinder.Root();
    std::cout << "Target statistics to exclude hypothesis at 3 sigma: " << targetStatistics << std::endl;
    return targetStatistics;
}