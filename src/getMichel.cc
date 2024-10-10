// this code is used to generated the michel electron distribution distribution

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <random>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TColor.h"
#include "TTreeFormula.h"
#include "TTreePlayer.h"
#include "TTreeFormulaManager.h"

#include "ROOT/RDataFrame.hxx"

#include "TRandom.h"

using namespace std;
///////////////////////////////////////////////////////////////////////////////////
// get the polarization distribution
///////////////////////////////////////////////////////////////////////////////////
// get the polarization distribution from the .root file

struct muonDataFormat // the data format for the muon data output
{
    float polarization;
    float muonEk;
    float zenith;
    float parentEk;
    int pdgID;
    int parentPDGID;
};

class getPmu
{
private:
public:
    std::map<int, std::map<int, muonDataFormat>> polarizationEnergyRelation(ROOT::RDataFrame &rawDataFrame, int muonID);  // find the muon data from the root file
    void FillTheTTree(TTree *tree, std::map<int, std::map<int, muonDataFormat>> muonDataForm, muonDataFormat &aMuonData); // fill the TTree with muon and parent particle data
    void findRelationPanM(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile);      // pick up all muon data from raw data file

    std::__1::pair<TH1D *, TH1D *> getPAndHeHist(ROOT::RDataFrame &p2MuonDataFrame, ROOT::RDataFrame &he2MuonDataFrame, float energyFloor, float energyCeil, std::string headerOfTH1D, int parentID); // get the proton and helium polarization distribution and mix together
    void TH1DDrawing(TH1D *hist, std::string xTitle, std::string yTitle);                                                                                                                             // TH1D drawing function
    void findCertainEnergyRangePolarization(std::string muonDataFile, std::string polarizationFile);                                                                                                  // generate a certain energy range polarization distribution and put it into a root file
    void getPolarizationFromRaw(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile, std::string polarizationFile);                                              // generate a certain energy range polarization distribution and put it into a root file
};

std::map<int, std::map<int, muonDataFormat>> getPmu::polarizationEnergyRelation(ROOT::RDataFrame &rawDataFrame, int muonID)
{
    // find the mother particle energy and the polarization
    std::vector<std::pair<int, std::vector<int>>> motherMap;
    std::map<int, std::map<int, muonDataFormat>> muonDataForm;

    std::pair<int, std::vector<int>> currentVector;
    currentVector.first = 0;

    rawDataFrame.Foreach([&motherMap, &currentVector, &muonDataForm, &muonID](int PDGID, std::string KillProc, int ParPDGID, float Ek, int ParTrkID, int EvtID, float Zenith, float Helicity)
                         {
                             if (PDGID == muonID && KillProc == "<0|" && Zenith < 1 - std::sqrt(3) / 2)
                             {
                                 if (EvtID != currentVector.first)
                                 {
                                    if (currentVector.second.size() > 0)
                                    {
                                        motherMap.emplace_back(currentVector);
                                    }
                                    currentVector.first = EvtID;
                                    std::vector<int> trkIDVec;
                                    trkIDVec.push_back(ParTrkID);
                                    currentVector.second = trkIDVec;
                                 }
                                 else
                                 {
                                    currentVector.second.emplace_back(ParTrkID);
                                 }
                                 muonDataForm[EvtID][ParTrkID] = {Helicity, Ek, Zenith, 0, PDGID, ParPDGID};
                             } },
                         {"PDGID", "KillProc", "ParPDGID", "Ek", "ParTrkID", "EvtID", "Zenith", "Helicity"});

    auto parentIter = motherMap.begin();
    // find the parent particle energy
    rawDataFrame.Foreach([&muonID, &parentIter, &motherMap, &muonDataForm](int PDGID, int EvtID, int TrkID, float Ek)
                         {
                             if (EvtID != parentIter->first && EvtID == std::next(parentIter)->first && parentIter != motherMap.end())
                             {
                                 parentIter++;
                             }
                             if (EvtID == parentIter->first && std::find(parentIter->second.begin(), parentIter->second.end(), TrkID) != parentIter->second.end())
                             {
                                 muonDataForm[EvtID][TrkID].parentEk = Ek;
                             } },
                         {"PDGID", "EvtID", "TrkID", "Ek"});

    return muonDataForm;
}

void getPmu::FillTheTTree(TTree *tree, std::map<int, std::map<int, muonDataFormat>> muonDataForm, muonDataFormat &aMuonData) // fill the TTree with muon and parent particle data
{
    for (auto &EvtPair : muonDataForm)
    {
        int EvtID = EvtPair.first;
        for (auto &muonPair : EvtPair.second)
        {
            int ParTrkID = muonPair.first;

            const muonDataFormat &muonData = muonPair.second;

            aMuonData.muonEk = muonData.muonEk;
            aMuonData.pdgID = muonData.pdgID;
            aMuonData.parentPDGID = muonData.parentPDGID;
            aMuonData.polarization = muonData.polarization;
            aMuonData.zenith = muonData.zenith;
            aMuonData.parentEk = muonData.parentEk;

            tree->Fill();
        }
    }
}

// pick up all muon data from raw data file
void getPmu::findRelationPanM(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile)
{

    // ROOT::RDataFrame PMuonDataFrame("G4Run0/ReactionChain(13)", "/Volumes/THEMATRIX/CRmuSim/2407muSim/musairs_240717_data/allROOTFile/MusAirS_AMS_polarization.root");
    // ROOT::RDataFrame PAntiMuonDataFrame("G4Run0/ReactionChain(-13)", "/Volumes/THEMATRIX/CRmuSim/2407muSim/musairs_240717_data/allROOTFile/MusAirS_AMS_polarization.root");
    // ROOT::RDataFrame HeMuonDataFrame("G4Run0/ReactionChain(13)", "/Volumes/THEMATRIX/CRmuSim/2407muSim/musairs_240717_data/allROOTFile/musAirSHelium.root");
    // ROOT::RDataFrame HeAntiMuonDataFrame("G4Run0/ReactionChain(-13)", "/Volumes/THEMATRIX/CRmuSim/2407muSim/musairs_240717_data/allROOTFile/musAirSHelium.root");

    ROOT::RDataFrame PMuonDataFrame("G4Run0/ReactionChain(13)", protonProducedData);
    ROOT::RDataFrame PAntiMuonDataFrame("G4Run0/ReactionChain(-13)", protonProducedData);
    ROOT::RDataFrame HeMuonDataFrame("G4Run0/ReactionChain(13)", heliumProducedData);
    ROOT::RDataFrame HeAntiMuonDataFrame("G4Run0/ReactionChain(-13)", heliumProducedData);

    auto PMuonData = polarizationEnergyRelation(PMuonDataFrame, 13);
    auto PAntiMuonData = polarizationEnergyRelation(PAntiMuonDataFrame, -13);
    auto HeMuonData = polarizationEnergyRelation(HeMuonDataFrame, 13);
    auto HeAntiMuonData = polarizationEnergyRelation(HeAntiMuonDataFrame, -13);

    // TFile *outputFile = new TFile("PmuEkRelation2.root", "RECREATE");
    TFile *outputFile = new TFile(muonDataFile.c_str(), "RECREATE");
    TTree *treeHe = new TTree("muonData(He)", "muonData(He)");
    TTree *treeP = new TTree("muonData(P)", "muonData(P)");

    muonDataFormat aMuonData;
    treeHe->Branch("Ek", &aMuonData.muonEk, "Ek/F");
    treeHe->Branch("PDGID", &aMuonData.pdgID, "PDGID/I");
    treeHe->Branch("ParPDGID", &aMuonData.parentPDGID, "ParPDGID/I");
    treeHe->Branch("Pmu", &aMuonData.polarization, "Pmu/F");
    treeHe->Branch("zenith", &aMuonData.zenith, "zenith/F");
    treeHe->Branch("ParEk", &aMuonData.parentEk, "ParEk/F");

    treeP->Branch("Ek", &aMuonData.muonEk, "Ek/F");
    treeP->Branch("PDGID", &aMuonData.pdgID, "PDGID/I");
    treeP->Branch("ParPDGID", &aMuonData.parentPDGID, "ParPDGID/I");
    treeP->Branch("Pmu", &aMuonData.polarization, "Pmu/F");
    treeP->Branch("zenith", &aMuonData.zenith, "zenith/F");
    treeP->Branch("ParEk", &aMuonData.parentEk, "ParEk/F");

    FillTheTTree(treeHe, HeMuonData, aMuonData);
    FillTheTTree(treeHe, HeAntiMuonData, aMuonData);
    treeHe->Write();

    FillTheTTree(treeP, PMuonData, aMuonData);
    FillTheTTree(treeP, PAntiMuonData, aMuonData);
    treeP->Write();

    outputFile->Close();
}

std::__1::pair<TH1D *, TH1D *> getPmu::getPAndHeHist(ROOT::RDataFrame &p2MuonDataFrame, ROOT::RDataFrame &he2MuonDataFrame, float energyFloor, float energyCeil, std::string headerOfTH1D, int parentID) // get the proton and helium polarization distribution and mix together
{
    if (gROOT->FindObject("MuonPmu") && gROOT->FindObject("AntiPmu") && gROOT->FindObject("MuonPmuHe") && gROOT->FindObject("AntiPmuHe"))
    {
        delete gROOT->FindObject("MuonPmu");
        delete gROOT->FindObject("AntiPmu");
        delete gROOT->FindObject("MuonPmuHe");
        delete gROOT->FindObject("AntiPmuHe");
    }

    TH1D *MuonPmuP = new TH1D("MuonPmu", "MuonPmu", 100, -1, 1);
    TH1D *AntiMuonPmuP = new TH1D("AntiPmu", "AntiPmu", 100, -1, 1);
    TH1D *MuonPmuHe = new TH1D("MuonPmuHe", "MuonPmuHe", 100, -1, 1);
    TH1D *AntiMuonPmuHe = new TH1D("AntiPmuHe", "AntiPmuHe", 100, -1, 1);

    p2MuonDataFrame.Foreach([&parentID, &energyFloor, &energyCeil, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu)
                            {
                                if (std::abs(ParPDGID) == parentID && Ek > energyFloor && Ek < energyCeil)
                                {
                                    if (PDGID == 13)
                                    {
                                        MuonPmuP->Fill(Pmu);
                                    }
                                    if (PDGID == -13)
                                    {
                                        AntiMuonPmuP->Fill(Pmu);
                                    }
                                } }, {"PDGID", "ParPDGID", "Ek", "ParEk", "Pmu"});

    he2MuonDataFrame.Foreach([&parentID, &energyFloor, &energyCeil, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu)
                             {
                                 if (std::abs(ParPDGID) == parentID && Ek > energyFloor && Ek < energyCeil)
                                 {
                                     if (PDGID == 13)
                                     {
                                         MuonPmuHe->Fill(Pmu);
                                     }
                                     if (PDGID == -13)
                                     {
                                         AntiMuonPmuHe->Fill(Pmu);
                                     }
                                 } },
                             {"PDGID", "ParPDGID", "Ek", "ParEk", "Pmu"});

    // combine two histograms into a new total histogram
    auto muonHistName = headerOfTH1D + "MuonPmu";
    auto antiMuonHistName = headerOfTH1D + "AntiMuonPmu";

    TH1D *totalMuonPmu = new TH1D(muonHistName.c_str(), muonHistName.c_str(), 100, -1, 1);
    TH1D *totalAntiMuonPmu = new TH1D(antiMuonHistName.c_str(), antiMuonHistName.c_str(), 100, -1, 1);

    MuonPmuP->Scale(9.704844263299673);
    AntiMuonPmuP->Scale(9.704844263299673);
    MuonPmuHe->Scale(100.0);
    AntiMuonPmuHe->Scale(100.0);

    totalMuonPmu->Add(MuonPmuP, MuonPmuHe);
    totalAntiMuonPmu->Add(AntiMuonPmuP, AntiMuonPmuHe);

    totalMuonPmu->Scale(1.0 / totalMuonPmu->Integral());
    totalAntiMuonPmu->Scale(1.0 / totalAntiMuonPmu->Integral());

    std::pair<TH1D *, TH1D *> histPair = {totalMuonPmu, totalAntiMuonPmu};
    return histPair;
}

void getPmu::TH1DDrawing(TH1D *hist, std::string xTitle, std::string yTitle) // TH1D drawing function
{
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->SetMarkerStyle(21);
    hist->SetMarkerSize(0.7);
}

// generate a certain energy range polarization distribution and put it into a root file
void getPmu::findCertainEnergyRangePolarization(std::string muonDataFile, std::string polarizationFile)
{
    // find the polarization of the pion

    // define the energy range for muon
    float energyFloor = 100;
    float energyCeil = 500;

    // find the related parent particles energy spectrum
    // ROOT::RDataFrame p2MuonDataFrame("muonData(P)", "PmuEkRelation2.root");
    // ROOT::RDataFrame he2MuonDataFrame("muonData(He)", "PmuEkRelation2.root");
    ROOT::RDataFrame p2MuonDataFrame("muonData(P)", muonDataFile.c_str());
    ROOT::RDataFrame he2MuonDataFrame("muonData(He)", muonDataFile.c_str());

    auto piHists = getPAndHeHist(p2MuonDataFrame, he2MuonDataFrame, energyFloor, energyCeil, "pi", 211);
    auto KHists = getPAndHeHist(p2MuonDataFrame, he2MuonDataFrame, energyFloor, energyCeil, "K", 321);
    auto KLongHists = getPAndHeHist(p2MuonDataFrame, he2MuonDataFrame, 0, 1e9, "KLong", 130);

    TH1DDrawing(piHists.first, "polarization (P_{#mu^{-}} (#pi^{-}))", "Normalized Frequency");
    TH1DDrawing(piHists.second, "polarization (P_{#mu^{+}} (#pi^{+}))", "Normalized Frequency");

    TH1DDrawing(KHists.first, "polarization (P_{#mu^{-}} (K^{-}))", "Normalized Frequency");
    TH1DDrawing(KHists.second, "polarization (P_{#mu^{+}} (K^{+}))", "Normalized Frequency");

    TH1DDrawing(KLongHists.first, "polarization (P_{#mu^{-}} (K_{L}^{0}))", "Normalized Frequency");
    TH1DDrawing(KLongHists.second, "polarization (P_{#mu^{+}} (K_{L}^{0}))", "Normalized Frequency");

    TFile *outputFile = new TFile(polarizationFile.c_str(), "RECREATE");
    piHists.first->Write();
    piHists.second->Write();

    KHists.first->Write();
    KHists.second->Write();

    KLongHists.first->Write();
    KLongHists.second->Write();

    outputFile->Close();
}

// total polarization calculation
void getPmu::getPolarizationFromRaw(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile, std::string polarizationFile)
{
    findRelationPanM(protonProducedData, heliumProducedData, muonDataFile);
    findCertainEnergyRangePolarization(muonDataFile, polarizationFile);
}

///////////////////////////////////////////////////////////////////////////////////
// MC for the Michel electron distribution
///////////////////////////////////////////////////////////////////////////////////

class GenMichel
{
private:
public:
    double getWeight(double alpha, double Pmu, double cosTheta, double phi);
    std::tuple<double, double, double> generator(TRandom &random, double Pmu);
    TH1D getSimulatedMichel(std::string inputFile);
    TH1D generatePolarization(std::string polarizationFile, double KPiRatio, bool charge);

    TH1D PmuMichelDis(double KPiRatio, bool charge, int numEvents, std::string polarizationFile); // generate the polarization distribution of the Michel electrons in a certain Pmu
    TH1D IsotropyDis(int numEvents);                                                              // generate the isotropy distribution of the Michel electrons
    TH2D PmuMichelDis2D(double KPiRatio, bool charge, int numEvents, std::string polarizationFile);
};

double GenMichel::getWeight(double alpha, double Pmu, double cosTheta, double phi)
{
    auto sinTheta = sqrt(1 - pow(cosTheta, 2));
    auto cosBeta = Pmu;
    auto sinBeta = sqrt(1 - pow(cosBeta, 2)); // beta < pi

    // auto son1 = pow(1 / cos(alpha) * sinBeta * sinTheta * pow(abs(cos(alpha)), 2) * (sin(phi) - tan(alpha) * cos(phi)) + cosTheta * cosTheta, 2);
    // // auto son1 = pow((sin(phi) - cos(phi) * tan(alpha)) * sinBeta * sinTheta * cos(alpha) + cosTheta * cosTheta, 2);

    // auto mum1 = pow(abs(cos(alpha)), 4) * pow((cosTheta * sinTheta * (sin(phi) - tan(alpha) * cos(phi)) - 1 / cos(alpha) * sinBeta * cosTheta + sin(alpha) * sinTheta * cos(alpha - phi) * pow(abs(1 / cos(alpha)), 2)), 2) + pow((cos(alpha) * sinTheta * cos(alpha - phi) + tan(alpha) * 1 / cos(alpha) * pow(abs(cos(alpha)), 2) * (cosTheta * sinTheta * sin(alpha - phi) + sinBeta * cosTheta)), 2) + pow((1 / cos(alpha) * sinBeta * sinTheta * pow(abs(cos(alpha)), 2) * (sin(phi) - tan(alpha) * cos(phi)) + cosTheta * cosTheta), 2);
    // // auto mum1 = pow((sin(phi) - cos(phi) * tan(alpha)) * sinBeta * sinTheta * cos(alpha) + cosTheta * cosTheta, 2) + pow((sin(phi) - cos(phi) * tan(alpha)) * sinTheta * pow(cos(alpha), 2) * cosTheta + sin(alpha) * sinTheta * cos(alpha - phi) - sinBeta * cos(alpha) * cosTheta, 2);

    // auto son2 = 1 / cos(alpha) * sinBeta * sinTheta * pow(abs(cos(alpha)), 2) * (sin(phi) - tan(alpha) * cos(phi)) + cosTheta * cosTheta;
    // // auto son2 = pow((sin(phi) - cos(phi) * tan(alpha)) * sinBeta * sinTheta * cos(alpha) + cosTheta * cosTheta, 2);

    // auto mum2 = 3 * sqrt(pow(abs(cos(alpha)), 4) * pow((cosTheta * sinTheta * (sin(phi) - tan(alpha) * cos(phi)) - 1 / cos(alpha) * sinBeta * cosTheta + sin(alpha) * sinTheta * cos(alpha - phi) * pow(abs(1 / cos(alpha)), 2)), 2) + pow((cos(alpha) * sinTheta * cos(alpha - phi) + tan(alpha) / cos(alpha) * pow(abs(cos(alpha)), 2) * (cosTheta * sinTheta * sin(alpha - phi) + sinBeta * cosTheta)), 2) + pow((1 / cos(alpha) * sinBeta * sinTheta * pow(abs(cos(alpha)), 2) * (sin(phi) - tan(alpha) * cos(phi)) + cosTheta * cosTheta), 2));
    // // auto mum2 = 3 * sqrt(pow((sinBeta * cosTheta + sinTheta * sin(alpha - phi) * cosTheta) * cosTheta * tan(alpha) + sinTheta * cos(alpha) * cos(alpha - phi), 2) + pow((sin(phi) - cos(phi) * tan(alpha)) * sinBeta * sinTheta * cos(alpha) + cosTheta * cosTheta, 2) + pow((sin(phi) - cos(phi) * tan(alpha)) * sinTheta * pow(cos(alpha), 2) * cosTheta + sin(alpha) * sinTheta * cos(alpha - phi) - sinBeta * cos(alpha) * cosTheta, 2));

    // auto weight = 2 * sqrt((1 - son1 / mum1) * pow(1 + son2 / mum2, 2));

    auto polar = 1 + (cosBeta * cosTheta + pow(abs(cos(alpha)), 2) * 1 / cos(alpha) * sinBeta * sinTheta * (sin(phi) - cos(phi) * tan(alpha))) /
                         (3. * sqrt(pow(cos(alpha) * cos(alpha - phi) * sinTheta + pow(abs(cos(alpha)), 2) * 1 / cos(alpha) * (cosTheta * sinBeta + cosBeta * sinTheta * sin(alpha - phi)) * tan(alpha), 2) +
                                    pow(abs(cos(alpha)), 4) * pow(-(cosTheta * 1 / cos(alpha) * sinBeta) + pow(abs(1 / cos(alpha)), 2) * cos(alpha - phi) * sin(alpha) * sinTheta + cosBeta * sinTheta * (sin(phi) - cos(phi) * tan(alpha)), 2) + pow(cosBeta * cosTheta + pow(abs(cos(alpha)), 2) * 1 / cos(alpha) * sinBeta * sinTheta * (sin(phi) - cos(phi) * tan(alpha)), 2)));
    return polar;
}

std::tuple<double, double, double> GenMichel::generator(TRandom &random, double Pmu)
{

    double alpha = random.Uniform(0, 2 * TMath::Pi());
    double cosTheta = random.Uniform(-1, 1);
    double phi = random.Uniform(0, 2 * TMath::Pi());

    auto weight = getWeight(alpha, Pmu, cosTheta, phi);

    return std::make_tuple(cosTheta, phi, weight);
}

TH1D GenMichel::getSimulatedMichel(std::string inputFile)
{
    // ROOT::RDataFrame simulatedMichelData("Michel Momentum Direction", "/Users/sunmingchen/RawData/simulationMichel.root");
    ROOT::RDataFrame simulatedMichelData("Michel Momentum Direction", inputFile.c_str());

    TH1D simulatedMichel("simulatedMichel", "simulatedMichel", 100, -1, 1);

    simulatedMichelData.Foreach([&simulatedMichel](double momentumX, double momentumY, double momentumZ)
                                { 
                                    double cosTheta = momentumZ / sqrt(momentumX * momentumX + momentumY * momentumY + momentumZ * momentumZ);
                                    simulatedMichel.Fill(- cosTheta); }, // cosTheta 前面的符号是对应缪子的电性
                                {"momentumX", "momentumY", "momentumZ"});
    return simulatedMichel;
}

TH1D GenMichel::generatePolarization(std::string polarizationFile, double KPiRatio, bool charge)
{
    // TFile polarizationData("100to500MuonPolarization.root", "READ");
    TFile polarizationData(polarizationFile.c_str(), "READ");
    if (charge)
    {
        auto piPolarization = (TH1D *)polarizationData.Get("piAntiMuonPmu");
        auto kaPolarization = (TH1D *)polarizationData.Get("KAntiMuonPmu");

        std::string headerOfTH1D = "KPiRatio = " + std::to_string(KPiRatio) + " Anti-muon";

        TH1D outputHist(headerOfTH1D.c_str(), headerOfTH1D.c_str(), 100, -1, 1);
        piPolarization->Scale(1 / piPolarization->Integral());
        kaPolarization->Scale(KPiRatio / kaPolarization->Integral());

        outputHist.Add(piPolarization, kaPolarization);
        outputHist.Scale(1 / outputHist.Integral());
        return outputHist;
    }
    else
    {
        auto piPolarization = (TH1D *)polarizationData.Get("piMuonPmu");
        auto kaPolarization = (TH1D *)polarizationData.Get("KMuonPmu");

        std::string headerOfTH1D = "KPiRatio = " + std::to_string(KPiRatio) + " Muon";

        TH1D outputHist(headerOfTH1D.c_str(), headerOfTH1D.c_str(), 100, -1, 1);
        piPolarization->Scale(1 / piPolarization->Integral());
        kaPolarization->Scale(KPiRatio / kaPolarization->Integral());

        outputHist.Add(piPolarization, kaPolarization);
        outputHist.Scale(1 / outputHist.Integral());
        return outputHist;
    }
}

TH1D GenMichel::PmuMichelDis(double KPiRatio, bool charge, int numEvents, std::string polarizationFile)
{
    std::string HistID = "KPiRatio = " + std::to_string(KPiRatio) + " ";
    if (charge)
    {
        HistID += "Anti-muon";
    }
    else
    {
        HistID += "Muon";
    }

    TH1D pmuMichelDis(HistID.c_str(), HistID.c_str(), 100, -1, 1);

    TRandom random(time(0));
    auto polarizationHist = generatePolarization(polarizationFile, KPiRatio, charge);
    for (int i = 0; i < numEvents; i++)
    {
        double Pmu = polarizationHist.GetRandom();

        auto coordinateTuple = generator(random, Pmu);
        double cosTheta = std::get<0>(coordinateTuple);
        double phi = std::get<1>(coordinateTuple);
        double weight = std::get<2>(coordinateTuple);

        pmuMichelDis.Fill(cosTheta, weight);
    }

    return pmuMichelDis;
}

TH2D GenMichel::PmuMichelDis2D(double KPiRatio, bool charge, int numEvents, std::string polarizationFile)
{
    std::string HistID = "KPiRatio = " + std::to_string(KPiRatio) + " ";
    if (charge)
    {
        HistID += "Anti-muon";
    }
    else
    {
        HistID += "Muon";
    }

    TH2D pmuMichelDis(HistID.c_str(), HistID.c_str(), 100, 0, 2 * TMath::Pi(), 100, -1, 1);

    TRandom random(time(0));
    auto polarizationHist = generatePolarization(polarizationFile, KPiRatio, charge);
    for (int i = 0; i < numEvents; i++)
    {
        double Pmu = polarizationHist.GetRandom();

        auto coordinateTuple = generator(random, Pmu);
        double cosTheta = std::get<0>(coordinateTuple);
        double phi = std::get<1>(coordinateTuple);
        double weight = std::get<2>(coordinateTuple);

        pmuMichelDis.Fill(phi, cosTheta, weight);
    }

    return pmuMichelDis;
}

TH1D GenMichel::IsotropyDis(int numEvents)
{
    TH1D isotropyDis("IsotropyDis", "IsotropyDis", 100, -1, 1);
    TRandom random(time(0));

    for (int i = 0; i < numEvents; i++)
    {
        double cosTheta = random.Uniform(-1, 1);
        double weight = 1;

        isotropyDis.Fill(cosTheta, weight);
    }

    return isotropyDis;
}

///////////////////////////////////////////////////////////////////////////////////
// main function: example for generating the polarization distribution of the Michel electrons in a certain Pmu
///////////////////////////////////////////////////////////////////////////////////

// int main()
// {
//     GenMichel genMichel;
//     TH1D michelDistribution = genMichel.PmuMichelDis(0.1, true, 1e7, "100to500MuonPolarization.root");
//     TH1D purePionDistribution = genMichel.PmuMichelDis(0, true, 1e7, "100to500MuonPolarization.root");
//     // TH1D iso = genMichel.IsotropyDis(1e8);

//     TFile output("output.root", "RECREATE");
//     TCanvas c1("c1", "c1", 800, 600);

//     michelDistribution.SetStats(0);
//     purePionDistribution.SetStats(0);

//     michelDistribution.SetLineColor(kRed);
//     purePionDistribution.SetLineColor(kBlue);

//     michelDistribution.SetLineWidth(2);
//     purePionDistribution.SetLineWidth(2);

//     michelDistribution.Draw("hist");
//     purePionDistribution.Draw("hist same");

//     c1.Write();
//     output.Close();

//     return 0;
// }