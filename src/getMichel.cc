//  ██████╗██████╗ ███╗   ███╗██╗   ██╗███████╗██████╗                ██████╗ ███╗   ███╗███████╗
// ██╔════╝██╔══██╗████╗ ████║██║   ██║██╔════╝██╔══██╗              ██╔════╝ ████╗ ████║██╔════╝
// ██║     ██████╔╝██╔████╔██║██║   ██║███████╗██████╔╝    █████╗    ██║  ███╗██╔████╔██║█████╗
// ██║     ██╔══██╗██║╚██╔╝██║██║   ██║╚════██║██╔══██╗    ╚════╝    ██║   ██║██║╚██╔╝██║██╔══╝
// ╚██████╗██║  ██║██║ ╚═╝ ██║╚██████╔╝███████║██║  ██║              ╚██████╔╝██║ ╚═╝ ██║███████╗
//  ╚═════╝╚═╝  ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝               ╚═════╝ ╚═╝     ╚═╝╚══════╝
// this code is used to generated the michel electron distribution distribution
// ==============================================================================================================
// version: 0.0.2
// date: 2025-08-27
// description: update to .root file with weight

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <limits>

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

struct muonDataFormat // the data format for the muon data output
{
    float polarization;
    float muonEk;
    float zenith;
    float parentEk;
    int pdgID;
    int parentPDGID;
};

struct muonDataWithWeight
{
    float polarization;
    float muonEk;
    float zenith;
    float parentEk;
    int pdgID;
    int parentPDGID;
    float weight;
};

///////////////////////////////////////////////////////////////////////////////////
// get the polarization distribution
///////////////////////////////////////////////////////////////////////////////////
// get the polarization distribution from the .root file

class getPmu
{
private:
public:
    std::map<int, std::map<int, muonDataFormat>> polarizationEnergyRelation(ROOT::RDataFrame &rawDataFrame, int muonID);  // find the muon data from the root file
    void FillTheTTree(TTree *tree, std::map<int, std::map<int, muonDataFormat>> muonDataForm, muonDataFormat &aMuonData); // fill the TTree with muon and parent particle data
    void findRelationPanM(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile);      // pick up all muon data from raw data file
    void processMuonDataByEvent(ROOT::RDataFrame &rawDataFrame, int muonID, TTree *tree, muonDataFormat &aMuonData);      // stream one event at a time and write directly to tree

    std::__1::pair<TH1D *, TH1D *> getPAndHeHist(ROOT::RDataFrame &pr_mu_df, ROOT::RDataFrame &he_mu_df, float E_lower_bound, float E_upper_bound, std::string hist_header, int parentID, int binNumber); // get the proton and helium polarization distribution and mix together
    void TH1DDrawing(TH1D *hist, std::string xTitle, std::string yTitle);                                                                                                                                 // TH1D drawing function
    void findCertainEnergyRangePolarization(std::string muonDataFile, std::string polarizationFile, int binNumber);                                                                                       // generate a certain energy range polarization distribution and put it into a root file
    void getPolarizationFromRaw(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile, std::string polarizationFile, int binNumber);                                   // generate a certain energy range polarization distribution and put it into a root file

    // get the polarization with weight
    std::map<int, std::map<int, muonDataWithWeight>> getPolarizationWithWeight(ROOT::RDataFrame &rawDataFrame, int muonID);                                                                                     // find the muon data from the root file (with weight)
    void FillTheTTreeWithWeight(TTree *tree, std::map<int, std::map<int, muonDataWithWeight>> muonDataForm, muonDataWithWeight &aMuonData);                                                                     // fill the TTree with muon and parent particle data (with weight)
    void findRelationPolarizationWithWeight(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile);                                                                          // pick up all muon data from raw data file (with weight)
    void findRelationPolarizationWithWeightByEvent(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile);                                                                   // pick up all muon data from raw data file (with weight, event-streaming)
    void processMuonDataByEventWithWeight(ROOT::RDataFrame &rawDataFrame, int muonID, TTree *tree, muonDataWithWeight &aMuonData);                                                                              // stream one event at a time and write directly to tree (with weight)
    std::__1::pair<TH1D *, TH1D *> pr_he_weighted_hist(ROOT::RDataFrame &pr_mu_df, ROOT::RDataFrame &he_mu_df, float E_lower_bound, float E_upper_bound, std::string hist_header, int parentID, int binNumber); // get the proton and helium polarization distribution and mix together
    void findCertainEnergyRangePolarizationWithWeight(std::string muonDataFile, std::string polarizationFile, int binNumber);                                                                                   // generate a certain energy range polarization distribution and put it into a root file
    void gen_pol_with_wgt(std::string proton_produced_data, std::string helium_produced_data, std::string muon_data_file, std::string polarization_file, int bin_number);                                       // generate a certain energy range polarization distribution and put it into a root file
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

    processMuonDataByEvent(HeMuonDataFrame, 13, treeHe, aMuonData);
    processMuonDataByEvent(HeAntiMuonDataFrame, -13, treeHe, aMuonData);
    treeHe->Write();

    processMuonDataByEvent(PMuonDataFrame, 13, treeP, aMuonData);
    processMuonDataByEvent(PAntiMuonDataFrame, -13, treeP, aMuonData);
    treeP->Write();

    outputFile->Close();
}

void getPmu::processMuonDataByEvent(ROOT::RDataFrame &rawDataFrame, int muonID, TTree *tree, muonDataFormat &aMuonData)
{
    struct MuonCandidate
    {
        float polarization;
        float muonEk;
        float zenith;
        int pdgID;
        int parentPDGID;
        int parentTrkID;
    };

    std::vector<MuonCandidate> currentEventMuons;
    std::unordered_map<int, float> currentEventTrackEnergy;
    int currentEvtID = std::numeric_limits<int>::min();
    bool hasEvent = false;
    auto flushCurrentEvent = [&]()
    {
        for (const auto &muon : currentEventMuons)
        {
            aMuonData.polarization = muon.polarization;
            aMuonData.muonEk = muon.muonEk;
            aMuonData.zenith = muon.zenith;
            aMuonData.pdgID = muon.pdgID;
            aMuonData.parentPDGID = muon.parentPDGID;

            auto parentIter = currentEventTrackEnergy.find(muon.parentTrkID);
            aMuonData.parentEk = parentIter == currentEventTrackEnergy.end() ? 0.0f : parentIter->second;
            tree->Fill();
        }

        currentEventMuons.clear();
        currentEventTrackEnergy.clear();
    };

    rawDataFrame.Foreach([&](int PDGID, std::string KillProc, int ParPDGID, float Ek, int ParTrkID, int TrkID, int EvtID, float Zenith, float Helicity)
                         {
                             if (!hasEvent)
                             {
                                 currentEvtID = EvtID;
                                 hasEvent = true;
                             }

                             if (EvtID != currentEvtID)
                             {
                                 flushCurrentEvent();
                                 currentEvtID = EvtID;
                             }

                             currentEventTrackEnergy[TrkID] = Ek;

                             if (PDGID == muonID && KillProc == "<0|" && Zenith < 1 - std::sqrt(3) / 2)
                             {
                                 currentEventMuons.emplace_back(MuonCandidate{Helicity, Ek, Zenith, PDGID, ParPDGID, ParTrkID});
                             } },
                         {"PDGID", "KillProc", "ParPDGID", "Ek", "ParTrkID", "TrkID", "EvtID", "Zenith", "Helicity"});

    if (hasEvent)
    {
        flushCurrentEvent();
    }
}

std::__1::pair<TH1D *, TH1D *> getPmu::getPAndHeHist(ROOT::RDataFrame &pr_mu_df, ROOT::RDataFrame &he_mu_df, float E_lower_bound, float E_upper_bound, std::string hist_header, int parentID, int binNumber) // get the proton and helium polarization distribution and mix together
{
    if (gROOT->FindObject("MuonPmu") && gROOT->FindObject("AntiPmu") && gROOT->FindObject("MuonPmuHe") && gROOT->FindObject("AntiPmuHe"))
    {
        delete gROOT->FindObject("MuonPmu");
        delete gROOT->FindObject("AntiPmu");
        delete gROOT->FindObject("MuonPmuHe");
        delete gROOT->FindObject("AntiPmuHe");
    }

    TH1D *MuonPmuP = new TH1D("MuonPmu", "MuonPmu", binNumber, -1, 1);
    TH1D *AntiMuonPmuP = new TH1D("AntiPmu", "AntiPmu", binNumber, -1, 1);
    TH1D *MuonPmuHe = new TH1D("MuonPmuHe", "MuonPmuHe", binNumber, -1, 1);
    TH1D *AntiMuonPmuHe = new TH1D("AntiPmuHe", "AntiPmuHe", binNumber, -1, 1);

    pr_mu_df.Foreach([&parentID, &E_lower_bound, &E_upper_bound, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu)
                     {
                                if (std::abs(ParPDGID) == parentID && Ek > E_lower_bound && Ek < E_upper_bound)
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

    he_mu_df.Foreach([&parentID, &E_lower_bound, &E_upper_bound, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu)
                     {
                                 if (std::abs(ParPDGID) == parentID && Ek > E_lower_bound && Ek < E_upper_bound)
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
    auto muonHistName = hist_header + "MuonPmu";
    auto antiMuonHistName = hist_header + "AntiMuonPmu";

    TH1D *totalMuonPmu = new TH1D(muonHistName.c_str(), muonHistName.c_str(), binNumber, -1, 1);
    TH1D *totalAntiMuonPmu = new TH1D(antiMuonHistName.c_str(), antiMuonHistName.c_str(), binNumber, -1, 1);

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
void getPmu::findCertainEnergyRangePolarization(std::string muonDataFile, std::string polarizationFile, int binNumber)
{
    // find the polarization of the pion

    // define the energy range for muon
    float E_lower_bound = 100;
    float E_upper_bound = 500;

    // find the related parent particles energy spectrum
    // ROOT::RDataFrame pr_mu_df("muonData(P)", "PmuEkRelation2.root");
    // ROOT::RDataFrame he_mu_df("muonData(He)", "PmuEkRelation2.root");
    ROOT::RDataFrame pr_mu_df("muonData(P)", muonDataFile.c_str());
    ROOT::RDataFrame he_mu_df("muonData(He)", muonDataFile.c_str());

    auto piHists = getPAndHeHist(pr_mu_df, he_mu_df, E_lower_bound, E_upper_bound, "pi", 211, binNumber);
    auto KHists = getPAndHeHist(pr_mu_df, he_mu_df, E_lower_bound, E_upper_bound, "K", 321, binNumber);
    auto KLongHists = getPAndHeHist(pr_mu_df, he_mu_df, 0, 1e9, "KLong", 130, binNumber);

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
void getPmu::getPolarizationFromRaw(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile, std::string polarizationFile, int binNumber)
{
    findRelationPanM(protonProducedData, heliumProducedData, muonDataFile);
    findCertainEnergyRangePolarization(muonDataFile, polarizationFile, binNumber);
}

std::map<int, std::map<int, muonDataWithWeight>> getPmu::getPolarizationWithWeight(ROOT::RDataFrame &rawDataFrame, int muonID)
{
    // find the mother particle energy and the polarization
    std::vector<std::pair<int, std::vector<int>>> motherMap;
    std::map<int, std::map<int, muonDataWithWeight>> muonDataForm;

    std::pair<int, std::vector<int>> currentVector;
    currentVector.first = 0;

    rawDataFrame.Foreach([&motherMap, &currentVector, &muonDataForm, &muonID](int PDGID, std::string KillProc, int ParPDGID, float Ek, int ParTrkID, int EvtID, float Zenith, float Helicity, float Weight)
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
                                 muonDataForm[EvtID][ParTrkID] = {Helicity, Ek, Zenith, 0, PDGID, ParPDGID, Weight};
                             }
                             //  mission accomplished
                         },
                         {"PDGID", "KillProc", "ParPDGID", "Ek", "ParTrkID", "EvtID", "Zenith", "Helicity", "Weight"});

    auto parentIter = motherMap.begin();
    // find the parent particle energy
    rawDataFrame.Foreach([&muonID, &parentIter, &motherMap, &muonDataForm](int PDGID, int EvtID, int TrkID, float Ek)
                         {
                             if (EvtID == parentIter->first && std::find(parentIter->second.begin(), parentIter->second.end(), TrkID) != parentIter->second.end())
                             {
                                 muonDataForm[EvtID][TrkID].parentEk = Ek;
                             }
                             if (EvtID == parentIter->first && std::find(parentIter->second.begin(), parentIter->second.end(), TrkID) != parentIter->second.end())
                             {
                                 muonDataForm[EvtID][TrkID].parentEk = Ek;
                             }

                             // mission accomplished
                         },
                         {"PDGID", "EvtID", "TrkID", "Ek"});

    return muonDataForm;
}

void getPmu::FillTheTTreeWithWeight(TTree *tree, std::map<int, std::map<int, muonDataWithWeight>> muonDataForm, muonDataWithWeight &aMuonData) // fill the TTree with muon and parent particle data (with weight)
{
    for (auto &EvtPair : muonDataForm)
    {
        int EvtID = EvtPair.first;
        for (auto &muonPair : EvtPair.second)
        {
            int ParTrkID = muonPair.first;

            const muonDataWithWeight &muonData = muonPair.second;

            aMuonData.muonEk = muonData.muonEk;
            aMuonData.pdgID = muonData.pdgID;
            aMuonData.parentPDGID = muonData.parentPDGID;
            aMuonData.polarization = muonData.polarization;
            aMuonData.zenith = muonData.zenith;
            aMuonData.parentEk = muonData.parentEk;
            aMuonData.weight = muonData.weight;

            tree->Fill();
        }
    }
}

void getPmu::findRelationPolarizationWithWeight(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile)
{
    ROOT::RDataFrame pr_mu_minus_df("G4Run0/ReactionChain(13)", protonProducedData);
    ROOT::RDataFrame pr_mu_plus_df("G4Run0/ReactionChain(-13)", protonProducedData);
    ROOT::RDataFrame he_mu_minus_df("G4Run0/ReactionChain(13)", heliumProducedData);
    ROOT::RDataFrame he_mu_plus_df("G4Run0/ReactionChain(-13)", heliumProducedData);

    auto pr_mu_minus = getPolarizationWithWeight(pr_mu_minus_df, 13);
    auto pr_mu_plus = getPolarizationWithWeight(pr_mu_plus_df, -13);
    auto he_mu_minus = getPolarizationWithWeight(he_mu_minus_df, 13);
    auto he_mu_plus = getPolarizationWithWeight(he_mu_plus_df, -13);

    TFile *output_file = new TFile(muonDataFile.c_str(), "RECREATE");
    TTree *tree_he = new TTree("muonData(He)", "muonData(He)");
    TTree *tree_pr = new TTree("muonData(P)", "muonData(P)");

    muonDataWithWeight a_muon_data;
    tree_he->Branch("Ek", &a_muon_data.muonEk, "Ek/F");
    tree_he->Branch("PDGID", &a_muon_data.pdgID, "PDGID/I");
    tree_he->Branch("ParPDGID", &a_muon_data.parentPDGID, "ParPDGID/I");
    tree_he->Branch("Pmu", &a_muon_data.polarization, "Pmu/F");
    tree_he->Branch("Zenith", &a_muon_data.zenith, "Zenith/F");
    tree_he->Branch("ParEk", &a_muon_data.parentEk, "ParEk/F");
    tree_he->Branch("Weight", &a_muon_data.weight, "Weight/F");

    tree_pr->Branch("Ek", &a_muon_data.muonEk, "Ek/F");
    tree_pr->Branch("PDGID", &a_muon_data.pdgID, "PDGID/I");
    tree_pr->Branch("ParPDGID", &a_muon_data.parentPDGID, "ParPDGID/I");
    tree_pr->Branch("Pmu", &a_muon_data.polarization, "Pmu/F");
    tree_pr->Branch("Zenith", &a_muon_data.zenith, "Zenith/F");
    tree_pr->Branch("ParEk", &a_muon_data.parentEk, "ParEk/F");
    tree_pr->Branch("Weight", &a_muon_data.weight, "Weight/F");

    FillTheTTreeWithWeight(tree_he, he_mu_minus, a_muon_data);
    FillTheTTreeWithWeight(tree_he, he_mu_plus, a_muon_data);
    tree_he->Write();

    FillTheTTreeWithWeight(tree_pr, pr_mu_minus, a_muon_data);
    FillTheTTreeWithWeight(tree_pr, pr_mu_plus, a_muon_data);
    tree_pr->Write();

    output_file->Close();
}

void getPmu::findRelationPolarizationWithWeightByEvent(std::string protonProducedData, std::string heliumProducedData, std::string muonDataFile)
{
    ROOT::RDataFrame pr_mu_minus_df("G4Run0/ReactionChain(13)", protonProducedData);
    ROOT::RDataFrame pr_mu_plus_df("G4Run0/ReactionChain(-13)", protonProducedData);
    ROOT::RDataFrame he_mu_minus_df("G4Run0/ReactionChain(13)", heliumProducedData);
    ROOT::RDataFrame he_mu_plus_df("G4Run0/ReactionChain(-13)", heliumProducedData);

    TFile *output_file = new TFile(muonDataFile.c_str(), "RECREATE");
    TTree *tree_he = new TTree("muonData(He)", "muonData(He)");
    TTree *tree_pr = new TTree("muonData(P)", "muonData(P)");

    muonDataWithWeight a_muon_data;
    tree_he->Branch("Ek", &a_muon_data.muonEk, "Ek/F");
    tree_he->Branch("PDGID", &a_muon_data.pdgID, "PDGID/I");
    tree_he->Branch("ParPDGID", &a_muon_data.parentPDGID, "ParPDGID/I");
    tree_he->Branch("Pmu", &a_muon_data.polarization, "Pmu/F");
    tree_he->Branch("Zenith", &a_muon_data.zenith, "Zenith/F");
    tree_he->Branch("ParEk", &a_muon_data.parentEk, "ParEk/F");
    tree_he->Branch("Weight", &a_muon_data.weight, "Weight/F");

    tree_pr->Branch("Ek", &a_muon_data.muonEk, "Ek/F");
    tree_pr->Branch("PDGID", &a_muon_data.pdgID, "PDGID/I");
    tree_pr->Branch("ParPDGID", &a_muon_data.parentPDGID, "ParPDGID/I");
    tree_pr->Branch("Pmu", &a_muon_data.polarization, "Pmu/F");
    tree_pr->Branch("Zenith", &a_muon_data.zenith, "Zenith/F");
    tree_pr->Branch("ParEk", &a_muon_data.parentEk, "ParEk/F");
    tree_pr->Branch("Weight", &a_muon_data.weight, "Weight/F");

    processMuonDataByEventWithWeight(he_mu_minus_df, 13, tree_he, a_muon_data);
    processMuonDataByEventWithWeight(he_mu_plus_df, -13, tree_he, a_muon_data);
    tree_he->Write();

    processMuonDataByEventWithWeight(pr_mu_minus_df, 13, tree_pr, a_muon_data);
    processMuonDataByEventWithWeight(pr_mu_plus_df, -13, tree_pr, a_muon_data);
    tree_pr->Write();

    output_file->Close();
}

void getPmu::processMuonDataByEventWithWeight(ROOT::RDataFrame &rawDataFrame, int muonID, TTree *tree, muonDataWithWeight &aMuonData)
{
    struct MuonCandidateWithWeight
    {
        float polarization;
        float muonEk;
        float zenith;
        int pdgID;
        int parentPDGID;
        int parentTrkID;
        float weight;
    };

    std::vector<MuonCandidateWithWeight> currentEventMuons;
    std::unordered_map<int, float> currentEventTrackEnergy;
    int currentEvtID = std::numeric_limits<int>::min();
    bool hasEvent = false;
    auto flushCurrentEvent = [&]()
    {
        for (const auto &muon : currentEventMuons)
        {
            aMuonData.polarization = muon.polarization;
            aMuonData.muonEk = muon.muonEk;
            aMuonData.zenith = muon.zenith;
            aMuonData.pdgID = muon.pdgID;
            aMuonData.parentPDGID = muon.parentPDGID;
            aMuonData.weight = muon.weight;

            auto parentIter = currentEventTrackEnergy.find(muon.parentTrkID);
            aMuonData.parentEk = parentIter == currentEventTrackEnergy.end() ? 0.0f : parentIter->second;
            tree->Fill();
        }

        currentEventMuons.clear();
        currentEventTrackEnergy.clear();
    };

    rawDataFrame.Foreach([&](int PDGID, std::string KillProc, int ParPDGID, float Ek, int ParTrkID, int TrkID, int EvtID, float Zenith, float Helicity, float Weight)
                         {
                             if (!hasEvent)
                             {
                                 currentEvtID = EvtID;
                                 hasEvent = true;
                             }

                             if (EvtID != currentEvtID)
                             {
                                 flushCurrentEvent();
                                 currentEvtID = EvtID;
                             }

                             currentEventTrackEnergy[TrkID] = Ek;

                             if (PDGID == muonID && KillProc == "<0|" && Zenith < 1 - std::sqrt(3) / 2)
                             {
                                 currentEventMuons.emplace_back(MuonCandidateWithWeight{Helicity, Ek, Zenith, PDGID, ParPDGID, ParTrkID, Weight});
                             } },
                         {"PDGID", "KillProc", "ParPDGID", "Ek", "ParTrkID", "TrkID", "EvtID", "Zenith", "Helicity", "Weight"});

    if (hasEvent)
    {
        flushCurrentEvent();
    }
}

std::__1::pair<TH1D *, TH1D *> getPmu::pr_he_weighted_hist(ROOT::RDataFrame &pr_mu_df, ROOT::RDataFrame &he_mu_df, float E_lower_bound, float E_upper_bound, std::string hist_header, int parentID, int binNumber)
{
    if (gROOT->FindObject("MuonPmu") && gROOT->FindObject("AntiPmu") && gROOT->FindObject("MuonPmuHe") && gROOT->FindObject("AntiPmuHe"))
    {
        delete gROOT->FindObject("MuonPmu");
        delete gROOT->FindObject("AntiPmu");
        delete gROOT->FindObject("MuonPmuHe");
        delete gROOT->FindObject("AntiPmuHe");
    }

    TH1D *MuonPmuP = new TH1D("MuonPmu", "MuonPmu", binNumber, -1, 1);
    TH1D *AntiMuonPmuP = new TH1D("AntiPmu", "AntiPmu", binNumber, -1, 1);
    TH1D *MuonPmuHe = new TH1D("MuonPmuHe", "MuonPmuHe", binNumber, -1, 1);
    TH1D *AntiMuonPmuHe = new TH1D("AntiPmuHe", "AntiPmuHe", binNumber, -1, 1);

    pr_mu_df.Foreach([&parentID, &E_lower_bound, &E_upper_bound, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu, float Weight)
                     {
                                if (std::abs(ParPDGID) == parentID && Ek > E_lower_bound && Ek < E_upper_bound)
                                {
                                    if (PDGID == 13)
                                    {
                                        MuonPmuP->Fill(Pmu, Weight);
                                    }
                                    if (PDGID == -13)
                                    {
                                        AntiMuonPmuP->Fill(Pmu, Weight);
                                    }
                                } }, {"PDGID", "ParPDGID", "Ek", "ParEk", "Pmu", "Weight"});

    he_mu_df.Foreach([&parentID, &E_lower_bound, &E_upper_bound, &MuonPmuP, &AntiMuonPmuP, &MuonPmuHe, &AntiMuonPmuHe](int PDGID, int ParPDGID, float Ek, float ParEk, float Pmu, float Weight)
                     {
                                 if (std::abs(ParPDGID) == parentID && Ek > E_lower_bound && Ek < E_upper_bound)
                                 {
                                     if (PDGID == 13)
                                     {
                                         MuonPmuHe->Fill(Pmu, Weight);
                                     }
                                     if (PDGID == -13)
                                     {
                                         AntiMuonPmuHe->Fill(Pmu, Weight);
                                     }
                                 } },
                     {"PDGID", "ParPDGID", "Ek", "ParEk", "Pmu", "Weight"});

    // combine two histograms into a new total histogram
    auto muonHistName = hist_header + "MuonPmu";
    auto antiMuonHistName = hist_header + "AntiMuonPmu";

    TH1D *totalMuonPmu = new TH1D(muonHistName.c_str(), muonHistName.c_str(), binNumber, -1, 1);
    TH1D *totalAntiMuonPmu = new TH1D(antiMuonHistName.c_str(), antiMuonHistName.c_str(), binNumber, -1, 1);

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

void getPmu::findCertainEnergyRangePolarizationWithWeight(std::string muonDataFile, std::string polarizationFile, int binNumber)
{
    // set the energy range for muon
    float E_lower_bound{100};
    float energy_upper_bound{500};

    ROOT::RDataFrame pr_mu_data("muonData(P)", muonDataFile.c_str());
    ROOT::RDataFrame he_mu_data("muonData(He)", muonDataFile.c_str());

    auto pi_hist = pr_he_weighted_hist(pr_mu_data, he_mu_data, E_lower_bound, energy_upper_bound, "pi", 211, binNumber);
    auto ka_hist = pr_he_weighted_hist(pr_mu_data, he_mu_data, E_lower_bound, energy_upper_bound, "K", 321, binNumber);
    auto kl_hist = pr_he_weighted_hist(pr_mu_data, he_mu_data, 0, 1e9, "KLong", 130, binNumber);

    TH1DDrawing(pi_hist.first, "polarization (P_{#mu^{-}} (#pi^{-}))", "Normalized Frequency");
    TH1DDrawing(pi_hist.second, "polarization (P_{#mu^{+}} (#pi^{+}))", "Normalized Frequency");

    TH1DDrawing(ka_hist.first, "polarization (P_{#mu^{-}} (K^{-}))", "Normalized Frequency");
    TH1DDrawing(ka_hist.second, "polarization (P_{#mu^{+}} (K^{+}))", "Normalized Frequency");

    TH1DDrawing(kl_hist.first, "polarization (P_{#mu^{-}} (K_{L}^{0}))", "Normalized Frequency");
    TH1DDrawing(kl_hist.second, "polarization (P_{#mu^{+}} (K_{L}^{0}))", "Normalized Frequency");

    TFile *outputFile = new TFile(polarizationFile.c_str(), "RECREATE");
    pi_hist.first->Write();
    pi_hist.second->Write();

    ka_hist.first->Write();
    ka_hist.second->Write();

    kl_hist.first->Write();
    kl_hist.second->Write();

    outputFile->Close();
}

void getPmu::gen_pol_with_wgt(std::string proton_produced_data, std::string helium_produced_data, std::string muon_data_file, std::string polarization_file, int bin_number)
{
    // findRelationPolarizationWithWeight(proton_produced_data, helium_produced_data, muon_data_file); // old map-based version
    findRelationPolarizationWithWeightByEvent(proton_produced_data, helium_produced_data, muon_data_file); // new event-streaming version
    findCertainEnergyRangePolarizationWithWeight(muon_data_file, polarization_file, bin_number);
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
    TH1D generateSinglePolarization(string polarizationFileName, int parentOption, bool charge);

    TH1D PmuMichelDis(double KPiRatio, bool charge, int numEvents, std::string polarizationFile); // generate the polarization distribution of the Michel electrons in a certain Pmu
    TH1D PmuEigenMichelDis(double KPiRatio, bool charge, int numEvents, TH1D polarizationHist);
    TH1D IsotropyDis(int numEvents); // generate the isotropy distribution of the Michel electrons
    TH2D PmuMichelDis2D(double KPiRatio, bool charge, int numEvents, std::string polarizationFile);
    TH1D getSingleParentMuon(int parentOption, int division, bool charge, std::string inputFile);
    void generateMichelMega(double KPiRatio, std::string inputFile, std::string outputFile);
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

        std::string hist_header = "KPiRatio = " + std::to_string(KPiRatio) + " Anti-muon";

        TH1D outputHist(hist_header.c_str(), hist_header.c_str(), piPolarization->GetNbinsX(), -1, 1);
        piPolarization->Scale(1 / piPolarization->Integral());
        kaPolarization->Scale(KPiRatio / kaPolarization->Integral());

        outputHist.Add(piPolarization, kaPolarization);
        outputHist.Scale(1 / outputHist.Integral());
        outputHist.SetDirectory(0);
        return outputHist;
    }
    else
    {
        auto piPolarization = (TH1D *)polarizationData.Get("piMuonPmu");
        auto kaPolarization = (TH1D *)polarizationData.Get("KMuonPmu");

        std::string hist_header = "KPiRatio = " + std::to_string(KPiRatio) + " Muon";

        TH1D outputHist(hist_header.c_str(), hist_header.c_str(), piPolarization->GetNbinsX(), -1, 1);
        piPolarization->Scale(1 / piPolarization->Integral());
        kaPolarization->Scale(KPiRatio / kaPolarization->Integral());

        outputHist.Add(piPolarization, kaPolarization);
        outputHist.Scale(1 / outputHist.Integral());
        outputHist.SetDirectory(0);
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

TH1D GenMichel::PmuEigenMichelDis(double KPiRatio, bool charge, int numEvents, TH1D polarizationHist)
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

TH1D GenMichel::generateSinglePolarization(string polarizationFileName, int parentOption, bool charge)
{
    TFile polarizationFile(polarizationFileName.c_str(), "READ");
    string parent, son;
    if (charge)
    {
        son = "AntiMuon";
    }
    else
    {
        son = "Muon";
    }

    if (parentOption == 0)
    {
        parent = "pi";
    }
    else if (parentOption == 1)
    {
        parent = "K";
    }
    else if (parentOption == 2)
    {
        parent = "KLong";
    }

    string histName = parent + son + "Pmu";
    TH1D outputHist(parent.c_str(), parent.c_str(), 100, -1, 1);

    auto polarizationHist = (TH1D *)polarizationFile.Get(histName.c_str());

    outputHist.Add(polarizationHist);

    outputHist.SetDirectory(0);
    return outputHist;
}

TH1D GenMichel::getSingleParentMuon(int parentOption, int division, bool charge, std::string inputFile)
{
    string parent, son;
    if (charge)
    {
        son = "AntiMuon";
    }
    else
    {
        son = "Muon";
    }

    if (parentOption == 0)
    {
        parent = "pi";
    }
    else if (parentOption == 1)
    {
        parent = "K";
    }
    else if (parentOption == 2)
    {
        parent = "KLong";
    }

    string HistID = parent + son;

    TH1D MichelDis(HistID.c_str(), HistID.c_str(), 100, -1, 1);

    auto polarizationHist = generateSinglePolarization(inputFile, parentOption, charge);

    TRandom random(time(0));
    for (int i = 0; i < division; i++)
    {
        double Pmu = polarizationHist.GetRandom();

        auto coordinateTuple = generator(random, Pmu);
        double cosTheta = std::get<0>(coordinateTuple);
        double weight = std::get<2>(coordinateTuple);

        MichelDis.Fill(cosTheta, weight);
    }

    return MichelDis;
}

void GenMichel::generateMichelMega(double KPiRatio, std::string polarizationFile, std::string outputFile)
{
    TFile polarizationData(polarizationFile.c_str(), "READ");

    auto piPolarization = (TH1D *)polarizationData.Get("piMuonPmu");
    auto kaPolarization = (TH1D *)polarizationData.Get("KMuonPmu");
    auto piPolarizationAnti = (TH1D *)polarizationData.Get("piAntiMuonPmu");
    auto kaPolarizationAnti = (TH1D *)polarizationData.Get("KAntiMuonPmu");

    TH1D piMichelDis("piMichelDis", "piMichelDis", 100, -1, 1);
    TH1D kaMichelDis("kaMichelDis", "kaMichelDis", 100, -1, 1);
    TH1D piMichelDisAnti("piMichelDisAnti", "piMichelDisAnti", 100, -1, 1);
    TH1D kaMichelDisAnti("kaMichelDisAnti", "kaMichelDisAnti", 100, -1, 1);

    auto tempHist = PmuEigenMichelDis(KPiRatio, false, 1e9, *piPolarization);
    piMichelDis.Add(&tempHist);
    tempHist.Clear();

    tempHist = PmuEigenMichelDis(KPiRatio, true, 1e9, *piPolarizationAnti);
    piMichelDisAnti.Add(&tempHist);
    tempHist.Clear();

    tempHist = PmuEigenMichelDis(KPiRatio, false, 1e9, *kaPolarization);
    kaMichelDis.Add(&tempHist);
    tempHist.Clear();

    tempHist = PmuEigenMichelDis(KPiRatio, true, 1e9, *kaPolarizationAnti);
    kaMichelDisAnti.Add(&tempHist);
    tempHist.Clear();

    TFile output(outputFile.c_str(), "RECREATE");

    piMichelDis.Write();
    piMichelDisAnti.Write();
    kaMichelDis.Write();
    kaMichelDisAnti.Write();

    output.Close();
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