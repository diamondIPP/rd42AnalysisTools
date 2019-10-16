#include <iostream>
#include <cstdlib>
#include <sstream>
#include <numeric>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <deque>
#include <TMath.h>
#include <TLeaf.h>
#include <TNamed.h>
#include "time.h"

int main(int argc, char* argv[]) {
    float_t time0 = clock();
    if (argc < 2) {
        std::cout << "Hello, World!. You have fucked up XD. give the path of the raw root file 'testData.root'"
                  << std::endl;
        return 0;
    }
    // get file with raw data
    TString InFileName = TString(argv[1]);
    TString StemName = InFileName.Copy();
    float_t threshold = (argc > 2)? strtof(argv[2], nullptr) : 10;
    StemName.ReplaceAll(".root", "");
    auto *blaf = new TFile(InFileName, "READ");
    auto *tree = (TTree*)blaf->Get("rawTree");
    uint32_t buff = 500;
    auto dets = uint8_t(tree->GetLeaf("adc")->GetLen());
    auto maxEvents = uint32_t(tree->GetEntries());
    uint16_t saturatedValue = 4095;

    // extra info in input file
    float_t inLandauPos = 0;
    float_t inLandauSc = 0;
    uint32_t inBuff = 0;
    uint16_t inPedMean = 0;
    float_t inPedSpread = 0;
    float_t inNoise = 0;
    float_t inCM = 0;

    if(blaf->Get("landau_pos"))
        inLandauPos = strtof(((TNamed *)blaf->Get("landau_pos"))->GetTitle(), nullptr);
    if(blaf->Get("landau_sc"))
        inLandauSc = strtof(((TNamed *)blaf->Get("landau_sc"))->GetTitle(), nullptr);
    if(blaf->Get("buff"))
        inBuff = uint32_t(strtol(((TNamed *)blaf->Get("buff"))->GetTitle(), nullptr, 10));
    if(blaf->Get("ped_mean"))
        inPedMean = uint16_t(strtoul(((TNamed *)blaf->Get("ped_mean"))->GetTitle(), nullptr, 10));
    if(blaf->Get("ped_spread"))
        inPedSpread = strtof(((TNamed *)blaf->Get("ped_spread"))->GetTitle(), nullptr);
    if(blaf->Get("noise"))
        inNoise = strtof(((TNamed *)blaf->Get("noise"))->GetTitle(), nullptr);
    if(blaf->Get("commonMode"))
        inCM = strtof(((TNamed *)blaf->Get("commonMode"))->GetTitle(), nullptr);

    // define variables to read rawTree
    uint16_t adcs[dets];
    uint32_t event;

    // define variables in pedTree
    uint8_t channels[dets];
    for(uint8_t deti=0; deti < dets; deti++){
        channels[deti] = deti;
    }
    float_t cm;
    float_t ped[dets];
    float_t pedCMC[dets];
    float_t sigm[dets];
    float_t sigmCMC[dets];
    float_t signal[dets];
    float_t signalCMC[dets];
    bool isPed[dets];
    bool isPedCMC[dets];
    bool isSaturated[dets];


    // define variables to esetimate noise and pedestals
    std::deque <uint16_t > detADC[dets];
    std::deque <float_t > detADCCMC[dets];
    std::deque <bool> evtIsPed[dets];
    std::deque <bool> evtIsPedCMC[dets];
    double_t detSum[dets];
    double_t detSumCMC[dets];
    double_t detSum2[dets];
    double_t detSum2CMC[dets];
    uint16_t numDetSum[dets];
    uint16_t numDetSumCMC[dets];

    // clear deques
    for(uint8_t det=0; det < dets; det++){
        detADC[det].clear();
        detADCCMC[det].clear();
        evtIsPed[det].clear();
        evtIsPedCMC[det].clear();
        detSum[det] = 0;
        detSum2[det] = 0;
        detSumCMC[det] = 0;
        detSum2CMC[det] = 0;
        numDetSum[det] = 0;
        numDetSumCMC[det] = 0;
    }

    // set variables to input tree
    if(tree->FindBranch("event")){
        tree->SetBranchAddress("event", &event);
    }
    if(tree->FindBranch("adc")){
        tree->SetBranchAddress("adc", &adcs);
    }

    // create pedestal file and tree
    auto *pedFile = new TFile(TString::Format("%sPed.root", StemName.Data()), "recreate");
    pedFile->cd();
    auto *pedTree = new TTree("pedTree", "pedTree");

    // assign branches to pedestal tree
    pedTree->Branch("channels", &channels, TString::Format("channels[%i]/b", dets));
    pedTree->Branch("cm", &cm, "cm/F");
    pedTree->Branch("ped", &ped, TString::Format("ped[%i]/F", dets));
    pedTree->Branch("pedCMC", &pedCMC, TString::Format("pedCMC[%i]/F", dets));
    pedTree->Branch("sigma", &sigm, TString::Format("sigma[%i]/F", dets));
    pedTree->Branch("sigmaCMC", &sigmCMC, TString::Format("sigmaCMC[%i]/F", dets));
    pedTree->Branch("signal", &signal, TString::Format("signal[%i]/F", dets));
    pedTree->Branch("signalCMC", &signalCMC, TString::Format("signalCMC[%i]/F", dets));
    pedTree->Branch("isPed", &isPed, TString::Format("isPed[%i]/O", dets));
    pedTree->Branch("isPedCMC", &isPedCMC, TString::Format("isPedCMC[%i]/O", dets));
    pedTree->Branch("isSaturated", &isSaturated, TString::Format("isSaturated[%i]/O", dets));

    // Loop over first buffered events:
    std::cout << "Buffer estimation:";std::cout << std::flush;
    for(uint32_t eventi = 0; eventi < buff; eventi++){
        tree->GetEvent(eventi);
        for(uint8_t det = 0; det < dets; det ++){
            detADC[det].push_back(adcs[det]);
            evtIsPed[det].push_back(true);
            detSum[det] += adcs[det];
            detSum2[det] += adcs[det] * adcs[det];
            if(eventi == 0){
                numDetSum[det] = buff;
            }
        }
        if(eventi % int(buff / 10) ==  0)
            std::cout << int(eventi * 10 / buff); std::cout << std::flush;
    }
    std::cout << std::endl; std::cout<<std::flush;
    for(uint8_t det = 0; det < dets; det++){
        ped[det] = float(detSum[det]) / float(numDetSum[det]);
        sigm[det] = float(TMath::Sqrt(float(detSum2[det]) / float(numDetSum[det]) - ped[det] * ped[det]));
    }
    std::cout << "Filling tree... "; std::cout << std::flush;
    for(uint32_t eventi = 0; eventi < buff; eventi++){
        tree->GetEvent(eventi);
        uint32_t cmevts = 0;
        cm = 0;
        float_t tempAdc;
        for(uint8_t det = 0; det < dets; det++){
            tempAdc = adcs[det];
            signal[det] = tempAdc - ped[det];
            cm += signal[det];
            cmevts ++;
            isSaturated[det] = bool(adcs[det] >= saturatedValue);
        }
        cm = cmevts > 0? float(cm) / float(cmevts) : 0;

        for(uint8_t det = 0; det < dets; det++){
            tempAdc = adcs[det] - cm;
            detADCCMC[det].push_back(tempAdc);
            evtIsPedCMC[det].push_back(true);
            detSumCMC[det] += tempAdc;
            detSum2CMC[det] += tempAdc * tempAdc;
            numDetSumCMC[det] ++;

            if(numDetSumCMC[det] > 0)
                pedCMC[det] = float(detSumCMC[det]) / float(numDetSumCMC[det]);
            else
                pedCMC[det] = ped[det];
            if(numDetSumCMC[det] > 1 && (detSum2CMC[det] > pedCMC[det] * pedCMC[det] * numDetSumCMC[det]))
                sigmCMC[det] = float(TMath::Sqrt(float(detSum2CMC[det]) / float(numDetSumCMC[det]) - pedCMC[det] * pedCMC[det]));
            else
                sigmCMC[det] = sigm[det];
            isPed[det] = true;
            isPedCMC[det] = true;
            signalCMC[det] = tempAdc - pedCMC[det];
        }
        pedTree->Fill();
        if(eventi % int(maxEvents / 10) == 0)
            std::cout << int(eventi * 10 / maxEvents); std::cout << std::flush;
    }

    // Loop over rest of events:
    for(uint32_t eventi = buff; eventi < maxEvents; eventi++){
        tree->GetEvent(eventi);
        // estimate cm
        uint32_t cmevts = 0;
        cm = 0;
        float_t tempMean, tempMeanCMC, tempSigma, tempSigmaCMC, tempSignal;
        for(int det = 0; det < dets; det++){
            tempSignal = adcs[det] - pedCMC[det];
            tempSigma = sigmCMC[det];
            if(TMath::Abs(tempSignal) >= threshold * tempSigma)
                continue;
            cm += tempSignal;
            cmevts++;
        }
        cm = cmevts > 0? float(cm) / float(cmevts) : 0;

        for(uint8_t det = 0; det < dets; det++){
            // no cmc
            detADC[det].push_back(adcs[det]);
            tempMean = float(detSum[det] / float(numDetSum[det]));
            tempSigma = float(TMath::Sqrt(float(detSum2[det] / float(numDetSum[det]) - tempMean * tempMean)));
            isPed[det] = false;


            if(evtIsPed[det].front()){
                detSum[det] -= detADC[det].front();
                detSum2[det] -= detADC[det].front() * detADC[det].front();
                numDetSum[det]--;
            }
            if(TMath::Abs(detADC[det].back() - tempMean) < threshold * tempSigma){
                detSum[det] += detADC[det].back();
                detSum2[det] += detADC[det].back() * detADC[det].back();
                numDetSum[det]++;
                evtIsPed[det].push_back(true);
                isPed[det] = true;
            }
            else
                evtIsPed[det].push_back(false);
            ped[det] = float(detSum[det] / float(numDetSum[det]));
            sigm[det] = float(TMath::Sqrt(float(detSum2[det] / float(numDetSum[det]) - ped[det] * ped[det])));
            signal[det] = adcs[det] - ped[det];
            isSaturated[det] = bool(adcs[det] >= saturatedValue);

            // with cmc
            detADCCMC[det].push_back(adcs[det] - cm);
            tempMeanCMC = float(detSumCMC[det] / float(numDetSumCMC[det]));
            tempSigmaCMC = float(TMath::Sqrt(float(detSum2CMC[det] / float(numDetSumCMC[det]) - tempMeanCMC * tempMeanCMC)));
            isPedCMC[det] = false;

            if(evtIsPedCMC[det].front()){
                detSumCMC[det] -= detADCCMC[det].front();
                detSum2CMC[det] -= detADCCMC[det].front() * detADCCMC[det].front();
                numDetSumCMC[det]--;
            }
            if(TMath::Abs(detADCCMC[det].back() - tempMeanCMC) < threshold * tempSigmaCMC){
                detSumCMC[det] += detADCCMC[det].back();
                detSum2CMC[det] += detADCCMC[det].back() * detADCCMC[det].back();
                numDetSumCMC[det]++;
                evtIsPedCMC[det].push_back(true);
                isPedCMC[det] = true;
            }
            else
                evtIsPedCMC[det].push_back(false);
            pedCMC[det] = float(detSumCMC[det] / float(numDetSumCMC[det]));
            sigmCMC[det] = float(TMath::Sqrt(float(detSum2CMC[det] / float(numDetSumCMC[det]) - pedCMC[det] * pedCMC[det])));
            signalCMC[det] = adcs[det] - pedCMC[det] - cm;

            if(!detADC[det].empty()) detADC[det].pop_front();
            if(!detADCCMC[det].empty()) detADCCMC[det].pop_front();
            if(!evtIsPed[det].empty()) evtIsPed[det].pop_front();
            if(!evtIsPedCMC[det].empty()) evtIsPedCMC[det].pop_front();
        }
        pedTree->Fill();
        if(eventi % int(maxEvents / 10) == 0) {
            std::cout << int(eventi * 10 / maxEvents);
            std::cout << std::flush;
        }
    }
    std::cout << std::endl;

    pedFile->cd();
    pedTree->AddFriend("rawTree", InFileName);
    pedFile->cd();
    pedTree->Write();

    // add info to new file
    auto landauPosInfo = new TNamed("landau_pos", TString::Format("%f", inLandauPos).Data());
    auto landauScInfo = new TNamed("landau_sc", TString::Format("%f", inLandauSc).Data());
    auto buffInfo = new TNamed("buff", TString::Format("%d", inBuff).Data());
    auto pedMeanInfo = new TNamed("ped_mean", TString::Format("%d", inPedMean).Data());
    auto pedSpreadInfo = new TNamed("ped_spread", TString::Format("%f", inPedSpread).Data());
    auto noiseInfo = new TNamed("noise", TString::Format("%f", inNoise).Data());
    auto cmInfo = new TNamed("commonMode", TString::Format("%f", inCM).Data());

    landauPosInfo->Write();
    landauScInfo->Write();
    buffInfo->Write();
    pedMeanInfo->Write();
    pedSpreadInfo->Write();
    noiseInfo->Write();
    cmInfo->Write();

//    pedFile->Write();
//    delete blaf;
//    delete pedFile;
//    pedTree->Delete();
//    tree->Delete();
    blaf->Close();
//    delete blaf;
    pedFile->Close();
//    delete pedFile;

    std::cout << "FINISHED in " << (clock() - time0)/CLOCKS_PER_SEC << " seconds" << std::endl;
    return 0;

}
