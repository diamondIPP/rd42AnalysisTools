#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <TUnixSystem.h>

using namespace std;

void ReadInputs(int argc, char **argv);
void PrintHelp();
void GetValuesFromTextFile();
void SetBranches(TTree *out_tree);
void LoopSourceTree(TTree *tree0, TTree *out_tree);
//void LoopSourceTree(TTree *tree0, TTree *out_tree, UChar_t (*pDet_ADCin)[8][256], UShort_t (*pDiaADCin)[128], UInt_t *pEventNumberin);
void SetBranchAddresses(TTree *in_tree);
//void SetBranchAddresses(TTree *in_tree, UChar_t (*pDet_ADCin)[8][256], UShort_t (*pDiaADCin)[128], UInt_t *pEventNumberin);
void ProgressBar(Long64_t i, Long64_t nEntries, int refresh, bool show, bool newLine);
void UpdateSilicon();
void UpdateDiamond();
void ResetArrays();

string run_dir = ".";
UShort_t run=0;
string outfile = "default";
Float_t sil = 0;
Float_t sil_each[8];
bool sil_each_pos[8];
UChar_t max_sil = 255;
UShort_t max_dia = 4095;
UChar_t measured_sil_adc = 0;
UShort_t measured_dia_adc = 0;
Float_t dia = 0;
UChar_t Det_ADC[8][256];
UChar_t Det_ADCin[8][256];
UShort_t DiaADC[128];
UShort_t DiaADCin[128];
UInt_t EventNumber = 0;
UInt_t EventNumberin = 0;
TBranch *bDet_ADCin[8];
TBranch *bDiaADCin;
TBranch *bEventNumberin;
Double_t real_sil_adc = 0;
UChar_t newSilADC = 0;
UChar_t prev_sil_adc = 0;
Double_t real_dia_adc = 0;
UShort_t newDiaADC = 0;
UShort_t prev_dia_adc = 0;
Long64_t entries_old = 0;

int main(int argc, char **argv) {

    ResetArrays();

    for(int det = 0; det < 8; det++){
        sil_each[det] = 0;
        sil_each_pos[det] = false;
    }
//    UChar_t (*pDet_ADCin)[8][256] = &Det_ADCin;
//    UShort_t (*pDiaADCin)[128] = &DiaADCin;
//    UInt_t *pEventNumberin = &EventNumberin;

    std::cout << "Starting CrossTalk corrections" << std::endl;
    ReadInputs(argc, argv);
    TString file_path = TString::Format("%s/rawData.%d.root", run_dir.c_str(), run);
    cout << "source file: " << file_path << endl;
    TFile *file0 = new TFile(file_path, "READ");
    TTree *tree0 = (TTree*)file0->Get("rawTree");
    tree0->SetName("rawTreeOld");
    entries_old = tree0->GetEntries();
    SetBranchAddresses(tree0);
    if(sil==0 || dia==0)
        GetValuesFromTextFile();
    TString out_name = outfile == "default"? TString::Format("rawData.%d-%05d-%05d.root", run, int(sil*10000), int(dia*10000)) : TString::Format("%s", outfile.c_str());
    TString out_file_path = TString::Format("%s/%s", run_dir.c_str(), out_name.Data());
    cout << "Corrected file will be: " << out_name << endl;
    TFile *file1 = new TFile(out_file_path, "RECREATE");
    TTree *out_tree = new TTree("rawTree", "rawTree");
    out_tree->Reset();
    SetBranches(out_tree);
//    LoopSourceTree(tree0, out_tree, pDet_ADCin, pDiaADCin, pEventNumberin);
    LoopSourceTree(tree0, out_tree);
    out_tree->Print();
    file1->Write();
    file0->Close();
//    file0->Close();
//    out_tree->SetNameTitle("rawTree", "rawTree");
//    out_tree->Write("rawTree");
    file1->Close();

    TUnixSystem bla;
    bla.Symlink(out_file_path, TString::Format("%s/rawData.%d0.root", run_dir.c_str(), run));
    return 0;
}

void ReadInputs(int argc, char **argv){
    for(int i=1; i< argc; i++){
        string string_arg(argv[i]);
        if(string_arg == "-h" || string_arg == "--help") PrintHelp();
    }
    for(int i=1; i < argc; i++){
        string string_arg(argv[i]);
        if((string_arg == "-o") && i+1 < argc){
            i++;
            run_dir = string(argv[i]);
        }
        if((string_arg == "-r") && i+1 < argc){
            i++;
            run = UShort_t(stoi(string(argv[i])));
        }
        if((string_arg == "-s") && i+1 < argc){
            i++;
            sil = Float_t(stof(string(argv[i])));
        }
        if((string_arg == "-s0") && i+1 < argc){
            i++;
            sil_each[0] = Float_t(stof(string(argv[i])));
            sil_each_pos[0] = true;
        }
        if((string_arg == "-s1") && i+1 < argc){
            i++;
            sil_each[1] = Float_t(stof(string(argv[i])));
            sil_each_pos[1] = true;
        }
        if((string_arg == "-s2") && i+1 < argc){
            i++;
            sil_each[2] = Float_t(stof(string(argv[i])));
            sil_each_pos[2] = true;
        }
        if((string_arg == "-s3") && i+1 < argc){
            i++;
            sil_each[3] = Float_t(stof(string(argv[i])));
            sil_each_pos[3] = true;
        }
        if((string_arg == "-s4") && i+1 < argc){
            i++;
            sil_each[4] = Float_t(stof(string(argv[i])));
            sil_each_pos[4] = true;
        }
        if((string_arg == "-s5") && i+1 < argc){
            i++;
            sil_each[5] = Float_t(stof(string(argv[i])));
            sil_each_pos[5] = true;
        }
        if((string_arg == "-s6") && i+1 < argc){
            i++;
            sil_each[6] = Float_t(stof(string(argv[i])));
            sil_each_pos[6] = true;
        }
        if((string_arg == "-s7") && i+1 < argc){
            i++;
            sil_each[7] = Float_t(stof(string(argv[i])));
            sil_each_pos[7] = true;
        }
        if((string_arg == "-d") && i+1 < argc){
            i++;
            dia = Float_t(stof(string(argv[i])));
        }
        if((string_arg == "-f") && i+1 < argc){
            i++;
            outfile = string(argv[i]);
        }
        if((string_arg == "-ss") && i+1 < argc){
            i++;
            max_sil = UChar_t(stoi(string(argv[i])));
        }
        if((string_arg == "-ds") && i+1 < argc){
            i++;
            max_dia = UShort_t(stoi(string(argv[i])));
        }
    }
    if(run==0){
        cout << "run number was not given. exit" << endl;
        PrintHelp();
    }
    for(int det = 0; det < 8; det ++){
        if(!sil_each_pos[det])
            sil_each[det] = sil;
    }
    if(!sil_each[0]&!sil_each[1]&!sil_each[2]&!sil_each[3]&!sil_each[4]&!sil_each[5]&!sil_each[6]&!sil_each[7]&!dia){
        cout << "no parameter was given for silicon or diamond detectors. Won't do anything. Exiting";
        exit(0);
    }
    cout << "will use a factor of " << dia << "% for diamond and " << sil << "% in general for silicon and " << sil_each[0] << "%, "<< sil_each[1] << "%, "<< sil_each[2] << "%, ";
    cout << sil_each[3] << "%, "<< sil_each[4] << "%, "<< sil_each[5] << "%, " << sil_each[6] << "%, " << sil_each[7] << "% for silicon planes: ";
    cout << "0, 1, 2, 3, 4, 5, 6 and 7 respectively."<<endl;
    cout << "Using saturation on silicon of " << int(max_sil) << " and on diamond of " << int(max_dia) << "." << endl;
}

void PrintHelp(){
    cout << "Showing help (this message is shown when -h flag is used)" << endl;
    cout << "Usage:" << endl;
    cout << "crossTalkCorrection -o RUNDIRECTORY -r RUNNUMBER -s SILALPHA -d DIAALPHA -ss SILSAT -ds DIASAT"  << endl;
    cout << "The values for SILALPHA and DIAALPHA are percentages of charge sharing estimated in the crossTalkCorrectionFactors text file"  << endl;
    cout << "Separate values can be given to silicon in addition to SILALPHA by setting -s0, -s1, -s2, ..., -s7" << endl;
    cout << "If one of them is not given, the program will use the text file." << endl;
    cout << "The parameters SILSAT and DIASAT are optional and they set the saturation value for silicon and diamond respectively." << endl;
    cout << "If desired, can specify the name for the new file with -f OUTNAME"<<endl;
    cout << "If not specified, the output will be of the type 'rawData.<RUNNUMBER>0.root'" << endl;
    exit(0);
}

void GetValuesFromTextFile(){
    //TODO
    cout << "Haven't implemented yet auto get values in this C script... sorry. But rd42Analysis.py does it automatically :D. Exiting" << endl;
    exit(0);
}

void SetBranches(TTree *out_tree){
    out_tree->Branch("D0X_ADC", &Det_ADC[0], "D0X_ADC[256]/b");
    out_tree->Branch("D0Y_ADC", &Det_ADC[1], "D0Y_ADC[256]/b");
    out_tree->Branch("D1X_ADC", &Det_ADC[2], "D1X_ADC[256]/b");
    out_tree->Branch("D1Y_ADC", &Det_ADC[3], "D1Y_ADC[256]/b");
    out_tree->Branch("D2X_ADC", &Det_ADC[4], "D2X_ADC[256]/b");
    out_tree->Branch("D2Y_ADC", &Det_ADC[5], "D2Y_ADC[256]/b");
    out_tree->Branch("D3X_ADC", &Det_ADC[6], "D3X_ADC[256]/b");
    out_tree->Branch("D3Y_ADC", &Det_ADC[7], "D3Y_ADC[256]/b");
    out_tree->Branch("DiaADC", DiaADC, "DiaADC[128]/s");
    out_tree->Branch("EventNumber", &EventNumber, "EventNumber/i");
}

//void LoopSourceTree(TTree *tree0, TTree *out_tree, UChar_t (*pDet_ADCin)[8][256], UShort_t (*pDiaADCin)[128], UInt_t *pEventNumberin){
void LoopSourceTree(TTree *tree0, TTree *out_tree){
//    SetBranchAddresses(tree0, pDet_ADCin, pDiaADCin, pEventNumberin);
//    SetBranchAddresses(tree0);
//    tree0->SetBranchStatus("*", 1);
//    UInt_t nEntries = UInt_t(tree0->GetEntries());
    for (Long64_t i=0; i<entries_old; i++){
        ProgressBar(i, entries_old, 100, false, false);
        tree0->GetEntry(i);
        EventNumber = EventNumberin;
        UpdateSilicon();
        UpdateDiamond();
        out_tree->Fill();
        ResetArrays();
    }
}

//void SetBranchAddresses(TTree *in_tree, UChar_t (*pDet_ADCin)[8][256], UShort_t (*pDiaADCin)[128], UInt_t *pEventNumberin){
void SetBranchAddresses(TTree *in_tree){
    if(in_tree->FindBranch("EventNumber"))
        in_tree->SetBranchAddress("EventNumber", &EventNumberin);
    if(in_tree->FindBranch("D0X_ADC"))
        in_tree->SetBranchAddress("D0X_ADC", &Det_ADCin[0]);
    if(in_tree->FindBranch("D0Y_ADC"))
        in_tree->SetBranchAddress("D0Y_ADC", &Det_ADCin[1]);
    if(in_tree->FindBranch("D1X_ADC"))
        in_tree->SetBranchAddress("D1X_ADC", &Det_ADCin[2]);
    if(in_tree->FindBranch("D1Y_ADC"))
        in_tree->SetBranchAddress("D1Y_ADC", &Det_ADCin[3]);
    if(in_tree->FindBranch("D2X_ADC"))
        in_tree->SetBranchAddress("D2X_ADC", &Det_ADCin[4]);
    if(in_tree->FindBranch("D2Y_ADC"))
        in_tree->SetBranchAddress("D2Y_ADC", &Det_ADCin[5]);
    if(in_tree->FindBranch("D3X_ADC"))
        in_tree->SetBranchAddress("D3X_ADC", &Det_ADCin[6]);
    if(in_tree->FindBranch("D3Y_ADC"))
        in_tree->SetBranchAddress("D3Y_ADC", &Det_ADCin[7]);
    if(in_tree->FindBranch("DiaADC"))
        in_tree->SetBranchAddress("DiaADC", DiaADCin);
}

void ProgressBar(Long64_t i, Long64_t nEntries, int refresh, bool show, bool newLine){
    if(i+1 >= nEntries) i++;
    cout.precision(3);
    int percentageLength = 50;
    if(i%(int)refresh == 0 || i >= nEntries - 1 || show){
        float percentage = (float)i / float(nEntries) * float(100);
        cout<<"\rFinished with " << std::setw(8) << i << " of " << std::setw(10) << nEntries << ": " << std::setw(6)<<fixed<<percentage<<"%\t\tSTATUS:\t\t";
        for(int j = 0; j < percentageLength; j++){
            if(j*10 < percentage * float(percentageLength) / float(10)) cout << "%";
            else cout << "-";
        }
        cout << " " << flush;
    }
    if(newLine && i + 1 >= nEntries) cout << endl;
}

void UpdateSilicon(){
    int det_max = 8;
    Int_t signalpha = 0;
    Double_t alpha = 0;
    UShort_t startCh = 0;
    UShort_t endCh = 0;
    for (int det = 0; det < det_max; det++){
        if(sil_each_pos[det]){
            signalpha = sil_each[det] >= 0 ? 1 : -1;
            alpha = fabs(sil_each[det] / 100.0);
        }
        else{
            signalpha = sil >= 0 ? 1 : -1;
            alpha = fabs(sil) / 100.0;
        }
        startCh = signalpha > 0? 0 : 255;
//        endCh = 127;
        endCh = signalpha > 0? 255 : 0;
//        if(det == 2 || det == 6){
////            startCh = 127;
//            startCh = signalpha > 0? 255 : 0;
//            endCh = signalpha > 0? 0 : 255;
//        }
        bool finished = false;
//        Double_t alpha = sil_each[det] / 100.0;
//        UChar_t adc = 0;
        prev_sil_adc = 0;
        for (UShort_t ch = startCh; !finished;){
            newSilADC = 0;
            real_sil_adc = 0;
            measured_sil_adc = Det_ADCin[det][ch];
            real_sil_adc = (Double_t(measured_sil_adc) - Double_t(prev_sil_adc) * alpha) / (1.0 - alpha);
            if(real_sil_adc < 0){
                newSilADC = 0;
            }
            else if(real_sil_adc >= float(max_sil)){
                newSilADC = max_sil;
            }
            else{
                newSilADC = UShort_t(std::floor(real_sil_adc + 0.5));
            }
            Det_ADC[det][ch] = newSilADC;
            prev_sil_adc = newSilADC;
            finished = (ch == endCh);

//            if(det == 2 || det == 6){
//                ch = signalpha > 0? ch - 1 : ch + 1;
//            }
//            else{
                ch = signalpha > 0? ch + 1 : ch - 1;
//            }
        }

//        for (UShort_t ch = startCh; !finished;){
//            UChar_t measured_adc = Det_ADCin[det][ch];
//  //          if((measured_adc < max_sil) && (adc < max_sil)){
//                Double_t real_adc = Double_t(measured_adc) - Double_t(adc) * alpha / (1.0 - alpha);
//                UChar_t newADC = real_adc >= max_sil? max_sil : UChar_t(std::floor(real_adc + 0.5));
//                Det_ADC[det][ch] = newADC <= 0? 0 : newADC;
//                adc = Det_ADC[det][ch];
//                finished = (ch==endCh);
//          //  }
//          //  else{
//          //      Det_ADC[det][ch] = measured_adc;
//          //      adc = Det_ADC[det][ch];
//          //  }
//            if(det==2 || det==6)
//                ch--;
//            else
//                ch++;
//        }
    }
}

void UpdateDiamond(){
//    UShort_t prev_adc = 0;
    prev_dia_adc = 0;
    TString out = TString::Format("%d_%5.3f: ", EventNumber, dia);
//    bool bPrint = false;
    Int_t signalpha = dia >= 0? 1 : -1;
    Double_t alpha = fabs(dia / 100.0);

    if (signalpha > 0){
        for (int ch = 0; ch < 128; ch++){
            newDiaADC = 0;
            real_dia_adc = 0;
            measured_dia_adc = DiaADCin[ch];
            real_dia_adc = (Double_t(measured_dia_adc) - Double_t(prev_dia_adc) * alpha) / (1.0 - alpha);
            if(real_dia_adc < 0){
                newDiaADC = 0;
            }
            else if(real_dia_adc >= float(max_dia)){
                newDiaADC = max_dia;
            }
            else{
                newDiaADC = UShort_t(std::floor(real_dia_adc + 0.5));
            }
            DiaADC[ch] = newDiaADC;
            prev_dia_adc = newDiaADC;
        }
    }
    else{
        for(int ch = 127; ch >= 0; ch--){
            newDiaADC = 0;
            real_dia_adc = 0;
            measured_dia_adc = DiaADCin[ch];
            real_dia_adc = (Double_t(measured_dia_adc) - Double_t(prev_dia_adc) * alpha) / (1.0 - alpha);
            if(real_dia_adc <= 0){
                newDiaADC = 0;
            }
            else if(real_dia_adc >= float(max_dia)){
                newDiaADC = max_dia;
            }
            else{
                newDiaADC = UShort_t(std::floor(real_dia_adc + 0.5));
            }
            DiaADC[ch] = newDiaADC;
            prev_dia_adc = newDiaADC;
        }
    }
//    for (int ch=0; ch<128; ch++){
//        measured_dia_adc = DiaADCin[ch];
//        if((measured_dia_adc < max_dia) && (prev_adc < max_dia) ){
//            Double_t real_adc = (Double_t(measured_dia_adc) - Double_t(prev_adc) * alpha) / (1.0 - alpha);
//            UShort_t newADC = real_adc >= max_dia? max_dia : UShort_t(std::floor(real_adc + 0.5));
//            newADC = newADC <= 0? 0 : newADC;
//            newADC = newADC > max_dia? max_dia : newADC;
//            DiaADC[ch] = newADC;
//            prev_adc = DiaADC[ch];
//        }
//        else{
//            DiaADC[ch] = measured_dia_adc;
//            prev_adc = DiaADC[ch];
//        }
//    }
}

void ResetArrays(){
    EventNumberin = 0;
    EventNumber = 0;
    for(int det = 0; det <= 8; det++){
        if(det < 8){
//            sil_each[det] = 0;
//            sil_each_pos[det] = false;
            for(int ch = 0; ch < 256; ch++){
                Det_ADC[det][ch] = 0;
                Det_ADCin[det][ch] = 0;
            }
        }
        else{
            for(int ch = 0; ch < 128; ch++){
                DiaADC[ch] = 0;
                DiaADCin[ch] = 0;
            }
        }
    }
}
