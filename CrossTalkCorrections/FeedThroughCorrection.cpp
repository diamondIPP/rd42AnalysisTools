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
void SetBranchAddresses(TTree *in_tree);
void ProgressBar(int i, UInt_t nEntries, int refresh, bool show, bool newLine);
void UpdateSilicon();
void UpdateDiamond();

string run_dir = ".";
UShort_t run=0;
string outfile = "default";
Float_t sil = 0;
Float_t dia = 0;
UChar_t Det_ADC[8][256];
UChar_t Det_ADCin[8][256];
UShort_t DiaADC[128];
UShort_t DiaADCin[128];
UInt_t EventNumber;
UInt_t EventNumberin;
TBranch *bDet_ADCin;
TBranch *bDiaADCin;
TBranch *bEventNumberin;

int main(int argc, char **argv) {
    std::cout << "Starting CrossTalk corrections" << std::endl;
    ReadInputs(argc, argv);
    TString file_path = TString::Format("%s/rawData.%d.root", run_dir.c_str(), run);
    cout << "source file: " << file_path << endl;
    TFile *file0 = new TFile(file_path, "READ");
    TTree *tree0 = (TTree*)file0->Get("rawTree");

    if(sil==0 || dia==0)
        GetValuesFromTextFile();
    TString out_name = outfile == "default"? TString::Format("rawData.%d-%05d-%05d.root", run, int(sil*10000), int(dia*10000)) : TString::Format("%s", outfile.c_str());
    TString out_file_path = TString::Format("%s/%s", run_dir.c_str(), out_name.Data());
    cout << "Corrected file will be: " << out_name << endl;
    TFile *file1 = new TFile(out_file_path, "RECREATE");
    TTree *out_tree = new TTree("rawTree2", "rawTree2");
    SetBranches(out_tree);
    LoopSourceTree(tree0, out_tree);
    file0->Close();
    out_tree->SetNameTitle("rawTree", "rawTree");
    out_tree->Write("rawTree");
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
        if((string_arg == "-d") && i+1 < argc){
            i++;
            dia = Float_t(stof(string(argv[i])));
        }
        if((string_arg == "-f") && i+1 < argc){
            i++;
            outfile = string(argv[i]);
        }
    }
    if(run==0){
        cout << "run number was not given. exit" << endl;
        exit(0);
    }
}

void PrintHelp(){
    cout << "Usage:" << endl;
    cout<<"crossTalkCorrection -o RUNDIRECTORY -r RUNNUMBER -s SILALPHA -d DIAALPHA"  << endl;
    cout<<"The values for SILALPHA and DIAALPHA are percentages of charge sharing estimated in the crossTalkCorrectionFactors text file"  << endl;
    cout<<"If one of them is not given, the program will use the text file." << endl;
    cout << "If desired, can specify the name for the new file with -f OUTNAME"<<endl;
    cout << "If not specified, the output will be of the type 'rawData.<RUNNUMBER>0.root'" << endl;
    exit(0);
}

void GetValuesFromTextFile(){
    //TODO
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
    out_tree->Branch("DiaADC", &DiaADC, "DiaADC[128]/s");
    out_tree->Branch("EventNumber", &EventNumber, "EventNumber/i");
}

void LoopSourceTree(TTree *tree0, TTree *out_tree){
    SetBranchAddresses(tree0);
    tree0->SetBranchStatus("*", 1);
    UInt_t nEntries = UInt_t(tree0->GetEntries());
    for (int i=0; i<nEntries; i++){
        ProgressBar(i, nEntries, 100, false, false);
        tree0->GetEntry(i);
        EventNumber = EventNumberin;
        UpdateSilicon();
        UpdateDiamond();
        out_tree->Fill();
    }
}

void SetBranchAddresses(TTree *in_tree){
    if(in_tree->FindBranch("EventNumber")) {
        bEventNumberin = in_tree->GetBranch("EventNumber");
        bEventNumberin->SetAddress(&EventNumberin);
        bEventNumberin->SetAutoDelete(true);
    }
    if(in_tree->FindBranch("D0X_ADC")) {
        bEventNumberin = in_tree->GetBranch("D0X_ADC");
        bEventNumberin->SetAddress(&EventNumberin);
        bEventNumberin->SetAutoDelete(true);
    }
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
        in_tree->SetBranchAddress("DiaADC", &DiaADCin);
}

void ProgressBar(int i, UInt_t nEntries, int refresh, bool show, bool newLine){
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
    UChar_t adcmax = 255;
    UShort_t det_max = 8;
    for (UShort_t det = 0; det < det_max; det++){
        UShort_t startCh = 0;
        UShort_t endCh = 127;
        if(det == 2 || det == 6){
            startCh = 127;
            endCh = 0;
        }
        bool finished = false;
        Double_t alpha = sil / 100.0;
        UChar_t adc = 0;
        for (UShort_t ch = startCh; !finished;){
            UChar_t measured_adc = Det_ADCin[det][ch];
            Double_t real_adc = Double_t(measured_adc) - Double_t(adc) * alpha / (1.0 - alpha);
            UChar_t newADC = real_adc >= adcmax? adcmax : UChar_t(std::floor(real_adc + 0.5));
            Det_ADC[det][ch] = newADC;
            finished = (ch==endCh);
            if(det==2 || det==6)
                ch--;
            else
                ch++;
        }
    }
}

void UpdateDiamond(){
    UShort_t  adcmax = 4095;
    UShort_t adc = 0;
    UChar_t startCh = 0;
    TString out = TString::Format("%d_%5.3f: ", EventNumber, dia);
//    bool bPrint = false;
    Double_t alpha = dia / 100.0;
    for (UChar_t ch=0; ch<128; ch++){
        UShort_t measured_adc = DiaADCin[ch];
        Double_t real_adc = Double_t(measured_adc) - Double_t(adc) * alpha / (1.0 - alpha);
        UShort_t newADC = real_adc >= adcmax? adcmax : UShort_t(std::floor(real_adc + 0.5));
        DiaADC[ch] = newADC;
    }
}