#include "../cut/SetCuts.C"
void TrackingEff()
{
    auto minE = 0.;
    auto maxE = 30.;
    auto nBin = 300;
    int detmap[40] = {0};
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }
    double SiThres[40] = {1.615, 100, 100, 100, 100, 100, 1.423, 1.455, 1.767, 1.511, 
                          1.427, 1.405, 1.453, 1.499, 1.451, 1.449, 1.443, 
                          100, 100, 100, 100, 100, 100, 100, 
                          1.511, 1.407, 1.509, 1.533, 1.887, 100, 1.483, 1.495, 
                          1.575, 1.537, 1.465, 1.525, 1.471, 1.423, 1.491, 1.563};
    TCutG *cuts[40];
    SetCutP(cuts);

    TCutG *cutCsI = new TCutG("CUTG", 7);
    //cutCsI->SetPoint(0, 11.1103, 8.48373);
    //cutCsI->SetPoint(1, 5.84527, 3.56509);
    //cutCsI->SetPoint(2, 9.53438, 1.93787);
    //cutCsI->SetPoint(3, 19.0258, 2.0858);
    //cutCsI->SetPoint(4, 19.6705, 4.1568);
    //cutCsI->SetPoint(5, 12.6146, 8.48373);
    //cutCsI->SetPoint(6, 11.1103, 8.48373);
    cutCsI->SetPoint(0,11.1772,8.65157);
    cutCsI->SetPoint(1,7.83429,5.6332);
    cutCsI->SetPoint(2,12.0046,3.271);
    cutCsI->SetPoint(3,24.2509,1.95866);
    cutCsI->SetPoint(4,23.5227,3.36942);
    cutCsI->SetPoint(5,11.5744,7.66732);
    cutCsI->SetPoint(6,11.1772,8.65157);

    //auto *file = new TFile("../ana_14Oap.root");
    //auto *file = new TFile("../ana_14OCO2p.root");
    auto *file = new TFile("../ana/alpha/ana_14Oaa.root");
    auto *tree = (TTree*) file->Get("event");
    auto nEvts = tree->GetEntries();

    TClonesArray *header = nullptr;
    TClonesArray *ender = nullptr;
    TClonesArray *track = nullptr;

    tree->SetBranchAddress("EventHeader",&header);
    tree->SetBranchAddress("EventEnder",&ender);
    tree->SetBranchAddress("Track",&track);

    //TH2D *his[2] = { new TH2D("hisF","hisF",nBin,minE,maxE,2,0,2), 
    //                 new TH2D("hisX","hisX",nBin,minE,maxE,2,0,2) };
    TH2D *his[5] = { new TH2D("his0","his0",nBin,minE,maxE,2,0,2),
                     new TH2D("his6","his6",nBin,minE,maxE,2,0,2),
                     new TH2D("his7","his7",nBin,minE,maxE,2,0,2),
                     new TH2D("his8","his8",nBin,minE,maxE,2,0,2),
                     new TH2D("his9","his9",nBin,minE,maxE,2,0,2) };
    for (int iEvt=0; iEvt<nEvts; iEvt++)
    {
        tree->GetEntry(iEvt);
        if (ender->GetEntries() != 1) continue;
        auto head = (TTEventHeader*) header->At(0);
        auto det = head->GetFiredDet();
        if (det>10) continue;
        int iDet = -1;
        for (int i = 0; i < 40; i++)
            if (detmap[i] == det)
                iDet = i;
        if (cuts[iDet]==nullptr) continue;

        auto end = (TTEventEnder*) ender->At(0);
        auto SiE = end->GetSiE();
        auto CsIE = end->GetCsIE();
        if (SiE<SiThres[iDet]) continue;

        auto Edet = end->GetEdet();
        auto dEdx = end->GetdEdx();
        //if (!cuts[iDet]->IsInside(Edet,dEdx) && dEdx>100) continue;
        //if (!cutCsI->IsInside(Edet,SiE)) continue;
        if (Edet<10 && CsIE != 0) continue;

        auto nTrack = track->GetEntries();
        //if (det<10) his[0]->Fill(Edet,nTrack);
        //else        his[1]->Fill(Edet,nTrack);
        if (det==0) his[0]->Fill(Edet,nTrack);
        else if (det<10) his[det-5]->Fill(Edet,nTrack);
    }
    //TH1D *his_ratio[2] = { new TH1D("his_ratioF","his_ratioF",nBin,minE,maxE),
    //                       new TH1D("his_ratioX","his_ratioX",nBin,minE,maxE) };
    TH1D *his_ratio[5] = { new TH1D("his_ratio0","his_ratio0",nBin,minE,maxE),
                           new TH1D("his_ratio6","his_ratio6",nBin,minE,maxE),
                           new TH1D("his_ratio7","his_ratio7",nBin,minE,maxE),
                           new TH1D("his_ratio8","his_ratio8",nBin,minE,maxE),
                           new TH1D("his_ratio9","his_ratio9",nBin,minE,maxE) };
    for (int i=0; i<5; i++)
        for (int iBin=1; iBin<=nBin; iBin++)
        {
            auto w = his[i]->GetBinContent(iBin,2);
            auto wo = his[i]->GetBinContent(iBin,1);
            auto ratio = w/(w+wo);
            if (w==0 && wo==0) ratio = 0;
            his_ratio[i]->SetBinContent(iBin,ratio);
        }
    auto *cvs = new TCanvas("cvs","cvs",3000,1200);
    cvs->Divide(5,2);
    for (int i=0; i<5; i++)
    {
        cvs->cd(i + 1); his[i]->Draw("colz");
        cvs->cd(i + 6); his_ratio[i]->GetYaxis()->SetRangeUser(0,1); his_ratio[i]->Draw();
    }
}
