#include "../macro_ana.h"
#include "../cut/SetCuts.C"
void macro_simple()
{
    gStyle->SetOptStat(0);
    TFile *fin = new TFile("ana_14Oap.root");
    TTree *tree = (TTree *)fin->Get("event");

    TClonesArray *EnderArray = nullptr;
    TClonesArray *HeaderArray = nullptr;
    tree->SetBranchAddress("EventEnder", &EnderArray);
    tree->SetBranchAddress("EventHeader", &HeaderArray);
    auto nEvts = tree->GetEntries();

    //int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    int iDet = 8;

    TH2D *his_dEdx = new TH2D("his_dEdx","his_dEdx",250,0,25,300,0,3000);

    SetCutP(cuts);
    ReadCut(cutEA,"A","14Oap");
    ReadCut(cutEZ,"Z","14Oap");
    TCutG *cutdEdx = cuts[0];
    TCutG *cutEvsA = cutEA[0];
    TCutG *cutEvsZ = cutEZ[0];

    TCutG *cut10MeV = new TCutG("cut10MeV", 10);
    cut10MeV->SetPoint(0, 8.75363, 510.244);
    cut10MeV->SetPoint(1, 9.64494, 501.108);
    cut10MeV->SetPoint(2, 9.96733, 430.305);
    cut10MeV->SetPoint(3, 9.73976, 334.377);
    cut10MeV->SetPoint(4, 9.05706, 309.253);
    cut10MeV->SetPoint(5, 8.58295, 343.513);
    cut10MeV->SetPoint(6, 8.2416, 386.909);
    cut10MeV->SetPoint(7, 8.2416, 428.021);
    cut10MeV->SetPoint(8, 8.4502, 464.564);
    cut10MeV->SetPoint(9, 8.75363, 510.244);

    TCutG *cutCsI = new TCutG("cutCsI", 6);
    cutCsI->SetPoint(0,7.67192,4.80042);
    cutCsI->SetPoint(1,5.88109,3.01471);
    cutCsI->SetPoint(2,9.71347,2.01681);
    cutCsI->SetPoint(3,14.4413,1.85924);
    cutCsI->SetPoint(4,10.1074,3.64496);
    cutCsI->SetPoint(5,7.67192,4.80042);

    for (int iEvt = 0; iEvt < nEvts; iEvt++)
    {
        //PrintProgress(iEvt,nEvts);
        //if (iEvt%10000==0) cout << iEvt << " / " << nEvts << endl;
        tree->GetEntry(iEvt);
        auto header = (TTEventHeader *)HeaderArray->At(0);
        auto ender = (TTEventEnder *)EnderArray->At(0);

        auto firedDet = header->GetFiredDet();
        if (firedDet != iDet) continue;
        //int iDet = -1;
        //for (int i = 0; i < 40; i++)
        //    if (detmap[i] == firedDet)
        //        iDet = i;
        auto dEdx = ender->GetdEdx();
        auto Edet = ender->GetEdet();
        auto SiE = ender->GetSiE();
        auto Alab = ender->GetAlab();
        auto Z = (ender->GetVertex()).Z();

        if (!cutdEdx->IsInside(Edet, dEdx)) continue;
        //if(!cutEvsA->IsInside(Alab,Edet)) continue;
        //if(!cutEvsA->IsInside(Z,Edet)) continue;
        //if(!cutCsI->IsInside(Edet,SiE)) continue;
        //if (cutg->IsInside(Edet, dEdx))
        {
            his_dEdx->Fill(Edet, dEdx);
            //cout << iEvt << " " << firedDet << " " << header->GetEventNumber() << endl;
        }
    }

    his_dEdx->Draw("colz");
}
