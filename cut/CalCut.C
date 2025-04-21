#include "../cut/SetCuts.C"
auto *ana = new TTAnalysisTask();
TF1 *fZE, *fEZ[3];
int iProton = 0;
int iAlpha = 1;
int i14O = 2;
int i14N = 3;
double Mass[4] = {1.00782, 4.00260, 14.00860, 14.00307};
void ThickTargetAnalysis(int mode, TVector3 SiHit, double Edet, int kLight, int kBeam, double &th, double &Z);

void CalCut()
{
    auto *run = new LKRun();
    auto *tt = new TexAT2();
    run->AddPar("/home/cens-alpha-00/lilak/texat_ana/temp_reco/config/config_common.mac");
    run->AddPar("/home/cens-alpha-00/lilak/texat_ana/temp_reco/config/config_reco.mac");
    run->AddPar("/home/cens-alpha-00/lilak/texat_ana/temp_reco/config/config_ana.mac");
    run->AddDetector(tt);
    run->Add(ana);
    run->Init();

    fZE = new TF1("fZE","pol7",0,375);
    fZE->SetParameters(8.98833, -0.00603845, -0.000362325, 4.78713e-06, -3.47837e-08, 1.36886e-10, -2.76809e-13, 2.24699e-16);
    fEZ[0] = new TF1("fEZcenter","pol7(0)+pol1(8)",0,10);
    fEZ[1] = new TF1("fEZmin",   "pol7(0)+pol1(8)",0,10); 
    fEZ[2] = new TF1("fEZmax",   "pol7(0)+pol1(8)",0,10);
    fEZ[0]->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465, 0, 0);
    fEZ[1]->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465, -33.3987, 0.558451);
    fEZ[2]->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465, 28.8964, 0.645681);

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
    TCutG *cuts[40]; SetCutP(cuts);

    auto *fin = new TFile("ana/p0/ana_14Oap.root");
    auto *tree = (TTree*) fin->Get("event");
    auto nEvts = tree->GetEntries();

    TClonesArray *header = nullptr;
    TClonesArray *ender = nullptr;
    TClonesArray *track = nullptr;

    tree->SetBranchAddress("EventHeader",&header);
    tree->SetBranchAddress("EventEnder",&ender);
    tree->SetBranchAddress("Track",&track);

    TH2D *his_EdvsSA[40], *his_EdvsSZ[40], *his_EdvsTA[40], *his_EdvsTZ[40]; 
    for (int i=0; i<40; i++)
    {
        his_EdvsTA[i] = new TH2D(Form("his_EdvsTA_%d", i), Form("his_EdvsTA_%d;TPC Ang;Edet", detmap[i]), 180, 0, 180, 30, 0, 30);
        his_EdvsTZ[i] = new TH2D(Form("his_EdvsTZ_%d", i), Form("his_EdvsTZ_%d;TPC Z;Edet", detmap[i]), 300, 0, 600, 300, 0, 30);
        his_EdvsSA[i] = new TH2D(Form("his_EdvsSA_%d", i), Form("his_EdvsSA_%d;Si Ang;Edet", detmap[i]), 180, 0, 180, 300, 0, 30);
        his_EdvsSZ[i] = new TH2D(Form("his_EdvsSZ_%d", i), Form("his_EdvsSZ_%d;Si Z;Edet", detmap[i]), 300, 0, 600, 300, 0, 30);
    }
    for (int iEvt = 0; iEvt < nEvts; iEvt++)
    {
        tree->GetEntry(iEvt);
        auto *head = (TTEventHeader*) header->At(0);
        auto det = head->GetFiredDet();
        int iDet = -1;
        for (int i=0; i<40; i++) if (detmap[i] == det) iDet = i;
        if (iDet == -1) continue;

        auto *end = (TTEventEnder*) ender->At(0);
        auto dEdx = end->GetdEdx();
        auto Edet = end->GetSiE();
        if (Edet < SiThres[iDet]) continue;
        if (cuts[iDet] != nullptr && !cuts[iDet]->IsInside(Edet, dEdx)) continue;
        auto SiHit = end->GetSiHit();
        double SiA = -1, SiZ = -1;
        ThickTargetAnalysis(0, SiHit, Edet, iProton, i14O, SiA, SiZ);
        his_EdvsSA[iDet]->Fill(SiA, Edet);
        his_EdvsSZ[iDet]->Fill(SiZ, Edet);
        ThickTargetAnalysis(1, SiHit, Edet, iProton, i14O, SiA, SiZ);
        his_EdvsSA[iDet]->Fill(SiA, Edet);
        his_EdvsSZ[iDet]->Fill(SiZ, Edet);
        ThickTargetAnalysis(2, SiHit, Edet, iProton, i14O, SiA, SiZ);
        his_EdvsSA[iDet]->Fill(SiA, Edet);
        his_EdvsSZ[iDet]->Fill(SiZ, Edet);

        auto TPCAlab = end->GetAlab();
        auto TPCvert = end->GetVertex();
        his_EdvsTA[iDet]->Fill(TPCAlab, Edet);
        his_EdvsTZ[iDet]->Fill(TPCvert.Z(), Edet);
    }

    gStyle->SetPalette(kBird);;
    TCanvas *cvs[2] = { new TCanvas("cvsA","cvsA",3000,2000), new TCanvas("cvsZ","cvsZ",3000,2000) };
    for (int iCvs=0; iCvs<2; iCvs++)
    {
        cvs[iCvs]->Divide(8, 5);
        for (int i = 0; i < 40; i++)
        {
            cvs[iCvs]->cd(i + 1);
            if (iCvs == 0) { his_EdvsSA[i]->Draw("colz"); }
            else           { his_EdvsSZ[i]->Draw("colz"); }
        }
    }

    TFile *fout = new TFile("FindCut_14Oap.root","recreate");
    fout->cd();
    for (int i = 0; i < 40; i++)
    {
        his_EdvsSZ[i]->Write();
        his_EdvsTZ[i]->Write();
        his_EdvsSA[i]->Write();
        his_EdvsTA[i]->Write();
    }
    fout->Close();
}

void ThickTargetAnalysis(int mode, TVector3 SiHit, double Edet, int kLight, int kBeam, double &th, double &Z)
{
    TVector3 BeamAtBM(0,0,24);
    TVector3 BeamDir(0,0,1);
    double IterLen = 327.68; // 2^15/100
    double BeamLen = IterLen/2;
    TVector3 Vertex(-1,-1,-1);
    TVector3 LightPtl(-1,-1,-1);
    double angle   = -1;
    double newElab = -1;
    double newEcm  = -1;
    int iter = 0;
    do
    {
        iter++;
        IterLen /= 2;
        Vertex.SetXYZ(BeamAtBM.X() + BeamDir.X() * BeamLen, 
                      BeamAtBM.Y() + BeamDir.Y() * BeamLen, 
                      BeamAtBM.Z() + BeamDir.Z() * BeamLen);
        LightPtl = SiHit - Vertex;
        angle = TMath::ACos(LightPtl.Dot(BeamDir) / LightPtl.Mag());

        newElab       = ana->RevertEnergyLoss(Edet, LightPtl.Mag(), kLight);
        newEcm        = ana->LabToCM(newElab, angle, kLight, kBeam);
        auto newZ = fEZ[mode]->Eval(newEcm);

        if (newZ > BeamLen + 24) BeamLen += IterLen;
        else                     BeamLen -= IterLen;
    }
    while(IterLen>0.01);

    th = angle * TMath::RadToDeg();
    Z  = BeamLen + 24;
    return;
}
