#include "macro_ana.h"
#include "cut/SetCuts.C"

TH1D *his_EAll, *his_E[40], *his_EGood[40];
TH1D *his_Eraw[40], *his_EGoodraw[40];
TH1D *his_EError, *his_weightedEError;
TH1D *his_EfromZAll, *his_EfromZ[40];
TH2D *his_EvsZAll, *his_EvsZ[40], *his_SEvsZ[40];
TH1D *his_Z[40];
TH2D *his_Ecm[40];
TH2D *his_dEdx[40];
TH2D *his_TEvsTA[40], *his_TEvsTL[40], *his_TEvsSA[40], *his_TEvsSL[40];
TH2D *his_SEvsTA[40], *his_SEvsTL[40], *his_SEvsSA[40], *his_SEvsSL[40];
TH2D *his_EdvsSA[40], *his_EdvsSZ[40], *his_EdvsTA[40], *his_EdvsTZ[40]; 
TH2D *his_EdvsTE[40], *his_EdvsSE[40];

void macro_ana(TString type = "14Oap")
{
    cout << minEcm << " - " << maxEcm << ", bin size = " << binSize << " -> nBins = " << nBins << endl;

    beam = type[2];
    reac = type[4];
    if (beam == 'N' || type(type.Length()-1) == 'N') { Is14N = true; cout << "14N RUN" << endl; }

    if (type(3,3) == "CO2")
    {
        IsBackground = true;
        reac = type[6];
        SetCutBG(cuts,reac,type(0,6));
        if (type == "14OCO2p_14N") SetCutBG(cuts,reac,"14NCO2");
        ReadCut(cutEA, "A", (TString)"14"+beam+"a"+reac);
        ReadCut(cutEZ, "Z", (TString)"14"+beam+"a"+reac);
    }
    else
    {
        if (reac == 'p') SetCutP(cuts);
        if (reac == 'a') SetCutA(cuts);
        ReadCut(cutEA, "A", type);
        ReadCut(cutEZ, "Z", type);
        if (type == "14Oap_14N" || type == "14OCO2p_14N") MoveCut(cutEZ,60);
    }
    if (type == "14Oap")            GetGoodEvent("selectmain.txt");
    else if (type == "14Oap_14N")   GetGoodEvent("select14N.txt");
    else if (type == "14OCO2p")     GetGoodEvent("selectbg.txt");
    else if (type == "14OCO2p_14N") GetGoodEvent("select14Nbg.txt");

    SetOthers();
    GetLiveTime();
    GetResolution();

    TFile *fin = new TFile(Form("ana/ana_%s.root", type.Data()));
    TTree *tree = (TTree *)fin->Get("event");

    TClonesArray *EnderArray = nullptr;
    TClonesArray *HeaderArray = nullptr;
    tree->SetBranchAddress("EventEnder", &EnderArray);
    tree->SetBranchAddress("EventHeader", &HeaderArray);
    auto nEvts = tree->GetEntries();

    TH1D *his[3] = { new TH1D("his0","his0",50,0,20),
                     new TH1D("his6","his6",50,0,20),
                     new TH1D("his8","his8",50,0,20) };
    // 8 TCutG *cutg = new TCutG("CUTG",6);
    // 8 cutg->SetPoint(0,0.661544,6.36474);
    // 8 cutg->SetPoint(1,0.761361,4.48102);
    // 8 cutg->SetPoint(2,3.92224,2.63034);
    // 8 cutg->SetPoint(3,4.72078,3.95225);
    // 8 cutg->SetPoint(4,1.26045,6.2656);
    // 8 cutg->SetPoint(5,0.661544,6.36474);
    // 0 TCutG *cutg = new TCutG("CUTG",5);
    // 0 cutg->SetPoint(0,0.561727,8.71114);
    // 0 cutg->SetPoint(1,8.21438,4.6793);
    // 0 cutg->SetPoint(2,8.04802,3.68787);
    // 0 cutg->SetPoint(3,0.46191,6.36474);
    // 0 cutg->SetPoint(4,0.561727,8.71114);

    int maxdEdx = 2000, maxEdet = 20;
    if (reac == 'a' || type(6, 3) == "all" || type(3, 3) == "CO2")
        maxdEdx *= 2, maxEdet *= 2;
    his_EAll = new TH1D("his_EAll","his_EAll",200,0,20);
    his_EfromZAll = new TH1D("his_EfromZAll","his_EfromZAll",200,0,20);
    his_EvsZAll = new TH2D("his_EvsZAll","his_EvsZAll",100,0,600,200,0,20);
    his_EError = new TH1D("his_EError", "his_EError", nBins, minEcm, maxEcm);
    his_weightedEError = new TH1D("his_weightedEError", "his_weightedEError", nBins, minEcm, maxEcm);
    for (int i = 0; i < 40; i++)
    {
        his_E[i] = new TH1D(Form("his_E_%d", i), Form("his_E_%d", i), nBins, minEcm, maxEcm);
        his_EGood[i] = new TH1D(Form("his_EGood_%d", i), Form("his_EGood_%d", i), nBins, minEcm, maxEcm);
        his_Eraw[i] = new TH1D(Form("his_Eraw_%d", i), Form("his_Eraw_%d", i), nBins, minEcm, maxEcm);
        his_EGoodraw[i] = new TH1D(Form("his_EGoodraw_%d", i), Form("his_EGoodraw_%d", i), nBins, minEcm, maxEcm);
        his_EfromZ[i] = new TH1D(Form("his_EfromZ_%d", i), Form("his_EfromZ_%d", i), nBins, minEcm, maxEcm);
        his_Z[i] = new TH1D(Form("his_Z_%d", i), Form("his_Z_%d", i), 600, 0, 600);
        his_Ecm[i] = new TH2D(Form("his_Ecm_%d", i), Form("his_Ecm_%d;TPC Ecm;Si Ecm", detmap[i]), nBins, minEcm, maxEcm, nBins, minEcm, maxEcm);
        his_dEdx[i] = new TH2D(Form("his_dEdx_%d", i), Form("his_dEdx_%d;Edet;dEdx", detmap[i]), 200, 0, maxEdet, 200, 0, maxdEdx);
        his_EvsZ[i] = new TH2D(Form("his_EvsZ_%d", i), Form("his_EvsZ_%d;Z;Ecm", detmap[i]), 100, 0, 600, nBins, minEcm, maxEcm);
        his_SEvsZ[i] = new TH2D(Form("his_SEvsZ_%d", i), Form("his_SEvsZ_%d;Z;Ecm", detmap[i]), 200, 0, 600, nBins, minEcm, maxEcm);

        his_TEvsTA[i] = new TH2D(Form("his_TEvsTA_%d", i), Form("his_TEvsTA_%d;TPC Ang;TPC Ecm", detmap[i]), 180, 0, 180, 200, 0, 20);
        his_TEvsTL[i] = new TH2D(Form("his_TEvsTL_%d", i), Form("his_TEvsTL_%d;TPC Len;TPC Ecm", detmap[i]), 300, 0, 600, 200, 0, 20);
        his_TEvsSA[i] = new TH2D(Form("his_TEvsSA_%d", i), Form("his_TEvsSA_%d;Si Ang ;TPC Ecm", detmap[i]), 180, 0, 180, 200, 0, 20);
        his_TEvsSL[i] = new TH2D(Form("his_TEvsSL_%d", i), Form("his_TEvsSL_%d;Si Len ;TPC Ecm", detmap[i]), 300, 0, 600, 200, 0, 20);
        his_SEvsTA[i] = new TH2D(Form("his_SEvsTA_%d", i), Form("his_SEvsTA_%d;TPC Ang;Si Ecm ", detmap[i]), 180, 0, 180, 200, 0, 20);
        his_SEvsTL[i] = new TH2D(Form("his_SEvsTL_%d", i), Form("his_SEvsTL_%d;TPC Len;Si Ecm ", detmap[i]), 300, 0, 600, 200, 0, 20);
        his_SEvsSA[i] = new TH2D(Form("his_SEvsSA_%d", i), Form("his_SEvsSA_%d;Si Ang;Si Ecm ", detmap[i]), 180, 0, 180, 200, 0, 20);
        his_SEvsSL[i] = new TH2D(Form("his_SEvsSL_%d", i), Form("his_SEvsSL_%d;Si Len;Si Ecm ", detmap[i]), 300, 0, 600, 200, 0, 20);
        his_EdvsTA[i] = new TH2D(Form("his_EdvsTA_%d", i), Form("his_EdvsTA_%d;TPC Ang;Edet", detmap[i]), 180, 0, 180, 100, 0, 30);
        his_EdvsTZ[i] = new TH2D(Form("his_EdvsTZ_%d", i), Form("his_EdvsTZ_%d;TPC Z;Edet", detmap[i]), 100, 0, 600, 100, 0, 30);
        his_EdvsSA[i] = new TH2D(Form("his_EdvsSA_%d", i), Form("his_EdvsSA_%d;Si Ang;Edet", detmap[i]), 180, 0, 180, 100, 0, 30);
        his_EdvsSZ[i] = new TH2D(Form("his_EdvsSZ_%d", i), Form("his_EdvsSZ_%d;Si Z;Edet", detmap[i]), 100, 0, 600, 100, 0, 30);
        his_EdvsTE[i] = new TH2D(Form("his_EdvsTE_%d", i), Form("his_EdvsTE_%d;TPC Ecm;Edet", detmap[i]), 200, 0, 20, 200, 0, 20);
        his_EdvsSE[i] = new TH2D(Form("his_EdvsSE_%d", i), Form("his_EdvsSE_%d;Si Ecm;Edet", detmap[i]), 200, 0, 20, 200, 0, 20);
    }
    TH2D *his_forw[5];
    for (int i = 0; i < 5; i++) his_forw[i] = new TH2D(Form("his_forw%d", drawmap[i]), Form("his_forw%d;CsIE;SiE", drawmap[i]), 100, 0, 20, 100, 0, 20);
    TH2D *his_CsI[5];
    for (int i = 0; i < 5; i++) his_CsI[i] = new TH2D(Form("his_CsI%d", drawmap[i]), Form("his_CsI%d", drawmap[i]), 200, 0, 20, 512, 0, 512);
    TH2D *his_Error[6] = { new TH2D("his_ErrorEcmF","his_ErrorEcmF;error[MeV];Ecm",200,0,0.5,100,0,10),
                           new TH2D("his_ErrorEcmX","his_ErrorEcmX;error[MeV];Ecm",200,0,0.5,100,0,10),
                           new TH2D("his_ErrorEprF","his_ErrorEprF;error[MeV];Ecm",200,0,0.5,100,0,10),
                           new TH2D("his_ErrorEprX","his_ErrorEprX;error[MeV];Ecm",200,0,0.5,100,0,10),
                           new TH2D("his_ErrorTheF","his_ErrorTheF;error[MeV];Ecm",200,0,0.5,100,0,10),
                           new TH2D("his_ErrorTheX","his_ErrorTheX;error[MeV];Ecm",200,0,0.5,100,0,10) };

    for (int iEvt = 0; iEvt < nEvts; iEvt++)
    {
        //PrintProgress(iEvt, nEvts);
        tree->GetEntry(iEvt);
        auto header = (TTEventHeader *)HeaderArray->At(0);
        auto ender = (TTEventEnder *)EnderArray->At(0);

        auto eventNo = header->GetEventNumber();
        auto firedDet = header->GetFiredDet();
        auto firedStrip = header->GetFiredStrip();
        if (!CutBadDetectors(firedDet, firedStrip)) continue;

        auto SiEcm = ender->GettEcm();
        auto TPCEcm = ender->GetpEcm();
        if (reac == 'a') TPCEcm = ender->GetaEcm();
        auto SiAlab = ender->GettAlab();
        auto TPCAlab = ender->GetAlab();
        auto Sivert = TVector3(0, 0, ender->GettZ());
        auto TPCvert = ender->GetVertex();
        auto SiPos = ender->GetSiHit();
        auto SiLen = (SiPos - Sivert).Mag();
        auto TPCLen = (SiPos - TPCvert).Mag();

        auto runNo = ender->GetRun();
        auto splitNo = ender->GetSplit();
        auto LT = liveTime[runNo];
        if (LT < 0.5) continue;

        int iDet = -1;
        for (int i = 0; i < 40; i++)
            if (detmap[i] == firedDet)
                iDet = i;
        if (ender->GetSiE()<SiThres[iDet]) continue;

        auto dEdx = ender->GetdEdx();
        auto Edet = ender->GetEdet();
        if (dEdx>0 && Edet>0) his_dEdx[iDet]->Fill(Edet,dEdx);
        if (TPCEcm == -1) continue;
        if (cuts[iDet] == nullptr) continue;
        if (!cuts[iDet]->IsInside(Edet, dEdx)) continue;

        auto SiE = ender->GetSiE();
        auto CsIE = ender->GetCsIE();
        if(iDet==0) his_forw[iDet]->Fill(CsIE,SiE);
        if(iDet>5 && iDet<10) his_forw[iDet-5]->Fill(CsIE,SiE);
        if (CsIE > 0.5)
        {
            if (iDet==0) his[0]->Fill(SiE+CsIE);
            if (iDet==6) his[1]->Fill(SiE+CsIE);
            if (iDet==8) his[2]->Fill(SiE+CsIE);
        }

        his_EdvsSA[iDet]->Fill(SiAlab, Edet);
        his_EdvsTA[iDet]->Fill(TPCAlab, Edet);
        his_EdvsSZ[iDet]->Fill(Sivert.Z(), Edet);
        his_EdvsTZ[iDet]->Fill(TPCvert.Z(), Edet);
        if (!cutEA[iDet]->IsInside(TPCAlab, Edet)) continue;
        if (!cutEZ[iDet]->IsInside(TPCvert.Z(), Edet)) continue;

        his_Z[iDet]->Fill(TPCvert.Z());
        his_Ecm[iDet]->Fill(TPCEcm, SiEcm);
        his_EvsZ[iDet]->Fill(TPCvert.Z(), TPCEcm);
        his_SEvsZ[iDet]->Fill(TPCvert.Z(), SiEcm);
        his_TEvsTA[iDet]->Fill(TPCAlab, TPCEcm);
        his_TEvsTL[iDet]->Fill(TPCLen, TPCEcm);
        his_TEvsSA[iDet]->Fill(SiAlab, TPCEcm);
        his_TEvsSL[iDet]->Fill(SiLen, TPCEcm);
        his_SEvsTA[iDet]->Fill(TPCAlab, SiEcm);
        his_SEvsTL[iDet]->Fill(TPCLen, SiEcm);
        his_SEvsSA[iDet]->Fill(SiAlab, SiEcm);
        his_SEvsSL[iDet]->Fill(SiLen, SiEcm);
        his_EdvsTE[iDet]->Fill(TPCEcm, Edet);
        his_EdvsSE[iDet]->Fill(SiEcm, Edet);

        //if (TPCEcm<10) cout << ender->GetRun() << "\t" << ender->GetSplit() << "\t" << eventNo << "\t" << TPCEcm << endl;
        
        //auto EcmFill = ZtoE(TPCvert.Z());
        auto EcmFill = TPCEcm;
        auto check = IsGoodEvent(runNo, splitNo, eventNo);
        if (check == 'Z' && EcmFill<=3.5) cout << "?" << "\t" << runNo << "\t" << splitNo << "\t" << eventNo << "\t" << EcmFill << "\t" << TPCvert.Z() << endl;
        if (check != 'A' && check != 'E') continue; // good or short
        //if (EcmFill >= 1.4 && EcmFill < 1.5) cout << runNo << "\t" << splitNo << "\t" << eventNo << "\t" << EcmFill << "\t" << ZtoE(TPCvert.Z()) << "\t" << TPCvert.Z() << "\t" << check << endl;
        //if (check == 'C')
        //{
        //    if (!GetManualEvent(makeKey(runNo, splitNo, eventNo), dEdx, EcmFill, TPCAlab, TPCvert))
        //        continue;
        //}
        
        auto trackEff= fTrackingEff->Eval(Edet);
        if (Edet>19) trackEff = fTrackingEff->Eval(19);
        //if (IsBackground && Edet<3.5) trackEff = 1;

        his_EvsZAll->Fill(TPCvert.Z(), EcmFill);
        his_EAll->Fill(EcmFill, 1 / LT/trackEff);
        his_E[iDet]->Fill(EcmFill, 1 / LT/trackEff);
        his_Eraw[iDet]->Fill(EcmFill);
        if (check == 'A')
        {
            his_EGood[iDet]->Fill(EcmFill, 1 / LT/trackEff);
            his_EGoodraw[iDet]->Fill(EcmFill);
        }

        if (TPCvert.Z() < 350)
        {
            his_EfromZAll->Fill(ZtoE(TPCvert.Z()), 1 / LT);
            his_EfromZ[iDet]->Fill(ZtoE(TPCvert.Z()), 1 / LT);
        }

        //auto errEp = CalculateEpError(Edet, firedDet, firedStrip);
        auto Elab = ender->GetpElab();
        auto straggling = Straggling(Elab, TPCLen);
        auto resolution = CalculateEpError(Edet, firedDet, firedStrip);
        auto errEp = sqrt(pow(resolution,2)+pow(straggling,2));
        auto errTh = CalculateThetaError(TPCAlab, SiPos, ender->GetSiHitError(), TPCvert);
        auto errEcm = CalculateEcmError(Edet, TPCAlab, errEp, errTh);
        //cout << TPCEcm << "\t" << straggling << "\t" << resolution << "\t" << errTh << "\t" << errEcm << endl;
        his_EError->Fill(EcmFill, errEcm*errEcm);
        his_weightedEError->Fill(EcmFill, 1/errEcm/errEcm);
        if (firedDet < 10)
        {
            his_Error[0]->Fill(errEcm, EcmFill);
            his_Error[2]->Fill(errEp , EcmFill);
            his_Error[4]->Fill(errTh , EcmFill);
        }
        else
        {
            his_Error[1]->Fill(errEcm, EcmFill);
            his_Error[3]->Fill(errEp , EcmFill);
            his_Error[5]->Fill(errTh , EcmFill);
        }
    }

    //FindCut();
    //Simulation();

    //DrawCuts();

    TCanvas *cvs_EvsZ = new TCanvas("cvs_EvsZ","cvs_EvsZ",800,600);
    //FindZtoE(his_EvsZAll,10,40);
    cvs_EvsZ->cd();
    his_EvsZAll->Draw("colz");
    if (beam=='O') Ecm_14O->Draw("PC");
    if (beam=='N') Ecm_14N->Draw("PC");
    //g_ZtoE->Draw("P"); f_ZtoE->Draw("same");
    //fZE->Draw("same");
    SaveBatch(cvs_EvsZ);

    //TCanvas *cvsforw = new TCanvas("cvsforw","cvsforw",2400,1600);
    //cvsforw->Divide(3,2);
    //for (int i=0; i<5; i++) { cvsforw->cd(i+1); his_forw[i]->Draw("colz"); }
    //cvsforw->cd(6); his->Draw("hist");
    auto *cvsforw = new TCanvas("cvsforw","cvsforw",1600,2000);
    cvsforw->Divide(2,3);
    cvsforw->cd(1); his_forw[0]->Draw("colz"); cvsforw->cd(2); his[0]->Draw("hist");
    cvsforw->cd(3); his_forw[1]->Draw("colz"); cvsforw->cd(4); his[1]->Draw("hist");
    cvsforw->cd(5); his_forw[3]->Draw("colz"); cvsforw->cd(6); his[2]->Draw("hist");

    TCanvas *cvs_error = new TCanvas("cvs_error","cvs_error",1500,1200);
    cvs_error->Divide(3,2);
    for (int i=0; i<6; i++) { cvs_error->cd(i+1); his_Error[i%3*2 + i/3]->Draw("colz"); }

    //if (type(3,3) == "CO2")
    //{
    //    auto scaleFactor = (double)1/63*30/nBeam14OC*nBeam14Oa;
    //    if (beam == 'N') scaleFactor = (double)1/63*30/nBeam14NC*nBeam14Na;
    //    cout << "Background run. Scale factor: " << scaleFactor << endl;
    //    for (int i=0; i<40; i++) {
    //        his_E[i]->Scale(scaleFactor);
    //        his_EfromZ[i]->Scale(scaleFactor);
    //    }
    //    his_EAll->Scale(scaleFactor);
    //    his_EfromZAll->Scale(scaleFactor);
    //}
    //if (type == "14Oap_14N")
    //{
    //    auto scaleFactor = (double)1/nBeam14Na*nBeam14Oa;
    //    cout << "14N run. Scale factor: " << scaleFactor << endl;
    //    for (int i=0; i<40; i++) {
    //        his_E[i]->Scale(scaleFactor);
    //        his_EfromZ[i]->Scale(scaleFactor);
    //    }
    //    his_EAll->Scale(scaleFactor);
    //    his_EfromZAll->Scale(scaleFactor);
    //}

    TFile *fout = new TFile(Form("results/Yield_%s_%d.root",type.Data(),(int)(binSize*1000)),"recreate");
    fout->cd();
    for (int i=0; i<40; i++)
    {
        his_E[i]->Write(); his_EfromZ[i]->Write();
        his_EGood[i]->Write();
        his_Eraw[i]->Write();
        his_EGoodraw[i]->Write();
    }
    his_EError->Write();
    his_weightedEError->Write();
    his_EAll->Write(); his_EfromZAll->Write();
    fout->Clear();
    fout->Close();
}

void FindZtoE(TH2D *his, int thres_entry, int thres_sigma)
{
    gROOT->SetBatch(true);
    if (g_ZtoE != nullptr) g_ZtoE->Clear();
    if (f_ZtoE != nullptr) f_ZtoE->Clear();

    // 14O kinematics caluclation
    auto calEcm = new TF1("calEcm", "pol6", 0, 10);
    calEcm->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465);

    auto gFit = new TGraph();
    TCanvas *cvs = new TCanvas("cvs","cvs",1000,800);

    double E[nBins], dE[nBins];
    double Z[nBins], dZ[nBins];
    for (int iE = 0; iE<nBins; iE++)
    {
        E[iE] = (double)iE * binSize;
        dE[iE] = (double)binSize / 2;

        auto proj = his->ProjectionX("proj", iE, iE);
        if (proj->GetEntries() > thres_entry)
        {
            auto xMin = calEcm->Eval(E[iE]) - 50;
            auto xMax = calEcm->Eval(E[iE]) + 100;
            auto *gaus = new TF1("gaus", "gaus", xMin, xMax);
            proj->Fit(gaus, "QR", "", xMin, xMax);
            //proj->Draw(); gaus->Draw("same");
            //gROOT->SetBatch(false); cvs->Update(); cvs->WaitPrimitive();
            Z[iE] = gaus->GetParameter(1);
            dZ[iE] = gaus->GetParameter(2);
            if (dZ[iE] < thres_sigma && Z[iE] < 600) gFit->SetPoint(gFit->GetN(), Z[iE], E[iE]);
        }
        else { Z[iE] = -1; dZ[iE] = 0; }
    }
    g_ZtoE = new TGraphErrors(nBins, Z, E, dZ, dE);
    g_ZtoE->SetMarkerColor(kViolet);
    g_ZtoE->SetMarkerStyle(8);

    f_ZtoE = new TF1("f_ZtoE", "pol7", 0, 375);
    f_ZtoE->SetLineColor(kViolet);
    f_ZtoE->SetLineWidth(3);
    f_ZtoE->FixParameter(2, -0.000362325);
    f_ZtoE->FixParameter(3,  4.78713e-06);
    f_ZtoE->FixParameter(4, -3.47837e-08);
    f_ZtoE->FixParameter(5,  1.36886e-10);
    f_ZtoE->FixParameter(6, -2.76809e-13);
    f_ZtoE->FixParameter(7,  2.24699e-16);

    gFit->Fit(f_ZtoE,"RB");

    gROOT->SetBatch(false);
}

void FindCut()
{
    TFile *fhout = new TFile("FindCut.root","recreate");
    fhout->cd();
    for (int i=0; i<40; i++) his_EdvsSZ[i]->Write();
    for (int i=0; i<40; i++) his_EdvsTZ[i]->Write();
    for (int i=0; i<40; i++) his_EdvsSA[i]->Write();
    for (int i=0; i<40; i++) his_EdvsTA[i]->Write();
    fhout->Close();
    return;
}
void Simulation()
{
    TH2D *hAna0 = (TH2D*) his_EvsZ[0]->Clone();
    for (int i=1; i<10; i++) hAna0->Add(his_EvsZ[i]);
    TH2D *hAna1 = (TH2D*) his_EvsZ[10]->Clone();
    for (int i=11; i<17; i++) hAna1->Add(his_EvsZ[i]);
    TH2D *hAna2 = (TH2D*) his_EvsZ[24]->Clone();
    for (int i=25; i<32; i++) hAna2->Add(his_EvsZ[i]);
    TH2D *hAna3 = (TH2D*) his_EvsZ[32]->Clone();
    for (int i=33; i<40; i++) hAna3->Add(his_EvsZ[i]);
    TH2D *hAna4 = (TH2D*) his_EvsZ[0]->Clone();
    for (int i=1; i<40; i++) hAna4->Add(his_EvsZ[i]);

    TFile *fsim = new TFile("saveAna.root","recreate");
    fsim->cd();
    hAna0->SetName("hAna0");hAna0->Write();
    hAna1->SetName("hAna1");hAna1->Write();
    hAna2->SetName("hAna2");hAna2->Write();
    hAna3->SetName("hAna3");hAna3->Write();
    hAna4->SetName("hAna4");hAna4->Write();
    fsim->Close();
    return;
}

bool CutBadDetectors(int det, int strip)
{
    if (det < 10) //forward
    {
        if (det >= 1 && det <= 5) return false;
        if (det == 8)
            if (strip == 1 || strip == 4) return false;
    }
    else
    {
        if (det > 20 && det < 30) return false;
        if (det == 106) return false;
        if (det == 105)
            if (strip >= 5) return false;
        if (det == 208)
            if (strip <= 4) return false;
        if (det == 207 && strip == 6) return false;
    }
    return true;
}

void GetGoodEvent(TString name)
{
    ifstream fin(Form("inputs/%s",name.Data()));
    string line;
    while (getline(fin,line))
    {
        istringstream iss(line);
        string token;
        int runNo, splitNo, eventNo;
        TString flag;
        while (iss >> token)
        {
            if (token == "Event")      iss >> eventNo;
            else if (token == "Run")   iss >> runNo;
            else if (token == "Split") iss >> splitNo;
        }
        flag = line.back();

        goodEvent[makeKey(runNo, splitNo, eventNo)] = flag;
    }
    fin.close();
}
TString IsGoodEvent(int runno, int splitno, int evtno)
{
    uint64_t key = makeKey(runno, splitno, evtno);
    auto it = goodEvent.find(key);
    return (it != goodEvent.end()) ? it->second : 'Z';
}
bool GetManualEvent(uint64_t key, double &dEdx, double &pEcm, double &Alab, TVector3 &vertex)
{
    ifstream fin("inputs/manualEvents.txt");
    int runNo, splitNo, eventNo;
    double cdEdx, cpEcm, cAlab, x, y, z;
    while (fin >> runNo >> splitNo >> eventNo >> cdEdx >> cpEcm >> cAlab >> x >> y >> z)
    {
        if (makeKey(runNo, splitNo, eventNo) == key)
        {
            dEdx = cdEdx;
            pEcm = cpEcm;
            Alab = cAlab;
            vertex = TVector3(x,y,z);
            return true;
        }
    }
    return false;
}

void GetResolution()
{
    ifstream fin("inputs/resolution.txt");
    int det, strip;
    double res;
    while (fin >> det >> strip >> res)
        resolution[det*10+strip] = res;
    fin.close();
}
double CalculateEpError(double energy, int det, int strip)
{
    double error = 0;
    auto res = resolution[det*10+strip];
    error = 5.462 * res / 2.35;
    return error;
}
double CalculateThetaError(double theta, TVector3 sihit, TVector3 errSihit, TVector3 vertex)
{
    double error = 0;
    TVector3 zAxis = TVector3(0, 0, 1);
    auto centerTrack = sihit - vertex;
    auto centerTheta = zAxis.Angle(centerTrack);
    double beamSize = 5;
    auto cos = TMath::Cos(centerTheta);
    auto z = sqrt(cos*cos*beamSize/(1-cos*cos));
    TVector3 errVertex = TVector3(beamSize, beamSize, z);

    double stdv = 0;
    int dirX[8] = { -1, -1, -1, -1, +1, +1, +1, +1 };
    int dirY[8] = { -1, -1, +1, +1, -1, -1, +1, +1 };
    int dirZ[8] = { -1, +1, -1, +1, -1, +1, -1, +1 };
    for (int iSi=0; iSi<8; iSi++)
        for (int iVer=0; iVer<8; iVer++)
        {
            TVector3 sierr = TVector3(errSihit.X()*dirX[iSi], errSihit.Y()*dirY[iSi], errSihit.Z()*dirZ[iSi]);
            TVector3 vererr = TVector3(errVertex.X()*dirX[iVer], errVertex.Y()*dirY[iVer], errVertex.Z()*dirZ[iVer]);
            TVector3 si = sihit + sierr;
            TVector3 ver = vertex + vererr;
            auto track = si - ver;
            auto th = zAxis.Angle(track);
            stdv += pow(th - centerTheta, 2);
        }
    error = sqrt(stdv/64);
    return error;
}
double CalculateEcmError(double energy, double theta, double errEne, double errThe)
{
    double error = 0;
    double cosTh = TMath::Cos(theta * TMath::DegToRad());
    double sinTh = TMath::Sin(theta * TMath::DegToRad());
    double amu = 931.5016;
    double Mbeam = 14.00860, Mtarget = 4.00260, Mlight = 1.00782, Mheavy = 17.00210;
    if (beam == 'N') { Mbeam = 14.00307; Mheavy = 16.99913; }

    if (reac == 'a') // Scattering
    {
        //Ecm = ene * ( Mbeam + Mtarget ) / ( 4 * Mbeam * cos * cos ); TODO
    }
    else // (a,p) reaction
    {
        double m1 = Mbeam   * amu; //MeV 
        double m2 = Mtarget * amu; //MeV 
        double m3 = Mlight  * amu; //MeV 
        double m4 = Mheavy  * amu; //MeV 
        double Q  = m1+m2-m3-m4;   //MeV 

        double f2   = sqrt(cosTh*cosTh * m1*m3*energy - (m4-m1)*(m4*Q - (m3+m4)*energy));
        double f    = -cosTh * sqrt(m1*m3*energy) + f2;
        double dfde = -cosTh * sqrt(m1*m3) / 2 / sqrt(energy) +
                      (cosTh*cosTh * m1 * m3 + (m4-m1) * (m3+m4)) / 2 / f2;
        double dfdt = sinTh * sqrt(m1*m3*energy) -
                      2*cosTh*sinTh * m1*m3*energy / 2 / f2;
        error = 2*f / ( pow(m4-m1,2)*(m1+m2)/m2 ) * sqrt( pow(errEne*dfde,2) + pow(errThe*dfdt,2) );
    }
    return error;
}

void GetLiveTime()
{
    ifstream fin("inputs/livetime.txt");
    int run, CRIBrun;
    double live;
    while (fin >> run >> CRIBrun >> live)
    {
        liveTime[run] = live;
        if (live < 0.5) cout << "bad LT: " << run << endl;
    }
    fin.close();
}

void DrawAllDetector(TH1D **his)
{
    TCanvas *cvs = new TCanvas("cvs_all", "cvs_all", 2000, 1200);
    cvs->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        if (his[drawmap[i]] == nullptr) continue;
        cvs->cd(i + 1);
        his[drawmap[i]]->Draw("hist");
    }
    SaveBatch(cvs);
}
void DrawAllDetector(TH2D **his)
{
    TLine *l = new TLine(0, 0, 20, 20);
    l->SetLineColor(kRed);
    l->SetLineStyle(kDashed);
    TCanvas *cvs = new TCanvas("cvs_all", "cvs_all", 2000, 1200);
    cvs->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        if (his[drawmap[i]] == nullptr) continue;
        cvs->cd(i + 1);
        his[drawmap[i]]->Draw("colz");
        if (strcmp(his[drawmap[i]]->GetName(), Form("his_Ecm_%d", drawmap[i])) == 0)
            l->Draw();
        if (strcmp(his[drawmap[i]]->GetName(), Form("his_EvsZ_%d", drawmap[i])) ==0)
        {
            if (beam=='O') { Ecm_14O->Draw("C"); 
                             if (reac=='a') Ecm_14N_Oaa->Draw("C"); }
            if (beam=='N') { Ecm_14N->Draw("C"); }
        }
        if (strcmp(his[drawmap[i]]->GetName(), Form("his_EdvsTA_%d", drawmap[i])) == 0 ||
            strcmp(his[drawmap[i]]->GetName(), Form("his_EdvsSA_%d", drawmap[i])) == 0)
            if (cutEA[drawmap[i]] != nullptr || cutEA[drawmap[i]]->GetN() > 0)
                cutEA[drawmap[i]]->Draw("same");
        if (strcmp(his[drawmap[i]]->GetName(), Form("his_EdvsTZ_%d", drawmap[i])) == 0 ||
            strcmp(his[drawmap[i]]->GetName(), Form("his_EdvsSZ_%d", drawmap[i])) == 0)
            if (cutEZ[drawmap[i]] != nullptr || cutEZ[drawmap[i]]->GetN() > 0)
                cutEZ[drawmap[i]]->Draw("same");
    }
    SaveBatch(cvs);
}
void DrawAllDetector(TGraphErrors **g)
{
    TCanvas *cvs = new TCanvas("cvs_all", "cvs_all", 2000, 1200);
    cvs->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        if (g[drawmap[i]] == nullptr) continue;
        cvs->cd(i + 1);
        g[drawmap[i]]->Draw("APL");
    }
    SaveBatch(cvs);
}
void DrawOneDetector(int Det, bool DrawAngle = true)
{
    int iDet = -1;
    for (int i = 0; i < 40; i++) if (detmap[i] == Det) iDet = i;
    if (iDet == -1) { cout << "Bad detector number" << endl; return; }

    TCanvas *cvs = new TCanvas(Form("cvs%d", Det), Form("cvs%d", Det), 1800, 1200);
    cvs -> Divide(3,3);
    cvs -> cd(1); his_Ecm[iDet] -> Draw("colz");
    cvs -> cd(2); his_dEdx[iDet] -> Draw("colz"); 
                  if (cuts[iDet]!=nullptr) cuts[iDet] -> Draw("same");
    cvs -> cd(3); his_EvsZ[iDet] -> Draw("colz");
                  FindZtoE(his_EvsZ[iDet], 8, 50);
                  cvs->cd(3);
                  g_ZtoE->Draw("PC");
                  f_ZtoE->Draw("same");
                  if (beam=='O') { Ecm_14O->Draw("PC"); 
                                   if (reac=='a') Ecm_14N_Oaa->Draw("PC"); }
                  if (beam=='N') { Ecm_14N->Draw("PC"); }
    cvs -> cd(4); if (DrawAngle) his_TEvsTA[iDet] -> Draw("colz");
                  else           his_TEvsTL[iDet] -> Draw("colz");
    cvs -> cd(7); if (DrawAngle) his_TEvsSA[iDet] -> Draw("colz");
                  else           his_TEvsSL[iDet] -> Draw("colz");
    cvs -> cd(5); if (DrawAngle) his_SEvsTA[iDet] -> Draw("colz");
                  else           his_SEvsTL[iDet] -> Draw("colz");
    cvs -> cd(8); if (DrawAngle) his_SEvsSA[iDet] -> Draw("colz");
                  else           his_SEvsSL[iDet] -> Draw("colz");
    cvs -> cd(6); if (DrawAngle) his_EdvsSA[iDet] -> Draw("colz");
                  else           his_EdvsSZ[iDet] -> Draw("colz");
    cvs -> cd(9); if (DrawAngle) his_EdvsTA[iDet] -> Draw("colz");
                  else           his_EdvsTZ[iDet] -> Draw("colz");
    SaveBatch(cvs);
}
void DrawCuts()
{
    TCanvas *cvs_dEdx = new TCanvas("cvs_dEdx", "cvs_dEdx", 2000, 1200);
    cvs_dEdx->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        cvs_dEdx->cd(i + 1);
        his_dEdx[drawmap[i]]->Draw("colz");
        if (cuts[drawmap[i]] != nullptr) cuts[drawmap[i]]->Draw("same");
    }
    SaveBatch(cvs_dEdx);

    TCanvas *cvs_EdvsA = new TCanvas("cvs_EdvsA", "cvs_EdvsA", 2000, 1200);
    cvs_EdvsA->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        cvs_EdvsA->cd(i + 1);
        his_EdvsTA[drawmap[i]]->Draw("colz");
        if (cutEA[drawmap[i]]!=nullptr && cutEA[drawmap[i]]->GetN()>0)
            cutEA[drawmap[i]]->Draw("same");
    }
    SaveBatch(cvs_EdvsA);
    TCanvas *cvs_EdvsZ = new TCanvas("cvs_EdvsZ", "cvs_EdvsZ", 2000, 1200);
    cvs_EdvsZ->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if (drawmap[i] == -1) continue;
        cvs_EdvsZ->cd(i + 1);
        his_EdvsTZ[drawmap[i]]->Draw("colz");
        if (cutEZ[drawmap[i]]!=nullptr && cutEZ[drawmap[i]]->GetN()>0)
            cutEZ[drawmap[i]]->Draw("same");
    }
    SaveBatch(cvs_EdvsZ);
}

void SetOthers()
{
    if (!IsBackground) {
        fTrackingEff= new TF1("fTrackingEff","pol7",1,19);
        fTrackingEff->SetParameters(1.20979, -0.218745, -0.00158357, 0.00664389, -0.0010288, 7.20742e-05, -2.46967e-06, 3.35301e-08);
    }
    else {
        fTrackingEff= new TF1("fTrackingEff","pol7",1,19);
        fTrackingEff->SetParameters(0.916231, 0.167327, -0.104535, 0.0214351, -0.00230307, 0.000138251, -4.3827e-06, 5.71362e-08);
    }
    // 5 um scintillator * 45 deg
    double Ecm_Z[15] = { 24, 49, 74, 99, 124, 149, 174, 199, 224, 249, 274, 299, 324, 349, 374 };
    double Ecm_O[15] = { 8.69179683752332, 8.21467248156703, 7.72021431109532, 7.20886678289065, 6.67662978591099, 6.12061435107045, 5.53793150928311, 4.92591451985431, 4.28011881495958, 3.59521091320956, 2.86585733321489, 2.08628013680377, 1.26203503375677, 0.455790430398863, 0.0233339810784401 };
    Ecm_14O = new TGraph(15, Ecm_Z, Ecm_O);
    Ecm_14O->SetMarkerStyle(22);
    Ecm_14O->SetMarkerColor(kRed);
    Ecm_14O->SetLineColor(kRed);
    Ecm_14O->SetLineWidth(2);
    double Ecm_N[15] = { 6.60554475340268, 6.15539405087398, 5.6870150236009, 5.19685092529187, 4.68067811972562, 4.13671823375637, 3.55741318151449, 2.93520487713037, 2.26031226830215, 1.52495497251699, 0.78137270093254, 0.170723821996071, 0, 0, 0 };
    Ecm_14N = new TGraph(15, Ecm_Z, Ecm_N);
    Ecm_14N->SetMarkerStyle(23);
    Ecm_14N->SetMarkerColor(kMagenta);
    Ecm_14N->SetLineColor(kMagenta);
    Ecm_14N->SetLineWidth(2);
    double Ecm_N_Oaa[15] = { 5.51251215867002, 5.01434538120492, 4.48905841326649, 3.93176072870379, 3.3384509879388, 2.6989035453832, 2.00178127223258, 1.25353077114042, 0.526398451154553, 0.0664666963239913, 0, 0, 0, 0, 0 };
    Ecm_14N_Oaa = new TGraph(15, Ecm_Z, Ecm_N_Oaa);
    Ecm_14N_Oaa->SetMarkerStyle(23);
    Ecm_14N_Oaa->SetMarkerColor(kOrange + 1);
    Ecm_14N_Oaa->SetLineColor(kOrange + 1);
    Ecm_14N_Oaa->SetLineWidth(2);

    //fZE = new TF1("fZE","pol7",0,375);
    //fZE->SetLineColor(kViolet+3);
    //fZE->SetParameters(8.98833, -0.00603845, -0.000362325, 4.78713e-06, -3.47837e-08, 1.36886e-10, -2.76809e-13, 2.24699e-16);
    //if (Is14N) fZE->SetParameters(6.44151, 0.129428, -0.00232846, 1.93606e-05, -9.24949e-08, 2.53522e-10, -3.71179e-13, 2.24697e-16);
    //fEZ = new TF1("fEZ","pol7",0,9);
    //fEZ->SetLineColor(kCyan+1);
    //fEZ->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465);
    for (int i=0; i<4; i++)
    {
        if (!Is14N)
        {
            fZE[i] = new TF1(Form("fZE_%d",i), "0.22222839*([0]+[1]*(x-10.413)+[2]*(x-10.413)^2+[3]*(x-10.413)^3+[4]*(x-10.413)^4)", 0, 375);
            fEZ[i] = new TF1(Form("fEZ_%d",i), "[0]+[1]*(4.4998751*x)+[2]*(4.4998751*x)^2+[3]*(4.4998751*x)^3+[4]*(4.4998751*x)^4 - 10.413", 0, 9);
        }
        else
        {
            fZE[i] = new TF1(Form("fZE_%d",i), "0.22229664*([0]+[1]*(x-2.86574)+[2]*(x-2.86574)^2+[3]*(x-2.86574)^3+[4]*(x-2.86574)^4)", 0, 375);
            fEZ[i] = new TF1(Form("fEZ_%d",i), "[0]+[1]*(4.4984935*x)+[2]*(4.4984935*x)^2+[3]*(4.4984935*x)^3+[4]*(4.4984935*x)^4 - 2.86574", 0, 9);
        }
    }
    if (!Is14N)
    {
        fZE[0]->SetParameters(39.9981, -0.0843282, -6.08828e-05, 3.27748e-08, -3.45099e-10);
        fZE[1]->SetParameters(93.5903, -0.861374, 0.00415962, -1.01406e-05, 8.83304e-09);
        fZE[2]->SetParameters(3452.16, -44.5551, 0.217312, -0.000472261, 3.84527e-07);
        fZE[3]->SetParameters(-1142.37, 18.9009, -0.103531, 0.000235856, -1.93395e-07);
        fEZ[0]->SetParameters(365.486, -35.0153, 21.5374, -8.47836, 1.37061);
        fEZ[1]->SetParameters(358.264, -16.9515, 4.05397, -0.733182, 0.0515715);
        fEZ[2]->SetParameters(349.372, -7.6718, 0.274643, -0.0260368, 0.000736422);
        fEZ[3]->SetParameters(344.024, -5.50135, -0.0752176, -4.97182e-05, -1.72572e-07);
    }
    else
    {
        fZE[0]->SetParameters(30, -0.0821361, -5.19621e-05, -1.29567e-07, -4.08549e-11);
        fZE[1]->SetParameters(31.0524, -0.107965, 0.000182799, -1.0634e-06, 1.32497e-09);
        fZE[2]->SetParameters(658.956, -10.6729, 0.0668874, -0.000188385, 1.98762e-07);
        fZE[3]->SetParameters(-5113.94, 74.558, -0.404784, 0.000971192, -8.69766e-07);
        fEZ[0]->SetParameters(298.599, -36.74, 23.9467, -10.3304, 1.81655);
        fEZ[1]->SetParameters(291.262, -17.2617, 3.73715, -0.691054, 0.050257);
        fEZ[2]->SetParameters(283.774, -9.15657, 0.322428, -0.0294115, 0.000814409);
        fEZ[3]->SetParameters(279.257, -7.09869, -0.0334126, -0.00193614, 1.98111e-05);
        OcutZ[0] = 0      ;
        OcutZ[1] = 157.068-2.86574;
        OcutZ[2] = 229.519-2.86574;
        OcutZ[3] = 261.217-2.86574;
        OcutZ[4] = 291    -2.86574;
        OcutE[0] = 0.5    *0.22229664;
        OcutE[1] = 1.59147*0.22229664;
        OcutE[2] = 3.57331*0.22229664;
        OcutE[3] = 10.9892*0.22229664;
        OcutE[4] = 30     *0.22229664;
    }

    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }
}

double EtoZ(double energy) 
{
    if (energy < OcutE[0] || energy >= OcutE[4]) return -1;
    for (int i=0; i<4; i++)
        if (energy < OcutE[i+1]) return fEZ[i]->Eval(energy);
    return -1;
}
double ZtoE(double dist)
{
    if (dist < OcutZ[0] || dist >= OcutZ[4]) return -1;
    for (int i=0; i<4; i++)
        if (dist < OcutZ[i+1]) return fZE[i]->Eval(dist);
    return -1;
}
