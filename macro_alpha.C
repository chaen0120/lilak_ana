#include "macro_ana.h"

void macro_alpha(TString type = "14Oaa")
{
    gStyle->SetPalette(kRainbow);
    cout << minEcm << " - " << maxEcm << ", bin size = " << binSize << " -> nBins = " << nBins << endl;

    beam = type[2];
    reac = type[4];
    if (beam == 'N') { Is14N = true; cout << "14N RUN" << endl; }

    if (type(3,3) == "CO2")
    {
        IsBackground = true;
        reac = type[6];
        SetCutBG(cuts,reac,type(0,6));
        if (type == "14OCO2a_14N") SetCutBG(cuts,reac,"14NCO2");
        ReadCut(cutEA, "A", (TString)"14"+beam+"a"+reac);
        ReadCut(cutEZ, "Z", (TString)"14"+beam+"a"+reac);
        ReadCutEcm(cutEcm, (TString)"14"+beam+"a"+reac);
    }
    else
    {
        if (reac == 'p') SetCutP(cuts);
        if (reac == 'a') SetCutA(cuts);
        ReadCut(cutEA, "A", type);
        ReadCut(cutEZ, "Z", type);
        ReadCutEcm(cutEcm, type);
        if (type == "14Oaa_14N" || type == "14OCO2a_14N") MoveCut(cutEZ,60);
    }
    if (type == "14Oap")            GetGoodEvent("selectmain.txt");
    else if (type == "14Oap_14N")   GetGoodEvent("select14N.txt");
    else if (type == "14OCO2p")     GetGoodEvent("selectbg.txt");
    else if (type == "14OCO2p_14N") GetGoodEvent("select14Nbg.txt");

    SetOthers();
    GetLiveTime();
    GetResolution();

    TFile *fin = new TFile(Form("ana/alpha/ana_%s.root", type.Data()));
    TTree *tree = (TTree *)fin->Get("event");

    TClonesArray *EnderArray = nullptr;
    TClonesArray *HeaderArray = nullptr;
    tree->SetBranchAddress("EventEnder", &EnderArray);
    tree->SetBranchAddress("EventHeader", &HeaderArray);
    auto nEvts = tree->GetEntries();

    TH1D *his[3] = { new TH1D("his0","his0",50,0,20),
                     new TH1D("his6","his6",50,0,20),
                     new TH1D("his8","his8",50,0,20) };

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
        his_EdvsTZ[i] = new TH2D(Form("his_EdvsTZ_%d", i), Form("his_EdvsTZ_%d;TPC Z;Edet", detmap[i]),   100, 0, 600, 100, 0, 30);
        his_EdvsSA[i] = new TH2D(Form("his_EdvsSA_%d", i), Form("his_EdvsSA_%d;Si Ang;Edet", detmap[i]),  180, 0, 180, 100, 0, 30);
        his_EdvsSZ[i] = new TH2D(Form("his_EdvsSZ_%d", i), Form("his_EdvsSZ_%d;Si Z;Edet", detmap[i]),    100, 0, 600, 100, 0, 30);
        his_EdvsTE[i] = new TH2D(Form("his_EdvsTE_%d", i), Form("his_EdvsTE_%d;TPC Ecm;Edet", detmap[i]), 200, 0,  20, 200, 0, 20);
        his_EdvsSE[i] = new TH2D(Form("his_EdvsSE_%d", i), Form("his_EdvsSE_%d;Si Ecm;Edet", detmap[i]),  200, 0,  20, 200, 0, 20);
        his_EdvsTZEcm[i] = new TH2D(Form("his_EdvsTZEcm_%d", i), Form("his_EdvsTZEcm_%d;Ecm from Z;Edet", detmap[i]), nBins, minEcm, maxEcm, 200, 0, 20);
        his_EdvsSZEcm[i] = new TH2D(Form("his_EdvsSZEcm_%d", i), Form("his_EdvsZEcm_%d;Ecm from Z;Edet", detmap[i]), nBins, minEcm, maxEcm, 200, 0, 20);
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

        auto SiE = ender->GetSiE();
        auto CsIE = ender->GetCsIE();
        auto SiPos = ender->GetSiHit();
        int iDet = -1;
        for (int i = 0; i < 40; i++)
            if (detmap[i] == firedDet)
                iDet = i;
        if (SiE<SiThres[iDet]) continue;

        auto runNo = ender->GetRun();
        auto splitNo = ender->GetSplit();
        auto LT = liveTime[runNo];
        if (LT < 0.5) continue;

        auto SiEcm = ender->GettEcm();
        auto SiAlab = ender->GettAlab();
        auto Sivert = TVector3(0, 0, ender->GettZ());
        auto SiLen = (SiPos - Sivert).Mag();
        auto TPCEcm = ender->GetpEcm();
        if (reac == 'a') TPCEcm = ender->GetaEcm();
        auto TPCAlab = ender->GetAlab();
        auto TPCvert = ender->GetVertex();
        auto TPCLen = (SiPos - TPCvert).Mag();

        auto dEdx = ender->GetdEdx();
        auto Edet = ender->GetEdet();
        if (dEdx>0 && Edet>0) his_dEdx[iDet]->Fill(Edet,dEdx);
        if (TPCEcm == -1) continue;
        if (cuts[iDet] == nullptr) continue;
        if (!cuts[iDet]->IsInside(Edet, dEdx)) continue;

        if(iDet==0) his_forw[iDet]->Fill(CsIE,SiE);
        if(iDet>5 && iDet<10) his_forw[iDet-5]->Fill(CsIE,SiE);
        if (CsIE > 0.5)
        {
            if (iDet==0) his[0]->Fill(SiE+CsIE);
            if (iDet==6) his[1]->Fill(SiE+CsIE);
            if (iDet==8) his[2]->Fill(SiE+CsIE);
        }

        his_EdvsTZEcm[iDet]->Fill(ZtoE(TPCvert.Z()), Edet);
        his_EdvsSZEcm[iDet]->Fill(ZtoE(Sivert.Z()), Edet);
        if (!cutEcm[iDet]->IsInside(ZtoE(TPCvert.Z()), Edet)) continue;
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
        
        auto EcmFill = TPCEcm;
        if (IsVertexEcm)
        {
            EcmFill = ZtoE(TPCvert.Z());
            if (type(type.Length()-1) == 'N') EcmFill = ZtoE(TPCvert.Z() - 60);
        }
        //auto check = IsGoodEvent(runNo, splitNo, eventNo);
        TString check = "A";
        if (check == 'Z' && EcmFill<=3.5) cout << "?" << "\t" << runNo << "\t" << splitNo << "\t" << eventNo << "\t" << EcmFill << "\t" << TPCvert.Z() << endl;
        if (check != 'A' && check != 'E') continue; // good or short
        
        //auto trackEff= fTrackingEff->Eval(Edet);
        //if (Edet>19) trackEff = fTrackingEff->Eval(19);
        double trackEff = 1.0;

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
        //XXX
        auto Elab = ender->GetaElab();
        auto straggling = PStraggling(Elab, TPCLen);
        auto resolution = CalculateEpError(Edet, firedDet, firedStrip);
        auto errEp = sqrt(pow(resolution,2)+pow(straggling,2));
        auto errTh = CalculateThetaError(TPCAlab, SiPos, ender->GetSiHitError(), TPCvert);
        auto errEcm = CalculateEcmError(Edet, TPCAlab, errEp, errTh);
        if (IsVertexEcm) errEcm = OStraggling(TPCvert.Z());
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
    auto *foutfile = new TFile("EpVSEcmZ_nocut.root","recreate");
    for (int i=0; i<40; i++) his_EdvsTZEcm[i]->Write();
    for (int i=0; i<40; i++) his_EdvsSZEcm[i]->Write();
    foutfile->Close();

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