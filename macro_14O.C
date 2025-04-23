#include "macro_ana.C"
#include "totxs/draw_cs.C"

TH1D *hisMain[40], *hisBG[40], *hisN[40], *hisNB[40], *hisGood[40];
TH1D *hisMainraw[40], *hisBGraw[40], *hisNraw[40], *hisNBraw[40], *hisGoodraw[40];
TH1D *hisMainError, *hisMainwEError;
TGraphErrors *g_xs[40], *g_xsGood[40], *g_sysError[40];

void macro_14O(bool isGoodEvents = false, bool isSubtractBG = true)
{
    // Init
    cout << "isGoodEvents = " << isGoodEvents << ", isSubtractBG = " << isSubtractBG << endl;

    beam = 'O';
    reac = 'p';
    nBeam = nBeam14Oa;

    SetOthers();
    GetEfficiency(effMain, nEcmMain, "proton");
    GetEfficiency(effBG, nEcmBG, "CO2_proton");
    GetEfficiency(effN, nEcmN, "Nproton");
    GetEfficiency(effNB, nEcmNB, "Nproton"); //XXX

    // Init histograms
    TH1D *hisMainAll     = new TH1D("hisMainAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisMainAllGood = new TH1D("hisMainAllGood", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisBGAll       = new TH1D("hisBGAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisNAll        = new TH1D("hisNAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisNBAll       = new TH1D("hisNBAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisAll         = new TH1D("hisAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);

    auto fnameMain = Form("results/Yield_14Oap_%d.root",(int)(binSize*1000));
    cout << "Main: " << fnameMain << endl;
    auto finMain = new TFile(fnameMain);
    hisMainError = (TH1D*) finMain->Get("his_EError");
    hisMainwEError = (TH1D*) finMain->Get("his_weightedEError");
    for (int i=0; i<40; i++)
    {
        hisMain[i] = (TH1D*) finMain->Get(Form("his_E_%d",i));
        hisMainraw[i] = (TH1D*) finMain->Get(Form("his_Eraw_%d",i));
        hisMainAll->Add(hisMain[i]);
        
        hisGood[i] = (TH1D*) finMain->Get(Form("his_EGood_%d",i));
        hisGoodraw[i] = (TH1D*) finMain->Get(Form("his_EGoodraw_%d",i));
        hisMainAllGood->Add(hisGood[i]);
    }

    auto fnameBG = Form("results/Yield_14OCO2p_%d.root",(int)(binSize*1000));
    cout << "CO2 : " << fnameBG << endl;
    auto finBG = new TFile(fnameBG);
    for (int i=0; i<40; i++)
    {
        hisBG[i] = (TH1D*) finBG->Get(Form("his_E_%d",i));
        hisBGraw[i] = (TH1D*) finBG->Get(Form("his_Eraw_%d",i));
        hisBGAll->Add(hisBG[i]);
    }

    auto fnameN = Form("results/Yield_14Oap_14N_%d.root",(int)(binSize*1000));
    cout << "14N : " << fnameN << endl;
    auto finN = new TFile(fnameN);
    for (int i=0; i<40; i++)
    {
        hisN[i] = (TH1D*) finN->Get(Form("his_E_%d",i));
        hisNraw[i] = (TH1D*) finN->Get(Form("his_Eraw_%d",i));
        hisNAll->Add(hisN[i]);
    }

    auto fnameNB = Form("results/Yield_14OCO2p_14N_%d.root",(int)(binSize*1000));
    cout << "14NB : " << fnameNB << endl;
    auto finNB = new TFile(fnameNB);
    for (int i=0; i<40; i++)
    {
        hisNB[i] = (TH1D*) finNB->Get(Form("his_E_%d",i));
        hisNBraw[i] = (TH1D*) finNB->Get(Form("his_Eraw_%d",i));
        hisNBAll->Add(hisNB[i]);
    }

    hisBGAll->Scale((double)(1./nBeam14OC*nBeam14Oa/2.1));
    hisNAll ->Scale((double)(1./nBeam14Na*nBeam14Oa));
    hisNBAll->Scale((double)(1./nBeam14NC*nBeam14Oa/2.1));
    hisAll = (TH1D *)hisMainAll->Clone();
    hisAll->Add(hisBGAll, -1);
    hisAll->Add(hisNAll, -1);
    hisAll->Add(hisNBAll, 1);
    
    // Get XS
    for (int iDet=0; iDet<40; iDet++)
    {
        g_xs[iDet] = GetCrossSection(false, iDet);
        g_xsGood[iDet] = GetCrossSection(true, iDet);
    }
    auto g_xsAll = GetAverage(g_xs, hisMainAll, hisMainError);
    g_xsAll->SetTitle("Total Cross Section of 14Oap;E_{cm} [MeV];XS [mb]");
    auto g_weightedxsAll = GetWeightedAverage(g_xs, hisMainwEError);
    g_weightedxsAll->SetTitle("Total Cross Section (weighted average) of 14Oap;E_{cm} [MeV];XS [mb]");

    auto g_xsGoodAll = GetAverage(g_xsGood, hisMainAll, hisMainError);
    auto g_sysErr = new TGraphErrors();
    for (int iPoint = 0; iPoint<nBins; iPoint++)
    {
        auto minAll = g_xsAll->GetPointY(iPoint) - g_xsAll->GetErrorY(iPoint);
        auto minGood = g_xsGoodAll->GetPointY(iPoint) - g_xsGoodAll->GetErrorY(iPoint);
        g_sysErr->SetPoint(iPoint, g_xsAll->GetPointX(iPoint), (minAll+minGood)/2);
        g_sysErr->SetPointError(iPoint, 0, (minAll-minGood)/2);
    }

    gStyle->SetOptStat(0);
    TCanvas *cvs_xs = new TCanvas("cvs_xs","cvs_xs",1500,1000);
    cvs_xs->Divide(2,2);
    cvs_xs->cd(4); g_weightedxsAll->Draw("APL");
    cvs_xs->cd(3); g_xsAll->Draw("APL"); g_sysErr->Draw("P");
    cvs_xs->cd(2); hisAll->Draw("HIST");
    cvs_xs->cd(1);
    hisMainAll->GetXaxis()->SetRangeUser(0,3.5);
    hisMainAll->SetLineColor(kBlack);   hisMainAll->Draw("HIST");
    hisBGAll  ->SetLineColor(kRed);     hisBGAll  ->Draw("HIST SAME");
    hisNAll   ->SetLineColor(kBlue);    hisNAll   ->Draw("HIST SAME");
    hisNBAll  ->SetLineColor(kGreen+1); hisNBAll  ->Draw("HIST SAME");
    TH1D *hisBAll = (TH1D*) hisBGAll->Clone();
    hisBAll->Add(hisNAll);
    hisBAll->Add(hisNBAll, -1);
    hisBAll->SetLineColor(kMagenta); hisBAll->Draw("HIST SAME");

    auto *leg = new TLegend(0.1,0.6,0.4,0.9);
    leg->AddEntry(hisMainAll, "14Oap", "l");
    leg->AddEntry(hisBGAll, "14OCO2", "l");
    leg->AddEntry(hisNAll, "14Nap", "l");
    leg->AddEntry(hisNBAll, "14NCO2", "l");
    leg->AddEntry(hisBAll, "14OCO2+14Nap-14NCO2", "l");
    leg->Draw();
    SaveBatch(cvs_xs);

    TCanvas *cvs_yield = new TCanvas("cvs_yield","cvs_yield",800,1200);
    cvs_yield->Divide(1,2);
    cvs_yield->cd(1);
    hisMainAll->Draw("HIST");
    hisBAll->Draw("HIST same");
    hisAll->SetFillColor(9); hisAll->SetFillStyle(3003); hisAll->Draw("HIST same");
    auto *leg2 = new TLegend(0.15,0.7,0.3,0.88);
    leg2->AddEntry(hisMainAll, "^{14}O(#alpha,p)", "l");
    leg2->AddEntry(hisBAll, "Background", "l");
    leg2->AddEntry(hisAll, "Net", "l");
    leg2->Draw();
    cvs_yield->cd(2);
    hisBAll->GetXaxis()->SetRangeUser(0,3.5);
    hisBAll->Draw("HIST");
    hisBGAll->SetFillColor(kRed); hisBGAll->SetFillStyle(3004); hisBGAll->Draw("HIST same");
    hisNAll ->SetFillColor(kBlue); hisNAll ->SetFillStyle(3005); hisNAll ->Draw("HIST same");
    hisNBAll->SetFillColor(kGreen+1); hisNBAll->SetFillStyle(3006); hisNBAll->Draw("HIST same");
    auto *leg3 = new TLegend(0.15,0.6,0.3,0.88);
    leg3->AddEntry(hisBAll, "Background", "l");
    leg3->AddEntry(hisBGAll, "^{14}O+CO_{2}", "l");
    leg3->AddEntry(hisNAll, "^{14}N(#alpha,p)", "l");
    leg3->AddEntry(hisNBAll, "^{14}N+CO_{2}", "l");
    leg3->Draw();

    ofstream outtxt("results/totxs_14Oap.txt");
    for (int i=0; i<g_xsAll->GetN(); i++)
        outtxt << g_xsAll->GetPointX(i) << "\t" << g_xsAll->GetPointY(i) << "\t"
               << g_xsAll->GetErrorX(i) << "\t" << g_xsAll->GetErrorY(i) << endl;
    outtxt.close();

    auto g = g_xsAll;
    outtxt.open("totxs/totxs_14Oap.txt");
    for (int i=0; i<g->GetN(); i++)
        outtxt << g->GetPointX(i) << "\t" << g->GetPointY(i) << "\t"
            << g->GetErrorX(i) << "\t" << g->GetErrorY(i) << endl;
    outtxt.close();
    outtxt.open("totxs/totxs_14Oap_sysErr.txt");
    for (int i=0; i<g_sysErr->GetN(); i++)
        outtxt << g_sysErr->GetPointX(i) << "\t" << g_sysErr->GetPointY(i) << "\t"
            << g_sysErr->GetErrorX(i) << "\t" << g_sysErr->GetErrorY(i) << endl;
    outtxt.close();
    draw_cs();
}

void GetEfficiency(vector<vector<double>> &eff, vector<double> &nEcm, TString fileName)
{
    for (int i = 0; i < 40; i++) eff[i].resize(nBins);
    nEcm.resize(nBins);
    TString effname = Form("%s_%d.root",fileName.Data(),(int)(binSize*1000));
    effname = "eff/" + effname;
    cout << "eff: " << effname << endl;
    TFile *fhis = new TFile(effname);
    auto his_eff = (TH2D*) fhis->Get("his_Eff");
    auto his_det = (TH2D*) fhis->Get("his_Det");
    auto his_nEcm = (TH1D*) fhis->Get("his_nEcm");

    auto thres = 5;
    for (int iE = 1; iE <= nBins; iE++)
    {
        nEcm[iE-1] = his_nEcm->GetBinContent(iE);
        for (int iDet = 1; iDet <= 40; iDet++)
            if (his_det->GetBinContent(iDet,iE)>thres)
                eff[iDet-1][iE-1] = his_eff->GetBinContent(iDet,iE);
    }
}

TH1D *ApplyEfficiency(TH1D *his, int iDet, bool isBackground)
{
    auto hisEff = (TH1D*) his->Clone();
    hisEff->Reset("ICES");
    hisEff->SetName(Form("eff_%s", his->GetName()));
    hisEff->SetTitle(Form("eff_%s", his->GetName()));
    for (int iE = 1; iE <= nBins; iE++)
    {
        auto yield = his->GetBinContent(iE);
        double solid = 0;
        if (iDet == 40) for (int i=0; i<40; i++) solid += effMain[i][iE-1];
        else solid = effMain[iDet][iE-1];
        if (isBackground)
        {
            if (iDet == 40) for (int i=0; i<40; i++) solid += effBG[i][iE-1];
            else solid = effBG[iDet][iE-1];
        }
        if (solid>0) hisEff->SetBinContent(iE, yield/solid);
        else hisEff->SetBinContent(iE,0);
    }
    return hisEff;
}

//TGraphErrors *GetCrossSection(TH1D *hisMain, int iDet)
//{
//    auto hisBG = (TH1D*) hisMain->Clone();
//    hisBG->Reset("ICES");
//    return GetCrossSection(hisMain, hisBG, iDet);
//}
//TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, int iDet)
//{
//    auto hisN = (TH1D*) hisMain->Clone();
//    hisN->Reset("ICES");
//    auto hisNBG = (TH1D*) hisN->Clone();
//    return GetCrossSection(hisMain, hisBG, hisN, hisNBG, iDet);
//}
//TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, TH1D *hisN, TH1D *hisNBG, int iDet)
TGraphErrors *GetCrossSection(bool IsGoodAna, int iDet)
{
    TH1D *hMain = hisMain[iDet];
    if (IsGoodAna) hMain = hisGood[iDet];
    TH1D *hBG   = hisBG  [iDet];
    TH1D *hN    = hisN   [iDet];
    TH1D *hNB   = hisNB  [iDet];
    TH1D *rMain = hisMainraw[iDet];
    if (IsGoodAna) rMain = hisGoodraw[iDet];
    TH1D *rBG   = hisBGraw  [iDet];
    TH1D *rN    = hisNraw   [iDet];
    TH1D *rNB   = hisNBraw  [iDet];
    TGraphErrors *gXS = new TGraphErrors();
    gXS->SetName(Form("xs_%d", iDet));
    gXS->SetTitle(Form("xs_%d", iDet));

    //hBG->Reset("ICES"); rBG->Reset("ICES");
    //hN ->Reset("ICES"); rN ->Reset("ICES");
    //hNB->Reset("ICES"); rNB->Reset("ICES");

    if (hMain->GetEntries() == 0) 
    {
        for (int iE = 1; iE <= nBins; iE++)
        {
            gXS->SetPoint(gXS->GetN(), hMain->GetBinCenter(iE), 0);
            gXS->SetPointError(gXS->GetN()-1, -1, -1);
        }
        return gXS;
    }

        auto *g = new TGraphErrors();
    for (int iE = 1; iE <= nBins; iE++)
    {
        // eff. target thickness
        double iEmin = hMain->GetBinCenter(iE) - binSize / 2;
        double iEmax = hMain->GetBinCenter(iE) + binSize / 2;
        double len = (EtoZ(iEmin) - EtoZ(iEmax))/10.; //[cm]
        double target = len * gasDensity / HeMolarMass * avo;

        if (hMain->GetBinCenter(iE)<9)
        {
            if (len<0.001 || len>2) continue;
            g->SetPoint(g->GetN(), hMain->GetBinCenter(iE), len);
            g->SetPointError(g->GetN() - 1, binSize / 2, 0);
        }

        // yield
        double yieldMain = hMain->GetBinContent(iE);
        double yieldBG   = hBG  ->GetBinContent(iE);
        double yieldN    = hN   ->GetBinContent(iE);
        double yieldNB   = hNB  ->GetBinContent(iE);

        double rawMain = rMain->GetBinContent(iE);
        double rawBG   = rBG  ->GetBinContent(iE);
        double rawN    = rN   ->GetBinContent(iE);
        double rawNB   = rNB  ->GetBinContent(iE);

        //XXX
        double scaleMain = yieldMain / rawMain;
        double scaleBG   = yieldBG   / rawBG  ;
        double scaleN    = yieldN    / rawN   ;
        double scaleNB   = yieldNB   / rawNB  ;
        if (isnan(scaleMain)) scaleMain = 0;
        if (isnan(scaleBG  )) scaleBG   = 0;
        if (isnan(scaleN   )) scaleN    = 0;
        if (isnan(scaleNB  )) scaleNB   = 0;

        // xs
        double xsMain = yieldMain / nBeam14Oa / target     / effMain[iDet][iE-1]* 1E+27; //[mb]
        double xsBG   = yieldBG   / nBeam14OC / target/2.1 / effBG  [iDet][iE-1]* 1E+27; //[mb]
        double xsN    = yieldN    / nBeam14Na / target     / effMain[iDet][iE-1]* 1E+27; //[mb]
        double xsNB   = yieldNB   / nBeam14NC / target/2.1 / effBG  [iDet][iE-1]* 1E+27; //[mb]
        double xs = xsMain - xsBG - xsN + xsNB;

        // stat. error
        double   yErrMain = sqrt(rawMain) * scaleMain;
        double effErrMain = effMain[iDet][iE-1] * sqrt(1/nEcmMain[iE-1] + 1/nEcmMain[iE-1]/effMain[iDet][iE-1]);
        double    ErrMain = xsMain * sqrt(pow(yErrMain/yieldMain, 2) + pow(effErrMain/effMain[iDet][iE-1], 2));
        double   yErrBG = sqrt(rawBG) * scaleBG;
        double effErrBG = effBG[iDet][iE-1] * sqrt(1/nEcmBG[iE-1] + 1/nEcmBG[iE-1]/effBG[iDet][iE-1]);
        double    ErrBG = xsBG * sqrt(pow(yErrBG/yieldBG, 2) + pow(effErrBG/effBG[iDet][iE-1], 2));
        double   yErrN = sqrt(rawN) * scaleN;
        double effErrN = effMain[iDet][iE-1] * sqrt(1/nEcmMain[iE-1] + 1/nEcmMain[iE-1]/effMain[iDet][iE-1]);
        double    ErrN = xsN * sqrt(pow(yErrN/yieldN, 2) + pow(effErrMain/effMain[iDet][iE-1], 2));
        double   yErrNB = sqrt(rawNB) * scaleNB;
        double effErrNB = effBG[iDet][iE-1] * sqrt(1/nEcmBG[iE-1] + 1/nEcmBG[iE-1]/effBG[iDet][iE-1]);
        double    ErrNB = xsNB * sqrt(pow(yErrNB/yieldNB, 2) + pow(effErrBG/effBG[iDet][iE-1], 2));

        double statErr = sqrt(pow(ErrMain,2) + pow(ErrBG,2) + pow(ErrN,2) + pow(ErrNB,2));

        double ecmerror = binSize / 2; // will be updated at averaging functions

        if (len <= 0 || isnan(xs) || isinf(xs)) xs = 0;
        if (isnan(statErr)) statErr = 0;
        gXS->SetPoint(gXS->GetN(), hMain->GetBinCenter(iE), xs);
        gXS->SetPointError(gXS->GetN()-1, ecmerror, statErr);
    }
    if (iDet==0) {
        auto *cvs = new TCanvas("cvs","cvs",1500,1200);
        cvs->cd();
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);
        g->Draw("AP");
        //cvs->SaveAs("TargetThickness.png");
    }
    return gXS;
}

TGraphErrors *GetAverage(TGraphErrors **gXS, TH1D *hisYield, TH1D *hisErr)
{
    TGraphErrors *gA = new TGraphErrors();
    for (int iPoint = 0; iPoint < nBins; iPoint++)
    {
        double sumXS = 0, sumXSError = 0, sumEError = 0;
        int nValidDet = 0;
        for (int iDet = 0; iDet < 40; iDet++)
        {
            if (gXS[iDet]->GetErrorX(0) == -1) continue;
            sumXS += gXS[iDet]->GetPointY(iPoint);
            sumXSError += pow(gXS[iDet]->GetPointY(iPoint), 2);
            nValidDet++;
        }
        if (nValidDet == 0) continue;
        sumEError = sqrt(hisErr->GetBinContent(iPoint+1) / hisYield->GetBinContent(iPoint+1));
        if (isnan(sumEError) || isinf(sumEError)) sumEError = 0;

        gA->SetPoint(gA->GetN(), gXS[0]->GetPointX(iPoint), sumXS / nValidDet);
        gA->SetPointError(gA->GetN() - 1, sumEError, sqrt(sumXSError) / nValidDet);
    }
    return gA;
}
TGraphErrors *GetWeightedAverage(TGraphErrors **gXS, TH1D *hisErr)
{
    TGraphErrors *gWA = new TGraphErrors();
    for (int iPoint = 0; iPoint < nBins; iPoint++)
    {
        double sumXSWeight = 0, sumWeight = 0;
        for (int iDet = 0; iDet < 40; iDet++)
        {
            double xs = gXS[iDet]->GetPointY(iPoint);
            double errXS = gXS[iDet]->GetErrorY(iPoint);
            if (xs > 0 && errXS > 0)
            {
                double weight = 1 / (errXS * errXS);
                sumXSWeight += weight * xs;
                sumWeight += weight;
            }
        }
        double errEcm = sqrt(1/hisErr->GetBinContent(iPoint+1));

        if (sumWeight > 0)
        {
            double weightedXS = sumXSWeight / sumWeight;
            double errorXS = sqrt(1 / sumWeight);
            gWA->SetPoint(gWA->GetN(), gXS[0]->GetPointX(iPoint), weightedXS);
            gWA->SetPointError(gWA->GetN() - 1, errEcm, errorXS);
        }
    }

    return gWA;
}
