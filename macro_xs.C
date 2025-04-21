#include "macro_ana.C"
#include "totxs/draw_cs.C"

TH1D *hisM_E[40], *hisB_E[40], *hisN_E[40], *hisM_EGood[40];
TH1D *hisM_EError, *hisM_wEError;
TGraphErrors *g_xs[40], *g_xsGood[40], *g_sysError[40];

void macro_xs(TString type = "14Oap", bool isGoodEvents = false, bool isSubtractBG = true)
{
    // Init
    bool runMain = false;
    bool runBG = false;
    cout << "isGoodEvents = " << isGoodEvents << ", isSubtractBG = " << isSubtractBG << endl;

    beam = type[2];
    reac = type[4];
    if (runMain) macro_ana(type);
    if (runBG)   macro_ana(type(0,3)+"CO2"+reac);
    if      (type(0, 5) == "14Oap" || type(0, 5) == "14Oaa") nBeam = nBeam14Oa;
    else if (type(0, 5) == "14Nap" || type(0, 5) == "14Naa") nBeam = nBeam14Na;
    else if (type(0, 6) == "14OCO2") nBeam = nBeam14OC;
    else if (type(0, 6) == "14NCO2") nBeam = nBeam14NC;
    else if (type(0, 6) == "14Ovac") nBeam = nBeam14Ov;
    else nBeam = 200000000000;

    SetOthers();
    GetEfficiency(effMain, nEcmMain, false);
    GetEfficiency(effBG, nEcmBG, true);

    // Init histograms
    TH1D *hisM_EAll     = new TH1D("hisM_EAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisM_EAllGood = new TH1D("hisM_EAllGood", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisB_EAll     = new TH1D("hisB_EAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *hisN_EAll     = new TH1D("hisN_EAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);
    TH1D *his_EAll      = new TH1D("his_EAll", Form("Yield;E_{cm}[MeV]; Counts / %d",(int)(binSize*1000)), nBins, minEcm, maxEcm);

    auto fnameMain = Form("results/Yield_%s_%d.root",type.Data(),(int)(binSize*1000));
    cout << "Main: " << fnameMain << endl;
    auto finMain = new TFile(fnameMain);
    hisM_EError = (TH1D*) finMain->Get("his_EError");
    hisM_wEError = (TH1D*) finMain->Get("his_weightedEError");
    for (int i=0; i<40; i++)
    {
        hisM_E[i] = (TH1D*) finMain->Get(Form("his_E_%d",i));
        hisM_EGood[i] = (TH1D*) finMain->Get(Form("his_EGood_%d",i));

        hisM_EAll->Add(hisM_E[i]);
        hisM_EAllGood->Add(hisM_EGood[i]);
    }
    if (isSubtractBG)
    {
        auto fnameBG = Form("results/Yield_%s_%d.root",(type(0,3)+"CO2"+reac).Data(),(int)(binSize*1000));
        cout << "CO2 : " << fnameBG << endl;
        auto finBG = new TFile(fnameBG);
        for (int i=0; i<40; i++)
        {
            hisB_E[i] = (TH1D*) finBG->Get(Form("his_E_%d",i));
            hisB_EAll->Add(hisB_E[i]);
        }
        auto fnameN = Form("results/Yield_14Oap_14N_%d.root",(int)(binSize*1000));
        cout << "14N : " << fnameN << endl;
        auto finN = new TFile(fnameN);
        for (int i=0; i<40; i++)
        {
            hisN_E[i] = (TH1D*) finN->Get(Form("his_E_%d",i));
            hisN_EAll->Add(hisN_E[i]);
        }
    }
    else
        for (int i=0; i<40; i++)
        {
            hisB_E[i] = (TH1D*) hisM_E[i]->Clone();
            hisB_E[i]->Reset("ICES");
            hisN_E[i] = (TH1D*) hisM_E[i]->Clone();
            hisN_E[i]->Reset("ICES");
        }

    his_EAll = (TH1D *)hisM_EAll->Clone();
    his_EAll->Add(hisB_EAll, -1);
    his_EAll->Add(hisN_EAll, -1);
    if (isGoodEvents)
    {
        his_EAll->Reset("ICES");
        his_EAll->Add(hisM_EAllGood);
        his_EAll->Add(hisB_EAll, -1);
        his_EAll->Add(hisN_EAll, -1);
    }
    
    // Get XS
    for (int iDet=0; iDet<40; iDet++)
    {
        g_xs[iDet] = GetCrossSection(hisM_E[iDet], hisB_E[iDet], hisN_E[iDet], iDet);
        g_xsGood[iDet] = GetCrossSection(hisM_EGood[iDet], hisB_E[iDet], hisN_E[iDet], iDet);
    }
    auto g_xsAll = GetAverage(g_xs, hisM_EAll, hisM_EError);
    g_xsAll->SetTitle(Form("Total Cross Section of %s;E_{cm} [MeV];XS [mb]",type.Data()));
    auto g_weightedxsAll = GetWeightedAverage(g_xs, hisM_wEError);
    g_weightedxsAll->SetTitle(Form("Total Cross Section (weighted average) of %s;E_{cm} [MeV];XS [mb]",type.Data()));

    auto g_xsGoodAll = GetAverage(g_xsGood, hisM_EAll, hisM_EError);
    auto g_sysErr = new TGraphErrors();
    for (int iPoint = 0; iPoint<nBins; iPoint++)
    {
        auto minAll = g_xsAll->GetPointY(iPoint) - g_xsAll->GetErrorY(iPoint);
        auto minGood = g_xsGoodAll->GetPointY(iPoint) - g_xsGoodAll->GetErrorY(iPoint);
        g_sysErr->SetPoint(iPoint, g_xsAll->GetPointX(iPoint), (minAll+minGood)/2);
        g_sysErr->SetPointError(iPoint, 0, (minAll-minGood)/2);
    }

    gStyle->SetOptStat(0);
    TCanvas *cvs_xs = new TCanvas("cvs_xs","cvs_xs",1500,1200);
    cvs_xs->Divide(2,2);
    cvs_xs->cd(1);
    hisM_EAll->SetLineColor(kBlack);
    hisM_EAll->Draw("HIST");
    hisM_EAllGood->SetLineColor(kRed);
    hisM_EAllGood->Draw("HISTSAME");
    hisB_EAll->Add(hisN_EAll);
    hisB_EAll->SetFillColor(kGreen + 1);
    hisB_EAll->SetLineColor(kGreen + 1);
    hisB_EAll->SetFillStyle(3005);
    hisB_EAll->Draw("HISTSAME");
    cvs_xs->cd(2); his_EAll->Draw("HIST");
    cvs_xs->cd(3); g_xsAll->Draw("APL"); g_sysErr->Draw("P");
    cvs_xs->cd(4); g_weightedxsAll->Draw("APL");
    SaveBatch(cvs_xs);

    ofstream outtxt(Form("results/totxs_%s.txt",type.Data()));
    for (int i=0; i<g_xsAll->GetN(); i++)
        outtxt << g_xsAll->GetPointX(i) << "\t" << g_xsAll->GetPointY(i) << "\t"
               << g_xsAll->GetErrorX(i) << "\t" << g_xsAll->GetErrorY(i) << endl;
    outtxt.close();
    if (type(0,5) == "14Oap" || type=="simap") 
    {
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
}

void GetEfficiency(vector<vector<double>> &eff, vector<double> &nEcm, bool isBackground)
{
    for (int i = 0; i < 40; i++) eff[i].resize(nBins);
    nEcm.resize(nBins);
    TString effname = Form("proton_%d.root",(int)(binSize*1000));
    if (reac=='a') effname = Form("alpha_%d.root",(int)(binSize*1000));
    if (isBackground) effname = "CO2_" + effname;
    effname = "eff/" + effname;
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

TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, TH1D *hisN, int iDet)
{
    TGraphErrors *gXS = new TGraphErrors();
    gXS->SetName(Form("xs_%s", hisMain->GetName()));
    gXS->SetTitle(Form("xs_%s", hisMain->GetName()));
    if (hisMain->GetEntries() == 0) 
    {
        for (int iE = 1; iE <= nBins; iE++)
        {
            gXS->SetPoint(gXS->GetN(), hisMain->GetBinCenter(iE), 0);
            gXS->SetPointError(gXS->GetN()-1, -1, -1);
        }
        return gXS;
    }

    auto hiseffMain = ApplyEfficiency(hisMain, iDet, false);
    auto hiseffBG   = ApplyEfficiency(hisBG  , iDet, true);
    auto hiseffN    = ApplyEfficiency(hisN   , iDet, false);

    auto gasDensity = 0.1313*0.9;   // mg/cm^3
    auto HeMolarMass = 4.0026*1000; // mg/mol
    auto avo = 6.022E+23;           // #/mol

    for (int iE = 1; iE <= nBins; iE++)
    {
        double iEmin = hisMain->GetBinCenter(iE) - binSize / 2;
        double iEmax = hisMain->GetBinCenter(iE) + binSize / 2;
        double len = (EtoZ(iEmin) - EtoZ(iEmax))/10; //[cm]

        double yieldMain = hisMain->GetBinContent(iE);
        double yieldBG   = hisBG  ->GetBinContent(iE);
        double yieldN    = hisN   ->GetBinContent(iE);
        double yield = hiseffMain->GetBinContent(iE) - hiseffBG->GetBinContent(iE) - hiseffN->GetBinContent(iE);
        double target = len * gasDensity / HeMolarMass * avo;
        double xs = yield / nBeam / target * 1E+27; //[mb]

        double ecmerror = binSize / 2; // will be updated at averaging functions

        auto effErrMain = effMain[iDet][iE-1] * sqrt(1/nEcmMain[iE-1] + 1/nEcmMain[iE-1]/effMain[iDet][iE-1]);
        auto  xsErrMain = sqrt(pow(sqrt(yieldMain)/effMain[iDet][iE-1],2) +
                               pow(yieldMain*effErrMain/effMain[iDet][iE-1]/effMain[iDet][iE-1],2));
        auto effErrBG  = effBG[iDet][iE-1] * sqrt(1/nEcmBG[iE-1] + 1/nEcmBG[iE-1]/effBG[iDet][iE-1]);
        auto  xsErrBG  = sqrt(pow(sqrt(yieldBG)/effBG[iDet][iE-1],2) +
                              pow(yieldBG*effErrBG/effBG[iDet][iE-1]/effBG[iDet][iE-1],2));
        auto effErrN   = effMain[iDet][iE-1] * sqrt(1/nEcmMain[iE-1] + 1/nEcmMain[iE-1]/effMain[iDet][iE-1]);
        auto  xsErrN   = sqrt(pow(sqrt(yieldN)/effMain[iDet][iE-1],2) +
                             pow(yieldN*effErrN/effMain[iDet][iE-1]/effMain[iDet][iE-1],2));
        auto xserror = sqrt(xsErrMain*xsErrMain + xsErrBG*xsErrBG + xsErrN*xsErrN) / nBeam / target * 1E+27;

        if (len <= 0 || isnan(xs) || isinf(xs)) xs = 0;
        if (isnan(xserror)) xserror = 0;
        gXS->SetPoint(gXS->GetN(), hisMain->GetBinCenter(iE), xs);
        gXS->SetPointError(gXS->GetN()-1, ecmerror, xserror);
    }
    //if (!gROOT->IsBatch())
    //{
    //    cout << "Draw " << iDet << endl;
    //    TCanvas *cvs = new TCanvas("cvs", "cvs", 1200, 600);
    //    cvs->Divide(2, 1);
    //    cvs->cd(1); hisMain->Draw();
    //    cvs->cd(2); gXS->Draw();
    //    cvs->WaitPrimitive();
    //}
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
        sumEError = sqrt(hisErr->GetBinContent(iPoint+1)) / hisYield->GetBinContent(iPoint+1);
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
