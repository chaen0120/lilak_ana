vector<double> rEcm, rXS; // experiment
vector<double> tEcm, tXS;
vector<double> t01Ecm, t01XS;
vector<double> rmb, tmb, t01mb;

const double Na = 6.022E23; // mol^-1
const double k = 1/11.6045; // MeV/GK
const double mu = 2898.5 * 9E-20; // MeV

const double Tmin = 0.2; // GK
const double Tmax = 10; // GK
const double Tstep = 0.01; // GK

double Rate(double T, TString mode)
{
    double prefactor = sqrt(8.0 / TMath::Pi() / mu) / pow(k * T, 1.5);
    double integral = 0;

    vector<double> *Ecm;
    vector<double> *XS;

    if      (mode == "r")   { Ecm = &rEcm;   XS = &rXS; }
    else if (mode == "t")   { Ecm = &tEcm;   XS = &tXS; }
    else if (mode == "t01") { Ecm = &t01Ecm; XS = &t01XS; }
    else { cerr << "Invalid mode: " << mode << endl; return 0; }

    for (size_t i = 0; i < Ecm->size() - 1; ++i)
    {
        double E1 = (*Ecm)[i];
        double E2 = (*Ecm)[i + 1];
        double XS1 = (*XS)[i];
        double XS2 = (*XS)[i + 1];
        double dE = E2 - E1;

        double integrand1 = E1 * XS1 * exp(-E1 / (k * T));
        double integrand2 = E2 * XS2 * exp(-E2 / (k * T));
        integral += 0.5 * (integrand1 + integrand2) * dE;
    }

    return Na * prefactor * integral; // cm^3/mol/s
}

double AsymGaus(TRandom3 &rng, double Ymean, double dy1, double dy2)
{
    double sigmaL = abs(dy1);
    double sigmaR = abs(dy2);
    double total = sigmaL + sigmaR;

    bool left = (rng.Uniform(0, total) < sigmaL);

    double val;
    do {
        val = rng.Gaus(0, left ? sigmaL : sigmaR);
    } while ((left && val > 0) || (!left && val < 0)); // 부호 조건

    return Ymean + val;
}

void ReacRate(bool isDraw = true)
{
    double scale = 1;
    ifstream fin;
    string line;
    double t1, t2, t3, t4, t5;

    /* XS */

    // Expriment
    vector<double> expEcm, expXS, expXSerr1, expXSerr2, expmb;
    fin.open("../totxs/totxs_14Oap_a.txt"); // XS(Ecm = 3.45) = 3.62684 mb : last point
    while (fin >> t1 >> t2 >> t3 >> t4)
    {
        if (t1 < 0.5 || t1 > 3.5) continue;
        expEcm.push_back(t1);
        expXS.push_back(t2 * 1E-27); //cm^2
        expmb.push_back(t2);
        expXSerr1.push_back(t4 * 1E-27); //cm^2
        if (abs(t1-3.45)<0.05) scale *= t2 * 1E-27;
    }
    cout << "scale: " << scale << endl;
    fin.close();
    fin.open("../totxs/totxs_14Oap_sysErr_a.txt");
    while (fin >> t1 >> t2 >> t3 >> t4)
    {
        if (t1 < 0.5 || t1 > 3.5) continue;
        expXSerr2.push_back(t4 * 1E-27); //cm^2
    }
    fin.close();
    //auto *gRes = new TGraph(expEcm.size(), expEcm.data(), expmb.data());
    //gRes->SetLineColor(kGreen+2);
    //gRes->SetLineWidth(3);
    auto *gRes = new TGraphErrors();
    gRes->SetLineColor(kGreen+2);
    gRes->SetLineWidth(3);
    gRes->SetFillColor(kGreen+1);
    gRes->SetFillStyle(3001);
    for (int i=0; i<expEcm.size(); i++)
    {
        gRes->SetPoint(i, expEcm[i], expmb[i]);
        gRes->SetPointError(i, 0.05, (expXSerr1[i] + expXSerr2[i]) * 1E27);
    }

    // TALYS; total
    fin.open("talys/ap.tot");
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        iss >> t1 >> t2 >> t3 >> t4 >> t5;
        tEcm.push_back(t1 * 14 / 18);
        tXS.push_back(t2 * 1E-27); //cm^2
        tmb.push_back(t2);
    }
    fin.close();
    auto *gXStot = new TGraph(tEcm.size(), tEcm.data(), tmb.data());
    gXStot->SetLineColor(kMagenta);
    gXStot->SetLineWidth(3);

    // TALYS; p0 + p1
    double scaletemp = 0;
    vector<double> t0XS, t1XS;
    fin.open("talys/ap.L00");
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        iss >> t1 >> t2 >> t3 >> t4;
        t01Ecm.push_back(t1 * 14 / 18);
        t0XS.push_back(t2 * 1E-27); //cm^2
    }
    fin.close();
    fin.open("talys/ap.L01");
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        iss >> t1 >> t2 >> t3 >> t4;
        t1XS.push_back(t2 * 1E-27); //cm^2
    }
    fin.close();
    for (int i=0; i<t01Ecm.size(); i++)
    {
        t01XS.push_back(t0XS[i] + t1XS[i]);
        t01mb.push_back((t0XS[i]+t1XS[i])*1E27);
        if (t01Ecm[i]*18/14 >= 4.4 && t01Ecm[i]*18/14 <= 4.5) scaletemp += t01XS[i];
    }
    auto *gXS = new TGraph(t01Ecm.size(), t01Ecm.data(), t01mb.data());
    gXS->SetLineColor(kBlack);
    gXS->SetLineWidth(3);

    scale /= scaletemp/2;
    cout << "scale: " << scale << endl;
    for (int i = 0; i < gXS->GetN(); i++)
    {
        double x, y;
        gXS->GetPoint(i, x, y);
        if (x >0.5) break;
        y *= scale*1E-27;
        rEcm.push_back(x);
        rXS.push_back(y);
        rmb.push_back(y*1E27);
    }
    for (int i=0; i < expEcm.size(); i++)
    {
        rEcm.push_back(expEcm[i]);
        rXS.push_back(expXS[i]);
        rmb.push_back(expmb[i]);
    }
    for (int i = 0; i < gXS->GetN(); i++)
    {
        double x, y;
        gXS->GetPoint(i, x, y);
        if (x <3.5) continue;
        y *= scale*1E-27;
        rEcm.push_back(x);
        rXS.push_back(y);
        rmb.push_back(y*1E27);
    }
    auto *gXSscale = new TGraph(rEcm.size(), rEcm.data(), rmb.data());
    gXSscale->SetLineColor(kGreen+1);
    gXSscale->SetLineWidth(3);
    //gXSscale->SetLineStyle(kDashed);

    /* Reaction Rate */

    // AZURE
    auto *gAZURE = new TGraph("AZURE.out");
    gAZURE->SetLineColor(kBlue);
    gAZURE->SetLineWidth(3);
    //gAZURE->SetLineStyle(kDotted);

    // Reaclib
    auto *gReaclib = new TGraph("Reaclib.dat");
    gReaclib->SetLineColor(kRed);
    gReaclib->SetLineWidth(3);
    //gReaclib->SetLineStyle(kDotted);

    // TALYS
    auto *gTALYS = new TGraph();
    gTALYS->SetLineColor(kMagenta+1);
    gTALYS->SetLineWidth(3);
    fin.open("talys/astrorate.p");
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        iss >> t1 >> t2 >> t3 >> t4;
        gTALYS->SetPoint(gTALYS->GetN(), t1, t2);
    }
    fin.close();

    // Calculate
    auto *gRate = new TGraph();
    gRate->SetLineColor(kBlack);
    gRate->SetLineWidth(3);
    auto *gRateTot = new TGraph();
    gRateTot->SetLineColor(kMagenta);
    gRateTot->SetLineWidth(3);
    auto *gExpRate = new TGraph();
    gExpRate->SetFillColor(kGreen+1);
    gExpRate->SetFillStyle(3001);
    gExpRate->SetLineColor(kGreen+1);
    gExpRate->SetLineWidth(3);
    for (double T = Tmin; T < Tmax; T += Tstep) //GK
    {
        auto rate = Rate(T, "t01");
        gRate->SetPoint(gRate->GetN(), T, rate);

        auto rateTot = Rate(T, "t");
        gRateTot->SetPoint(gRateTot->GetN(), T, rateTot);

        auto rateRes = Rate(T, "r");
        gExpRate->SetPoint(gExpRate->GetN(), T, rateRes);
    }

    auto *gMExpRate = new TGraphErrors();
    gMExpRate->SetFillColor(kGreen + 1);
    gMExpRate->SetFillStyle(3001);
    if (isDraw)
    {
        fin.open("ExpRate.dat");
        while (fin >> t1 >> t2 >> t3 >> t4)
        {
            gMExpRate->SetPoint(gMExpRate->GetN(), t1, t2);
            gMExpRate->SetPointError(gMExpRate->GetN() - 1, 0.05, t4);
        }
    }
    else
    {
        TRandom3 rng(0);
        int iter = 6E5;
        auto *hMExpRate = new TH2D("hMExpRate", "hMExpRate", (Tmax - Tmin) / Tstep + 1, Tmin, Tmax, 1000, -15, 6);
        for (int i = 0; i < iter; i++)
        {
            cout << "\r" << i << "/" << iter << flush;
            bool first = true;
            int start = 0;
            for (int iE = 0; iE < rEcm.size(); iE++)
            {
                if (rEcm[iE] <= 0.5) continue;
                if (rEcm[iE] >= 3.5) break;
                if (first) { start = iE; first = false; }
                auto newXS = AsymGaus(rng, expXS[iE - start], expXSerr1[iE - start] + expXSerr2[iE - start], expXSerr1[iE - start]);
                rXS[iE] = newXS;
            }

            for (double T = Tmin; T < Tmax; T += Tstep) // GK
            {
                auto rateRes = Rate(T, "r");
                hMExpRate->Fill(T, log10(rateRes));
            }
        }
        cout << endl;

        auto *gMExpRatelog = new TGraphErrors();
        gMExpRatelog->SetFillColor(kGreen + 1);
        gMExpRatelog->SetFillStyle(3001);
        TH1D *proj;
        for (int iT = 0; iT < hMExpRate->GetNbinsX(); iT++)
        {
            proj = hMExpRate->ProjectionY("proj", iT + 1, iT + 1);
            double x = Tmin + iT * Tstep;
            double y = proj->GetMean();
            double err = proj->GetRMS();
            gMExpRatelog->SetPoint(iT, x, y);
            gMExpRatelog->SetPointError(iT, 0.05, err);
            proj->Reset("ICES");
        }

        const int window = 3; // 양쪽 1포인트씩 평균 (총 2N+1 포인트 사용)
        int nPoints = gMExpRatelog->GetN();

        auto *gSmoothed = new TGraphErrors();
        gSmoothed->SetFillColor(kGreen + 1);
        gSmoothed->SetFillStyle(3001);

        for (int i = 0; i < nPoints; ++i)
        {
            double x, y, ex, ey;
            gMExpRatelog->GetPoint(i, x, y);
            ex = gMExpRatelog->GetErrorX(i);
            ey = gMExpRatelog->GetErrorY(i);

            // 에러 스무딩
            double sumErr = 0.0;
            int count = 0;
            for (int j = i - window; j <= i + window; ++j)
            {
                if (j < 0 || j >= nPoints)
                    continue;
                sumErr += gMExpRatelog->GetErrorY(j);
                count++;
            }
            double smoothErr = (count > 0) ? sumErr / count : ey;

            gMExpRate->SetPoint(i, x, pow(10, y));
            gMExpRate->SetPointError(i, ex, pow(10, y) * log(10) * smoothErr);
            cout << x << "\t" << pow(10, y) << "\t" << ey << "\t" << pow(10, y) * log(10) * smoothErr << endl;
        }
    }

    /* Draw */

    auto *cvs = new TCanvas("cvs","cvs", 1000,1500);
    cvs->Divide(1,2);

    cvs->cd(1); gPad->SetLogy();
    gXS->GetXaxis()->SetTitle("E_{cm} [MeV]");
    gXS->GetXaxis()->CenterTitle();
    gXS->GetXaxis()->SetTitleOffset(0.9);
    gXS->GetXaxis()->SetLabelSize(0.04);
    gXS->GetXaxis()->SetTitleSize(0.05);
    gXS->GetYaxis()->SetTitle("Cross Section [mb]");
    gXS->GetYaxis()->CenterTitle();
    gXS->GetYaxis()->SetTitleOffset(0.8);
    gXS->GetYaxis()->SetLabelSize(0.04);
    gXS->GetYaxis()->SetTitleSize(0.05);
    gXS->GetXaxis()->SetLimits(0,10);
    gXS->GetYaxis()->SetRangeUser(0.01, 1E+3);
    gXS->Draw("LA");
    gXStot->Draw("Lsame");
    gXSscale->Draw("Lsame");
    gRes->Draw("xLsame");
    gRes->Draw("e3same");

    auto *leg1 = new TLegend(0.7, 0.15, 0.88, 0.3);
    leg1->AddEntry(gXStot, "TALYS (total)", "l");
    leg1->AddEntry(gXS, "TALYS (p0+p1)", "l");
    leg1->AddEntry(gRes, "Exp. + TALYS", "lf");
    leg1->Draw();

    cvs->cd(2); gPad->SetLogy(); gPad->SetLogx();
    gRate->GetXaxis()->SetTitle("Temperature [GK]");
    gRate->GetXaxis()->CenterTitle();
    gRate->GetXaxis()->SetTitleOffset(1);
    gRate->GetXaxis()->SetLabelSize(0.04);
    gRate->GetXaxis()->SetTitleSize(0.05);
    gRate->GetYaxis()->SetTitle("Reaction Rate [cm^{3}/s/mol]");
    gRate->GetYaxis()->CenterTitle();
    gRate->GetYaxis()->SetTitleOffset(1);
    gRate->GetYaxis()->SetLabelSize(0.04);
    gRate->GetYaxis()->SetTitleSize(0.05);
    gRate->GetXaxis()->SetLimits(0.2,10);
    gRate->GetYaxis()->SetRangeUser(1E-16, 1E+8);
    gRate->Draw("CA"); 
    gRateTot->Draw("Csame");
    gReaclib->Draw("Csame");
    gAZURE->Draw("Csame");
    //gTALYS->Draw("Csame");
    gExpRate->Draw("Csame");
    gMExpRate->Draw("e3same");

    auto *leg2 = new TLegend(0.7, 0.15, 0.88, 0.4);
    leg2->AddEntry(gRateTot, "TALYS (total)", "l");
    leg2->AddEntry(gRate, "TALYS (p0+p1)", "l");
    //leg2->AddEntry(gTALYS, "TALYS (rate)", "l");
    leg2->AddEntry(gReaclib, "REACLIB", "l");
    leg2->AddEntry(gAZURE, "AZURE", "l");
    leg2->AddEntry(gExpRate, "Exp. + TALYS", "lf");
    leg2->Draw();

    //cvs->cd(3); gPad->SetLogx();
    //hMExpRate->Draw("colz");
    //gMExpRatelog->Draw("e3same");

    cvs->SaveAs("cvs.root");
}
