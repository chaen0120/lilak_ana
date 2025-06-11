int iProton = 0;
int iAlpha = 1;
int i14O = 2;
int i14N = 3;
double Mass[4] = {1.00782, 4.00260, 14.00860, 14.00307};
double PEtoD[2][5] = {{5.22557432373778E+04, -6.00914623289058E+01, -1.34161139073868E+02, 9.16540747672604E+00, -6.89539665679287E-01},  // P energy: 0.1 ~ 4   MeV
                      {5.23942866946498E+04, -1.82265845381518E+02, -8.91299661752110E+01, 6.54622882848444E-01, -6.04320727892540E-03}}; // P energy: 3   ~ 25  MeV
double PDtoE[4][5] = {{2.49988393317393E+01, -2.55403523023795E-04, -2.55098349329561E-09, 6.02647152629786E-14, -1.31563325443243E-18},  // P range : 25  ~ 10  MeV
                      {-4.55315614010784E+02, 4.46346723563283E-02, -1.57870032529034E-06, 2.47102850343635E-11, -1.46257450716399E-16},  // P range : 11  ~ 5   MeV
                      {-9.75646431885352E+04, 7.89021796723329E+00, -2.39267276275527E-04, 3.22506827360814E-09, -1.63051663902410E-14},  // P range : 6   ~ 1.5 MeV
                      {2.83261922161464E+05, -1.09457228987476E+01, -1.55106404277375E-09, 4.08688822452100E-09, -3.94889582762636E-14}}; // P range : 2   ~ 0.5 MeV
double AEtoD[2][5] = {{7.82670812310234E+03, -2.61215704458080E+01, -4.26209335992390E+00, -5.48396317089452E-01, 4.23840487550776E-02},  // A energy: 0.5 ~ 6   MeV
                      {7.83067095649276E+03, -2.31238465294816E+01, -6.91057278066634E+00, 4.81956819578733E-02, -4.15494339786851E-04}}; // A energy: 5   ~ 35  MeV
double ADtoE[3][5] = {{3.49948641392046E+01, -2.42085902074895E-03, -1.62791655473277E-07, 2.63644908159824E-11, -3.84386471183837E-15},  // A range : 35  ~ 12  MeV
                      {-1.88048698664164E+03, 1.15059460028287E+00, -2.60622547470514E-04, 2.62018030116540E-08, -9.91549323278050E-13},  // A range : 14  ~ 5   MeV
                      {-8.99004284147101E+05, 4.74858901136244E+02, -9.40600597752470E-02, 8.28100911278438E-06, -2.73415678999220E-10}}; // A range : 6   ~ 0.5 MeV
double OEtoD[3][5] = {{4.90820164348495E+02, -3.35284909595191E+01, 1.92121659126051E+01, -6.95182844235064E+00, 1.01350570608564E+00},   // O energy: 0.5 ~ 2   MeV
                      {4.82046153674622E+02, -1.40051689818663E+01, 2.41085083960792E+00, -3.43478505539321E-01, 1.82290221194483E-02},   // O energy: 1.5 ~ 6   MeV
                      {4.70839842824462E+02, -5.71929910986907E+00, -6.09882744374245E-02, -4.33206818870420E-04, 3.48137103177498E-06}}; // O energy: 5   ~ 50  MeV
double ODtoE[3][5] = {{4.99930173144693E+01, -7.35087975311014E-02, -5.43002600832763E-05, 8.10653271596185E-08, -2.38126239549678E-10},  // O range : 50  ~ 5   MeV
                      {1.99655747200540E+04, -1.82175276588626E+02, 6.24489884331181E-01, -9.52119612429623E-04, 5.44297113930965E-07},   // O range : 6   ~ 1.5 MeV
                      {1.04797874811022E+02, 2.46808178093234E+00, -1.56300323586775E-02, 2.93889916819240E-05, -1.76717558851694E-08}};  // O range : 3   ~ 0.5 MeV
double NEtoD[3][5] = {{5.84850663725171E+02, -3.37657122503956E+01, 1.92672813676005E+01, -7.23658364307753E+00, 1.08690987249433E+00},   // N energy: 0.5 ~ 2   MeV
                      {5.76465628537150E+02, -1.45311682515657E+01, 2.17686274683726E+00, -3.09886654049272E-01, 1.64871681961712E-02},   // N energy: 1.5 ~ 6   MeV
                      {5.65760176310850E+02, -6.92487337357859E+00, -5.31767634784621E-02, -1.10640745668768E-03, 8.28226330637759E-06}}; // N energy: 5   ~ 50  MeV
double NDtoE[3][5] = {{4.99971402056198E+01, -6.00781330573650E-02, -3.33644717466442E-05, 2.96371334518131E-08, -1.04379359597430E-10},  // N range : 50  ~ 5   MeV
                      {2.12919329493845E+04, -1.60800182917995E+02, 4.56233442298669E-01, -5.75770288657990E-04, 2.72481714178871E-07},   // N range : 6   ~ 1.5 MeV
                      {-8.21525151356121E+02, 3.39714947471399E+00, -3.08648036253048E-04, -1.23527581451279E-05, 1.20659208024749E-08}}; // N range : 2   ~ 0.5 MeV
// [range]
double PcutE[3] = {0.1, 3.5, 25};
double PcutD[5] = {0, 40380, 48800, 51804.73, 52193.52}; // 25, 11, 5.5, 1.7, 0.5 MeV
double AcutE[3] = {0.5, 5.5, 35};
double AcutD[4] = {0, 6460, 7502, 7812.86}; // 35, 13, 5.5, 0.5      MeV
double OcutE[4] = {0.5, 1.7, 5.5, 50};
double OcutD[4] = {0, 437.50, 463.66, 490.31}; // 50, 5.5, 1.7, 0.5     MeV
double NcutE[4] = {0.5, 1.7, 5.5, 50};
double NcutD[4] = {0, 525.94, 556.66, 584.72}; // 50, 5.5, 1.7, 0.5     MeV
double CMtoLab(double ene, double angle, int light, int beam);
double ApplyEnergyLoss(double ene, double len, int ptl);

void SpreaddEdx(double spread = 1/1200., double thres = 0.0012)
//void SpreaddEdx(double spread = 1/1500., double thres = 0.00065)
{
    //gStyle->SetOptStat(0);
    gRandom->SetSeed(0);

    auto *cvs = new TCanvas("cvs","cvs",2000,1600);
    cvs->Divide(2,2);
    cvs->cd(1); gPad->SetMargin(0.15,0.1,0.15,0.1);
    cvs->cd(2); gPad->SetMargin(0.15,0.1,0.15,0.1);
    cvs->cd(3); gPad->SetMargin(0.15,0.1,0.15,0.1);
    cvs->cd(4); gPad->SetMargin(0.15,0.1,0.15,0.1);
    
    auto *fdEdx = new TF1("fdEdx","[0]/([1]+x)",1,19);
    fdEdx->SetParameters(0.00617753,1.28108);
    auto *hisdEdx = new TH2D("hisdEdx","hisdEdx;E_{det} [MeV];dE/dx [MeV/mm]",201,0,20,200,0,0.005);
    for (double iEcm = 1; iEcm < 19; iEcm+=0.1)
    {
        auto mean = fdEdx->Eval(iEcm);
        auto sigma = spread;
        for (int i=0; i<1000; i++) 
        {
            auto random = gRandom->Gaus(mean,sigma);
            hisdEdx->Fill(iEcm, random);
        }
    }
    cvs->cd(1); hisdEdx->Draw("colz");
    fdEdx->Draw("same");

    auto *g = new TGraph();
    g->SetPoint(g->GetN(), 2.26, thres);
    g->SetPoint(g->GetN(), 2.26, thres*2.62);
    g->SetPoint(g->GetN(), 10.7 ,thres);
    g->SetPoint(g->GetN(), 10.7 ,thres*1.83);
    //g->SetPoint(g->GetN(), 2.62, thres);
    //g->SetPoint(g->GetN(), 2.62, thres*1405/327);
    //g->SetPoint(g->GetN(), 11.2 ,thres);
    //g->SetPoint(g->GetN(), 11.2 ,thres*720/327);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(kRed);
    g->Draw("P");

    auto *line = new TLine(0, thres, 20, thres);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("same");

    auto *hisRatio = new TH1D("hisRatio","hisRatio",200,0,20);
    for (int iE = 1; iE<=200; iE++)
    {
        double nTot = 0, nAbove = 0;
        for (int idEdx = 1; idEdx<=200; idEdx++)
        {
            auto cont = hisdEdx->GetBinContent(iE,idEdx);
            nTot += cont;
            if (hisdEdx->GetYaxis()->GetBinCenter(idEdx)>thres) nAbove += cont;
        }
        if (nTot>0) hisRatio->SetBinContent(iE,nAbove/nTot);
    }
    hisRatio->GetYaxis()->SetRangeUser(0,1);
    //hisRatio->Fit("pol7","R","",1,19);

    auto *fiin = new TFile("his.root");
    auto *hiss = (TH2D*) fiin->Get("h");
    auto *fin = new TFile("eff.root");
    auto *his = (TH1D*) fin->Get("his_ratio0");
    his->GetXaxis()->SetRangeUser(0,20);
    his->SetLineColor(kRed);

    hisRatio->GetXaxis()->SetTitle("Tracking Eff.");
    hisRatio->GetYaxis()->SetTitle("E_{det} [MeV]");
    hiss->GetXaxis()->SetTitle("E_{det} [MeV]");
    hiss->GetYaxis()->SetTitle("dE/dx [arb.]");

    hisdEdx->GetXaxis()->SetTitleSize(0.05);
    hisdEdx->GetXaxis()->CenterTitle();
    hisdEdx->GetXaxis()->SetLabelSize(0.04);
    hisdEdx->GetYaxis()->SetTitleSize(0.05);
    hisdEdx->GetYaxis()->CenterTitle();
    hisdEdx->GetYaxis()->SetLabelSize(0.04);
    hisdEdx->GetYaxis()->SetTitleOffset(1.6);
    hiss->GetXaxis()->SetTitleSize(0.05);
    hiss->GetXaxis()->CenterTitle();
    hiss->GetXaxis()->SetLabelSize(0.04);
    hiss->GetYaxis()->SetTitleSize(0.05);
    hiss->GetYaxis()->CenterTitle();
    hiss->GetYaxis()->SetLabelSize(0.04);
    hiss->GetYaxis()->SetTitleOffset(1.1);
    hisRatio->GetXaxis()->SetTitleSize(0.05);
    hisRatio->GetXaxis()->CenterTitle();
    hisRatio->GetXaxis()->SetLabelSize(0.04);
    hisRatio->GetYaxis()->SetTitleSize(0.05);
    hisRatio->GetYaxis()->CenterTitle();
    hisRatio->GetYaxis()->SetLabelSize(0.04);

    cvs->cd(2); hiss->Draw("colz");
    cvs->cd(3); hisRatio->Draw();
    cvs->cd(4); hisRatio->Draw("hist"); his->Draw("same hist");

    auto *leg = new TLegend(0.70,0.75,0.88,0.88);
    leg->AddEntry(hisRatio, "Calculated", "l");
    leg->AddEntry(his, "Experiment", "l");
    leg->Draw();
    cvs->SaveAs("draw.png");
}

void CaldEdx()
{
    auto *run = new LKRun();
    auto *tt = new TexAT2();
    auto *ana = new TTAnalysisTask();
    run->AddPar("../config/config_common.mac");
    run->AddPar("../config/config_reco.mac");
    run->AddPar("../config/config_ana.mac");
    run->AddDetector(tt);
    run->Add(ana);
    run->Init();

    double lenA = 40.3;
    double lenB = 75;
    double lenF = 50;
    int sectionA[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    int sectionB[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};

    auto fZE = new TF1("fZE", "pol7", 0, 375);
    fZE->SetParameters(8.98833, -0.00603845, -0.000362325, 4.78713e-06, -3.47837e-08, 1.36886e-10, -2.76809e-13, 2.24699e-16);
    auto fEZ = new TF1("fEZ", "pol7", 0, 9);
    fEZ->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465);

    TH2D *his_dEdx = new TH2D("his_dEdx", "his_dEdx;Edet;dEdx", 190, 1, 20, 200, 0, 0.05);
    for (int iDet : {0, 6, 7, 8, 9,
                     11, 12, 13, 14, 15, 16, 17,
                     101, 102, 103, 104, 105, 106, 107, 108,
                     201, 202, 203, 204, 205, 206, 207, 208})
    {
        double detX = tt->GetX6posx(iDet, 4);
        double detY = tt->GetX6posy(iDet, 4);
        double detZ = tt->GetX6posz(iDet, 4);
        if (iDet < 10)
        {
            detZ = 468.29;
            if (iDet == 0) { detX = 124;   detY = 35.8;  }
            if (iDet == 6) { detX = -62.2; detY = 30.8;  }
            if (iDet == 7) { detX = -62.2; detY = -30.8; }
            if (iDet == 8) { detX = -124;  detY = 35.8;  }
            if (iDet == 9) { detX = -124;  detY = -35.8; }
        }

        for (int section = 0; section < 9; section++)
        {
            double x, y, z;
            if (iDet < 10)
            {
                x = detX + sectionA[section] * lenF / 2;
                y = detY + sectionB[section] * lenF / 2;
                z = detZ;
            }
            else if (iDet < 20)
            {
                x = detX;
                y = detY + sectionB[section] * lenB / 2;
                z = detZ + sectionA[section] * lenA / 2;
            }
            else if (iDet == 101 || iDet == 204)
            {
                x = detX + sectionA[section] * lenA / 2;
                y = detY;
                z = detZ + sectionB[section] * lenB / 2;
            }
            else
            {
                x = detX + sectionB[section] * lenB / 2;
                y = detY;
                z = detZ + sectionA[section] * lenA / 2;
            }
            auto si = TVector3(x, y, z);

            for (double Ecm = 8.6; Ecm > 0; Ecm -= 0.1)
            {
                auto vertex = TVector3(0, 0, fEZ->Eval(Ecm));
                auto angle = si.Angle(vertex);
                auto length = (si - vertex).Mag();

                auto Elab = CMtoLab(Ecm, angle, iAlpha, i14O);
                cout << Ecm << " " << Elab << endl;
                auto Edet = ApplyEnergyLoss(Elab, length, iAlpha);
                auto dEdx = (Elab - Edet) / length;
                his_dEdx->Fill(Edet, dEdx);
            }
        }
    }

    gStyle->SetPalette(kBird);
    his_dEdx->Draw("colz");
    auto *f = new TF1("f","[0]/([1]+x)",1,19);
    his_dEdx->Fit(f);
}

double CMtoLab(double ene, double angle, int light, int beam)
{
    double Edet = -1;
    double Mbeam, Mtarget, Mlight, Mheavy;
    double cos = TMath::Cos(angle);
    double amu = 931.5016;

    Mbeam = Mass[beam];
    Mtarget = Mass[iAlpha];
    if (light == iAlpha) // Scattering
        Edet = ene / Mtarget * (Mbeam + Mtarget) * (Mbeam + Mtarget) / (4 * Mbeam * cos * cos);
    else // (a,p) reaction
    {
        Mlight = 1.00782;
        if (beam == i14O)
            Mheavy = 17.00210;
        if (beam == i14N)
            Mheavy = 16.99913;
        double m1 = Mbeam * amu;      // MeV
        double m2 = Mtarget * amu;    // MeV
        double m3 = Mlight * amu;     // MeV
        double m4 = Mheavy * amu;     // MeV
        double Q = m1 + m2 - m3 - m4; // MeV

        double E1 = ene / Mtarget * (Mbeam + Mtarget);
        double a = -m3 - m4;
        double b = sqrt(m1 * m3 * E1) * cos;
        double c = m4 * Q + E1 * (m4 - m1);

        double Esq = (-b - sqrt(b * b - a * c)) / a;
        Edet = Esq * Esq;
    }
    return Edet;
}
double ApplyEnergyLoss(double ene, double len, int ptl)
{
    double Eadd = 0;
    int nRangeE, nRangeD;
    double cutE[5], cutD[5];       // [range]
    double EtoD[4][5], DtoE[4][5]; // [range][constant]
    if (ptl == iProton)
    {
        nRangeE = 2;
        nRangeD = 4;
        for (int i = 0; i < nRangeE; i++)
        {
            cutE[i] = PcutE[i];
            for (int j = 0; j < 5; j++)
                EtoD[i][j] = PEtoD[i][j];
        }
        for (int i = 0; i < nRangeD; i++)
        {
            cutD[i] = PcutD[i];
            for (int j = 0; j < 5; j++)
                DtoE[i][j] = PDtoE[i][j];
        }
        cutE[nRangeE] = PcutE[nRangeE];
        cutD[nRangeD] = PcutD[nRangeD];
    }
    else if (ptl == iAlpha)
    {
        nRangeE = 2;
        nRangeD = 3;
        for (int i = 0; i < nRangeE; i++)
        {
            cutE[i] = AcutE[i];
            for (int j = 0; j < 5; j++)
                EtoD[i][j] = AEtoD[i][j];
        }
        for (int i = 0; i < nRangeD; i++)
        {
            cutD[i] = AcutD[i];
            for (int j = 0; j < 5; j++)
                DtoE[i][j] = ADtoE[i][j];
        }
        cutE[nRangeE] = AcutE[nRangeE];
        cutD[nRangeD] = AcutD[nRangeD];
    }
    else if (ptl == i14O)
    {
        nRangeE = 3;
        nRangeD = 3;
        for (int i = 0; i < nRangeE; i++)
        {
            cutE[i] = OcutE[i];
            for (int j = 0; j < 5; j++)
                EtoD[i][j] = OEtoD[i][j];
        }
        for (int i = 0; i < nRangeD; i++)
        {
            cutD[i] = OcutD[i];
            for (int j = 0; j < 5; j++)
                DtoE[i][j] = ODtoE[i][j];
        }
        cutE[nRangeE] = OcutE[nRangeE];
        cutD[nRangeD] = OcutD[nRangeD];
    }
    else if (ptl == i14N)
    {
        nRangeE = 3;
        nRangeD = 3;
        for (int i = 0; i < nRangeE; i++)
        {
            cutE[i] = NcutE[i];
            for (int j = 0; j < 5; j++)
                EtoD[i][j] = NEtoD[i][j];
        }
        for (int i = 0; i < nRangeD; i++)
        {
            cutD[i] = NcutD[i];
            for (int j = 0; j < 5; j++)
                DtoE[i][j] = NDtoE[i][j];
        }
        cutE[nRangeE] = NcutE[nRangeE];
        cutD[nRangeD] = NcutD[nRangeD];
    }

    int range = -1;
    for (int i = 0; i < nRangeE + 1; ++i)
    {
        if (ene < cutE[i])
        {
            if (i == 0)
                return cutE[0];
            range = i - 1;
            break;
        }
    }
    if (range == -1)
        return cutE[nRangeE];

    double Dcal = 0;
    for (Int_t i = 0; i < 5; i++)
        Dcal += EtoD[range][i] * TMath::Power(ene, i);
    Dcal += len;

    range = -1;
    for (int i = 0; i < nRangeD + 1; ++i)
    {
        if (Dcal < cutD[i])
        {
            if (i == 0)
                return cutE[nRangeE];
            range = i - 1;
            break;
        }
    }
    if (range == -1)
        return cutE[0];

    for (Int_t i = 0; i < 5; i++)
        Eadd += DtoE[range][i] * TMath::Power(Dcal, i);

    return Eadd;
}
