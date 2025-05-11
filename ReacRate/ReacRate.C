vector<double> rEcm, rXS; // experiment
vector<double> tEcm, tXS;
vector<double> t01Ecm, t01XS;

const double Na = 6.022E23; // mol^-1
const double k = 1/11.6045; // MeV/GK
const double mu = 2898.5 * 9E-20; // MeV

double Rate(double T, TString mode)
{
    double prefactor = sqrt(8.0 / TMath::Pi() / mu) / pow(k * T, 1.5);
    double integral = 0;

    vector<double> *Ecm;
    vector<double> *XS;

    if (mode == "r")        { Ecm = &rEcm;   XS = &rXS; }
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

void ReacRate()
{
    double scale = 1;
    ifstream fin;
    string line;
    double t1, t2, t3, t4, t5;

    /* XS */

    // Expriment
    vector<double> expEcm, expXS;
    fin.open("../totxs/totxs_14Oap_a.txt"); // XS(Ecm = 3.45) = 3.62684 mb : last point
    while (fin >> t1 >> t2 >> t3 >> t4)
    {
        if (t1 < 0.5 || t1 > 3.5) continue;
        expEcm.push_back(t1);
        expXS.push_back(t2 * 1E-27); //cm^2
        if (t1 == 3.45) scale *= t2 * 1E-27;
    }
    fin.close();
    auto *gRes = new TGraph(expEcm.size(), expEcm.data(), expXS.data());
    gRes->SetLineColor(kGreen+2);
    gRes->SetLineWidth(3);

    // TALYS; total
    fin.open("talys/ap.tot");
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        iss >> t1 >> t2 >> t3 >> t4 >> t5;
        tEcm.push_back(t1 * 14 / 18);
        tXS.push_back(t2 * 1E-27); //cm^2
    }
    fin.close();
    auto *gXStot = new TGraph(tEcm.size(), tEcm.data(), tXS.data());
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
        if (t01Ecm[i]*18/14 >= 4.4 && t01Ecm[i]*18/14 <= 4.5) scaletemp += t01XS[i];
    }
    auto *gXS = new TGraph(t01Ecm.size(), t01Ecm.data(), t01XS.data());
    gXS->SetLineColor(kBlack);
    gXS->SetLineWidth(3);

    scale /= scaletemp/2;
    for (int i = 0; i < gXS->GetN(); i++)
    {
        double x, y;
        gXS->GetPoint(i, x, y);
        if (x >0.5) break;
        y *= scale;
        rEcm.push_back(x);
        rXS.push_back(y);
    }
    for (int i=0; i < expEcm.size(); i++)
    {
        rEcm.push_back(expEcm[i]);
        rXS.push_back(expXS[i]);
    }
    for (int i = 0; i < gXS->GetN(); i++)
    {
        double x, y;
        gXS->GetPoint(i, x, y);
        if (x <3.5) continue;
        y *= scale;
        rEcm.push_back(x);
        rXS.push_back(y);
    }
    auto *gXSscale = new TGraph(rEcm.size(), rEcm.data(), rXS.data());
    gXSscale->SetLineColor(kGreen+1);
    gXSscale->SetLineWidth(3);
    gXSscale->SetLineStyle(kDashed);

    /* Reaction Rate */

    // AZURE
    auto *gAZURE = new TGraph("AZURE.out");
    gAZURE->SetLineColor(kBlue);
    gAZURE->SetLineWidth(3);
    gAZURE->SetLineWidth(kDotted);

    // Reaclib
    auto *gReaclib = new TGraph("Reaclib.dat");
    gReaclib->SetLineColor(kRed);
    gReaclib->SetLineWidth(3);
    gReaclib->SetLineWidth(kDotted);

    // Calculate
    auto *gRate = new TGraph();
    gRate->SetLineColor(kBlack);
    gRate->SetLineWidth(3);
    auto *gRateTot = new TGraph();
    gRateTot->SetLineColor(kMagenta);
    gRateTot->SetLineWidth(3);
    auto *gExpRate = new TGraph();
    gExpRate->SetLineColor(kGreen+1);
    gExpRate->SetLineWidth(3);
    for (double T = 0.2; T < 10; T += 0.1) //GK
    {
        auto rate = Rate(T, "t01");
        gRate->SetPoint(gRate->GetN(), T, rate);

        auto rateTot = Rate(T, "t");
        gRateTot->SetPoint(gRateTot->GetN(), T, rateTot);

        auto rateRes = Rate(T, "r");
        gExpRate->SetPoint(gExpRate->GetN(), T, rateRes);
    }

    /* Draw*/

    auto *cvs = new TCanvas("cvs","cvs", 1500,800);
    cvs->Divide(2,1);
    cvs->cd(1); gPad->SetLogy();
    gXS->GetXaxis()->SetTitle("E_{cm} [MeV]");
    gXS->GetYaxis()->SetTitle("Cross Section [mb]");
    gXS->Draw("CA");
    gXStot->Draw("Csame");
    gXSscale->Draw("Csame");
    gRes->Draw("Csame");
    cvs->cd(2); gPad->SetLogy(); gPad->SetLogx();
    gRate->GetXaxis()->SetTitle("Temperature [GK]");
    gRate->GetYaxis()->SetTitle("Reaction Rate [cm^3/s/mol]");
    gRate->Draw("CA"); 
    gRateTot->Draw("Csame");
    gExpRate->Draw("Csame");
    gReaclib->Draw("Csame");
    gAZURE->Draw("Csame");
}