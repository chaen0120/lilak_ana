const double e2 = 1.44; // MeV * fm
const double hbarc = 197.327053; // MeV * fm
const double k = 8.61733326E-11; // MeV/K
const double T = 1.0E9; // K

double Sommerfeld(double E, double Z1, double Z2, double mu) {
    double v = sqrt(2 * E / mu);
    double eta = (Z1 * Z2 * e2) / (hbarc * v);
    return eta;
}
double yfunc(double *x, double *par) {
    double E = x[0]; // MeV
    double scale = par[0];
    double EkT = E / (k * T);
    return scale * exp(-EkT);
}
double yfunc2(double *x, double *par) {
    double E = x[0]; // MeV
    double scale = par[0];
    double Z1 = par[1];
    double Z2 = par[2];
    double mu = par[3]; // MeV/c^2

    double eta = Sommerfeld(E, Z1, Z2, mu);
    return scale * exp(-2 * TMath::Pi() * eta);
}
double yfunc3(double *x, double *par) {
    double E = x[0]; // MeV
    double EkT = E / (k * T);
    double scale = par[0];
    double Z1 = par[1];
    double Z2 = par[2];
    double mu = par[3]; // MeV/c^2

    double eta = Sommerfeld(E, Z1, Z2, mu);
    return scale * exp(-EkT -2 * TMath::Pi() * eta);
}

void drawGamow()
{
    auto Emin = 0.01;
    auto Emax = 5;

    auto Z1 = 8;
    auto Z2 = 2;
    auto mu = 2897.55; // MeV/c^2

    auto *f1 = new TF1("f1", yfunc, Emin, Emax, 1);
    f1->SetParameter(0,1);
    f1->SetLineColor(kBlue);
    f1->SetLineWidth(2);
    auto *f2 = new TF1("f2", yfunc2, Emin, Emax, 4);
    f2->SetParameters(1, Z1, Z2, mu);
    f2->SetLineColor(kRed);
    f2->SetLineWidth(2);

    auto *fall = new TF1("fall", yfunc3, Emin, Emax, 4);
    fall->SetParameters(1, Z1, Z2, mu);
    fall->SetLineColor(kBlack);
    fall->SetLineWidth(2);

    auto *cvs = new TCanvas("Gamow","Gamow",1500,2000);
    cvs->Divide(1,2);
    cvs->cd(1); gPad->SetLogy();
    f1->Draw();
    f2->Draw("same");
    fall->Draw("same");

    cvs->cd(2);
    auto *g = new TGraph();
    auto nBins = 10000;
    for (int i=0; i<nBins; i++)
    {
        double x = Emin + (Emax-Emin)/nBins*i;
        double y = fall->Eval(x);
        g->SetPoint(i, x, y);
    }
    auto *fg = new TF1("fg","gaus",Emin,Emax);
    fg->SetLineColor(kGreen+1);
    fg->SetLineWidth(2);
    g->Fit(fg);
    fg->Draw();
    g->Draw("L");
    fall->Draw("same");
    cvs->SaveAs("gamow.png");
    cvs->SaveAs("gamow.eps");
    cvs->SaveAs("gamow.root");
}