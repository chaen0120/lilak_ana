double Mass[4] = { 1.00782, 4.00260, 14.00860, 14.00307 };
double OEtoD[3][5] = { {  4.90820164348495E+02, -3.35284909595191E+01,  1.92121659126051E+01, -6.95182844235064E+00,  1.01350570608564E+00 },   // O energy: 0.5 ~ 2   MeV
                       {  4.82046153674622E+02, -1.40051689818663E+01,  2.41085083960792E+00, -3.43478505539321E-01,  1.82290221194483E-02 },   // O energy: 1.5 ~ 6   MeV
                       {  4.70839842824462E+02, -5.71929910986907E+00, -6.09882744374245E-02, -4.33206818870420E-04,  3.48137103177498E-06 } }; // O energy: 5   ~ 50  MeV
double ODtoE[3][5] = { {  4.99930173144693E+01, -7.35087975311014E-02, -5.43002600832763E-05,  8.10653271596185E-08, -2.38126239549678E-10 },   // O range : 50  ~ 5   MeV
                       {  1.99655747200540E+04, -1.82175276588626E+02,  6.24489884331181E-01, -9.52119612429623E-04,  5.44297113930965E-07 },   // O range : 6   ~ 1.5 MeV
                       {  1.04797874811022E+02,  2.46808178093234E+00, -1.56300323586775E-02,  2.93889916819240E-05, -1.76717558851694E-08 } }; // O range : 3   ~ 0.5 MeV
double NEtoD[3][5] = { {  5.84850663725171E+02, -3.37657122503956E+01,  1.92672813676005E+01, -7.23658364307753E+00,  1.08690987249433E+00 },   // N energy: 0.5 ~ 2   MeV
                       {  5.76465628537150E+02, -1.45311682515657E+01,  2.17686274683726E+00, -3.09886654049272E-01,  1.64871681961712E-02 },   // N energy: 1.5 ~ 6   MeV
                       {  5.65760176310850E+02, -6.92487337357859E+00, -5.31767634784621E-02, -1.10640745668768E-03,  8.28226330637759E-06 } }; // N energy: 5   ~ 50  MeV
double NDtoE[3][5] = { {  4.99971402056198E+01, -6.00781330573650E-02, -3.33644717466442E-05,  2.96371334518131E-08, -1.04379359597430E-10 },   // N range : 50  ~ 5   MeV
                       {  2.12919329493845E+04, -1.60800182917995E+02,  4.56233442298669E-01, -5.75770288657990E-04,  2.72481714178871E-07 },   // N range : 6   ~ 1.5 MeV
                       { -8.21525151356121E+02,  3.39714947471399E+00, -3.08648036253048E-04, -1.23527581451279E-05,  1.20659208024749E-08 } }; // N range : 2   ~ 0.5 MeV
double OcutE[4] = { 0.5, 1.7, 5.5, 50 };
double OcutD[4] = { 0, 437.50, 463.66, 490.31 };           // 50, 5.5, 1.7, 0.5     MeV
double NcutE[4] = { 0.5, 1.7, 5.5, 50 };
double NcutD[4] = { 0, 525.94, 556.66, 584.72 };           // 50, 5.5, 1.7, 0.5     MeV
TF1 *f[3];
double EtoD(double energy) 
{
    if (energy < OcutE[0]) return -1;
    else if (energy < OcutE[1]) return f[0]->Eval(energy);
    else if (energy < OcutE[2]) return f[1]->Eval(energy);
    else if (energy < OcutE[3]) return f[2]->Eval(energy);
    else return -1;
}

void drawThickTarget()
{
    auto *cvs = new TCanvas("cvs","cvs",1500,1200);

    f[0] = new TF1("f1", "pol4 - 136.007", OcutE[0], OcutE[1]);
    f[1] = new TF1("f2", "pol4 - 136.007", OcutE[1], OcutE[2]);
    f[2] = new TF1("f3", "pol4 - 136.007", OcutE[2], OcutE[3]);
    for (int i=0; i<3; i++)
        for (int j=0; j<5; j++)
            f[i]->SetParameter(j, OEtoD[i][j]);

    auto *g = new TGraphErrors();
    //for (double i=0; i<40; i+=0.01)
    //{
    //    auto dist = EtoD(i);
    //    if (fabs(dist) < 0.1) cout << i << endl;
    //    if (dist>0) g->SetPoint(g->GetN(), i*Mass[1]/(Mass[1]+Mass[2]), dist);
    //}
    for (int i=0; i<1000; i++)
    {
        double Ecm1 = i/100.;
        double Ebeam1 = Ecm1 * (Mass[1]+Mass[2]) / Mass[1];
        double dist1 = EtoD(Ebeam1);
        double Ecm2 = (i+1)/100.;
        double Ebeam2 = Ecm2 * (Mass[1]+Mass[2]) / Mass[1];
        double dist2 = EtoD(Ebeam2);
        if (dist1>0 && dist2>0)
        {
            g->SetPoint(g->GetN(), (Ebeam1+Ebeam2)/2, dist1 - dist2);
            g->SetPointError(g->GetN()-1, 0.05, 0);
        }
    }
    g->SetMarkerStyle(8);
    g->Draw("AP");

    //auto *fEZ = new TF1("fEZ","pol7",0,9);
    //fEZ->SetLineColor(kCyan+1);
    //fEZ->SetParameters(375.415, -72.9445, 45.3516, -22.2208, 5.69877, -0.805735, 0.0589968, -0.0017465);
    //fEZ->Draw("same");
}
