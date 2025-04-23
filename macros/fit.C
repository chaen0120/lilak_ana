double Mass[4] = { 1.00782, 4.00260, 14.00860, 14.00307 };

/* 14O
*/
const int nF = 4;
const int nP = 4; //5
double maxE = 40;
double OcutEmin[nF] = { 0.5, 1.5, 3.6, 10.0 };
double OcutEmax[nF] = { OcutEmin[1]+0.1, OcutEmin[2]+0.5, OcutEmin[3]+1, maxE };
double OcutDmax[nF] = { 270, 310, 348, 367};
double OcutDmin[nF] = { 0, OcutDmax[0]-10, OcutDmax[1]-10, OcutDmax[2]-10};
/* 14N
const int nF = 4;
const int nP = 4; //5
double maxE = 30;
double OcutEmin[nF] = { 0.5, 1.5, 3.3, 10.0 };
double OcutEmax[nF] = { OcutEmin[1]+0.1, OcutEmin[2]+0.5, OcutEmin[3]+1, maxE };
double OcutDmax[nF] = { 160, 230, 270, 291};
double OcutDmin[nF] = { 0, OcutDmax[0]-10, OcutDmax[1]-10, OcutDmax[2]-10};
*/
/*alpha
const int nF = 4;
const int nP = 4; //5
double maxE = 30;
double OcutEmin[nF] = { 0.1, 1.0, 3.3, 5.5 };
double OcutEmax[nF] = { OcutEmin[1]+0.1, OcutEmin[2]+0.5, OcutEmin[3]+1, maxE };
double OcutDmax[nF] = { 4730, 5800, 5920, 5950};
double OcutDmin[nF] = { 0, OcutDmax[0]-450, OcutDmax[1]-290, OcutDmax[2]-70};
*/
/*proton
const int nF = 4;
const int nP = 4; //5
double maxE = 20;
double OcutEmin[nF] = { 0.1, 0.9, 1.7, 5. };
double OcutEmax[nF] = { OcutEmin[1]+0.2, OcutEmin[2]+0.5, OcutEmin[3]+1, maxE };
double OcutDmax[nF] = { 28500, 33000, 34400, 34850};
double OcutDmin[nF] = { 0, OcutDmax[0]-8000, OcutDmax[1]-4000, OcutDmax[2]-1080};
*/

double OcutE[nF+1], OcutD[nF+1];
TF1 *f[nF];
double EtoD(double energy) 
{
    if (energy < OcutE[0] || energy >= OcutE[nF]) return -1;
    for (int i=0; i<nF; i++)
        if (energy < OcutE[i+1]) return f[i]->Eval(energy);
    return -1;
}
double DtoE(double dist)
{
    if (dist < OcutD[0] || dist >= OcutD[nF]) return -1;
    for (int i=0; i<nF; i++)
        if (dist < OcutD[i+1]) return f[i]->Eval(dist);
    return -1;
}

/* DtoE
void fit()
{
    auto *cvs = new TCanvas("cvs","cvs",2000,1200);
    cvs->Divide(2,1);

    ifstream fin("14N.txt");
    map<double, double> Eloss;
    double tE, tR;
    while ( fin >> tE >> tR) Eloss[tE] = tR;

    double maxR = Eloss[maxE];
    cout << maxR << endl;
    auto *gFit = new TGraph();
    for (auto it = Eloss.begin(); it != Eloss.end(); ++it)
    {
        if (it->first > maxE) continue;
        auto D = maxR - it->second;
        gFit->SetPoint(gFit->GetN(), D, it->first);
    }

    cvs->cd(1); gFit->Draw("AP*");

    for (int i=0; i<nF; i++)
    {
        f[i] = new TF1(Form("f%d",i),Form("pol%d",nP));
        f[i]->SetLineColor(i+2);
        gFit->Fit(f[i],"RQ+","",OcutDmin[i],OcutDmax[i]);
    }

    OcutD[0] = OcutDmin[0];
    OcutD[nF] = OcutDmax[nF-1];
    for (int i=1; i<nF; i++)
    {
        auto *ftemp = new TF1("ftemp", 
                              [=](double *x, double*) {return f[i-1]->Eval(x[0]) - f[i]->Eval(x[0]);}, 
                              OcutDmin[i], OcutDmax[i],0);
        OcutD[i] = ftemp->GetX(0,OcutDmin[i], OcutDmax[i]);
        cout << OcutD[i] << endl;
    }

    auto *gDif = new TGraphErrors();
    for (int i=0; i<1000; i++)
    {
        double dist1 = i/1000.*maxR;
        double Ebeam1 = DtoE(dist1);
        double dist2 = (i+1)/1000.*maxR;
        double Ebeam2 = DtoE(dist2);
        if (dist1>0 && dist2>0)
        {
            if (Ebeam1 - Ebeam2 > 1) continue;
            gDif->SetPoint(gDif->GetN(), (dist1+dist2)/2, Ebeam1 - Ebeam2);
            //gDif->SetPoint(gDif->GetN(), dist1, Ebeam1);
        }
    }
    gDif->SetMarkerStyle(8);
    cvs->cd(2); gDif->Draw("AP");
    cvs->SaveAs("cvs.png");

    cout << "{ ";
    for (int iF=0; iF<=nF; iF++) cout << OcutD[iF] << ", ";
    cout << "};" << endl;
    for (int iF=0; iF<nF; iF++)
    {
        if (iF==0) cout << "{ " << endl;
        for (int iP=0; iP<nP+1; iP++)
        {
            if (iP==0) cout << "{ ";
            cout << f[iF]->GetParameter(iP);
            if (iP!=nP) cout << ", ";
            else 
            {
                if (iF==nF-1) cout << " }";
                else cout << " },";
            }
        }
        cout << endl;
        if (iF==nF-1) cout << " };" << endl;
    }
}
*/
/* EtoD */
void fit()
{
    auto *cvs = new TCanvas("cvs","cvs",2000,1200);
    cvs->Divide(2,1);

    ifstream fin("14O.txt");
    map<double, double> Eloss;
    double tE, tR;
    while ( fin >> tE >> tR) Eloss[tE] = tR;

    double maxR = Eloss[maxE];
    auto *gFit = new TGraph();
    //for (auto& [E, R] : Eloss)
    //{
    //    if (E>maxE) continue;
    //    auto D = maxR - R;
    //    gFit->SetPoint(gFit->GetN(), E, D);
    //}
    for (auto it = Eloss.begin(); it != Eloss.end(); ++it)
    {
        if (it->first > maxE) continue;
        auto D = maxR - it->second;
        gFit->SetPoint(gFit->GetN(), it->first, D);
    }

    cvs->cd(1); gFit->Draw("AP*");

    for (int i=0; i<nF; i++)
    {
        f[i] = new TF1(Form("f%d",i),Form("pol%d",nP));
        f[i]->SetLineColor(i+2);
        gFit->Fit(f[i],"RQ+","",OcutEmin[i],OcutEmax[i]);
    }

    OcutE[0] = OcutEmin[0];
    OcutE[nF] = OcutEmax[nF-1];
    for (int i=1; i<nF; i++)
    {
        auto *ftemp = new TF1("ftemp", 
                              [=](double *x, double*) {return f[i-1]->Eval(x[0]) - f[i]->Eval(x[0]);}, 
                              OcutEmin[i], OcutEmax[i],0);
        OcutE[i] = ftemp->GetX(0,OcutEmin[i], OcutEmax[i]);
        cout << OcutE[i] << endl;
    }

    auto *gDif = new TGraphErrors();
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
            gDif->SetPoint(gDif->GetN(), (Ebeam1+Ebeam2)/2, dist1 - dist2);
            gDif->SetPointError(gDif->GetN()-1, 0.05, 0);
            //gDif->SetPoint(gDif->GetN(), Ebeam1, dist1);
        }
    }
    gDif->SetMarkerStyle(8);
    cvs->cd(2); gDif->Draw("AP");
    cvs->SaveAs("cvs.png");

    cout << "{ ";
    for (int iF=0; iF<=nF; iF++) cout << OcutE[iF] << ", ";
    cout << "};" << endl;
    for (int iF=0; iF<nF; iF++)
    {
        if (iF==0) cout << "{ " << endl;
        for (int iP=0; iP<nP+1; iP++)
        {
            if (iP==0) cout << "{ ";
            cout << f[iF]->GetParameter(iP);
            if (iP!=nP) cout << ", ";
            else 
            {
                if (iF==nF-1) cout << " }";
                else cout << " },";
            }
        }
        cout << endl;
        if (iF==nF-1) cout << " };" << endl;
    }
    cout << EtoD(39.112) << endl;
}
    //auto *cvs2 = new TCanvas("TargetThickness","TargetThickness",1500,1000);
    //gDif->SetMarkerColor(kBlue);
    //gDif->SetMarkerSize(1.3);
    //gDif->GetXaxis()->SetLabelSize(0.04);
    //gDif->GetXaxis()->SetTitle("E_{cm} [MeV]");
    //gDif->GetXaxis()->SetTitleSize(0.05);
    //gDif->GetXaxis()->SetTitleOffset(0.90);
    //gDif->GetXaxis()->CenterTitle();
    //gDif->GetYaxis()->SetLabelSize(0.04);
    //gDif->GetYaxis()->SetTitle("Effective Target Thickness [mm]");
    //gDif->GetYaxis()->SetTitleSize(0.05);
    //gDif->GetYaxis()->SetTitleOffset(0.75);
    //gDif->GetYaxis()->CenterTitle();
    //cvs2->cd(); gDif->Draw("AP");
    //cvs2->SaveAs("TargetThickness.root");
    //cvs2->SaveAs("TargetThickness.eps");


    //double dist = 0;
    //int nRange = 0;
    //for (int i=0; i<nF; i++)
    //    if (energy>=OcutEmin[i] && energy<OcutEmax[i])
    //    {
    //        nRange++;
    //        dist += f[i]->Eval(energy);
    //    }
    //if (nRange==0) return -1;
    //else return (double)(dist/nRange);
