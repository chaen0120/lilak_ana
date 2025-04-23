void BeamRes2()
{
    Double_t fx[4] = { 9.12603, 7.50399, 5.9015, 3.88722};
    //Double_t fy[4] = { -3.51768, 84.4275, 158.865, 238.117};
    Double_t fy[4] = { 0, 469/6, 469/3, 469/2 };
    Double_t fex[4] = { 0.444438, 0.644592, 0.71113, 0.891728};
    Double_t fey[4] = { 0, 0, 0, 0};
    TGraphErrors *gre = new TGraphErrors(4,fy,fx,fey,fex);
    gre->SetName("Graph");
    gre->SetTitle("Graph");
    //gre->Draw("AP*");

    TF1 *f[3];
    f[0] = new TF1("f0","pol7",24,375);
    f[1] = new TF1("f1","[10]*pol7(0)+pol1(8)",24,375);
    f[2] = new TF1("f2","[10]*pol7(0)+pol1(8)",24,400);
    for (int i=0; i<3; i++)
    {
        f[i]->FixParameter(0,     8.98833);
        f[i]->FixParameter(1, -0.00603845);
        f[i]->FixParameter(2,-0.000362325);
        f[i]->FixParameter(3, 4.78713e-06);
        f[i]->FixParameter(4,-3.47837e-08);
        f[i]->FixParameter(5, 1.36886e-10);
        f[i]->FixParameter(6,-2.76809e-13);
        f[i]->FixParameter(7, 2.24699e-16);
    }
    f[0]->SetLineColor(kRed);
    //f[0]->Draw("same");

    double min[4], max[4];
    for (int i=0; i<4; i++) { min[i] = fx[i]-fex[i]; max[i] = fx[i]+fex[i]; }
    TGraph *g1 = new TGraph(4,fy,min);
    TGraph *g2 = new TGraph(4,fy,max);
    //g1->Draw("P*same");
    //g2->Draw("P*same");
    g1->Fit(f[1]);
    g2->Fit(f[2]);

    auto *cvs = new TCanvas("cvs","cvs",800,1200);
    cvs->Divide(1,2);
    cvs->cd(1);
    f[0]->Draw("");
    f[1]->Draw("same");
    f[2]->Draw("same");
    gre->Draw("*same");

    //auto *his = new TH2D("his","his",400,0,400,100,0,10);
    //for (double iEcm = 0; iEcm < 10; iEcm += 0.1)
    //{
    //    auto mean = f[0]->Eval(iEcm);
    //    auto diff1 = mean - f[1]->Eval(iEcm);
    //    auto diff2 = f[2]->Eval(iEcm) - mean;
    //    auto sigma = diff1 > diff2 ? diff1 : diff2;
    //    for (int i=0; i<100; i++) his->Fill(gRandom->Gaus(mean,sigma),iEcm);
    //}
    //his->Draw("colz");
    cvs->cd(2);
    auto *fdiff = new TF1("fdiff", [=](double *x, double*) {return (f[2]->Eval(x[0]) - f[1]->Eval(x[0]))*4.00260/(4.00260+14.00860);}, 24,375,0);
    fdiff->Draw();
}
