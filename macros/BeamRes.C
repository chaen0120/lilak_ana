void BeamRes()
{
    Double_t fx[4] = { 9.12603, 7.50399, 5.9015, 3.88722};
    Double_t fy[4] = { -3.51768, 84.4275, 158.865, 238.117};
    //Double_t fy[4] = { 0, 469/6, 469/3, 469/2 };
    Double_t fex[4] = { 0.444438, 0.644592, 0.71113, 0.891728};
    Double_t fey[4] = { 0, 0, 0, 0};
    TGraphErrors *gre = new TGraphErrors(4,fx,fy,fex,fey);
    gre->SetName("Graph");
    gre->SetTitle("Graph");
    //gre->Draw("AP*");

    TF1 *f[3];
    f[0] = new TF1("f0","pol7",0,9);
    f[1] = new TF1("f1","pol7(0)+pol1(8)",0,9);
    f[2] = new TF1("f2","pol7(0)+pol1(8)",0,9);
    for (int i=0; i<3; i++)
    {
        f[i]->FixParameter(0,375.415    );
        f[i]->FixParameter(1,-72.9445   );
        f[i]->FixParameter(2,45.3516    );
        f[i]->FixParameter(3,-22.2208   );
        f[i]->FixParameter(4,5.69877    );
        f[i]->FixParameter(5,-0.805735  );
        f[i]->FixParameter(6,0.0589968  );
        f[i]->FixParameter(7,-0.0017465 );
    }
    f[0]->SetLineColor(kRed);
    //f[0]->Draw("same");

    double min[4], max[4];
    for (int i=0; i<4; i++) { min[i] = fx[i]-fex[i]; max[i] = fx[i]+fex[i]; }
    TGraph *g1 = new TGraph(4,min,fy);
    TGraph *g2 = new TGraph(4,max,fy);
    //g1->Draw("C*same");
    //g2->Draw("C*same");
    g1->Fit(f[1]);
    g2->Fit(f[2]);

    f[0]->Draw();
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
}
