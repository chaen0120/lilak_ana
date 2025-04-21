void Threshold()
{
    TFile *fin = new TFile("threshold.root");
    TH2D *his = (TH2D*) fin->Get("h");

    auto Emin = his->GetYaxis()->GetXmin();
    auto Emax = his->GetYaxis()->GetXmax();
    auto nE   = his->GetNbinsY();

    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    TCanvas *cvs = new TCanvas("cvs","cvs",1200,1000);
    cvs->cd();
    for (int iDet=0; iDet<40; iDet++)
    {
        int nEntry = 2;
        if (iDet>=10) nEntry = 3;
        if (iDet== 8 || iDet==28) nEntry = 1;
        if (iDet==11 || iDet==25 || iDet==26 || iDet==34) nEntry = 4;
        if (iDet==33 || iDet==35) nEntry = 5;
        if (iDet==24) nEntry = 6;

        double thres = 0;
        auto proj = his->ProjectionY(Form("proj_%d",iDet),detmap[iDet]+1,detmap[iDet]+1);
        if (proj->GetEntries()==0) { cout << iDet << " " << thres << endl; continue; }
        for (int iE=1; iE<=nE; iE++)
        {
            if (proj->GetBinContent(iE)>nEntry) { thres = proj->GetBinCenter(iE); break; }
        }
        cout << iDet << " " << thres << endl;
        proj->GetYaxis()->SetRangeUser(0,10);
        proj->Draw();
        TLine *l = new TLine(thres,0,thres,10);
        l->SetLineColor(kRed);
        l->Draw("same");
        cvs->Update();
        cvs->WaitPrimitive();
    }
}
/* result
0 1.615
1 0
2 0
3 0
4 0
5 0
6 1.423
7 1.455
8 1.767
9 1.511
10 1.427
11 1.405
12 1.453
13 1.499
14 1.451
15 1.449
16 1.443
17 0
18 0
19 0
20 0
21 0
22 0
23 0
24 1.511
25 1.407
26 1.509
27 1.533
28 1.887
29 0
30 1.483
31 1.495
32 1.575
33 1.537
34 1.465
35 1.525
36 1.471
37 1.423
38 1.491
39 1.563
*/