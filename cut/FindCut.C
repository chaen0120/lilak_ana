int drawmap[40] = {0, 6, 7, 8, 9, -1, -1, -1,
                   14, 15, 16, 10, 11, 12, 13, -1,
                   23, 22, 21, 20, 19, 18, 17, -1,
                   28, 29, 30, 31, 25, 26, 27, 24,
                   39, 38, 37, 36, 34, 33, 32, 35};

void PrintProgress(int iEvt, int nEvt) { if (iEvt % 100 == 0) { cout << ">> " << iEvt << " / " << nEvt << "\r"; cout.flush(); } if (iEvt+1 == nEvt) cout << endl; }
void DrawAllDetector(TH2D **his, TString CvsName = "cvs")
{
    TCanvas* cvs = new TCanvas(CvsName,CvsName,2000,1200);
    cvs -> Divide(8,5);
    for(int i=0; i<40; i++) 
    { 
        if(his[drawmap[i]]==nullptr) continue;
        cvs -> cd(i+1); 
        if(drawmap[i]!=-1) his[drawmap[i]] -> Draw("colz"); 
    }
}

void MakeCut();
void DrawCut();
void DrawCut2();
TH2D* his_angle[40];
TH2D* his_angle_real[40];
TH2D* his_angleXY[40];
TH2D* his_angleZY[40];

void FindCut()
{
    TFile *fin = new TFile("../ana_14Oap.root");
    TTree *tree = (TTree *)fin->Get("event");

    TClonesArray *EnderArray = nullptr;
    TClonesArray *HeaderArray = nullptr;
    tree->SetBranchAddress("EventEnder", &EnderArray);
    tree->SetBranchAddress("EventHeader", &HeaderArray);
    auto nEvts = tree->GetEntries();

    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }
    for (int i = 0; i < 40; i++)
    {
        his_angle[i] = new TH2D(Form("his_angle_%d",i),Form("his_angle_%d",detmap[i]),180,0,180,40,0,40);
        his_angle_real[i] = new TH2D(Form("his_angle_real_%d",i),Form("his_angle_real_%d",detmap[i]),180,0,180,40,0,40);
        his_angleXY[i] = new TH2D(Form("his_angleXY_%d",i),Form("his_angleXY_%d",detmap[i]),180,0,180,40,0,40);
        his_angleZY[i] = new TH2D(Form("his_angleZY_%d",i),Form("his_angleZY_%d",detmap[i]),180,0,180,40,0,40);
    }

    for (int iEvt = 0; iEvt < nEvts; iEvt++)
    {
        PrintProgress(iEvt, nEvts);
        tree->GetEntry(iEvt);
        auto header = (TTEventHeader *)HeaderArray->At(0);
        auto ender = (TTEventEnder *)EnderArray->At(0);

        auto firedDet = header->GetFiredDet();
        auto Edet = ender->GetEdet();
        auto Sivert = TVector3(0, 0, ender->GettZ());
        auto SiPos = ender->GetSiHit();
        auto realang = ender->GetAlab();

        auto track = SiPos - Sivert;
        auto ang = TMath::ACos(track.Z() / track.Mag()) * TMath::RadToDeg();
        auto trackXY = TVector2(track.X(), track.Y());
        auto trackZY = TVector2(track.Z(), track.Y());
        auto angXY = trackXY.Phi() * TMath::RadToDeg() + 90;
        auto angZY = trackZY.Phi() * TMath::RadToDeg() + 90;
        while (angXY >= 180) angXY -= 180;
        while (angZY >= 180) angZY -= 180;
        if(angXY<0 || angXY>=180) cout << angXY << endl;
        if(angZY<0 || angZY>=180) cout << angZY << endl;

        int iDet = -1;
        for (int i = 0; i < 40; i++)
            if (detmap[i] == firedDet)
            {
                his_angle[i]->Fill(ang, Edet);
                his_angle_real[i]->Fill(realang, Edet);
                his_angleXY[i]->Fill(angXY, Edet);
                his_angleZY[i]->Fill(angZY, Edet);
            }
    }
    TFile *fout = new TFile("FindCut.root","recreate");
    fout -> cd();
    for (int i=0; i<40; i++)
    {
        his_angle[i] -> Write();
        his_angle_real[i] -> Write();
        his_angleXY[i] -> Write();
        his_angleZY[i] -> Write();
    }
    fout -> Close();
}

void MakeCut()
{
    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    TH2D *his_angle[2][40];
    TFile *fin = new TFile("FindCut.root");
    for(int i=0; i<40; i++)
    {
        his_angle[0][i] = (TH2D*) fin -> Get(Form("his_angleXY_%d",i));
        his_angle[1][i] = (TH2D*) fin -> Get(Form("his_angleZY_%d",i));
    }

    double thres = 1;
    double AMinMax[40][25][2][2][2]; //[det][energy][XY/ZY][first/second cut][min/max]
    TGraph *g[2][40];
    for (int j = 0; j < 2; j++) for (int i = 0; i < 40; i++) g[j][i] = new TGraph();
    for (int iXZ: {0, 1}) {
        for (int iDet = 0; iDet < 40; iDet++) {
            for (int iHalf : {0, 1})
            {
                for (int i = 0; i < 25; i++)
                    for (int j = 0; j < 2; j++)
                        AMinMax[iDet][i][iXZ][iHalf][j] = 0;
                if((detmap[iDet]<100 || iXZ==0) && iHalf==1) continue;

                int Rmin = 1;
                int Rmax = 180;
                if (detmap[iDet] > 100 && iXZ==1)
                {
                    if (iHalf == 0) Rmax = 90;
                    if (iHalf == 1) Rmin = 90;
                }

                for (int Edet = 1; Edet <= 25; Edet++)
                {
                    double AMin = -1;
                    double AMax = -1;
                    for (int Alab = Rmin; Alab <= Rmax; Alab++) {
                        if (his_angle[iXZ][iDet]->GetBinContent(Alab, Edet) >= thres) {
                            if (AMin == -1)
                                AMin = his_angle[iXZ][iDet]->GetXaxis()->GetBinCenter(Alab);
                            AMax = his_angle[iXZ][iDet]->GetXaxis()->GetBinCenter(Alab);
                        }
                    }
                    if (AMin != -1 && AMax != -1)
                    {
                        AMinMax[iDet][Edet - 1][iXZ][iHalf][0] = AMin - 2;
                        AMinMax[iDet][Edet - 1][iXZ][iHalf][1] = AMax + 5;
                        //AMinMax[iDet][Edet - 1][iXZ][iHalf][1] = AMax + 2;
                    }
                    else
                        //cout << Edet << " no entries" << endl;

                    if (AMinMax[iDet][Edet - 1][iXZ][iHalf][0] == 0)
                    {
                        AMinMax[iDet][Edet - 1][iXZ][iHalf][0] = 0;
                        AMinMax[iDet][Edet - 1][iXZ][iHalf][1] = AMinMax[iDet][Edet - 2][iXZ][iHalf][1];
                    }
                }
            }
        }
    }

    ofstream outtxt("cut_EdetvsAlab.txt");
    for (int iDet = 0; iDet < 40; iDet++) {
        for (int iEne = 0; iEne < 25; iEne++) {
            outtxt << iDet << "\t" << iEne << "\t";
            for (int iXZ = 0; iXZ < 2; iXZ++) {
                for (int iHalf = 0; iHalf < 2; iHalf++) {
                    outtxt << AMinMax[iDet][iEne][iXZ][iHalf][0] << "\t" << AMinMax[iDet][iEne][iXZ][iHalf][1] << "\t";
                    if (detmap[iDet] < 100 || iXZ == 0) {
                        if (iHalf==1) continue;
                        g[iXZ][iDet]->SetPoint(iEne, AMinMax[iDet][iEne][iXZ][iHalf][0] - 2, iEne);
                        g[iXZ][iDet]->SetPoint(49 - iEne, AMinMax[iDet][iEne][iXZ][iHalf][1] + 2, iEne);
                    }
                    else {
                        g[iXZ][iDet]->SetPoint(iEne+iHalf*50, AMinMax[iDet][iEne][iXZ][iHalf][0] - 2, iEne);
                        g[iXZ][iDet]->SetPoint(49 - iEne+iHalf*50, AMinMax[iDet][iEne][iXZ][iHalf][1] + 2, iEne);
                    }
                }
            }
            outtxt << endl;
        }
    }

    TCanvas *cvs_make0 = new TCanvas("cvs_make0", "cvs_make0", 2000, 1200); cvs_make0->Divide(8, 5);
    TCanvas *cvs_make1 = new TCanvas("cvs_make1", "cvs_make1", 2000, 1200); cvs_make1->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if(his_angle[0][drawmap[i]]==nullptr) continue;
        if(drawmap[i]==-1) continue;
        cvs_make0->cd(i + 1);
        his_angle[0][drawmap[i]]->Draw("colz");
        g[0][drawmap[i]] -> SetLineColor(kRed);
        g[0][drawmap[i]] -> Draw("SAME");
        cvs_make1->cd(i + 1);
        his_angle[1][drawmap[i]]->Draw("colz");
        g[1][drawmap[i]] -> SetLineColor(kRed);
        g[1][drawmap[i]] -> Draw("SAME");
    }
}
void DrawCut()
{
    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    TH2D *his_angle[2][40];
    TFile *fin = new TFile("FindCut.root");
    for (int i = 0; i < 40; i++)
    {
        his_angle[0][i] = (TH2D *)fin->Get(Form("his_angleXY_%d", i));
        his_angle[1][i] = (TH2D *)fin->Get(Form("his_angleZY_%d", i));
    }

    ifstream intxt("cut_EdetvsAlab.txt");
    double cut[40][25][2][2][2]; //[det][energy][XY/ZY][cut1/cut2][min/max]
    int Det, Energy;
    double cutXY[4], cutZY[4];
    while (intxt >> Det >> Energy 
            >> cutXY[0] >> cutXY[1] >> cutXY[2] >> cutXY[3]
            >> cutZY[0] >> cutZY[1] >> cutZY[2] >> cutZY[3])
    {
        for (int i = 0; i < 4; i++)
        {
            cut[Det][Energy][0][(int)i / 2][i % 2] = cutXY[i];
            cut[Det][Energy][1][(int)i / 2][i % 2] = cutZY[i];
        }
        if (Det >= 24 && cutZY[2] + cutZY[3] == 0 && cutZY[0] > 90)
        {
            cut[Det][Energy][1][0][0] = cutZY[2];
            cut[Det][Energy][1][0][1] = cutZY[3];
            cut[Det][Energy][1][1][0] = cutZY[0];
            cut[Det][Energy][1][1][1] = cutZY[1];
        }
    }

    TGraph *g[2][40];
    for (int i = 0; i < 2; i++) for (int j = 0; j < 40; j++) g[i][j] = new TGraph();
    TCanvas *cvs_draw[2];
    for (int i = 0; i < 2; i++)
    {
        cvs_draw[i] = new TCanvas(Form("cvs_draw%d", i), Form("cvs_draw%d", i), 2000, 1200);
        cvs_draw[i]->Divide(8, 5);
    }
    for (int iDet = 0; iDet < 40; iDet++)
    {
        int mapdraw = -1;
        for (int i=0; i<40; i++) if(drawmap[i]==iDet) { mapdraw=i; break; }
        if (mapdraw == -1) continue;

        int nSecond= 0;
        for (int iEne = 0; iEne < 25; iEne++) if (cut[iDet][iEne][1][1][0] > 0) nSecond++;

        int iSecond = 0;
        for (int iEne = 0; iEne < 25; iEne++)
        {
            for (int iXZ : {0,1})
            {
                g[iXZ][iDet]->SetPoint(iEne,      cut[iDet][iEne][iXZ][0][0], iEne+0.5);
                g[iXZ][iDet]->SetPoint(49 - iEne, cut[iDet][iEne][iXZ][0][1], iEne+0.5);
            }
            if (nSecond > 0 && cut[iDet][iEne][1][1][0] > 0)
            {
                g[1][iDet]->SetPoint(50 + iSecond,               cut[iDet][iEne][1][1][0], iEne+0.5);
                g[1][iDet]->SetPoint(50 + 2*nSecond-1 - iSecond, cut[iDet][iEne][1][1][1], iEne+0.5);
                iSecond++;
            }
        }

        for (int iXZ : {0, 1})
        {
            cvs_draw[iXZ]->cd(mapdraw+1);
            his_angle[iXZ][iDet] -> Draw("colz");
            g[iXZ][iDet] -> SetLineColor(kRed);
            g[iXZ][iDet] -> Draw("SAME");
        }
    }
}

void MakeCut2()
{
    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    TH2D *his_angle[40];
    TFile *fin = new TFile("FindCut_14Oaa.root");
    for(int i=0; i<40; i++)
        his_angle[i] = (TH2D*) fin -> Get(Form("his_EdvsSZEcm_%d",i));

    double thres = 1;
    double AMinMax[40][25][2]; //[det][energy][min/max]
    TGraph *g[40];
    for (int i = 0; i < 40; i++) g[i] = new TGraph();
    for (int iDet = 0; iDet < 40; iDet++)
    {
        for (int i = 0; i < 25; i++)
            for (int j = 0; j < 2; j++)
                AMinMax[iDet][i][j] = 0;

        int Rmin = 1;
        int Rmax = 180;
        for (int Edet = 1; Edet <= 25; Edet++)
        {
            double AMin = -1;
            double AMax = -1;
            for (int Alab = Rmin; Alab <= Rmax; Alab++)
            {
                if (his_angle[iDet]->GetBinContent(Alab, Edet) >= thres)
                {
                    if (AMin == -1)
                        AMin = his_angle[iDet]->GetXaxis()->GetBinCenter(Alab);
                    AMax = his_angle[iDet]->GetXaxis()->GetBinCenter(Alab);
                }
            }
            if (AMin != -1 && AMax != -1)
            {
                AMinMax[iDet][Edet - 1][0] = AMin - 0.1;
                AMinMax[iDet][Edet - 1][1] = AMax + 0.1;
            }

            if (AMinMax[iDet][Edet - 1][0] == 0)
            {
                AMinMax[iDet][Edet - 1][0] = 0;
                AMinMax[iDet][Edet - 1][1] = AMinMax[iDet][Edet - 2][1];
            }
        }
    }

    ofstream outtxt("cut.txt");
    for (int iDet = 0; iDet < 40; iDet++) {
        for (int iEne = 0; iEne < 25; iEne++) {
            outtxt << iDet << "\t" << iEne << "\t";
            outtxt << AMinMax[iDet][iEne][0] << "\t" << AMinMax[iDet][iEne][1] << "\t";
            g[iDet]->SetPoint(iEne, AMinMax[iDet][iEne][0] - 0.2, iEne);
            g[iDet]->SetPoint(49 - iEne, AMinMax[iDet][iEne][1] + 0.2, iEne);
            outtxt << endl;
        }
    }

    TCanvas *cvs_make0 = new TCanvas("cvs_make0", "cvs_make0", 2000, 1200); cvs_make0->Divide(8, 5);
    for (int i = 0; i < 40; i++)
    {
        if(his_angle[drawmap[i]]==nullptr) continue;
        if(drawmap[i]==-1) continue;
        cvs_make0->cd(i + 1);
        his_angle[drawmap[i]]->Draw("colz");
        g[drawmap[i]] -> SetLineColor(kRed);
        g[drawmap[i]] -> Draw("SAME");
    }
}
void DrawCut2()
{
    int detmap[40];
    for (int i = 0; i < 40; i++) switch (i) {
        case 0  ... 9 : detmap[i] = i;                  break;
        case 10 ... 16: detmap[i] = i - 10 + 11;        break;
        case 17 ... 23: detmap[i] = i - 10 + 21 - 7;    break;
        case 24 ... 31: detmap[i] = i - 10 + 101 - 14;  break;
        case 32 ... 39: detmap[i] = i - 10 + 201 - 22;  break;
        default: break;
    }

    TH2D *his_angle[2][40];
    TFile *fin = new TFile("FindCut_Ecm_14Oap.root");
    for (int i = 0; i < 40; i++)
    {
        his_angle[0][i] = (TH2D *)fin->Get(Form("his_EdvsSZEcm_%d", i));
        his_angle[1][i] = (TH2D *)fin->Get(Form("his_EdvsTZEcm_%d", i));
    }

    ifstream intxt("cutEcm_14Oap.txt");
    double cut[40][25][2]; //[det][energy][min/max]
    int Det, Energy;
    double tcut[4];
    while (intxt >> Det >> Energy >> tcut[0] >> tcut[1])
    {
        cut[Det][Energy][0] = tcut[0];
        cut[Det][Energy][1] = tcut[1];
    }

    TGraph *g[2][40];
    for (int i = 0; i < 2; i++) for (int j = 0; j < 40; j++) g[i][j] = new TGraph();
    TCanvas *cvs_draw[2];
    for (int i = 0; i < 2; i++)
    {
        cvs_draw[i] = new TCanvas(Form("cvs_draw%d", i), Form("cvs_draw%d", i), 2000, 1200);
        cvs_draw[i]->Divide(8, 5);
    }
    for (int iDet = 0; iDet < 40; iDet++)
    {
        int mapdraw = -1;
        for (int i=0; i<40; i++) if(drawmap[i]==iDet) { mapdraw=i; break; }
        if (mapdraw == -1) continue;

        for (int iEne = 0; iEne < 25; iEne++)
        {
            for (int iXZ : {0,1})
            {
                g[iXZ][iDet]->SetPoint(iEne,      cut[iDet][iEne][0], iEne+0.5);
                g[iXZ][iDet]->SetPoint(49 - iEne, cut[iDet][iEne][1], iEne+0.5);
            }
        }

        for (int iXZ : {0, 1})
        {
            cvs_draw[iXZ]->cd(mapdraw+1);
            his_angle[iXZ][iDet] -> Draw("colz");
            g[iXZ][iDet] -> SetLineColor(kRed);
            g[iXZ][iDet] -> Draw("SAME");
        }
    }
}