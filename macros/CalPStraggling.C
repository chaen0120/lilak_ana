int iProton = 0;
int iAlpha = 1;
int i14O = 2;
int i14N = 3;
double Mass[4] = {1.00782, 4.00260, 14.00860, 14.00307};
double CMtoLab(double ene, double angle, int light, int beam);
double straggling(double Ep, double dist) { return 0.000628384*sqrt(dist)/(1-0.253199/Ep); }

void CalPStraggling()
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

    //TH2D *his = new TH2D("his","his;distance;Ep",100,0,400,100,0,20);
    double min[100], max[100];
    for (int i=0; i<100; i++) { min[i] = 9999; max[i] = -9999; }
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

        //for (int section = 0; section < 9; section++)
        for (int section : {4} )
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

            for (double Ecm = 8.6; Ecm > 0; Ecm -= 0.01)
            {
                auto vertex = TVector3(0, 0, fEZ->Eval(Ecm));
                auto length = (si - vertex).Mag();
                auto angle = si.Angle(vertex);

                auto Elab = CMtoLab(Ecm, angle, iProton, i14O);
                //his->Fill(length, Elab);
                //if (Ecm==3) cout <<Ecm << "\t"<< length << "\t" << Elab << endl;
                auto st = straggling(Elab, length);
                if (st < min[(int)(Ecm*10)]) min[(int)(Ecm*10)] = st;
                if (st > max[(int)(Ecm*10)]) max[(int)(Ecm*10)] = st;
            }
        }
    }

    gStyle->SetPalette(kBird);
    //his->Draw("colz");
    auto *g = new TGraphErrors();
    for (int i=0; i<86; i++)
    {
        g->SetPoint(g->GetN(), (double)(i/10.), (min[i]+max[i])/2);
        g->SetPointError(g->GetN()-1, 0.05, (min[i]+max[i])/2 - min[i]);
    }
    g->SetFillColor(kRed);
    g->SetFillStyle(3004);
    g->Draw("Ae3");
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
        Edet = ene * (Mbeam + Mtarget) / (4 * Mbeam * cos * cos);
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
