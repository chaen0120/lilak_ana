auto tt = new TexAT2();
auto ta = new TTAnalysisTask();
double CMToLab(double beamene, double angle, TString beam, TString light);
void CalEff()
{
    auto nDet = 40;
    auto nEcm = 100;
    auto minEcm = 0;
    auto maxEcm = 10;
    auto Ebeam = 39.112;

    TGraphErrors *gFo = new TGraphErrors();
    TGraphErrors *gLS = new TGraphErrors();
    TGraphErrors *gRS = new TGraphErrors();
    TGraphErrors *gAB = new TGraphErrors();

    for (int iDet = 0; iDet < 40; iDet++)
    {
        int Det = -1;
        double x = 0, y = 0, z = 0;
        double dx = 0, dy = 0, dz = 0;
        switch (iDet) {
        case  0 ...  9: Det = iDet;            break;
        case 10 ... 16: Det = iDet - 10 + 11;  break;
        case 17 ... 23: Det = iDet - 17 + 21;  break;
        case 24 ... 31: Det = iDet - 24 + 101; break;
        case 32 ... 39: Det = iDet - 32 + 201; break;
        default: break; }

        if (Det<10)
            tt->CAACToGlobalPosition(10100+Det, x,y,z, dx,dy,dz);
        else
        {
            x = tt->GetX6posx(Det, 4);
            y = tt->GetX6posy(Det, 4);
            z = tt->GetX6posz(Det, 4);
            if (iDet<24) { dx = 3.7/2; dy = 75/2;  dz = 40.3/2; }
            else         { dx = 75/2;  dy = 3.7/2; dz = 40.3/2; }
            if (Det == 101 || Det == 204) { dx = 40.3/2; dz = 75/2; }
        }

        switch (iDet) {
        case 0 ... 9:
            gFo->SetPoint(gFo->GetN(), x, y); gFo->SetPointError(gFo->GetN() - 1, dx, dy); break;
        case 10 ... 16:
            gLS->SetPoint(gLS->GetN(), z, y); gLS->SetPointError(gLS->GetN() - 1, dz, dy); break;
        case 17 ... 23:
            gRS->SetPoint(gRS->GetN(), z, y); gRS->SetPointError(gRS->GetN() - 1, dz, dy); break;
        case 24 ... 39:
            gAB->SetPoint(gAB->GetN(), x, z); gAB->SetPointError(gAB->GetN() - 1, dx, dz); break;
        default: break; }

        for (int iEcm = 0; iEcm < nEcm; iEcm++)
        {
            double Ecm = minEcm + (double)(maxEcm - minEcm) / nEcm * iEcm;
        }
    }
    //TCanvas *cvs = new TCanvas("cvs","cvs",2000,2000);
    //cvs->Divide(2,2);
    //cvs->cd(1); gFo->SetFillColor(kBlue); gFo->Draw("A2");
    //cvs->cd(2); gLS->SetFillColor(kBlue); gLS->Draw("A2");
    //cvs->cd(3); gRS->SetFillColor(kBlue); gRS->Draw("A2");
    //cvs->cd(4); gAB->SetFillColor(kBlue); gAB->Draw("A2");
}

void GetEff()
{
    auto run = new LKRun();
    run->AddPar("../config/config_common.mac");
    run->AddPar("../config/config_reco.mac");
    run->AddPar("../config/config_ana.mac");
    run->AddDetector(tt);
    run->Add(ta);
    run->Init();

    CalEff();
}

double CMToLab(double beamene, double angle, TString beam, TString light)
{
    double Edet = -1;
    double Mbeam, Mtarget, Mlight, Mheavy;
    double cos = TMath::Cos(angle);
    double amu = 931.5016;

    //Mbeam = Mass[beam];
    //Mtarget = Mass[kAlpha];
    //if (light == kAlpha) // Scattering
    //    Ecm = beamene * ( Mbeam + Mtarget ) / ( 4 * Mbeam * cos * cos );
    //else // (a,p) reaction
    //{
    //    Mlight    =  1.00782;
    //    if (beam==k14O) Mheavy = 17.00210;
    //    if (beam==k14N) Mheavy = 16.99913;
    //    double m1 = Mbeam   * amu; //MeV 
    //    double m2 = Mtarget * amu; //MeV 
    //    double m3 = Mlight  * amu; //MeV 
    //    double m4 = Mheavy  * amu; //MeV 
    //    double Q  = m1+m2-m3-m4;   //MeV 

    //    double a = m4 - m1;
    //    double b = sqrt(m1*m3*beamene) * cos;
    //    double c = m4*Q - beamene * (m3+m4);

    //    double Esq = ( -b + sqrt(b*b-a*c) ) / a;
    //    double Eb  = Esq * Esq;
    //    Ecm = Eb * Mtarget / (Mbeam+Mtarget);
    //}
    return Edet;
}
