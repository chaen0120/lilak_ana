TF1 *fZE[4], *fEZ[4];
double OcutZ[5] = {0, 262.011-10.413+24, 305.331-10.413+24, 339.753-10.413+24, 367-10.413+24};
double OcutE[5] = {0.5*0.22222839, 1.51387*0.22222839, 3.80064*0.22222839, 10.5704*0.22222839, 40*0.22222839};

double ZtoE(double dist)
{
    if (dist < OcutZ[0] || dist >= OcutZ[4]) return -1;
    for (int i=0; i<4; i++)
        if (dist < OcutZ[i+1]) return fZE[i]->Eval(dist);
    return -1;
}
double EtoZ(double energy) 
{
    if (energy < OcutE[0] || energy >= OcutE[4]) return -1;
    for (int i=0; i<4; i++)
        if (energy < OcutE[i+1]) return fEZ[i]->Eval(energy);
    return -1;
}

void GetSolidAngle(bool Is14N = false)
{
    double detX = -124;
    double detY = 35;
    double detZ = 470;
    auto si = TVector3(detX, detY, detZ);

    for (int i=0; i<4; i++)
    {
        if (!Is14N)
        {
            fZE[i] = new TF1(Form("fZE_%d",i), "0.22222839*([0]+[1]*(x-10.413)+[2]*(x-10.413)^2+[3]*(x-10.413)^3+[4]*(x-10.413)^4)", 0, 375);
            fEZ[i] = new TF1(Form("fEZ_%d",i), "[0]+[1]*(4.4998751*x)+[2]*(4.4998751*x)^2+[3]*(4.4998751*x)^3+[4]*(4.4998751*x)^4 - 10.413+24", 0, 9);
        }
        else
        {
            fZE[i] = new TF1(Form("fZE_%d",i), "0.22229664*([0]+[1]*(x-2.86574-24)+[2]*(x-2.86574-24)^2+[3]*(x-2.86574-24)^3+[4]*(x-2.86574-24)^4)", 0, 375);
            fEZ[i] = new TF1(Form("fEZ_%d",i), "[0]+[1]*(4.4984935*x)+[2]*(4.4984935*x)^2+[3]*(4.4984935*x)^3+[4]*(4.4984935*x)^4 - 2.86574+24", 0, 9);
        }
    }
    if (!Is14N)
    {
        fZE[0]->SetParameters(39.9981, -0.0843282, -6.08828e-05, 3.27748e-08, -3.45099e-10);
        fZE[1]->SetParameters(93.5903, -0.861374, 0.00415962, -1.01406e-05, 8.83304e-09);
        fZE[2]->SetParameters(3452.16, -44.5551, 0.217312, -0.000472261, 3.84527e-07);
        fZE[3]->SetParameters(-1142.37, 18.9009, -0.103531, 0.000235856, -1.93395e-07);
        fEZ[0]->SetParameters(365.486, -35.0153, 21.5374, -8.47836, 1.37061);
        fEZ[1]->SetParameters(358.264, -16.9515, 4.05397, -0.733182, 0.0515715);
        fEZ[2]->SetParameters(349.372, -7.6718, 0.274643, -0.0260368, 0.000736422);
        fEZ[3]->SetParameters(344.024, -5.50135, -0.0752176, -4.97182e-05, -1.72572e-07);
    }
    else
    {
        fZE[0]->SetParameters(30, -0.0821361, -5.19621e-05, -1.29567e-07, -4.08549e-11);
        fZE[1]->SetParameters(31.0524, -0.107965, 0.000182799, -1.0634e-06, 1.32497e-09);
        fZE[2]->SetParameters(658.956, -10.6729, 0.0668874, -0.000188385, 1.98762e-07);
        fZE[3]->SetParameters(-5113.94, 74.558, -0.404784, 0.000971192, -8.69766e-07);
        fEZ[0]->SetParameters(298.599, -36.74, 23.9467, -10.3304, 1.81655);
        fEZ[1]->SetParameters(291.262, -17.2617, 3.73715, -0.691054, 0.050257);
        fEZ[2]->SetParameters(283.774, -9.15657, 0.322428, -0.0294115, 0.000814409);
        fEZ[3]->SetParameters(279.257, -7.09869, -0.0334126, -0.00193614, 1.98111e-05);
        OcutZ[0] = 0      ;
        OcutZ[1] = 157.068-2.86574+24;
        OcutZ[2] = 229.519-2.86574+24;
        OcutZ[3] = 261.217-2.86574+24;
        OcutZ[4] = 291    -2.86574+24;
        OcutE[0] = 0.5    *0.22229664;
        OcutE[1] = 1.59147*0.22229664;
        OcutE[2] = 3.57331*0.22229664;
        OcutE[3] = 10.9892*0.22229664;
        OcutE[4] = 30     *0.22229664;
    }

    auto *his = new TH1D("his", "Solid Angle", 101, 0, 10);
    for (double iE=0; iE<8; iE+=0.1)
    {
        double iZ = EtoZ(iE);
        auto beam = TVector3(0, 0, iZ);
        auto track = si-beam;
        auto dist = track.Mag();
        auto angle = track.Angle(beam);
        cout << iE << " "<< TMath::RadToDeg()*angle << " " << iZ << endl;
        auto SA = 2500/2 * 4 * pow(cos(angle),2)/(pow(dist, 2));
        his->Fill(iE, SA);
    }
    his->Draw("hist");
}
