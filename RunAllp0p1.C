#include "macro_p0p1.C"
#include "macros/drawGamow.C"
void DrawAngleVertex();
void DrawSFactor();

void RunAllp0p1(bool VertexEcm = false)
{
    if (VertexEcm) { // vertex mode
        IsVertexEcm = true;
        gROOT->SetBatch(true);
        macro_p0p1("14Oap");
        macro_p0p1("14OCO2p");
        macro_p0p1("14Oap_14N");
        macro_p0p1("14OCO2p_14N");
        gROOT->SetBatch(false);
        macro_14O();
        gSystem->Exec("mv totxs/totxs_14Oap.txt totxs/totxs_14Oap_v.txt");
        gSystem->Exec("mv totxs/totxs_14Oap_sysErr.txt totxs/totxs_14Oap_sysErr_v.txt");
    }
    else { // angle mode
        IsVertexEcm = false;
        gROOT->SetBatch(true);
        macro_p0p1("14Oap");
        macro_p0p1("14OCO2p");
        macro_p0p1("14Oap_14N");
        macro_p0p1("14OCO2p_14N");
        gROOT->SetBatch(false);
        macro_14O();
        gSystem->Exec("mv totxs/totxs_14Oap.txt totxs/totxs_14Oap_a.txt");
        gSystem->Exec("mv totxs/totxs_14Oap_sysErr.txt totxs/totxs_14Oap_sysErr_a.txt");
    }
}

void DrawAngleVertex()
{
    ifstream fin;
    double e, xs, de, dxs;
    fin.open("totxs/totxs_14Oap_a.txt");
    vector<double> ea, xsa, dea, dxsa;
    while (fin >> e >> xs >> de >> dxs)
    {
        if (e < 0.5 || e > 3.6) continue;
        ea.push_back(e);
        xsa.push_back(xs);
        dea.push_back(de);
        dxsa.push_back(dxs);
    }
    fin.close();
    fin.open("totxs/totxs_14Oap_sysErr_a.txt");
    vector<double> dsxsa;
    while (fin >> e >> xs >> de >> dxs)
    {
        if (e < 0.5 || e > 3.6) continue;
        dsxsa.push_back(dxs);
    }
    fin.close();
    fin.open("totxs/totxs_14Oap_v.txt");
    vector<double> ev, xsv, dev, dxsv;
    while (fin >> e >> xs >> de >> dxs)
    {
        if (e < 0.5 || e > 3.6) continue;
        ev.push_back(e);
        xsv.push_back(xs);
        dev.push_back(de);
        dxsv.push_back(dxs);
    }
    fin.close();
    fin.open("totxs/totxs_14Oap_sysErr_v.txt");
    vector<double> dsxsv;
    while (fin >> e >> xs >> de >> dxs)
    {
        if (e < 0.5 || e > 3.6) continue;
        dsxsv.push_back(dxs);
    }
    fin.close();

    auto *gAngle = new TGraphAsymmErrors();
    for (int i=0; i<ea.size(); i++)
    {
        gAngle->SetPoint(i,ea[i], xsa[i]);
        gAngle->SetPointError(i,dea[i],dea[i], dxsa[i]+dsxsa[i]*2,dxsa[i]);
    }
    auto *gVertex= new TGraphAsymmErrors();
    for (int i=0; i<ev.size(); i++)
    {
        gVertex->SetPoint(i,ev[i], xsv[i]);
        gVertex->SetPointError(i,dev[i],dev[i], dxsv[i]+dsxsv[i]*2,dxsv[i]);
    }
    //auto *gAngle = new TGraphErrors(ea.size(), ea.data(), xsa.data(), dea.data(), dxsa.data());
    //auto *gVertex = new TGraphErrors(ev.size(), ev.data(), xsv.data(), dev.data(), dxsv.data());

    gAngle->GetXaxis()->SetLimits(0.5, 3.6);
    gAngle->GetXaxis()->SetTitle("E_{cm} (MeV)");
    gAngle->GetYaxis()->SetTitle("Total Cross Section (mb)");
    gAngle->SetMarkerStyle(20);
    gAngle->SetMarkerSize(0.5);
    gAngle->SetMarkerColor(kRed);
    gAngle->SetLineColor(kRed);
    gAngle->SetFillColor(kRed);
    gAngle->SetFillStyle(3004);
    gAngle->SetLineWidth(2);
    gAngle->GetXaxis()->SetLimits(0.5, 3.6);
    gAngle->GetXaxis()->SetTitle("E_{cm} [MeV]");
    gAngle->GetYaxis()->SetTitle("Total Cross Section [mb]");
    gVertex->SetMarkerStyle(20);
    gVertex->SetMarkerSize(0.5);
    gVertex->SetMarkerColor(kBlue);
    gVertex->SetLineColor(kBlue);
    gVertex->SetFillColor(kBlue);
    gVertex->SetFillStyle(3005);
    gVertex->SetLineWidth(2);

    gVertex->GetXaxis()->SetLabelSize(0.04);
    gVertex->GetXaxis()->SetTitleSize(0.05);
    gVertex->GetXaxis()->SetTitleOffset(0.5);
    gVertex->GetXaxis()->CenterTitle();
    gVertex->GetYaxis()->SetLabelSize(0.04);
    gVertex->GetYaxis()->SetTitleSize(0.05);
    gVertex->GetYaxis()->SetTitleOffset(0.9);
    gVertex->GetYaxis()->CenterTitle();

    //gAngle->Draw("3 same");
    gVertex->Draw("AP");
    gAngle->Draw("Psame");
    gVertex->Draw("e3same");
    gAngle->Draw("e3same");
    //gVertex->Draw("3 same");

    auto *leg = new TLegend(0.15, 0.6, 0.3, 0.88);
    leg->AddEntry(gAngle, "Angle", "pf");
    leg->AddEntry(gVertex, "Vertex", "pf");
    leg->Draw();
    gPad->SetLogy();
    gPad->SetGrid();
}

void DrawSFactor()
{
    ifstream fin;
    double e, xs, de, dxs;
    fin.open("totxs/totxs_14Oap_a.txt");
    vector<double> ea, xsa, dea, dxsa;
    while (fin >> e >> xs >> de >> dxs)
    {
        if (e < 0.5 || e > 3.6) continue;
        ea.push_back(e);
        xsa.push_back(xs);
        dea.push_back(de);
        dxsa.push_back(dxs);
    }
    fin.close();

    auto *gSFactor = new TGraphErrors();
    auto Z1 = 8;
    auto Z2 = 2;
    auto mu = 2897.55; // MeV/c^2
    for (int i=0; i<ea.size(); i++)
    {
        auto Ecm = ea[i];
        auto XS = xsa[i];
        auto eta = Sommerfeld(Ecm, Z1, Z2, mu);
        auto sommer = exp(2 * TMath::Pi() * eta);
        auto Sfactor = XS * Ecm * sommer;

        gSFactor->SetPoint(i, Ecm, Sfactor);
    }
    gSFactor->Draw("APLE*");
}
