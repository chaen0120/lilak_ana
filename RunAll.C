#include "macro_14O.C"
void p0p1();

void RunAll(bool VertexEcm = false)
{
    if (VertexEcm) { // vertex mode
        IsVertexEcm = true;
        gROOT->SetBatch(true);
        macro_ana("14Oap");
        macro_ana("14OCO2p");
        macro_ana("14Oap_14N");
        macro_ana("14OCO2p_14N");
        gROOT->SetBatch(false);
        macro_14O();
    }
    else { // angle mode
        IsVertexEcm = false;
        gROOT->SetBatch(true);
        gSystem->Exec("cd ana; ./p0.sh; cd ..");
        macro_ana("14Oap");
        macro_ana("14OCO2p");
        macro_ana("14Oap_14N");
        macro_ana("14OCO2p_14N");
        macro_14O();
        gSystem->Exec("mv totxs/totxs_14Oap.txt temp_14Oap0.txt");
        gSystem->Exec("mv totxs/totxs_14Oap_sysErr.txt temp_14Oap0_sysErr.txt");

        gSystem->Exec("cd ana; ./p1.sh; cd ..");
        macro_ana("14Oap");
        macro_ana("14OCO2p");
        macro_ana("14Oap_14N");
        macro_ana("14OCO2p_14N");
        macro_14O();
        gSystem->Exec("mv totxs/totxs_14Oap.txt temp_14Oap1.txt");
        gSystem->Exec("mv totxs/totxs_14Oap_sysErr.txt temp_14Oap1_sysErr.txt");
        gROOT->SetBatch(false);

        p0p1();
    }
}
void p0p1()
{
    ifstream fin;
    double e, xs, de, dxs;
    fin.open("temp_14Oap0.txt");
    vector<double> e0, xs0, de0, dxs0;
    while (fin >> e >> xs >> de >> dxs)
    {
        e0.push_back(e);
        xs0.push_back(xs);
        de0.push_back(de);
        dxs0.push_back(dxs);
    }
    fin.close();
    fin.open("temp_14Oap1.txt");
    vector<double> e1, xs1, de1, dxs1;
    while (fin >> e >> xs >> de >> dxs)
    {
        e1.push_back(e);
        xs1.push_back(xs);
        de1.push_back(de);
        dxs1.push_back(dxs);
    }
    fin.close();
    fin.open("temp_14Oap0_sysErr.txt");
    vector<double> dsxs0;
    while (fin >> e >> xs >> de >> dxs)
    {
        dsxs0.push_back(dxs);
    }
    fin.close();
    fin.open("temp_14Oap1_sysErr.txt");
    vector<double> dsxs1;
    while (fin >> e >> xs >> de >> dxs)
    {
        dsxs1.push_back(dxs);
    }
    fin.close();
    gSystem->Exec("rm temp_14Oap*.txt");

    double Ecm[nBins], XS[nBins], dEcm[nBins], dXS[nBins], sXS[nBins], dsXS[nBins];
    for (int iE=0; iE<nBins; iE++)
    {
        Ecm[iE] = e0[iE];
        dEcm[iE] = (de0[iE] + de1[iE]) / 2;
        XS[iE] = (xs0[iE] + xs1[iE]) / 2;
        auto stat = sqrt(dxs0[iE]*dxs0[iE] + dxs1[iE]*dxs1[iE]) / 2;
        auto syst = fabs(xs0[iE] - xs1[iE]) / 2;
        dXS[iE] = sqrt(stat*stat + syst*syst);
        dsXS[iE] = sqrt(dsxs0[iE]*dsxs0[iE] + dsxs1[iE]*dsxs1[iE]) / 2;
        sXS[iE] = XS[iE] - dXS[iE] - dsXS[iE];
    }

    ofstream fout;
    fout.open("totxs/totxs_14Oap.txt");
    for (int iE=0; iE<nBins; iE++)
        fout << Ecm[iE] << "\t" << XS[iE] << "\t" << dEcm[iE] << "\t" << dXS[iE] << endl;
    fout.close();
    fout.open("totxs/totxs_14Oap_sysErr.txt");
    for (int iE=0; iE<nBins; iE++)
        fout << Ecm[iE] << "\t" << sXS[iE] << "\t" << 0 << "\t" << dsXS[iE] << endl;
    fout.close();
    draw_cs();
}