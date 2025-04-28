void draw_cs(bool sorryKubono = true)
{
	double xmin = 0.8, xmax = 3.5;
	double Ex = 5.114;

	auto readData = [](const string& filename, vector<Double_t>& x, vector<Double_t>& y, vector<Double_t>& err) 
    {
		ifstream file(filename);
        Double_t tempx, tempy, temperr;
		while (file >> tempx >> tempy >> temperr)
        {
            x.push_back(tempx);
            y.push_back(tempy);
            err.push_back(temperr);
        }
		file.close();
	};

	auto readData2 = [](const string& filename, vector<Double_t>& x, vector<Double_t>& y) 
    {
		ifstream file(filename);
        Double_t tempx, tempy;
		while (file >> tempx >> tempy)
        {
            x.push_back(tempx);
            y.push_back(tempy);
        }
		file.close();
	};

	auto readData4 = [](const string& filename, vector<Double_t>& x, vector<Double_t>& y, vector<Double_t>& xerr, vector<Double_t>& yerr) 
    {
		ifstream file(filename);
        Double_t tempx, tempy, tempxerr, tempyerr;
		while (file >> tempx >> tempy >> tempxerr >> tempyerr)
        {
            x.push_back(tempx);
            y.push_back(tempy);
            xerr.push_back(tempxerr);
            yerr.push_back(tempyerr);
        }
		file.close();
	};


	vector<Double_t> cmE, cs_b, cs_mb;
	readData2("totxs/AZUREOut_14Oap.txt.new", cmE, cs_b);
	for (int i = 0; i < cmE.size(); ++i) cs_mb.push_back(cs_b[i] * 1000);
	TGraph * gAZURE = new TGraph(cmE.size(),cmE.data(),cs_mb.data());
	gAZURE->SetLineColor(kRed);

	vector<Double_t> cmEr, cs_br, cs_mbr;
	readData2("totxs/AZUREOut_14Oap_res.txt.new", cmEr, cs_br);
	for (int i = 0; i < cmEr.size(); ++i) cs_mbr.push_back(cs_br[i] * 1000);
	TGraph * gAZUREr = new TGraph(cmEr.size(),cmEr.data(),cs_mbr.data());
	gAZUREr->SetLineColor(kBlue);

	TGraph * gTalys = new TGraph("totxs/TALYS_14Oap.txt");
	//gTalys->SetMarkerStyle(22);
	//gTalys->SetMarkerColor(6);
	//gTalys->SetMarkerSize(0.85);
	gTalys->SetLineColor(6);
	gTalys->SetLineWidth(2);
	gTalys->SetLineStyle(2);
	//gTalys->GetXaxis()->SetRangeUser(0.9,3);

	vector<Double_t> kcmE, kcs_mb, kcs_err;
	readData("totxs/Kubono_cs.txt", kcmE, kcs_mb, kcs_err);
    vector<Double_t> kcmE_err(kcmE.size(), 0);
	TGraphErrors * gKub = new TGraphErrors(kcmE.size(),kcmE.data(),kcs_mb.data(),kcmE_err.data(),kcs_err.data());
	gKub->SetMarkerStyle(4);
	gKub->SetMarkerColor(4);
	gKub->SetMarkerSize(0.7);

	vector<Double_t> acmE, acs_mb, acs_err;
	readData("totxs/Aram_cs.txt", acmE, acs_mb, acs_err);
    vector<Double_t> acmE_err(acmE.size(), 0);
	TGraphErrors * gAram = new TGraphErrors(acmE.size(),acmE.data(),acs_mb.data(),acmE_err.data(),acs_err.data());
	gAram->SetMarkerStyle(22);
	gAram->SetMarkerSize(0.95);
	gAram->SetMarkerColor(kYellow+3);

	const Int_t l=4;
	Double_t exp_cmE[l]={1.036,1.172,1.936,2.236};
	Double_t exp_cs[l]={0.2397,0.07174,3.427,114.6};
	Double_t exp_err[l]={0.03094,0.016823,0.113,0.644};
	Double_t exp_cmE_err[l]={0,0,0,0};
	TGraphErrors * gPoint = new TGraphErrors(l,exp_cmE,exp_cs,exp_cmE_err,exp_err);
	gPoint->SetMarkerStyle(29);
	gPoint->SetMarkerSize(1.5);

	vector<Double_t> bbcmE_keV, bbcmE_MeV, bbcs_mb, bbcs_err;
	readData("totxs/Blackmon_2001_cs.txt", bbcmE_keV, bbcs_mb, bbcs_err);
    vector<Double_t> bbcmE_err(bbcmE_keV.size(), 0);
	for (int i = 0; i < bbcmE_keV.size(); ++i) bbcmE_MeV.push_back(bbcmE_keV[i] - 1.2);
	TGraphErrors * gBlack = new TGraphErrors(bbcmE_keV.size(),bbcmE_MeV.data(),bbcs_mb.data(),bbcmE_err.data(),bbcs_err.data());
	gBlack->SetMarkerStyle(21);
	gBlack->SetMarkerSize(0.7);
	gBlack->SetMarkerColor(8);


	vector<Double_t> expcmE_MeV, expcs_mb, expcmE_err, expcs_err;
	readData4("totxs/totxs_14Oap.txt", expcmE_MeV, expcs_mb, expcmE_err, expcs_err);
	TGraphErrors * gExp = new TGraphErrors(expcmE_MeV.size(),expcmE_MeV.data(),expcs_mb.data(),expcmE_err.data(),expcs_err.data());
	gExp->SetMarkerStyle(20);
	gExp->SetMarkerSize(0.8);
	gExp->SetMarkerColor(kRed);
	gExp->SetLineColor(kRed);
	gExp->SetFillColor(kRed);
	gExp->SetFillStyle(3004);

	vector<Double_t> syserrcmE_MeV, syserrcs_mb, syserrcmE_err, syserrcs_err;
	readData4("totxs/totxs_14Oap_sysErr.txt", syserrcmE_MeV, syserrcs_mb, syserrcmE_err, syserrcs_err);
	TGraphErrors * gSys = new TGraphErrors(syserrcmE_MeV.size(),syserrcmE_MeV.data(),syserrcs_mb.data(),syserrcmE_err.data(),syserrcs_err.data());
	gSys->SetLineColor(kBlue);
	gSys->SetFillColor(kBlue);
	gSys->SetFillStyle(3004);

	//vector<Double_t> expcmE_MeV1, expcs_mb1, expcmE_err1, expcs_err1;
	//readData4("totxs/totxs_14Oap1.txt", expcmE_MeV1, expcs_mb1, expcmE_err1, expcs_err1);
	//TGraphErrors * gExp1 = new TGraphErrors(expcmE_MeV1.size(),expcmE_MeV1.data(),expcs_mb1.data(),expcmE_err1.data(),expcs_err1.data());
	//gExp1->SetMarkerStyle(20);
	//gExp1->SetMarkerSize(0.8);
	//gExp1->SetMarkerColor(kYellow+1);
	//gExp1->SetLineColor(kYellow+1);
	//gExp1->SetFillColor(kYellow+1);
	//gExp1->SetFillStyle(3004);

	//vector<Double_t> syserrcmE_MeV1, syserrcs_mb1, syserrcmE_err1, syserrcs_err1;
	//readData4("totxs/totxs_14Oap1_sysErr.txt", syserrcmE_MeV1, syserrcs_mb1, syserrcmE_err1, syserrcs_err1);
	//TGraphErrors * gSys1 = new TGraphErrors(syserrcmE_MeV1.size(),syserrcmE_MeV1.data(),syserrcs_mb1.data(),syserrcmE_err1.data(),syserrcs_err1.data());
	//gSys1->SetLineColor(kGreen+1);
	//gSys1->SetFillColor(kGreen+1);
	//gSys1->SetFillStyle(3004);


	TGaxis *axis[2];
	for (int i=0; i<2; i++)
	{
		if (i==0) axis[i] = new TGaxis(xmin,    3, xmax,    3, xmin + Ex, xmax + Ex, 510, "-");
		if (i==1) axis[i] = new TGaxis(xmin, 1000, xmax, 1000, xmin + Ex, xmax + Ex, 510, "-");
		axis[i]->SetName("axis1");
		axis[i]->CenterTitle();
		axis[i]->SetTitle("E_{x}(^{18}Ne) [MeV]");
		axis[i]->SetTitleOffset(0.9);
		axis[i]->SetLabelFont(42);
		axis[i]->SetLabelSize(0.05);
		axis[i]->SetTitleFont(42);
		axis[i]->SetTitleSize(0.05);
	}

	TLegend * led = new TLegend(0.63,0.12,0.86,0.35);
	led->SetTextSize(0.03);
	led->SetFillColor(0);
	//led->AddEntry(gPoint, "Expected result","ep");
	if (!sorryKubono) led->AddEntry(gKub, "S. Kubono et al.","ep");
	led->AddEntry(gAram, "A. Kim et al.","ep");
	led->AddEntry(gBlack, "Blackmon et al.","ep");
	led->AddEntry(gAZURE, "AZURE", "l");
	led->AddEntry(gAZUREr, "AZURE_res", "l");
	led->AddEntry(gTalys, "TALYS 1.95", "l");
	led->AddEntry(gExp, "Expt. Data", "ep");


	for (auto g : { gAZURE, gAZUREr })
	{
		g->GetXaxis()->SetTitle("E_{c.m.}(^{14}O+#alpha) [MeV]");
		g->GetXaxis()->CenterTitle();
		g->GetXaxis()->SetLabelSize(0.050);
		g->GetXaxis()->SetTitleSize(0.05);
		g->GetXaxis()->SetLimits(xmin, xmax);
		g->GetYaxis()->SetTitle("#sigma_{total} [mb]");
		g->GetYaxis()->CenterTitle();
		g->GetYaxis()->SetLabelSize(0.050);
		g->GetYaxis()->SetTitleSize(0.05);
		g->GetYaxis()->SetTitleOffset(1.3);
		if (g==gAZURE) g->GetYaxis()->SetRangeUser(0.0001, 1000);
		else 		   g->GetYaxis()->SetRangeUser(0.0001,  3);
	}
	//TCanvas * c1 = new TCanvas("c1","c1",600,500);
	//c1->cd();
	TCanvas * c1 = new TCanvas("c1","c1",1200,500);
	c1->Divide(2,1);
	for (int i=0; i<2; i++)
	{
		c1->cd(i+1);
		gStyle->SetOptTitle(kFALSE);
		gPad->SetGrid();
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.12);
		gPad->SetTopMargin(0.12);

		if (i==0) { gAZUREr->Draw("AL"); gAZURE->Draw("same"); }
		if (i==1) { gAZURE->Draw("AL"); gAZUREr->Draw("same"); }
		gTalys->Draw("L");
		if (!sorryKubono) gKub->Draw("P");
		gAram->Draw("P");
		// gPoint->Draw("P");
		gBlack->Draw("P");
		gExp->Draw("e3");
		gExp->Draw("P");
		gSys->Draw("e3");
		gSys->Draw("P");
		//gExp1->Draw("e3");
		//gExp1->Draw("P");
		//gSys1->Draw("e3");
		//gSys1->Draw("P");
		axis[i]->Draw();
		if (i == 1) { gPad->SetLogy(); led->Draw(); }
	}

	c1->SaveAs("png/all.png");
}
