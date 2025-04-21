void GaussianKernelSmoothing(TH1 *hist, double res) {
    int nBins = hist->GetNbinsX();  // 히스토그램의 빈 수

    // 새로운 스무딩된 히스토그램 생성
    TH1F *smoothedHist = (TH1F*) hist->Clone("smoothedHist");
    smoothedHist->Reset(); // 빈 값을 모두 0으로 리셋

    // 각 빈에 대해 가우시안 커널을 적용하여 스무딩
    for (int i = 1; i <= nBins; ++i) {
        double binCenter_i = hist->GetBinCenter(i);
        double binContent_i = hist->GetBinContent(i);
        double binWidth_i = hist->GetBinWidth(i);  // 빈 폭 고려
        auto sigma = res*binCenter_i;

        for (int j = 1; j <= nBins; ++j) {
            double binCenter_j = hist->GetBinCenter(j);
            double binWidth_j = hist->GetBinWidth(j);  // 빈 폭 고려
            double distance = binCenter_i - binCenter_j;

            // 가우스 커널 계산 (표준편차 sigma 사용)
            double kernel = TMath::Exp(-0.5 * (distance * distance) / (sigma * sigma)) / (sigma * TMath::Sqrt(2 * TMath::Pi()));

            // 커널 가중치를 적용하여 새로운 히스토그램에 누적
            smoothedHist->AddBinContent(j, binContent_i * kernel * binWidth_i); // 원본 빈 크기 고려
        }
    }

    // 스무딩된 히스토그램 그리기
	ofstream fout("AZURE_resolution.dat");
	for (int i=0; i<nBins; i++) fout << smoothedHist->GetBinLowEdge(i) << "\t" << smoothedHist->GetBinContent(i) << endl;
    smoothedHist->SetLineColor(kRed);  // 스무딩된 히스토그램 색상 설정
    smoothedHist->SetLineWidth(2);     // 선 두께 설정
    smoothedHist->Draw("SAME");        // 기존 히스토그램 위에 덧그리기
}

void applyResolution(double detRes=0.05)
{
    //detRet [MeV]
    vector<double> Ecm; //MeV
    vector<double> crossSection;
    ifstream inputData("AZUREOut_14Oap.txt");

    double temp[5];
    while(inputData >> temp[0] >> temp[1])
    {
        Ecm.push_back(temp[0]);
        crossSection.push_back(temp[1]);
    }

    int nBins = Ecm.size() - 1;
    Ecm.push_back(Ecm[nBins-1]*2-Ecm[nBins-2]);
    TH1D *his = new TH1D("his", "his", nBins, &Ecm[0]);

    for (int i = 0; i < nBins; i++) 
        his->Fill(Ecm[i],crossSection[i]);

    TCanvas *c1 = new TCanvas("c1", "Histogram with Varying Bin Widths and Gaussian Smoothing", 800, 600);
    // 원본 히스토그램 그리기
    his->SetLineColor(kBlue);  // 원본 히스토그램 색상 설정
    his->SetLineWidth(2);      // 선 두께 설정
    his->Draw("HIST");               // 원본 히스토그램 그리기

    // 가우시안 스무딩 적용 (원하는 표준편차 설정, 예: 0.5)
    GaussianKernelSmoothing(his, detRes);

    // 캔버스를 업데이트하고 저장
    c1->Update();
}