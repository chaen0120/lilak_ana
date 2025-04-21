int iProton = 0;
int iAlpha = 1;
int i14O = 2;
int i14N = 3;
double Mass[4] = {1.00782, 4.00260, 14.00860, 14.00307};

void test()
{
    auto *run = new LKRun();
    auto *tt = new TexAT2();
    auto *ana = new TTAnalysisTask();
    run->AddPar("/home/cens-alpha-00/lilak_old2/texat_ana/temp_reco/config/config_common.mac");
    run->AddPar("/home/cens-alpha-00/lilak_old2/texat_ana/temp_reco/config/config_reco.mac");
    run->AddPar("/home/cens-alpha-00/lilak_old2/texat_ana/temp_reco/config/config_ana.mac");
    run->AddDetector(tt);
    run->Add(ana);
    run->Init();

    cout << ana->RevertEnergyLoss(7.15,368.27,iProton) << endl;
    cout << ana->RevertEnergyLoss(6.543,355.19,iProton) << endl;
}
