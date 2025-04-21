LKRun *run;
TexAT2* tt;
void draw()
{
    TFile *fin = new TFile("/mnt/CRIBdisk/o14apf17/lilak_reco/proton_checkBox_newGeo/texat_main_highI3_reco.root");
    TTree *tree = (TTree*) fin->Get("event");

}

void DrawBeam()
{
    run = new LKRun();
    run->AddInputFile("/mnt/CRIBdisk/o14apf17/lilak_reco/proton_checkBox_newGeo/texat_main_highI3_reco.root");
    tt = new TexAT2();
    run->AddDetector(tt);

    run->Init();
    draw();
}