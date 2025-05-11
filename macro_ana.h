const double minEcm = 0;     // [MeV]
const double maxEcm = 10;    // [MeV]
const int nBins = 100;
const double binSize = (maxEcm - minEcm) / nBins; // [MeV]

bool IsVertexEcm = true;
bool IsBackground = false;
bool Is14N = false;

TString beam;
TString reac;
long nBeam;
long nBeam14Oa = 24133085012; //wo run 711, 817
//long nBeam14Oa = 24190716832; //wo run 711, 817
//long nBeam14Oa = 24485489726;
long nBeam14Na = 5420073573;
long nBeam14OC = 4537165179; //wo run 928
//long nBeam14OC = 4553559610;
long nBeam14NC = 2717193980;
long nBeam14Ov = 419158806;

int detmap[40] = {0};
int drawmap[40] = {0, 6, 7, 8, 9, -1, -1, -1,
                   14, 15, 16, 10, 11, 12, 13, -1,
                   23, 22, 21, 20, 19, 18, 17, -1,
                   28, 29, 30, 31, 25, 26, 27, 24,
                   39, 38, 37, 36, 34, 33, 32, 35};
double SiThres[40] = {1.615, 100, 100, 100, 100, 100, 1.423, 1.455, 1.767, 1.511, 
                      1.427, 1.405, 1.453, 1.499, 1.451, 1.449, 1.443, 
                      100, 100, 100, 100, 100, 100, 100, 
                      1.511, 1.407, 1.509, 1.533, 1.887, 100, 1.483, 1.495, 
                      1.575, 1.537, 1.465, 1.525, 1.471, 1.423, 1.491, 1.563};

TGraph *Ecm_14O;
TGraph *Ecm_14N;
TGraph *Ecm_14N_Oaa;

void SetOthers();     
                      
TGraphErrors *g_ZtoE; 
TF1 *f_ZtoE;
void FindZtoE(TH2D *his, int thres_entry, int thres_sigma);

TF1 *fTrackingEff = nullptr;
TF1 *fZE[4], *fEZ[4];
double OcutE[5] = {0.5*0.22222839, 1.51387*0.22222839, 3.80064*0.22222839, 10.5704*0.22222839, 40*0.22222839};
double OcutZ[5] = {24, 262.011-10.413+24, 305.331-10.413+24, 339.753-10.413+24, 367-10.413+24};
double ZtoE(double Z);
double EtoZ(double E);


// macro_xs
vector<vector<double>> effMain(40), effBG(40), effN(40), effNB(40);
vector<double> nEcmMain, nEcmBG, nEcmN, nEcmNB;
void GetEfficiency(vector<vector<double>> &eff, vector<double> &nEcm, TString fileName);
TH1D *ApplyEfficiency(TH1D *his, int iDet, bool isBackground);

double gasDensity = 0.1313*0.9;   // mg/cm^3
double HeMolarMass = 4.0026*1000; // mg/mol
double avo = 6.022E+23;           // #/mol
//TGraphErrors *GetCrossSection(TH1D *hisMain, int iDet);
//TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, int iDet);
//TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, TH1D *hisN, TH1D *hisNBG, int iDet);
TGraphErrors *GetCrossSection(bool IsGoodAna, int iDet);
TGraphErrors *GetCrossSection(TH1D *hisMain, TH1D *hisBG, TH1D *hisN, int iDet);

TGraphErrors *GetAverage(TGraphErrors **gXS, TH1D *hisYield, TH1D *hisErr);
TGraphErrors *GetWeightedAverage(TGraphErrors **gXS, TH1D *hisErr);


// SetCuts
TCutG *cuts[40];
TCutG *cutEA[40];
TCutG *cutEZ[40];
void SetCutP(TCutG **cutP);
void SetCutA(TCutG **cutA);
void SetCutBG(TCutG **cutBG, TString pora, TString type);
void ReadCut(TCutG **cutE, TString AorZ, TString type);
void MoveCut(TCutG **cut, double shift);


// macro_ana
void FindCut();
void Simulation();

bool CutBadDetectors(int det, int strip);

unordered_map<uint64_t, TString> goodEvent;
uint64_t makeKey(int runno, int splitno, int evtno) { return ((uint64_t)runno << 44) | ((uint64_t)splitno << 32) | evtno; }
void GetGoodEvent(TString name);
TString IsGoodEvent(int runno, int splitno, int evtno);
bool GetManualEvent(uint64_t key, double &dEdx, double &pEcm, double &Alab, TVector3 &vertex);

map<int, double> resolution;
void GetResolution();
TF1 *fOStr1, *fOStr2, *fOStr;
double OStraggling(double Z) { return fOStr->Eval(Z); }
double PStraggling(double Ep, double dist) { return 0.000628384*sqrt(dist)/(1-0.253199/Ep); }
double CalculateEpError(double energy, int det, int strip);
double CalculateThetaError(double theta, TVector3 sihit, TVector3 errSihit, TVector3 vertex);
double CalculateEcmError(double energy, double theta, double errEne, double errThe);

map<int, double> liveTime;
void GetLiveTime();

void DrawAllDetector(TH1D **his);
void DrawAllDetector(TH2D **his);
void DrawAllDetector(TGraphErrors **g);
void DrawOneDetector(int Det, bool DrawAngle = true);
void DrawCuts();

void PrintProgress(int iEvt, int nEvt)
{
    if (iEvt % 100 == 0)
    {
        cout << ">> " << iEvt << " / " << nEvt << "\r";
        cout.flush();
    }
    if (iEvt + 1 == nEvt) cout << endl;
}
void SaveBatch(TCanvas *cvs)
{
    if (gROOT->IsBatch())
        cvs->SaveAs(Form("png/%s.png", cvs->GetName()));
}
