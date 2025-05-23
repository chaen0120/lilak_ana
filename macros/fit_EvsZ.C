void fit_EvsZ()
{
    const int nPoint = 93;
    double E[nPoint] = {0.1, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11, 11.1, 11.2};
    double Z[nPoint] = {374, 313.226, 317.897, 241.955, 236.868, 301.041, 272.674, 270.178, 268.544, 286.842, 258.966, 263.964, 268.155, 265.21, 258.037, 253.037, 253.803, 243.676, 243.416, 239.25, 235.929, 232.976, 198.038, 214.398, 211.824, 218.901, 207.6, 209.949, 206.154, 205.158, 203.117, 190.95, 193.045, 194.207, 190.627, 182.472, 185.454, 171.276, 172.283, 162.863, 164.904, 167.125, 164.67, 159.918, 161.206, 155.375, 156.972, 150.22, 148.61, 142.331, 136.142, 143.098, 135.813, 128.175, 132.193, 126.238, 120.305, 126.845, 117.982, 120.423, 122.953, 134.831, 122.99, 136.09, 307.59, 1697.14, 2110.19, 1379.11, 186.819, 120.377, 1077.69, 1107.81, 665.835, 190.171, 617.883, 475.896, 1452.53, 87.2384, 106.938, 87.0014, 93.8322, 813.078, 1144.1, 185.89, 58.9556, 293.294, 321.659, 226.998, 122.473, 231.674, 627.697, 97.0353, 181.601};
    double dZ[nPoint] = {1, 62.8432, 80.5294, 87.4727, 70.8361, 35.1548, 38.9713, 33.7362, 41.581, 29.0942, 48.1029, 35.884, 37.7726, 37.1849, 39.1644, 40.1355, 45.577, 51.6166, 44.2842, 47.1702, 43.5157, 43.6921, 61.8863, 52.1195, 49.0811, 46.8619, 47.3037, 44.6961, 45.1092, 41.9162, 39.5373, 52.688, 46.6196, 41.035, 42.8048, 43.4743, 39.3215, 48.394, 42.2025, 46.4837, 44.6278, 35.3366, 42.1034, 38.674, 36.0031, 37.818, 33.7568, 43.1284, 35.8574, 43.2324, 40.34, 41.5766, 46.7387, 45.3099, 49.7274, 46.5251, 50.3283, 51.3476, 60.5772, 49.239, 51.9782, 61.49, 66.2109, 69.4432, 149.553, 491.537, 493.901, 411.089, 109.878, 83.8549, 416.557, 452.34, 460.762, 103.028, 438.091, 388.556, 338.402, 57.3138, 80.501, 46.5995, 66.3071, 361.646, 397.754, 203.265, 77.7492, 398.454, 359.215, 449.983, 210.444, 215.553, 499.781, 143.57, 423.79};

    int thres = 45;
    /* E(Z) */
    //TGraph *g = new TGraph();
    //for (int i=0; i<nPoint; i++)
    //    if (dZ[i]<thres && Z[i]<600) g->SetPoint(g->GetN(),Z[i],E[i]);
    //g->GetXaxis()->SetLimits(24,380);
    //g->GetYaxis()->SetRangeUser(0,10);
    //TF1 *f = new TF1("f", "pol7", 0, 380);
    //f->FixParameter(2, -0.000362325);
    //f->FixParameter(3,  4.78713e-06);
    //f->FixParameter(4, -3.47837e-08);
    //f->FixParameter(5,  1.36886e-10);
    //f->FixParameter(6, -2.76809e-13);
    //f->FixParameter(7,  2.24699e-16);

    /* Z(E) */
    TGraph *g = new TGraph();
    for (int i=0; i<nPoint; i++)
        if (dZ[i]<thres && Z[i]<600) g->SetPoint(g->GetN(),E[i],Z[i]);
    g->GetXaxis()->SetLimits(0,10);
    g->GetYaxis()->SetRangeUser(24,380);
    TF1 *f = new TF1("f", "pol7", 0, 10);
    f->FixParameter(2,    45.3516);
    f->FixParameter(3,   -22.2208);
    f->FixParameter(4,    5.69877);
    f->FixParameter(5,  -0.805735);
    f->FixParameter(6,  0.0589968);
    f->FixParameter(7, -0.0017465);

    g->Fit(f,"RB");
    g->Draw("AP*");
}
/* results : E(Z)
Chi2                      =      1.47335
NDf                       =           32
Edm                       =  7.65263e-19
NCalls                    =           39
p0                        =       9.9466   +/-   0.141782
p1                        =  -0.00934571   +/-   0.000647674
p2                        = -0.000362325                         (fixed)
p3                        =  4.78713e-06                         (fixed)
p4                        = -3.47837e-08                         (fixed)
p5                        =  1.36886e-10                         (fixed)
p6                        = -2.76809e-13                         (fixed)
p7                        =  2.24699e-16                         (fixed)
*/
/* results : Z(E)
Chi2                      =      2035.63
NDf                       =           32
Edm                       =  1.82761e-16
NCalls                    =           37
p0                        =      360.535   +/-   4.08026
p1                        =     -67.4762   +/-   0.802834
p2                        =      45.3516                         (fixed)
p3                        =     -22.2208                         (fixed)
p4                        =      5.69877                         (fixed)
p5                        =    -0.805735                         (fixed)
p6                        =    0.0589968                         (fixed)
p7                        =   -0.0017465                         (fixed)
*/