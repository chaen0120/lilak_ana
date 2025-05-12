void cvs()
{
    auto copyGraphs = [](TPad *srcPad, Int_t newPadNum, TCanvas *destCanvas) 
    {
        destCanvas->cd(newPadNum);

        TIter next(srcPad->GetListOfPrimitives());
        TObject *obj;
        TList *legendList = new TList();
        bool isFirst = true;

        while ((obj = next())) 
        {
            if (obj->InheritsFrom(TGraph::Class())) 
            {
                TObject *clone = obj->Clone();

                TString drawOpt = "L"; 
                if (obj->InheritsFrom(TGraphErrors::Class())) { drawOpt = "E3"; }

                if (isFirst) { clone->Draw("A" + drawOpt); isFirst = false; } 
                else { clone->Draw(drawOpt + " same"); }
            }
            else if (obj->InheritsFrom(TLegend::Class())) {
                legendList->Add(obj);
            }
        }
        TIter legIter(legendList);
        while ((obj = legIter())) { TLegend *leg = (TLegend*)obj->Clone(); leg->Draw(); }
    };
    auto *fin = new TFile("cvs.root");
    auto *cvs = (TCanvas*) fin->Get("cvs");

    auto *pad1 = (TPad*) cvs->GetPad(1);
    auto *pad2 = (TPad*) cvs->GetPad(2);

    auto *cvsnew = new TCanvas("cvsnew","cvsnew",1000,1500);
    cvsnew->Divide(1,2);
    cvsnew->cd(1); gPad->SetLogy();
    cvsnew->cd(2); gPad->SetLogx(); gPad->SetLogy();
    copyGraphs(pad1, 1, cvsnew);
    copyGraphs(pad2, 2, cvsnew);
}
