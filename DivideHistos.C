#include <stdio.h>
#include <stdlib.h>





TGraphAsymmErrors *  DivideHistos(TH1F * hnum, TH1F* hden , TString xtitle ="",  bool saveplot = true){


  TString fname = (TString) _file0->GetName();
  
  TGraphAsymmErrors * gph= new TGraphAsymmErrors ;
   
  TString hname = hnum->GetName();
  TString hnameden = hden->GetName();
  bool isbarrel = hname.Index("_barrel") > 0;
  bool isendcaps = hname.Index("_endcaps") > 0;
  
  gph->BayesDivide(hnum,hden);
  gph->GetXaxis()->SetTitle(xtitle);gph->GetXaxis()->SetTitleSize(0.045);gph->GetXaxis()->SetTitleOffset(1.);
  gph->GetYaxis()->SetTitle("Efficiency");gph->GetYaxis()->SetTitleSize(0.045);
  
  gph->SetMinimum(0);
  gph->SetMaximum(1.05);
    
  TCanvas *c = new TCanvas ;
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.2);
  c->SetGridx();
  c->SetGridy();
  
  
  gph->Draw("AP");
    
  
  TPaveLabel *labellumi = new TPaveLabel(0.6,0.75,0.9,0.9,"2017, 13 TeV","brNDC");
  labellumi->SetFillColor(0);
  labellumi->SetTextColor(kBlack);
  labellumi->SetFillStyle(0);
  labellumi->SetBorderSize(0);
  labellumi->SetTextSize(0.25);
  labellumi->SetTextAlign(12);
  labellumi->Draw("same");

  
  if(saveplot)c->SaveAs(hname+"_over_"+hnameden+".pdf");
  
  return gph;
  
 
  
}


