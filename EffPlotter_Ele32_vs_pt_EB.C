#include <iostream>
#include <algorithm>
#include "TLatex.h"
#include <iomanip>
#include <vector>
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPostScript.h"
#include <TPaveStats.h>
#include "TLegend.h"
#include <TProfile.h>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"

void EffPlotter_Ele32_vs_pt_EB() {

  TFile *file_sig = new TFile("out.root");
  
  //--Plotting Styles//
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);  
  gStyle->SetPadTopMargin(0.05);   
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.05);
  gStyle->SetOptStat();

  TH1D* den_sig = (TH1D*)file_sig->Get("den_ele_pt_EB");
  TH1D* num_sig = (TH1D*)file_sig->Get("num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EB");

  TEfficiency* pEff_sig = 0;

  TCanvas* my_canvas1 = new TCanvas("canvas","canvas",800,700);
  my_canvas1->cd();

  if(TEfficiency::CheckConsistency(*num_sig,*den_sig)) {
    pEff_sig = new TEfficiency(*num_sig,*den_sig);
    pEff_sig->SetLineColor(kCyan+1);
    pEff_sig->SetMarkerColor(kCyan+1);
    pEff_sig->SetMarkerStyle(21);
    pEff_sig->SetLineWidth(4);
    pEff_sig->SetTitle("Efficiency vs pT;HLT pT [GeV];signal efficiency");
    pEff_sig->Draw();
    gPad->Update();
    auto graph = pEff_sig->GetPaintedGraph(); 
    graph->SetMinimum(0.0);
    graph->SetMaximum(1.01); 
  }


  TLegend *leg_example2 = new TLegend(0.38,0.35,0.78,0.5);
  leg_example2->SetHeader("Ele32","C"); // option "C" allows to center the header
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  leg_example2->AddEntry(pEff_sig, "Barrel", "lp");
  leg_example2->Draw("same");
  
  my_canvas1->SetGrid();

  //  TLine *line = new TLine(32,0.001,32,1.1);
  //line->SetLineColor(kRed+1);
  // line->SetLineStyle(9);
  //line->SetLineWidth(3);
  //line->Draw("same");

  my_canvas1->Update();
  my_canvas1->SaveAs("Ele32_EB_eff_vs_pt.png");
}
