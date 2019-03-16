#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "/u/user/joonblee/Headers/tdrstyle.C"
#include "/u/user/joonblee/Headers/CMS_lumi.C"
using namespace std;

void setDataHist(TH1D* hist);
void setMCHist(TH1D* hist, const int& color);
void setFRStyle(TH1D* hist, const int& color, int MarkerStyle, double MarkerSize, double LineWidth);
TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator);

void estimateFR1(int muon_order = 0) {


  TString filename = "histB";

  TFile* f[14];

  f[0] = new TFile("root_box/"+filename+"0.root","READ"); // qq->ZZ
  f[1] = new TFile("root_box/"+filename+"1.root","READ"); // gg->ZZ
  f[2] = new TFile("root_box/"+filename+"2.root","READ"); // higgs
  f[3] = new TFile("root_box/"+filename+"3.root","READ"); // ttZ
  f[4] = new TFile("root_box/"+filename+"4.root","READ"); // WWZ
  f[5] = new TFile("root_box/"+filename+"5.root","READ"); // WZZ
  f[6] = new TFile("root_box/"+filename+"6.root","READ"); // ZZZ
  f[7] = new TFile("root_box/"+filename+"7.root","READ"); // WZ

  f[8] = new TFile("root_box/"+filename+"8.root","READ"); // tW (antitop)
  f[9] = new TFile("root_box/"+filename+"9.root","READ"); // tW (top)
  f[10] = new TFile("root_box/"+filename+"10_.root","READ"); // DY
  f[11] = new TFile("root_box/"+filename+"20_.root","READ"); // tt

  f[12] = new TFile("root_box/"+filename+"10_.root","READ");

  f[13] = new TFile("root_box/"+filename+"100_.root","READ"); // data


	TH1D* FR_endcap[3];
	TH1D* FR_barrel[3];

	int oo = 0;
	for(int iOrder = oo; iOrder<oo+1; iOrder++) {
	  TString MuOrder = ""; TString MuOrder_ = "";
	  if( oo == 1 ) {
	    MuOrder = "l";
	    MuOrder_ = "_lead";
	  }
	  else if( oo == 2 ) {
	    MuOrder = "s";
	    MuOrder_ = "_sub";
	  }

	  TH1D* denominator_pt_barrel[14];
	  TH1D* denominator_pt_endcap[14];
	
	  TH1D* numerator_pt_barrel[14];
	  TH1D* numerator_pt_endcap[14];
	

	  for(int i=0;i<14;i++) {
	    if( i==10 ) { // DY //
	      denominator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_DiMu_pt_barrel")->Clone("denominator_pt_barrel"+TString::Itoa(i,10));
	      numerator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_DiMu_pt_barrel")->Clone("numerator_pt_barrel"+TString::Itoa(i,10));
	      denominator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_DiMu_pt_endcap")->Clone("denominator_pt_endcap"+TString::Itoa(i,10));
	      numerator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_DiMu_pt_endcap")->Clone("numerator_pt_endcap"+TString::Itoa(i,10));
	    }
	    else if( i==12 ) {
	      denominator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_FourMu_pt_barrel")->Clone("denominator_pt_barrel"+TString::Itoa(i,10));
	      numerator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_FourMu_pt_barrel")->Clone("numerator_pt_barrel"+TString::Itoa(i,10));
	      denominator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_FourMu_pt_endcap")->Clone("denominator_pt_endcap"+TString::Itoa(i,10));
	      numerator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_FourMu_pt_endcap")->Clone("numerator_pt_endcap"+TString::Itoa(i,10));
	    }
	    else {
	      // barrel
	      denominator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_pt_barrel")->Clone("denominator_pt_barrel"+TString::Itoa(i,10));
	      numerator_pt_barrel[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_pt_barrel")->Clone("numerator_pt_barrel"+TString::Itoa(i,10));

	      // endcap
	      denominator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"denominator_pt_endcap")->Clone("denominator_pt_endcap"+TString::Itoa(i,10));
	      numerator_pt_endcap[i] = (TH1D*)f[i]->Get(MuOrder+"numerator_pt_endcap")->Clone("numerator_pt_endcap"+TString::Itoa(i,10));
	    }

			const int Nbins = 6;
      Double_t xbins[Nbins+1] = {5,10,20,40,70,100,500};

	    denominator_pt_barrel[i] = (TH1D*)denominator_pt_barrel[i]->Rebin(Nbins,"",xbins);
	    numerator_pt_barrel[i] = (TH1D*)numerator_pt_barrel[i]->Rebin(Nbins,"",xbins);
	    denominator_pt_endcap[i] = (TH1D*)denominator_pt_endcap[i]->Rebin(Nbins,"",xbins);
	    numerator_pt_endcap[i] = (TH1D*)numerator_pt_endcap[i]->Rebin(Nbins,"",xbins);

			denominator_pt_barrel[i]->SetBinContent(Nbins-1, denominator_pt_barrel[i]->GetBinContent(Nbins) + denominator_pt_barrel[i]->GetBinContent(Nbins-1));
     	numerator_pt_barrel[i]->SetBinContent(Nbins-1, numerator_pt_barrel[i]->GetBinContent(Nbins) + numerator_pt_barrel[i]->GetBinContent(Nbins-1));
      denominator_pt_endcap[i]->SetBinContent(Nbins-1, denominator_pt_endcap[i]->GetBinContent(Nbins) + denominator_pt_endcap[i]->GetBinContent(Nbins-1));
      numerator_pt_endcap[i]->SetBinContent(Nbins-1, numerator_pt_endcap[i]->GetBinContent(Nbins) + numerator_pt_endcap[i]->GetBinContent(Nbins-1));

	    // histogram shape
	    if(i==13) {
	      setDataHist( denominator_pt_barrel[i] );
	      setDataHist( numerator_pt_barrel[i] );
	      setDataHist( denominator_pt_endcap[i] );
	      setDataHist( numerator_pt_endcap[i] );
	    }
	    else {
	      setMCHist( denominator_pt_barrel[i], i );
	      setMCHist( numerator_pt_barrel[i], i );
	      setMCHist( denominator_pt_endcap[i], i );
	      setMCHist( numerator_pt_endcap[i], i );
	    }
	  }

	  // qq, gg, higgs, ttZ, WWZ, WZZ, ZZZ, WZ, tbarW, tW, DY, tt, data
	  double nEvts[14] = {6669990,995200,999800,1992440,250000,246800,249237,1000000,6933090,6852830,80780990.0,77061150,80780990.0,1.0}; // ws
	  double xsec[14] = {1.256,0.00159,0.01212,0.2529,0.1651,0.05565,0.01398,22.82,76.18,76.18,5765.4,815.9,5765.4,1.0}; // xsec
	  double norm_xsec[14];
	  double norm_fit_barrel[14];
	  double norm_fit_endcap[14];
	  double lumi = 35900.0;

	  for(int i=0; i<13; i++) {
	    norm_xsec[i] = 0.8*lumi*xsec[i]/nEvts[i];
	  }

	  for(int i=0; i<13; i++) {
		  // ---------------------------------------------------------- //
		  // num (tight) / den (loose)                  //
		  // num part : xsec normalized -> num_data * num_QCD / num_MC  //
		  // den part : xsec normalized -> den_data - den_4L - den_3L   //
		  // ---------------------------------------------------------- //

	    denominator_pt_barrel[i]->Scale(norm_xsec[i]);
	    numerator_pt_barrel[i]->Scale(norm_xsec[i]);

	    denominator_pt_endcap[i]->Scale(norm_xsec[i]);
	    numerator_pt_endcap[i]->Scale(norm_xsec[i]);
	  }
  
	  FR_barrel[iOrder] = (TH1D*)FRByTemplate(numerator_pt_barrel, denominator_pt_barrel); // -> central value
	  FR_endcap[iOrder] = (TH1D*)FRByTemplate(numerator_pt_endcap, denominator_pt_endcap);
	}



  int W = 1200;
  int H = 1200;

  int H_ref = 1200;
  int W_ref = 1200;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.15*H_ref; //0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  lumi_13TeV = "35.9 fb^{-1}";
  writeExtraText = true;
  extraText = "Preliminary";
  drawLogo = false;
  tdrGrid(true);
  lumiTextSize = 0.5;
  cmsTextSize = 0.75;

  TH1D* ptFrame = new TH1D("ptFrame","",495,5,500);
  ptFrame->SetStats(kFALSE);
  ptFrame->GetXaxis()->SetTitle("p_{T}[GeV]");
  ptFrame->GetYaxis()->SetTitle("Fake Rate");

  vector<double> v;
  for(int i=0; i<20; i++) {
    v.push_back(FR_barrel[oo]->GetBinContent(i));
    v.push_back(FR_endcap[oo]->GetBinContent(i));
  }
  double YRangeMax = *max_element(v.begin(),v.end()) * 1.1;
  double YRangeMin = *min_element(v.begin(),v.end()) * 0.9;
  YRangeMax = 1.0;

  ptFrame->SetMinimum(0);
  ptFrame->SetMaximum(YRangeMax/*1.0*/); 
  ptFrame->GetXaxis()->SetRangeUser(5,100);
  ptFrame->GetXaxis()->SetTitleOffset(1.1);
  ptFrame->GetYaxis()->SetTitleOffset(1.1);
  ptFrame->GetXaxis()->SetTitleSize(0.05);
  ptFrame->GetYaxis()->SetTitleSize(0.05);  
  ptFrame->GetXaxis()->SetLabelSize(0.045);
  ptFrame->GetYaxis()->SetLabelSize(0.045); 
    
  ptFrame->GetXaxis()->SetMoreLogLabels(); 

  TCanvas* canv = new TCanvas("canv","",900,900);
  //canv->SetFillColor(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );

	setFRStyle(FR_barrel[oo], kBlue, 20, 1, 1);
	setFRStyle(FR_endcap[oo], kRed, 21, 1, 1);

  canv->SetGrid();
  canv->cd();
  canv->SetGrid();

  double xLeg1 = .7;
  double yLeg1 = .75;
  double xLeg2 = xLeg1 + .19;
  double yLeg2 = yLeg1 + .14;

  TLegend* legend2 = new TLegend(xLeg1,yLeg1,xLeg2,yLeg2);

  legend2->AddEntry(FR_barrel[oo],"Barrel");
  legend2->AddEntry(FR_endcap[oo],"Endcap");
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(0);


  ptFrame->Draw();
  CMS_lumi(canv,4,11);
  canv->Update();
  //canv->RedrawAxis();
  canv->SetGrid();
  //canv->GetFrame()->Draw();
  FR_barrel[oo]->Draw("E1PSAME");

  //ptFrame->Draw();
  //CMS_lumi(canv,4,11);
  //canv->Update();
  //canv->RedrawAxis();
  //canv->GetFrame()->Draw();
  FR_endcap[oo]->Draw("E1PSAME");
  canv->RedrawAxis();
  canv->SetGrid();
  legend2->Draw("SAME");
  canv->Print("./FR1.png");

  TFile* g = new TFile("./FR1.root","RECREATE");
  FR_barrel[oo]->Write();
  FR_endcap[oo]->Write();
  g->Close();
}

void setDataHist(TH1D* hist) {
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(33);
  hist->SetMarkerSize(3);
  hist->SetStats(kFALSE);
}

void setMCHist(TH1D* hist, const int& color) {
  hist->SetFillColor(color+2);
  hist->SetStats(kFALSE);
}

void setFRStyle(TH1D* hist, const int& color, int MarkerStyle, double MarkerSize = 1.0, double LineWidth = 1.0) {
  hist->SetMarkerSize(MarkerSize);
  hist->SetLineColor(color);
  hist->SetLineWidth(LineWidth);
  hist->SetMarkerStyle(MarkerStyle);
  hist->SetMarkerColor(color);
}


TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator) {

  // What is 5? There is no 5 in selectDen~.cc.
  TString name = ( ((TString)(denominator[5]->GetName())).Contains("barrel") ) ? "FR_subtraction_barrel" : "FR_subtraction_endcap";

  TH1D* num = (TH1D*)numerator[13]->Clone(name);
  num->Add(numerator[0],-1.0);
  num->Add(numerator[1],-1.0);
  num->Add(numerator[2],-1.0);
  num->Add(numerator[3],-1.0);
  num->Add(numerator[4],-1.0);
  num->Add(numerator[5],-1.0);
  num->Add(numerator[6],-1.0);
  num->Add(numerator[7],-1.0);
  num->Add(numerator[12],-1.0);

  TH1D* den = (TH1D*)denominator[13]->Clone(name+"_den"); // den_data - den_4l -den_3l
  // qq, gg, higgs, ttZ, WWZ, WZZ, ZZZ, WZ, tbarW, tW, DY, tt, data
  den->Add(denominator[0],-1.0);
  den->Add(denominator[1],-1.0);
  den->Add(denominator[2],-1.0);
  den->Add(denominator[3],-1.0);
  den->Add(denominator[4],-1.0);
  den->Add(denominator[5],-1.0);
  den->Add(denominator[6],-1.0);
  den->Add(denominator[7],-1.0);
  den->Add(denominator[12],-1.0);

  vector<double> eNum, eDen;
  for(int i=0; i<num->GetNbinsX(); i++) {
    eNum.push_back(num->GetBinContent(i));
    eDen.push_back(den->GetBinContent(i));
  }

  num->Divide(den);
  for(int i=0; i<num->GetNbinsX(); i++) {
    num->SetBinError(i, sqrt( eNum[i]*(1-eNum[i]/eDen[i]) ) / eDen[i] );
  }

  delete den;
  return num;
}
