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
void setFRStyle(TH1D* hist, const int& color, int MarkerStyle, int LineStyle, double MarkerSize, double LineWidth);
TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator);

void estimateFR2() {

	TString filename = "hist";

  TFile* f[13];

  for(int i=0; i<13; i++) {
    string str = Form("%d",i);
    if( i == 10 ) str = "10_";
    else if( i == 11 ) str = "20_";
    else if( i == 12 ) str = "100_";
    f[i] = new TFile("root_box/"+filename+str+".root","READ");
  }

    TH1D* denominator_pt_barrel[13];
    TH1D* denominator_pt_endcap[13];

    TH1D* numerator_pt_barrel[13];
    TH1D* numerator_pt_endcap[13];

    for(int i=0;i<13;i++) {
        // barrel
        denominator_pt_barrel[i] = (TH1D*)f[i]->Get("denominator_pt_barrel")->Clone("denominator_pt_barrel"+TString::Itoa(i,10));
        numerator_pt_barrel[i] = (TH1D*)f[i]->Get("numerator_pt_barrel")->Clone("numerator_pt_barrel"+TString::Itoa(i,10));

        // endcap
        denominator_pt_endcap[i] = (TH1D*)f[i]->Get("denominator_pt_endcap")->Clone("denominator_pt_endcap"+TString::Itoa(i,10));
        numerator_pt_endcap[i] = (TH1D*)f[i]->Get("numerator_pt_endcap")->Clone("numerator_pt_endcap"+TString::Itoa(i,10));

        const int Nbins = 10;
        Double_t xbins[Nbins+1] = {5,7,10,15,20,30,40,50,70,100,500};
        denominator_pt_barrel[i] = (TH1D*)denominator_pt_barrel[i]->Rebin(Nbins,"",xbins);
        numerator_pt_barrel[i] = (TH1D*)numerator_pt_barrel[i]->Rebin(Nbins,"",xbins);
        denominator_pt_endcap[i] = (TH1D*)denominator_pt_endcap[i]->Rebin(Nbins,"",xbins);
        numerator_pt_endcap[i] = (TH1D*)numerator_pt_endcap[i]->Rebin(Nbins,"",xbins);

        // histogram shape
        if(i==12) {
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
    double nEvts[13] = {6669990,995200,999800,1992440,250000,246800,249237,1000000,6933090,6852830,80780990.0,77061150,1.0}; // ws
    double xsec[13] = {1.256,0.00159,0.01212,0.2529,0.1651,0.05565,0.01398,22.82,76.18,76.18,5765.4,815.9,1.0}; // xsec
    double norm_xsec[13];
    double norm_fit_barrel[13];
    double norm_fit_endcap[13];
    double lumi = 35900.0;

    for(int i=0; i<12; i++) {
        norm_xsec[i] = lumi*xsec[i]/nEvts[i];
    }

    for(int i=0; i<12; i++) {
        // ---------------------------------------------------------- //
        // num (tight) / den (loose)                                  //
        // num part : xsec normalized -> num_data * num_QCD / num_MC  //
        // den part : xsec normalized -> den_data - den_4L - den_3L   //
        // ---------------------------------------------------------- //

        denominator_pt_barrel[i]->Scale(norm_xsec[i]);
        numerator_pt_barrel[i]->Scale(norm_xsec[i]);

        denominator_pt_endcap[i]->Scale(norm_xsec[i]);
        numerator_pt_endcap[i]->Scale(norm_xsec[i]);
    }

		cout<<endl<<"Subtraction Barrel"<<endl;
    TH1D* FR_barrel = (TH1D*)FRByTemplate(numerator_pt_barrel, denominator_pt_barrel); // -> central value
		cout<<endl<<"Subtraction Endcap"<<endl;
    TH1D* FR_endcap = (TH1D*)FRByTemplate(numerator_pt_endcap, denominator_pt_endcap);

		vector<double> nSubtBar, nSubtEnd, eSubtBar, eSubtEnd, RelErrBar, RelErrEnd, RelDiffBar, RelDiffEnd;
		for(int i=1; i<FR_barrel->GetNbinsX()+1; i++) {
			nSubtBar.push_back( FR_barrel->GetBinContent(i) );
			nSubtEnd.push_back( FR_endcap->GetBinContent(i) );
			eSubtBar.push_back( FR_barrel->GetBinError(i) );
			eSubtEnd.push_back( FR_endcap->GetBinError(i) );
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

    TH1D* ptFrame = new TH1D("ptFrame","",/*8,10,500*/495,5,500);
    ptFrame->SetStats(kFALSE);
    ptFrame->GetXaxis()->SetTitle("p_{T}[GeV]");
    ptFrame->GetYaxis()->SetTitle("Fake Rate");

    vector<double> v;
    for(int i=0; i<20; i++) {
      v.push_back(FR_barrel->GetBinContent(i));
      v.push_back(FR_endcap->GetBinContent(i));
    }
    double YRangeMax = *max_element(v.begin(),v.end()) * 1.5;
    double YRangeMin = *min_element(v.begin(),v.end()) * 0.9;
  	YRangeMax = 0.9;

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

    setFRStyle(FR_barrel, kBlue, 20, 1, 1, 1);
    setFRStyle(FR_endcap, kRed, 21, 1, 1, 1);

		canv->SetGrid();
    canv->cd();

    TLegend* legend2 = new TLegend(.50,.75,.9,.9);
    legend2->AddEntry(FR_barrel,"Barrel");
    legend2->AddEntry(FR_endcap,"Endcap");
    legend2->SetBorderSize(0);
		legend2->SetFillStyle(0);

    ptFrame->Draw();
    CMS_lumi(canv,4,11);
    canv->Update();
		canv->SetGrid();
    FR_barrel->Draw("E1PSAME");
    FR_endcap->Draw("E1PSAME");
		canv->RedrawAxis();
		canv->SetGrid();
    legend2->Draw("SAME");
    canv->Print("./FR_2.png");

    TFile* g = new TFile("./FR2.root","RECREATE");
    FR_barrel->Write();
    FR_endcap->Write();
    g->Close();
}

void setDataHist(TH1D* hist) {
    hist->SetLineWidth(2);
    hist->SetMarkerStyle(33);
    hist->SetMarkerSize(3);
    hist->SetStats(kFALSE);
		hist->Sumw2();
}

void setMCHist(TH1D* hist, const int& color) {
    hist->SetFillColor(color+2);
    hist->SetStats(kFALSE);
		hist->Sumw2();
}

void setFRStyle(TH1D* hist, const int& color, int MarkerStyle, int LineStyle = 1, double MarkerSize = 1.0, double LineWidth = 1.0) {
  hist->SetMarkerSize(MarkerSize);
  hist->SetLineColor(color);
  hist->SetLineWidth(LineWidth);
  hist->SetLineStyle(LineStyle);
  hist->SetMarkerStyle(MarkerStyle);
  hist->SetMarkerColor(color);
	hist->Sumw2();
}

TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator) {

    // What is 5? There is no 5 in selectDen~.cc.
    TString name = ( ((TString)(denominator[5]->GetName())).Contains("barrel") ) ? "FR_subtraction_barrel" : "FR_subtraction_endcap";

    TH1D* num = (TH1D*)numerator[12]->Clone(name);
    num->Add(numerator[0],-1.0);
    num->Add(numerator[1],-1.0);
    num->Add(numerator[2],-1.0);
    num->Add(numerator[3],-1.0);
    num->Add(numerator[4],-1.0);
    num->Add(numerator[5],-1.0);
    num->Add(numerator[6],-1.0);
    num->Add(numerator[7],-1.0);

    TH1D* den = (TH1D*)denominator[12]->Clone(name+"_den"); // den_data - den_4l -den_3l
    // qq, gg, higgs, ttZ, WWZ, WZZ, ZZZ, WZ, tbarW, tW, DY, tt, data
    den->Add(denominator[0],-1.0);
    den->Add(denominator[1],-1.0);
    den->Add(denominator[2],-1.0);
    den->Add(denominator[3],-1.0);
    den->Add(denominator[4],-1.0);
    den->Add(denominator[5],-1.0);
    den->Add(denominator[6],-1.0);
    den->Add(denominator[7],-1.0);

    num->Divide(den);


    delete den;
    return num;
}

