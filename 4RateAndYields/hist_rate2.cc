#define higgs_cs 0.01212
#define higgs_ws 999800.0
#define dy_cs 5795.4 //6025.2 //5765.4
#define dy_ws 80780990.0 
//11442000+12629500+12152400+12664100+12255800+11492100+12438700+12285200+12474000+12221400
#define wz_cs 22.82
#define wz_ws 1000000.0
#define tw_cs 76.18
#define tw_ws 13885920.0 //6933090 + 6952830
#define tbarW_ws 1.0
#define tW_ws 1.0

#define tt_cs 815.9 //679.1 //815.9 //679.1
#define tt_ws 77061150.0 
//7447700+7913880+7988370+7611900+7750470+7526200+7914710+7484670+8085290+7337960
#define qq_cs 1.256
#define qq_ws 6669990.0
#define gg_cs 0.00159
#define gg_ws 995200.0
#define wwz_cs 0.1651
#define wwz_ws 250000.0
#define wzz_cs 0.05565
#define wzz_ws 246800.0
#define zzz_cs 0.01398
#define zzz_ws 249237.0
#define ttz_cs 0.2529
#define ttz_ws 1992440.0

#define zg_cs 
#define zg_ws 9635131

#define lumi 35900.0//35900.0

Double_t nPol(Double_t *x, Double_t *par) {
	// const x0 y0 n
	Double_t t = x[0];
	//return 0.00001*par[0]*t + par[1];
	return 0.00001*(par[0]*t + par[1])*(1.0-par[2]/4000000.0*pow(t,par[3]));
}

Double_t FitPol4(Double_t *x, Double_t *par) {
	Double_t t = x[0];
	return par[0]*pow(t,4) + par[1]*pow(t,3) + par[2]*pow(t,2) + par[3]*t + par[4];
}

Double_t pol4(double t, double A, double B, double C, double D, double E) {
	return A*pow(t,4) + B*pow(t,3) + C*pow(t,2) + D*t + E;
}



void hist_rate2(TString StrMA = "A50", TString Object = "rate", double bin = 1) {

  bool letLogX = false; bool letLogY = false;
	double XRangeMin = 0; double XRangeMax = 2200;
  double YRangeMin = 0; double YRangeMax = 0;

  TString xTitle = StrMA;
  TString yTitle = "events";
  if(letLogY) yTitle = "events (log)";

  TCanvas *c1 = new TCanvas("c1","histogram",800,800);
  c1->SetBottomMargin(0.2);
  c1->cd();

  auto pad1 = new TPad("pad1","",0.01,0.01,0.99,0.99,0);
  pad1->SetLeftMargin(0.1);
  //pad1->SetBottomMargin(0.3);
  if(letLogX) pad1 -> SetLogx();
  if(letLogY) pad1 -> SetLogy();
  pad1->Draw();

  /////////////////////////////////////////////////////////////////////////////
  pad1->cd();

	gStyle->SetOptFit(111111);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.98);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);
	//gStyle->SetTextSize(0.02);

  // ---------- Call HIST ---------- //
  TFile *f0 = new TFile("rate_"+StrMA+".root");

  TH1D *h0 = (TH1D*)f0 -> Get(Object);

  h0 -> SetMarkerStyle(21);
  h0 -> SetMarkerSize(1.0);
  h0 -> SetMarkerColor(kBlack);
  h0 -> SetLineColor(kBlack);
  h0 -> SetStats(1);

	double BinContent;
	for(int i=0; i<h0->GetNbinsX(); i++) {
		BinContent = h0->GetBinContent(i);
		if( BinContent == 0 ) continue;
		h0->SetBinError(i, 0.1 );
	}

	vector<double> v;
	for(int i=0; i<10; i++) {
		v.push_back( h0->GetBinContent(i) );
	}

	YRangeMin = 0; YRangeMax = *max_element(v.begin(), v.end()) * 1.2;

  h0 -> GetXaxis()->SetRangeUser(XRangeMin,XRangeMax);
  h0 -> GetYaxis()->SetRangeUser(YRangeMin,h0->GetMaximum()*1.2);
  if(letLogX) h0->GetXaxis()->SetMoreLogLabels();

  h0 -> Draw("P");

/*
	TF1 *fit = new TF1("fit",nPol,0,2200,4);
	fit->SetParameters(3.0,0.01,1.0,2.0);
	fit->SetParNames("n1","n0","Decay","DecayExp");
*/
  TF1 *fit = new TF1("fit",FitPol4,0,2200,5);
  fit->SetParameters(0.0,0.0,0.0,0.0,10.0);
  fit->SetParNames("A","B","C","D","E");

	h0->Fit(fit,"RM");

	fit->Draw("SAME");

	double A, B, C, D, E;
	A = fit->GetParameter(0);
	B = fit->GetParameter(1);
	C = fit->GetParameter(2);
	D = fit->GetParameter(3);
	E = fit->GetParameter(4);

  h0 -> GetYaxis()->SetTitle("Signal rate (accXeff)");
  h0 -> GetYaxis()->SetTitleSize(0.04);
  h0 -> GetYaxis()->SetLabelSize(0.025);
  h0 -> GetYaxis()->SetTitleOffset(1.1);

  h0 -> GetXaxis()->SetTitle("M_{4l} [GeV]");
  h0 -> GetXaxis()->SetTitleSize(0.04);
  h0 -> GetXaxis()->SetTitleOffset(0.7);
  h0 -> GetXaxis()->SetLabelSize(0.025);

  //TLegend* leg;
  //leg = new TLegend(0.80,0.85,1.0,1.0);
  //leg -> AddEntry(h0,"resonance","lep");
  //leg -> Draw();

  TLatex *tex = new TLatex(0.1,0.92,"#font[60]{CMS SIM}");
  tex->SetNDC();
  tex->SetLineWidth(2);
  tex->Draw();

  gPad->RedrawAxis();

  c1->SaveAs("fit_rate_"+StrMA+".png");

	/////////////////////////////////////////////////

	TH1D* rate = new TH1D("rate","",2500,0,2500); rate->Sumw2();

  vector<int> M_H;
  for(int i=0; i<20; i++) {
    M_H.push_back( 200 + 10 * i );
  }
  for(int i=0; i<10; i++) {
    M_H.push_back( 400 + 20 * i );
  }
  for(int i=0; i<10; i++) {
    M_H.push_back( 600 + 40 * i );
  }
  for(int i=0; i<11; i++) {
    M_H.push_back( 1000 + 100 * i );
  }

  cout<<setprecision(5);
  cout<<right;

  cout<<setw(6)<<"M_H"<<setw(10)<<"rate"<<endl;;
	double r;
  for(int i=0; i<M_H.size(); i++) {
		r = pol4(M_H[i],A,B,C,D,E);
    cout<<setw(6)<<M_H[i]<<setw(10)<<r<<endl;;
		rate->Fill( M_H[i] , r );
  }

	rate->SetMarkerStyle(21);

	TFile *f = new TFile("fit_rate_"+StrMA+".root","RECREATE");
	rate->Write();
	f->Close();

  cout<<"end"<<endl;
}
