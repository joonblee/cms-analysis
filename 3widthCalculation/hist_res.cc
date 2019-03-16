#define lumi 35900.0//35900.0

Double_t nPol(Double_t *x, Double_t *par) {
	Double_t t = x[0];
	return 0.00001*(par[0]*t + par[1])*(1.0-par[2]/4000000.0*pow(t,par[3]));
}

Double_t FitPol4(Double_t *x, Double_t *par) {
	Double_t t = x[0];
	return par[0]*pow(t,4) + par[1]*pow(t,3) + par[2]*pow(t,2) + par[3]*t + par[4];
}

Double_t pol4(double t, double A, double B, double C, double D, double E) {
	return A*pow(t,4) + B*pow(t,3) + C*pow(t,2) + D*t + E;
}

void hist_res(TString StrMA = "A50", TString Object = "sigma", double bin = 1) {

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

  pad1->cd();

  // ---------- Call HIST ---------- //
  TFile *f0 = new TFile("sigma_"+StrMA+".root");

  TH1D *h0 = (TH1D*)f0 -> Get(Object);

  h0 -> SetMarkerStyle(21);
  h0 -> SetMarkerSize(1.0);
  h0 -> SetMarkerColor(kBlack);
  h0 -> SetLineColor(kBlack);
  h0 -> SetStats(1);

	TH1D *h1 = (TH1D*)f0 -> Get("error");
	h1->SetStats(0);

	cout<<right;
	cout<<setprecision(4);
	cout<<"+++ ----- Resonance ----- +++"<<endl;
	cout<<setw(5)<<"Mass"<<setw(10)<<"sig"<<setw(10)<<"res"<<endl;
	double BinContent, LowEdge, BinError;
	for(int i=0; i<h0->GetNbinsX(); i++) {
		BinContent = h0->GetBinContent(i);
		if( BinContent == 0 ) continue;
		LowEdge  = h0->GetXaxis()->GetBinLowEdge(i);
		BinError = h1->GetBinContent(i);
		cout<<setw(5)<<LowEdge<<setw(10)<<BinContent<<setw(10)<<BinContent/LowEdge<<endl;
		h0->SetBinContent(i, BinContent / LowEdge );
		h0->SetBinError(i, BinError / LowEdge );
	}
	cout<<endl<<endl;

	vector<double> v;
	for(int i=0; i<10; i++) {
		v.push_back( h0->GetBinContent(i) );
	}

	YRangeMin = 0; YRangeMax = *max_element(v.begin(), v.end()) * 1.2;


  h0 -> GetXaxis()->SetRangeUser(XRangeMin,XRangeMax);
  h0 -> GetYaxis()->SetRangeUser(YRangeMin,h0->GetMaximum()*1.2);
  if(letLogX) h0->GetXaxis()->SetMoreLogLabels();

  h0 -> Draw("P");

  TF1 *fit = new TF1("fit",FitPol4,0,2200,5);
  fit->SetParameters(0.0,0.0,0.0,0.0,0.01);
  fit->SetParNames("A","B","C","D","E");

	h0->Fit(fit,"RM");

	fit->Draw("SAME");

	double A, B, C, D, E;
	A = fit->GetParameter(0);
	B = fit->GetParameter(1);
	C = fit->GetParameter(2);
	D = fit->GetParameter(3);
	E = fit->GetParameter(4);

  h0 -> GetYaxis()->SetTitle("M_"+StrMA+"  Mass Resolution");
  h0 -> GetYaxis()->SetTitleSize(0.03);
  h0 -> GetYaxis()->SetLabelSize(0.025);
  h0 -> GetYaxis()->SetTitleOffset(1.1);

  h0 -> GetXaxis()->SetTitle("M(H) [GeV]");
  h0 -> GetXaxis()->SetTitleSize(0.04);
  h0 -> GetXaxis()->SetTitleOffset(0.7);
  h0 -> GetXaxis()->SetLabelSize(0.025);

  TLegend* leg;
  leg = new TLegend(0.80,0.85,1.0,1.0);
  leg -> AddEntry(h0,"resonance","lep");

  leg -> Draw();

  TLatex *tex = new TLatex(0.1,0.92,"#font[60]{CMS SIM}");
  tex->SetNDC();
  tex->SetLineWidth(2);
  tex->Draw();

  gPad->RedrawAxis();

  c1->SaveAs("hist_res_"+StrMA+"_"+Object+".png");

	/////////////////////////////////////////////////

	TH1D* MWindow = new TH1D("MWindow","",2500,0,2500); MWindow->Sumw2();

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

  cout<<setw(6)<<"M_H"<<setw(10)<<"res"<<setw(10)<<"M window"<<endl;;
	double res;
  for(int i=0; i<M_H.size(); i++) {
		res = pol4(M_H[i],A,B,C,D,E);
    cout<<setw(6)<<M_H[i]<<setw(10)<<res<<setw(10)<<res*M_H[i]<<endl;;
		MWindow->Fill( M_H[i] , res * M_H[i] );
  }

	MWindow->SetMarkerStyle(21);

	TFile *f = new TFile("MassWindow"+StrMA+".root","RECREATE");
	MWindow->Write();
	f->Close();

  cout<<"end"<<endl;
}
