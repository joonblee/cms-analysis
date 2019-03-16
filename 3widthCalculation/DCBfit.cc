
#define lumi 35900.0

Double_t fitf(Double_t *x,Double_t *par) {
  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t DoubleCB(Double_t *x, Double_t *par) {
  // 0(const), 1(mean), 2(sigma), 3(alpha1), 4(n1), 5(alpha2), 6(n2)
  Double_t t = (x[0]-par[1])/par[2];
  if (t>=-1*par[3] && t<=par[5]) {
    return par[0]*exp(-0.5*t*t);
  }
  else if (t<-1*par[3]) {
    Double_t A1 = pow(par[4]/fabs(par[3]),par[4])*exp(-1*par[3]*par[3]/2);
    Double_t B1 = par[4]/fabs(par[3])-fabs(par[3]);
    return par[0]*A1*pow(B1-t,-par[4]);
  }
  else if (t>par[5]) {
    Double_t A2 = pow(par[6]/fabs(par[5]),par[6])*exp(-par[5]*par[5]/2);
    Double_t B2 = par[6]/fabs(par[5])-fabs(par[5]);
    return par[0]*A2*pow(B2+t,-par[6]);
  }
  else {
    //cout << "ERROR evaluating range..." << endl;
    return 99;
  }
}


Double_t funVoigt(Double_t *x, Double_t *par) {
  return par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3], 4);
}

void DCBfit(double M_H = 2000, double M_A = 50, double bin = 20) {

  TString Object = "PPPP_M";

  double XRangeMin = M_H * 0.9;
  double XRangeMax = M_H * 1.1;

	bin = 1;
	if( M_H > 1000 ) bin = 4;

  TString str_M_H = to_string( (int)M_H );
  TString str_M_A = to_string( (int)M_A );

  bool letLogX = false; bool letLogY = false;
  double YRangeMin = 0; double YRangeMax = 0;

	double DoubleBin = (double)bin / 4.0;
	stringstream stream;
	stream << fixed << setprecision(2) << DoubleBin;
	string StringBin = stream.str();

  TString xTitle = Object;
  TString yTitle = "events / "+StringBin+" (GeV)";
  if(letLogY) yTitle = "events (log)";

  TCanvas *c1 = new TCanvas("c1","histogram",800,800);
  c1->SetBottomMargin(0.2);
  c1->cd();

  auto pad1 = new TPad("pad1","",0,0.0,1,1,0);
  pad1->SetLeftMargin(0.1);
  if(letLogX) pad1 -> SetLogx();
  if(letLogY) pad1 -> SetLogy();
  pad1->Draw();
  pad1->cd();

  vector<double> v,v1,v2, ve;

  // ---------- Call HIST ---------- //
  TFile *f0 = new TFile("../final_result/zp"+str_M_H+"phi"+str_M_A+".root");
  TH1D *h0 = (TH1D*)f0 -> Get(Object);
  h0 -> SetMarkerStyle(21);
  h0 -> SetMarkerSize(0.5);
  h0 -> SetLineColor(kBlack);

	TH1D *h1 = (TH1D*)h0->Clone(Object+"_");

  h0->Rebin(bin);
  for(int i=0; i<h0->GetNbinsX(); i++) {
    v.push_back(h0->GetBinContent(i));
		ve.push_back(h0->GetBinError(i));
		if( h0->GetXaxis()->GetBinLowEdge(i) < M_H ) v1.push_back(h0->GetBinContent(i));
		else v2.push_back(h0->GetBinContent(i));
  }
  h0 -> SetStats(0);

	double MaxElement = *max_element(v.begin(),v.end());
  double MaxElement1 = *max_element(v1.begin(),v1.end());
	double MaxElement2 = *max_element(v2.begin(),v2.end());
	double MaxElementError;

  YRangeMax = MaxElement * 1.2;

  if(XRangeMax!=0) h0 -> GetXaxis()->SetRangeUser(XRangeMin,XRangeMax);
  if(YRangeMax!=0) h0 -> GetYaxis()->SetRangeUser(YRangeMin,YRangeMax);
  if(letLogX) h0->GetXaxis()->SetMoreLogLabels();

  h0->Draw("P");

  double m = h0->GetMean(); double rms = h0->GetRMS()/2.0;
	double u = 5.0;
	if( M_H > 1000 ) u = 7.0;

	// +++ ----- gaus fit ----- +++ //

	TF1 *gausfit = new TF1("gausfit","gaus",m-rms,m+rms);
	gausfit->SetParameters(h0->GetMaximum(),m,rms);
	gausfit->SetParNames("Constant","Mean","Sigma");
	gausfit->SetParLimits(0,0.7*h0->GetMaximum(),1.3*h0->GetMaximum());
	gausfit->SetParLimits(1,max(m-rms,M_H-u),M_H);

	h0->Fit(gausfit, "R");
	gausfit->SetLineColor(kOrange);
	gausfit->Draw("SAME");

  double gConst = gausfit->GetParameter(0);
  double gMean = gausfit->GetParameter(1);
  double gSigma = gausfit->GetParameter(2);
  double gChi2 = gausfit->GetChisquare();
  double gNdf = gausfit->GetNDF();


	// +++ ----- DCB fit ----- +++ //

	double PlotXLimit = 5;
	if( M_A == 1 ) {
		PlotXLimit = 8;
	}

	double unc1 = 0.1; double unc2 = 0.1; double unc3 = 0.1; double unc4 = 0.1; double unc5 = 0.1;
	double AlphaLimit = 0.5;
	double ggExp1 = 10.0; ggExp2 = 10.0;

  TF1 *func = new TF1("fit",DoubleCB,gMean-gSigma*PlotXLimit,gMean+gSigma*PlotXLimit,7);
  func->SetParameters(gConst, gMean, gSigma, 0.7, gExp1, 0.7, gExp2);
  func->SetParNames("Const","Mean","Sigma","alpha1","Exp1","alpha2","Exp2");


	func->SetParLimits(0, ggConst*(1.0-unc1), ggConst*(1.0+unc1));
	func->SetParLimits(1, ggMean*(1.0-unc2), ggMean*(1.0+unc3));
	func->SetParLimits(2, ggSigma*(1.0-unc4), ggSigma*(1.0+unc5));

  func->SetParLimits(3, AlphaLimit, 2.0);
	func->SetParLimits(4, 1.0, 9999.0);
	func->SetParLimits(5, AlphaLimit, 2.0);
	func->SetParLimits(6, 1.0, 9999.0);

  h0->Fit(func, "RL");

	func->SetLineColor(kBlue);
  func->Draw("SAME");

  double tConst = func->GetParameter(0);
  double tMean = func->GetParameter(1);
  double tSigma = func->GetParameter(2);
  double tAlpha1 = func->GetParameter(3);
  double tExp1 = func->GetParameter(4);
  double tAlpha2 = func->GetParameter(5);
  double tExp2 = func->GetParameter(6);

  double tChi2 = func->GetChisquare();
  double tNdf = func->GetNDF();


  ///////////////////////////////////////////////////////////////

  cout<<endl<<endl;
  cout<<"+ --- fit result --- +"<<endl;
  cout<<"const : "<<tConst<<"  mean : "<<tMean<<"  sigma : "<<tSigma<<endl;
  cout<<"alpha1 : "<<tAlpha1<<"  exp1 : "<<tExp1<<"  alpha2 : "<<tAlpha2<<"  exp2 : "<<tExp2<<endl;
  cout<<"chi2 : "<<tChi2<<"  ndf : "<<tNdf<<"  chi2/dof : "<<tChi2/tNdf<<endl<<endl;


	TH1D* sigma = new TH1D("sigma","",2500,0,2500); sigma->Sumw2();
	sigma->Fill( M_H , tSigma );
	sigma->SetMarkerStyle(21);
	sigma->SetMarkerSize(1.0);

	TH1D* Esigma = new TH1D("error","",2500,0,2500); Esigma->Sumw2();
	Esigma->Fill( M_H , tSigmaError );

	TH1D* RMS = new TH1D("rms","",2500,0,2500);	 RMS->Sumw2();
	RMS->Fill( M_H , rms );
	RMS->SetMarkerStyle(23);
	RMS->SetMarkerSize(1.0);


	TFile *f = new TFile("sigma_"+str_M_H+"_"+str_M_A+".root","RECREATE");
	sigma->Write();
	Esigma->Write();
	RMS->Write();
	f->Close();

  cout<<"end"<<endl;
}
