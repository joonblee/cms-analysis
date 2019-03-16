#define higgs_cs 0.01212
#define higgs_ws 999800.0
#define dy_cs 5765.4 //5765.4
#define dy_ws 80780990.0 //122055200.0 //11442000+12629500+12152400+12664100+12255800+11492100+12438700+12285200+12474000+12221400

#define dylow_cs 18810.0
#define dylow_ws 99831100.0

#define dy100to200_cs 226.6
#define dy100to200_ws 10310600.0

#define dy200to400_cs 7.77
#define dy200to400_ws 1760976.0

#define dy400to500_cs 0.4065 
#define dy400to500_ws 151190.0

#define zg_cs 117.864
#define zg_ws 9709530

#define wz_cs 22.82
#define wz_ws 1000000.0
#define tw_cs 76.18
#define tw_ws 13885920.0 //6933090 + 6952830
#define tbarW_ws 1.0
#define tW_ws 1.0

#define tt_cs 815.9 //679.1
#define tt_ws 77061150.0 //7447700+7913880+7988370+7611900+7750470+7526200+7914710+7484670+8085290+7337960
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

//#define lumi 35900.0

void setDataStyle(TH1D *h, int bin, int MarkerStyle, double MarkerSize , int MarkerColor);
void setMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double LineWidth, double lumi, double xsec, double ws);

void hist_final_mass(double XRangeMin = 50, double XRangeMax = 500, double bin = 50) {
  TString Object = "PPPP_M";

  bool letLogX = false; bool letLogY = true;
  double YRangeMin = 0; double YRangeMax = 0;

	double LW = 1.0;

  TString xTitle = Object;
  TString yTitle = "Events / 5 GeV";
  //if(letLogY) yTitle = "events (log)";

	double lumi = 35900;

  TCanvas *c1 = new TCanvas("c1","histogram",1000,1600);
  //c1->SetBottomMargin(0.2);
  c1->cd();

  auto pad1 = new TPad("pad1","",0,0.26,1,1,0);
  pad1->SetLeftMargin(0.1);
  //pad1->SetBottomMargin(0.3);
  if(letLogX) pad1 -> SetLogx();
  if(letLogY) pad1 -> SetLogy();
  pad1->Draw();

  auto pad2 = new TPad("pad2","",0,0,1,0.33,0);
  pad2->SetLeftMargin(0.1);
  pad2->SetBottomMargin(0.3);
  if(letLogX) pad2->SetLogx();
  pad2->Draw();

  /////////////////////////////////////////////////////////////////////////////
  pad1->cd();

  vector<double> v;

  // ---------- Call HIST ---------- //
  TFile *f0 = new TFile("final_result/FR_data.root");
  TH1D *h0 = (TH1D*)f0 -> Get(Object);
	setDataStyle(h0, bin, 20, 1.0, kBlack);


  // Irreducible Bg //
  TFile *f1 = new TFile("final_result/FR_GGToZZ.root");
  TFile *f2 = new TFile("final_result/FR_qqToZZ.root");
  TFile *f3 = new TFile("final_result/FR_DY.root");
  //TFile *f4 = new TFile("final_result/FR_TT.root");
  TFile *f5 = new TFile("final_result/FR_higgs.root");
  TFile *f6 = new TFile("final_result/FR_WWZ.root");
  TFile *f7 = new TFile("final_result/FR_WZZ.root");
  TFile *f8 = new TFile("final_result/FR_ZZZ.root");
  TFile *f9 = new TFile("final_result/FR_ttZ.root");
  TFile *f10 = new TFile("final_result/FR_WZ.root");
  TFile *f11 = new TFile("final_result/FR_tW.root");
	TFile *f16 = new TFile("final_result/FR_ZG.root");
	TFile *f17 = new TFile("final_result/FR_DY_VirtualGamma.root");
	//TFile *f18 = new TFile("final_result/FR_DY-M200to400_VirtualGamma.root");
	//TFile *f19 = new TFile("final_result/FR_DY-M400to500_VirtualGamma.root");

  TH1D *h1 = (TH1D*)f1 -> Get(Object); // gg
	setMCStyle(h1, bin, kYellow+1, kYellow+1, LW, lumi, gg_cs, gg_ws);

  TH1D *h2 = (TH1D*)f2 -> Get(Object); // qq
  setMCStyle(h2, bin, kBlack, kYellow+1, LW, lumi, qq_cs, qq_ws);

  //TH1D *h3 = (TH1D*)f3 -> Get("PPPP_FourMu");
  //setMCStyle(h3, bin, kViolet, kViolet, LW, lumi, dy_cs, dy_ws);

  TH1D *h5 = (TH1D*)f5 -> Get(Object); // higgs
  setMCStyle(h5, bin, kBlack, kCyan+2, LW, lumi, higgs_cs, higgs_ws);

  TH1D *h6 = (TH1D*)f6 -> Get(Object); // WWZ
  setMCStyle(h6, bin, kBlue, kBlue, LW, lumi, wwz_cs, wwz_ws);

  TH1D *h7 = (TH1D*)f7 -> Get(Object);
  setMCStyle(h7, bin, kBlue, kBlue, LW, lumi, wzz_cs, wzz_ws);

  TH1D *h8 = (TH1D*)f8 -> Get(Object);
  setMCStyle(h8, bin, kBlue, kBlue, LW, lumi, zzz_cs, zzz_ws);

  TH1D *h9 = (TH1D*)f9 -> Get(Object); // ttZ
	setMCStyle(h9, bin, kBlue, kBlue, LW, lumi, ttz_cs, ttz_ws);

  TH1D *h10 = (TH1D*)f10 -> Get(Object); // WZ
  setMCStyle(h10, bin, kBlue, kBlue, LW, lumi, wz_cs, wz_ws);

  TH1D *h11 = (TH1D*)f11 -> Get(Object); // tW
  setMCStyle(h11, bin, kBlack, kBlue, LW, lumi, tw_cs, tw_ws);


	TH1D *h16 = (TH1D*)f16->Get("PPPP_FourMu");
	setMCStyle(h16, bin, kViolet+5, kViolet+5, LW, lumi, zg_cs, zg_ws);

	TH1D *h17 = (TH1D*)f17->Get("PPPP_FourMu");
	setMCStyle(h17, bin, kBlack, kViolet+5, LW, lumi, dy_cs, dy_ws);

  //TH1D *h18 = (TH1D*)f18->Get("PPPP_FourMu");
  //setMCStyle(h18, bin, kRed, kRed, LW, lumi, dy200to400_cs, dy200to400_ws);

  //TH1D *h19 = (TH1D*)f19->Get("PPPP_FourMu");
  //setMCStyle(h19, bin, kRed+3, kRed+3, LW, lumi, dy400to500_cs, dy400to500_ws);

	h0->Sumw2(); h1->Sumw2(); h2->Sumw2(); h5->Sumw2(); h6->Sumw2(); h7->Sumw2(); h8->Sumw2(); h9->Sumw2(); h10->Sumw2(); h11->Sumw2(); h16->Sumw2(); h17->Sumw2(); 


  // Reducible Bg //
  TH1D *h_3P1F_data = (TH1D*)f0->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(data)");
  TH1D *h_2P2F_data = (TH1D*)f0->Get("PPFF_M_UR03_ReWt")->Clone("2P2F(data)");
	TH1D *h_1P3F_data = (TH1D*)f0->Get("PFFF_M_UR03_ReWt")->Clone("1P3F(data)");
	TH1D *h_4F_data = (TH1D*)f0->Get("FFFF_M_UR03_ReWt")->Clone("4F(data)");

  TH1D *h_3P1F_gg   = (TH1D*)f1->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(gg)");
  TH1D *h_3P1F_qq   = (TH1D*)f2->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(qq)");
  TH1D *h_3P1F_higgs= (TH1D*)f5->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(h0)");
  TH1D *h_3P1F_WWZ  = (TH1D*)f6->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(WWZ)");
  TH1D *h_3P1F_WZZ  = (TH1D*)f7->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(WZZ)");
  TH1D *h_3P1F_ZZZ  = (TH1D*)f8->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(ZZZ)");
  TH1D *h_3P1F_ttZ  = (TH1D*)f9->Get("PPPF_M_UR03_ReWt")->Clone("3P1F(ttZ)");
  TH1D *h_3P1F_zg    = (TH1D*)f16->Get("PPPF_FourMu_UR03_ReWt")->Clone("3P1F(ZG)");

  TH1D *h_2P2F_gg   = (TH1D*)f1->Get("PPFF_M_UR03_ReWt")->Clone("2P2F(gg)");
  TH1D *h_2P2F_qq   = (TH1D*)f2->Get("PPFF_M_UR03_ReWt")->Clone("2P2F(qq)");
  TH1D *h_2P2F_higgs= (TH1D*)f5->Get("PPFF_M_UR03_ReWt")->Clone("2P2F(h0)");
  TH1D *h_2P2F_zg    = (TH1D*)f16->Get("PPFF_FourMu_UR03_ReWt")->Clone("2P2F(ZG)");

  TH1D *h_3P1F_data_R03 = (TH1D*)f0->Get("PPPF_M_LR03_ReWt")->Clone("3P1F_R03(data)");
  TH1D *h_2P2F_data_R03 = (TH1D*)f0->Get("PPFF_M_LR03_ReWt")->Clone("2P2F_R03(data)");

  TH1D *h_3P1F_gg_R03    = (TH1D*)f1->Get("PPPF_M_LR03_ReWt")->Clone("3P1F_R03(gg)");
  TH1D *h_3P1F_qq_R03    = (TH1D*)f2->Get("PPPF_M_LR03_ReWt")->Clone("3P1F_R03(qq)");
  TH1D *h_3P1F_higgs_R03 = (TH1D*)f5->Get("PPPF_M_LR03_ReWt")->Clone("3P1F_R03(h0)");
  TH1D *h_3P1F_zg_R03    = (TH1D*)f16->Get("PPPF_FourMu_LR03_ReWt")->Clone("3P1F_R03(ZG)");

  TH1D *h_2P2F_gg_R03    = (TH1D*)f1->Get("PPFF_M_LR03_ReWt")->Clone("2P2F_R03(gg)");
  TH1D *h_2P2F_qq_R03    = (TH1D*)f2->Get("PPFF_M_LR03_ReWt")->Clone("2P2F_R03(qq)");
  TH1D *h_2P2F_higgs_R03 = (TH1D*)f5->Get("PPFF_M_LR03_ReWt")->Clone("2P2F_R03(h0)");
  TH1D *h_2P2F_zg_R03    = (TH1D*)f16->Get("PPFF_FourMu_LR03_ReWt")->Clone("2P2F_R03(ZG)");

	h_3P1F_data->Sumw2(); h_2P2F_data->Sumw2();  h_1P3F_data->Sumw2();  h_4F_data->Sumw2();  h_3P1F_gg->Sumw2(); h_3P1F_qq->Sumw2(); h_3P1F_higgs->Sumw2();  h_3P1F_WWZ->Sumw2();  h_3P1F_WZZ->Sumw2();  h_3P1F_ZZZ->Sumw2(); h_3P1F_ttZ->Sumw2(); h_3P1F_zg->Sumw2(); h_2P2F_gg->Sumw2(); h_2P2F_qq->Sumw2(); h_2P2F_higgs->Sumw2(); h_2P2F_zg->Sumw2(); h_3P1F_data_R03->Sumw2();  h_2P2F_data_R03->Sumw2();  h_3P1F_gg_R03->Sumw2();  h_3P1F_qq_R03->Sumw2(); h_3P1F_higgs_R03->Sumw2(); h_3P1F_zg_R03->Sumw2(); h_2P2F_gg_R03->Sumw2(); h_2P2F_qq_R03->Sumw2(); h_2P2F_higgs_R03->Sumw2(); h_2P2F_zg_R03->Sumw2();

  h_3P1F_gg->Scale(lumi*gg_cs/qq_ws);
  h_3P1F_qq->Scale(lumi*qq_cs/qq_ws);
  h_3P1F_higgs->Scale(lumi*higgs_cs/higgs_ws);
  h_3P1F_WWZ->Scale(lumi*wwz_cs/wwz_ws);
  h_3P1F_WZZ->Scale(lumi*wzz_cs/wzz_ws);
  h_3P1F_ZZZ->Scale(lumi*zzz_cs/zzz_ws);
  h_3P1F_ttZ->Scale(lumi*ttz_cs/ttz_ws);
  h_3P1F_zg->Scale(lumi*zg_cs/zg_ws);

  h_2P2F_gg->Scale(lumi*gg_cs/qq_ws);
  h_2P2F_qq->Scale(lumi*qq_cs/qq_ws);
  h_2P2F_higgs->Scale(lumi*higgs_cs/higgs_ws);
  h_2P2F_zg->Scale(lumi*zg_cs/zg_ws);

  h_3P1F_gg_R03->Scale(lumi*gg_cs/qq_ws);
  h_3P1F_qq_R03->Scale(lumi*qq_cs/qq_ws);
  h_3P1F_higgs_R03->Scale(lumi*higgs_cs/higgs_ws);
  h_3P1F_zg_R03->Scale(lumi*zg_cs/zg_ws);

  h_2P2F_gg_R03->Scale(lumi*gg_cs/qq_ws);
  h_2P2F_qq_R03->Scale(lumi*qq_cs/qq_ws);
  h_2P2F_higgs_R03->Scale(lumi*higgs_cs/higgs_ws);
  h_2P2F_zg_R03->Scale(lumi*zg_cs/zg_ws);


  h_3P1F_data->Add(h_3P1F_gg,-1.0);
  h_3P1F_data->Add(h_3P1F_qq,-1.0);
  h_3P1F_data->Add(h_3P1F_higgs,-1.0);
  h_3P1F_data->Add(h_3P1F_WWZ,-1.0);
  h_3P1F_data->Add(h_3P1F_WZZ,-1.0);
  h_3P1F_data->Add(h_3P1F_ZZZ,-1.0);
  h_3P1F_data->Add(h_3P1F_ttZ,-1.0);
  h_3P1F_data->Add(h_3P1F_zg,-1.0);

  h_2P2F_data->Add(h_2P2F_gg,-1.0);
  h_2P2F_data->Add(h_2P2F_qq,-1.0);
  h_2P2F_data->Add(h_2P2F_higgs,-1.0);
  h_2P2F_data->Add(h_2P2F_zg,-1.0);

  h_3P1F_data->Add(h_2P2F_data,-2.0);

  for(int i=0; i<h0->GetNbinsX(); i++) {
    if(h_3P1F_data->GetBinContent(i) < 0) {
      h_3P1F_data->SetBinContent(i,0.0);
      h_3P1F_data->SetBinError(i,0.0);
    }
  }

  h_3P1F_data->Add(h_2P2F_data);
	h_3P1F_data->Add(h_1P3F_data);


  h_3P1F_data_R03->Add(h_3P1F_gg_R03,-1.0);
  h_3P1F_data_R03->Add(h_3P1F_qq_R03,-1.0);
  h_3P1F_data_R03->Add(h_3P1F_higgs_R03,-1.0);
  h_3P1F_data_R03->Add(h_3P1F_zg_R03,-1.0);

  h_2P2F_data_R03->Add(h_2P2F_gg_R03,-1.0);
  h_2P2F_data_R03->Add(h_2P2F_qq_R03,-1.0);
  h_2P2F_data_R03->Add(h_2P2F_higgs_R03,-1.0);
	h_2P2F_data_R03->Add(h_2P2F_zg_R03,-1.0);

	h_3P1F_data_R03->Add(h_2P2F_data_R03);


	setMCStyle(h_3P1F_data, bin, kBlack, kRed, LW, 1.0, 1.0, 1.0);
	setMCStyle(h_3P1F_data_R03, bin, kBlack, kOrange+1, LW, 1.0, 1.0, 1.0);

  //h_3P1F_data->Draw("HIST SAME");

  THStack *hs = new THStack("hs","");
  hs->Add(h_3P1F_data); hs->Add(h_3P1F_data_R03); hs->Add(h6); hs->Add(h7); hs->Add(h8); hs->Add(h9); hs->Add(h10); hs->Add(h11); hs->Add(h16); hs->Add(h17); hs->Add(h1); hs->Add(h2); hs->Add(h5); 

  TH1D *h_mc = (TH1D*)h1->Clone("h_mc");
  h_mc->Add(h2); h_mc->Add(h5); h_mc->Add(h6); h_mc->Add(h7); h_mc->Add(h8); h_mc->Add(h9); h_mc->Add(h16); h_mc->Add(h17); h_mc->Add(h10); h_mc->Add(h11);
  h_mc->Add(h_3P1F_data); h_mc->Add(h_3P1F_data_R03);
	h_mc->Sumw2();

  TH1D *h_Others = (TH1D*)h6->Clone("h_Others"); h_Others->Sumw2();
  h_Others->Add(h7);
  h_Others->Add(h8);
  h_Others->Add(h9);
  h_Others->Add(h10);
  h_Others->Add(h11);

  TH1D *h_ZG = (TH1D*)h16->Clone("h_ZG"); h_ZG->Sumw2();
  h_ZG->Add(h17);

  TH1D *h_ZZ = (TH1D*)h1->Clone("h_ZZ"); h_ZZ->Sumw2();
  h_ZZ->Add(h2);

  TH1D *h_BKG = (TH1D*)h5->Clone("h_BKG"); h_BKG->Sumw2();
  h_BKG->Add(h_ZZ);
  h_BKG->Add(h_ZG);
  h_BKG->Add(h_Others);
  h_BKG->Add(h_3P1F_data_R03);
  h_BKG->Add(h_3P1F_data);

	/*
  const int Nbins = 6;
  Double_t xbins[Nbins+1] = {50,100,200,300,400,600,3000};
  h0 = (TH1D*)h0->Rebin(Nbins,"",xbins);
  h5 = (TH1D*)h5->Rebin(Nbins,"",xbins);
  h_ZZ = (TH1D*)h_ZZ->Rebin(Nbins,"",xbins);
  h_ZG = (TH1D*)h_ZG->Rebin(Nbins,"",xbins);
  h_Others = (TH1D*)h_Others->Rebin(Nbins,"",xbins);
  h_3P1F_data_R03 = (TH1D*)h_3P1F_data_R03->Rebin(Nbins,"",xbins);
  h_3P1F_data = (TH1D*)h_3P1F_data->Rebin(Nbins,"",xbins);

  h_BKG = (TH1D*)h_BKG->Rebin(Nbins,"",xbins);

  cout<<right;
  cout<<setprecision(4);

  cout<<setw(11)<<"data"<<setw(11)<<"higgs"<<setw(11)<<"ZZ"<<setw(11)<<"ZG"<<setw(11)<<"Others"<<setw(11)<<"R1"<<setw(11)<<"R2"<<endl;
  for(int j=1; j<h0->GetNbinsX()+1; j++) {
    cout<<setw(11)<<h0->GetBinContent(j)<<"+"<<h0->GetBinError(j);
    cout<<setw(11)<<h5->GetBinContent(j)<<"+"<<h5->GetBinError(j);
    cout<<setw(11)<<h_ZZ->GetBinContent(j)<<"+"<<h_ZZ->GetBinError(j);
    cout<<setw(11)<<h_ZG->GetBinContent(j)<<"+"<<h_ZG->GetBinError(j);
    cout<<setw(11)<<h_Others->GetBinContent(j)<<"+"<<h_Others->GetBinError(j);
    cout<<setw(11)<<h_3P1F_data_R03->GetBinContent(j)<<"+"<<h_3P1F_data_R03->GetBinError(j);
    cout<<setw(11)<<h_3P1F_data->GetBinContent(j)<<"+"<<h_3P1F_data->GetBinError(j);
    cout<<setw(11)<<h_BKG->GetBinContent(j)<<"+"<<h_BKG->GetBinError(j);
    cout<<endl;
  }
	*/


  for(int i=0; i<h0->GetNbinsX(); i++) {
		v.push_back(h0->GetBinContent(i));
    v.push_back(h_mc->GetBinContent(i));
  }

	double multiplier = 1.05;
	if( letLogY ) {
		multiplier = 3.0;
		YRangeMin = 0.1;
	}
  YRangeMax = *max_element(v.begin(),v.end()) * multiplier;
	//YRangeMax = 70;

  if(XRangeMax!=0) h0 -> GetXaxis()->SetRangeUser(XRangeMin,XRangeMax);
  if(YRangeMax!=0) h0 -> GetYaxis()->SetRangeUser(YRangeMin,YRangeMax);
  h0->GetXaxis()->SetMoreLogLabels();
  h0->GetYaxis()->SetMoreLogLabels();
	h0->GetYaxis()->SetNoExponent();

  h0 -> Draw("P");
  hs -> Draw("HIST SAME");
  h0 -> Draw("P SAME");

  h0 -> GetYaxis()->SetTitle(yTitle);
  h0 -> GetYaxis()->SetTitleSize(0.045);
  h0 -> GetYaxis()->SetLabelSize(0.040);
  h0 -> GetYaxis()->SetTitleOffset(1.0);

  h0 -> GetXaxis()->SetTitle(xTitle);
  h0 -> GetXaxis()->SetTitleSize(0/*0.04*/);
  h0 -> GetXaxis()->SetTitleOffset(0.7);
  h0 -> GetXaxis()->SetLabelSize(0.025);

  TLegend* leg;
  leg = new TLegend(0.56,0.52,0.90,0.89);
  leg -> AddEntry(h0,"DATA","lep");
  leg -> AddEntry(h5,"higgs","f");
  leg -> AddEntry(h2,"ZZ","f");
	leg -> AddEntry(h17,"ZG","f");
  leg -> AddEntry(h11,"Others","f");
	leg->AddEntry(h_3P1F_data_R03,"type I Reducible bg","f");
  leg->AddEntry(h_3P1F_data,"type II Reducible bg","f");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg -> Draw();


  TLatex *tex = new TLatex(0.1,0.92,"#font[60]{CMS Preliminary            35.9 fb^{-1} (13 TeV)}");
  tex->SetNDC();
  tex->SetLineWidth(2);
  tex->Draw();

  gPad->RedrawAxis();

  ////////////////////////////////////////////////////////////////////////
  pad2->cd();

  TH1D *h_ratio = (TH1D*)h0->Clone("h_ratio");
	h_ratio->Sumw2();

  h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  if(XRangeMax!=0) h_ratio->GetXaxis()->SetRangeUser(XRangeMin,XRangeMax);
  h_ratio->Divide(h_mc);


  h_ratio->GetXaxis()->SetTitle("M_{4#mu} [GeV]");
  h_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_ratio->GetXaxis()->SetTitleSize(0.1);
  h_ratio->GetXaxis()->SetLabelSize(0.09);
  h_ratio->GetXaxis()->SetLabelOffset(0.01);

  h_ratio->GetYaxis()->SetTitle("Data / MC");
  h_ratio->GetYaxis()->SetTitleSize(0.095);
  h_ratio->GetYaxis()->SetTitleOffset(0.40);
  h_ratio->GetYaxis()->SetLabelSize(0.075);

  h_ratio->SetStats(0);
  h_ratio->Draw("P");

  TF1 *hline = new TF1("","1",XRangeMin,XRangeMax);
  hline->SetLineColor(kRed);
  hline->Draw("SAME");

	TString logTag = "";
	if( letLogY ) logTag = "_logY";

  c1->SaveAs("hist_4m_lowMass.png");
	c1->SaveAs("hist_4m_lowMass.pdf");
  cout<<"end"<<endl;
}

void setDataStyle(TH1D *h, int bin, int MarkerStyle, double MarkerSize , int MarkerColor) {
  h->Rebin(bin);

  h->SetMarkerStyle(MarkerStyle);
  h->SetMarkerSize(MarkerSize);
  h->SetLineColor(MarkerColor);
  h->SetStats(0);
}


void setMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double LineWidth, double lumi, double xsec, double ws) {
  h->Rebin(bin);

  h->SetLineColor(LineColor);
  h->SetFillColor(FillColor);
	h->SetLineWidth(LineWidth);
  h->Scale(lumi*xsec/ws);
}

