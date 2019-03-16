#define higgs_cs 0.01212
#define higgs_ws 999800.0
#define dy_cs 5765.4 //5765.4
#define dy_ws 80780990.0 //122055200.0 //11442000+12629500+12152400+12664100+12255800+11492100+12438700+12285200+12474000+12221400

#define dylow_cs 18810.0
#define dylow_ws 99831100.0

#define dy100to200_cs 226.6
#define dy100to200_ws 10310600.0

#define dy200to400_cs 7.77
#define dy200to400_ws 169676.0

#define dy400to600_cs 0.4065 
#define dy400to600_ws 151190.0

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
void setMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double lumi, double xsec, double ws);
void print_contents(TH1D histo);

void nContents(int M_A = 50) {
	TString String_M_A = to_string( M_A );

  TString Object = "PPPP_M";

	double lumi = 35900;

  // ---------- Call HIST ---------- //
  TFile *f0 = new TFile("../final_result/FR_data.root");
  TH1D *h0 = (TH1D*)f0 -> Get(Object);
	setDataStyle(h0, bin, 21, 0.5, kBlack);

  // Irreducible Bg //
  TFile *f1 = new TFile("../final_result/FR_GGToZZ.root");
  TFile *f2 = new TFile("../final_result/FR_qqToZZ.root");
  TFile *f3 = new TFile("../final_result/FR_DY.root");
  //TFile *f4 = new TFile("../final_result/FR_TT.root");
  TFile *f5 = new TFile("../final_result/FR_higgs.root");
  TFile *f6 = new TFile("../final_result/FR_WWZ.root");
  TFile *f7 = new TFile("../final_result/FR_WZZ.root");
  TFile *f8 = new TFile("../final_result/FR_ZZZ.root");
  TFile *f9 = new TFile("../final_result/FR_ttZ.root");
  //TFile *f10 = new TFile("../final_result/FR_WZ.root");
  //TFile *f11 = new TFile("../final_result/FR_tW.root");
	TFile *f16 = new TFile("../final_result/FR_ZG.root");
	TFile *f17 = new TFile("../final_result/FR_DY_VirtualGamma.root");

  TH1D *h1 = (TH1D*)f1 -> Get(Object); // gg
	setMCStyle(h1, bin, kYellow+1, kYellow+1, lumi, gg_cs, gg_ws);

  TH1D *h2 = (TH1D*)f2 -> Get(Object); // qq
  setMCStyle(h2, bin, kYellow+1, kYellow+1, lumi, qq_cs, qq_ws);

  //TH1D *h3 = (TH1D*)f3 -> Get("PPPP_FourMu");
  //setMCStyle(h3, bin, kViolet, kViolet, lumi, dy_cs, dy_ws);

  TH1D *h5 = (TH1D*)f5 -> Get(Object); // higgs
  setMCStyle(h5, bin, kCyan+2, kCyan+2, lumi, higgs_cs, higgs_ws);

  TH1D *h6 = (TH1D*)f6 -> Get(Object); // WWZ
  setMCStyle(h6, bin, kGreen, kGreen, lumi, wwz_cs, wwz_ws);

  TH1D *h7 = (TH1D*)f7 -> Get(Object);
  setMCStyle(h7, bin, kGreen, kGreen, lumi, wzz_cs, wzz_ws);

  TH1D *h8 = (TH1D*)f8 -> Get(Object);
  setMCStyle(h8, bin, kGreen, kGreen, lumi, zzz_cs, zzz_ws);

  TH1D *h9 = (TH1D*)f9 -> Get(Object); // ttZ
	setMCStyle(h9, bin, kBlue, kBlue, lumi, ttz_cs, ttz_ws);

	TH1D *h16 = (TH1D*)f16->Get("PPPP_FourMu");
	setMCStyle(h16, bin, kViolet+5, kViolet+5, lumi, zg_cs, zg_ws);

	TH1D *h17 = (TH1D*)f17->Get("PPPP_FourMu");
	setMCStyle(h17, bin, kAzure+8, kAzure+8, lumi, dy_cs, dy_ws);

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

  for(int i=1; i<1000; i++) {
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


  h_3P1F_data->Rebin(bin);
  h_3P1F_data->SetLineColor(kRed);
  h_3P1F_data->SetFillColor(kRed);

	h_3P1F_data_R03->Rebin(bin);
	h_3P1F_data_R03->SetLineColor(kOrange+1);
	h_3P1F_data_R03->SetFillColor(kOrange+1);

  //h_3P1F_data->Draw("HIST SAME");

  THStack *hs = new THStack("hs","");
  hs->Add(h_3P1F_data); hs->Add(h_3P1F_data_R03);/*hs->Add(h4); hs->Add(h11); hs->Add(h10);*/ hs->Add(h6); hs->Add(h7); hs->Add(h8); hs->Add(h9); /*hs->Add(h12); hs->Add(h3); hs->Add(h13); hs->Add(h14);*/ hs->Add(h16); hs->Add(h17); hs->Add(h1); hs->Add(h2); hs->Add(h5); 

  TH1D *h_mc = (TH1D*)h1->Clone("h_mc");
  h_mc->Add(h2); /*h_mc->Add(h4);*/ h_mc->Add(h5); h_mc->Add(h6); h_mc->Add(h7); h_mc->Add(h8); h_mc->Add(h9); h_mc->Add(h16); h_mc->Add(h17); /*h_mc->Add(h12); h_mc->Add(h3); h_mc->Add(h13); h_mc->Add(h14);*/ /*h_mc->Add(h10); h_mc->Add(h11);*/
  h_mc->Add(h_3P1F_data); h_mc->Add(h_3P1F_data_R03);

	/////////////////////////////////////////////////////////

	vector<int> mH;
  for(int i=0; i<20; i++) {
    mH.push_back( 200 + 10 * i );
  }
  for(int i=0; i<10; i++) {
    mH.push_back( 400 + 20 * i );
  }
  for(int i=0; i<10; i++) {
    mH.push_back( 600 + 40 * i );
  }
  for(int i=0; i<11; i++) {
    mH.push_back( 1000 + 100 * i );
  }


	vector<double> wH;

	TFile *f_ = new TFile("../widthCalculation/MassWindow_phi"+String_M_A+".root");
	TH1D *hh = (TH1D*)f_->Get("MWindow");

	for(int i=0; i<mH.size(); i++) {
		for(int j=0; j<hh->GetNbinsX(); j++) {
			if( hh->GetXaxis()->GetBinLowEdge(j) == mH[i] ) {
				wH.push_back( hh->GetBinContent(j) );
				break;
			}
		}
	}

	TFile *frate = new TFile("fit_rate_phi"+String_M_A+".root");	
	TH1D *hr = (TH1D*)frate->Get("rate");

	cout<<endl<<"+++ ----- rate table for limit plot ----- +++"<<endl<<endl;
	string list_sample[] = {"data","sig","higgs","gg->ZZ","qq->ZZ","DY","ZG","RB_I","RB_II","Others"};
	vector<string> sample;
	sample.insert( sample.end() , begin(list_sample) , end(list_sample) );

	TH1D *h_Others = (TH1D*)h6->Clone("h_Others");
	h_Others->Add(h7);
	h_Others->Add(h8);
	h_Others->Add(h9);

	TH1D* Lhist[] = {h0,hr,h5,h1,h2,h17,h16,h_3P1F_data,h_3P1F_data_R03,h_Others};
	TH1D* histo;
	double binContent;

	cout<<setprecision(5);
	cout<<right;

	for(int i=0; i<mH.size(); i++) {
		cout<<"        ";
    if(i!=0) cout<<"el";
    cout<<"if theta_B == "<<setw(4)<<mH[i]<<":"<<endl;
    cout<<"          rates = {";
    for(int j=0; j<sample.size(); j++) {
      histo = Lhist[j];

      if(j!=0) cout<<", ";
      cout<<"\""<<sample[j]<<"\":";

			if( j == 1 ) {
				for(int k=0; k<histo->GetNbinsX(); k++) {
					double BinEdge = histo->GetXaxis()->GetBinLowEdge(k);
					if( BinEdge == mH[i] ) {
						cout<<setw(10)<<histo->GetBinContent(k);
						break;
					}
				}
			}
			else {
	      binContent = 0;
	      for(int k=0; k<histo->GetNbinsX(); k++) {
					double BinCenter = histo->GetXaxis()->GetBinCenter(k);
	        if( mH[i] - 5*wH[i] < BinCenter && BinCenter < mH[i] + 5*wH[i] ) 
						binContent += histo->GetBinContent(k);
					if( mH[i] + 5*wH[i] < BinCenter ) break;
	      }
				if( binContent >= 0 )
					cout<<setw(10)<<binContent;
				else
					cout<<setw(10)<<0;
			}
		}
		cout<<"}"<<endl;
	}
	cout<<endl<<endl;

	cout<<"Success"<<endl;
}

void setDataStyle(TH1D *h, int bin, int MarkerStyle, double MarkerSize , int MarkerColor) {
  h->Rebin(bin);
  h -> SetMarkerStyle(MarkerStyle);
  h -> SetMarkerSize(MarkerSize);
  h -> SetLineColor(MarkerColor);
  h -> SetStats(0);
}


void setMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double lumi, double xsec, double ws) {
  h->Rebin(bin);
  h -> SetLineColor(LineColor);
  h -> SetFillColor(FillColor);
  h -> Scale(lumi*xsec/ws);
}

void print_contents(TH1D histo) {
  double binContent = 0;
  int n = 0;
  cout<<"115   120   125   130   135   140"<<endl;
  for(int i = 100; i < 150; i++) {
    n = i-1;
    if( n==115 || n==120 || n==125 || n==130 || n==135 || n==140 ) {
      binContent = 0;
      for(int j=0; j<5; j++) {
        binContent += histo.GetBinContent(i+j-2);
      }
      cout<<binContent<<",";
    }
  }
}



