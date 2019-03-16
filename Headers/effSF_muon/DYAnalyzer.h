#pragma once

#include "/u/user/joonblee/Headers/Object.h"
#include "/u/user/joonblee/Headers/NtupleHandle.h"

#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

class DYAnalyzer {
public:

  Double_t Eff_ID_data_BtoF[4][6];
  Double_t Eff_ID_MC_BtoF[4][6];

  Double_t Eff_ID_data_GtoH[4][6];
  Double_t Eff_ID_MC_GtoH[4][6];

  Double_t Eff_LooseID_data_BtoF[4][6];
  Double_t Eff_LooseID_MC_BtoF[4][6];

  Double_t Eff_LooseID_data_GtoH[4][6];
  Double_t Eff_LooseID_MC_GtoH[4][6];


  DYAnalyzer();

  void SetupEfficiencyScaleFactor_BtoF();
  void SetupEfficiencyScaleFactor_GtoH();
  void SetupEfficiencyScaleFactorLoose_BtoF();
  void SetupEfficiencyScaleFactorLoose_GtoH();

  Double_t MuonSF_BtoF_ID( Muon mu );
  Double_t MuonSF_GtoH_ID( Muon mu );
  Double_t LooseMuonSF_BtoF_ID( Muon mu );
  Double_t LooseMuonSF_GtoH_ID( Muon mu );

  Double_t EfficiencySF_EventWeight_HLT_BtoF_ID( vector<Muon> muons );
  Double_t EfficiencySF_EventWeight_HLT_GtoH_ID( vector<Muon> muons );
  Double_t EfficiencySF_EventWeight_HLT_BtoF_LooseID( vector<Muon> muons );
  Double_t EfficiencySF_EventWeight_HLT_GtoH_LooseID( vector<Muon> muons );

	Int_t Find_muon_PtBin_ID(Double_t pt);
	Int_t Find_muon_EtaBin_ID(Double_t eta);

};

DYAnalyzer::DYAnalyzer() {}

void DYAnalyzer::SetupEfficiencyScaleFactor_BtoF() 
{
  TString Location = "/u/user/joonblee/Headers/effSF_muon/";
  cout<<endl<<"[Tag&Probe efficiency is from " << Location <<"]"<<endl;

  TFile *f1 = new TFile( Location + "EfficienciesAndSF_BCDEF.root" );

	TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
	TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins1 = h_ID_data->GetNbinsX();
	Int_t nPtBins1 = h_ID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_ID_data_BtoF[iter_x][iter_y] = ID_data;
			Eff_ID_MC_BtoF[iter_x][iter_y] = ID_MC;
		}
	}

	cout << "Setting for efficiency correction factors (BtoF) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactor_GtoH()
{
	TString Location = "/u/user/joonblee/Headers/effSF_muon/";
	cout << "[Tag&Probe efficiency is from " << Location << "]" << endl; 

	TFile *f1 = new TFile( Location+"EfficienciesAndSF_GH.root" );

	TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
	TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

	Int_t nEtaBins1 = h_ID_data->GetNbinsX();
	Int_t nPtBins1 = h_ID_data->GetNbinsY();

	for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
	{
		for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
		{
			Int_t i_etabin = iter_x + 1;
			Int_t i_ptbin = iter_y + 1;

			Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
			Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

			Eff_ID_data_GtoH[iter_x][iter_y] = ID_data;
			Eff_ID_MC_GtoH[iter_x][iter_y] = ID_MC;
		}
	}

	cout << "Setting for efficiency correction factors (GtoH) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactorLoose_BtoF()
{
  TString Location = "/u/user/joonblee/Headers/effSF_muon/";
  cout<<endl<<"[Tag&Probe efficiency is from " << Location <<"]"<<endl;

  TFile *f1 = new TFile( Location + "EfficienciesAndSF_BCDEF.root" );

  TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
  TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

  Int_t nEtaBins1 = h_ID_data->GetNbinsX();
  Int_t nPtBins1 = h_ID_data->GetNbinsY();

  for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
  {
    for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
    {
      Int_t i_etabin = iter_x + 1;
      Int_t i_ptbin = iter_y + 1;

      Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
      Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

      Eff_LooseID_data_BtoF[iter_x][iter_y] = ID_data;
      Eff_LooseID_MC_BtoF[iter_x][iter_y] = ID_MC;
    }
  }

  cout << "Setting for loose muon efficiency correction factors (BtoF) is completed" << endl;
}

void DYAnalyzer::SetupEfficiencyScaleFactorLoose_GtoH()
{
  TString Location = "/u/user/joonblee/Headers/effSF_muon/";
  cout << "[Tag&Probe efficiency is from " << Location << "]" << endl;

  TFile *f1 = new TFile( Location+"EfficienciesAndSF_GH.root" );

  TH2F *h_ID_data = (TH2F*)f1->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/abseta_pt_DATA"); //Tight
  TH2F *h_ID_MC = (TH2F*)f1->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

  Int_t nEtaBins1 = h_ID_data->GetNbinsX();
  Int_t nPtBins1 = h_ID_data->GetNbinsY();

  for(Int_t iter_x = 0; iter_x < nEtaBins1; iter_x++)
  {
    for(Int_t iter_y = 0; iter_y < nPtBins1; iter_y++)
    {
      Int_t i_etabin = iter_x + 1;
      Int_t i_ptbin = iter_y + 1;

      Double_t ID_data = h_ID_data->GetBinContent(i_etabin, i_ptbin);
      Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);

      Eff_LooseID_data_GtoH[iter_x][iter_y] = ID_data;
      Eff_LooseID_MC_GtoH[iter_x][iter_y] = ID_MC;
    }
  }

  cout << "Setting for loose muon efficiency correction factors (GtoH) is completed" << endl;
}


Double_t DYAnalyzer::MuonSF_BtoF_ID( Muon mu )
{
  Double_t wt = -999;

  Double_t pt = mu.Pt;
  Double_t eta = mu.eta;

  Int_t ptbin_ID = Find_muon_PtBin_ID( pt );
  Int_t etabin_ID = Find_muon_EtaBin_ID( eta );

  Double_t sf = Eff_ID_data_BtoF[etabin_ID][ptbin_ID] / Eff_ID_MC_BtoF[etabin_ID][ptbin_ID];

  return sf;
}

Double_t DYAnalyzer::LooseMuonSF_BtoF_ID( Muon mu )
{
  Double_t wt = -999;

  Double_t pt = mu.Pt;
  Double_t eta = mu.eta;

  Int_t ptbin_ID = Find_muon_PtBin_ID( pt );
  Int_t etabin_ID = Find_muon_EtaBin_ID( eta );

  Double_t sf = Eff_LooseID_data_BtoF[etabin_ID][ptbin_ID] / Eff_LooseID_MC_BtoF[etabin_ID][ptbin_ID];

  return sf;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_ID( vector<Muon> muons )
{
  Double_t sf = 1;

  for(int i=0; i<muons.size(); i++) {
    sf = sf * MuonSF_BtoF_ID( muons[i] );
  }

  return sf;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_BtoF_LooseID( vector<Muon> muons )
{
  Double_t sf = 1;

  for(int i=0; i<muons.size(); i++) {
    sf = sf * LooseMuonSF_BtoF_ID( muons[i] );
  }

  return sf;
}


Double_t DYAnalyzer::MuonSF_GtoH_ID( Muon mu )
{
  Double_t wt = -999;

  Double_t pt = mu.Pt;
  Double_t eta = mu.eta;

  Int_t ptbin_ID = Find_muon_PtBin_ID( pt );
  Int_t etabin_ID = Find_muon_EtaBin_ID( eta );

  Double_t sf = Eff_ID_data_GtoH[etabin_ID][ptbin_ID] / Eff_ID_MC_GtoH[etabin_ID][ptbin_ID];

  return sf;
}

Double_t DYAnalyzer::LooseMuonSF_GtoH_ID( Muon mu )
{
  Double_t wt = -999;

  Double_t pt = mu.Pt;
  Double_t eta = mu.eta;

  Int_t ptbin_ID = Find_muon_PtBin_ID( pt );
  Int_t etabin_ID = Find_muon_EtaBin_ID( eta );

  Double_t sf = Eff_LooseID_data_GtoH[etabin_ID][ptbin_ID] / Eff_LooseID_MC_GtoH[etabin_ID][ptbin_ID];

  return sf;
}


Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_ID( vector<Muon> muons )
{
  Double_t sf = 1;

  for(int i=0; i<muons.size(); i++) {
    sf = sf * MuonSF_GtoH_ID( muons[i] );
  }

  return sf;
}

Double_t DYAnalyzer::EfficiencySF_EventWeight_HLT_GtoH_LooseID( vector<Muon> muons )
{
  Double_t sf = 1;

  for(int i=0; i<muons.size(); i++) {
    sf = sf * LooseMuonSF_GtoH_ID( muons[i] );
  }

  return sf;
}


Int_t DYAnalyzer::Find_muon_PtBin_ID(Double_t pt)
{
  const Int_t nPtBins = 6;
	Double_t PtBinEdges[nPtBins+1] = {20, 25, 30, 40, 50, 60, 120};

	Int_t ptbin = 9999;

	if( pt > PtBinEdges[nPtBins] )
		ptbin = nPtBins-1;
	else if( pt < PtBinEdges[0] )
		ptbin = 0;
	else
	{
		for(Int_t i=0; i<nPtBins; i++)
		{
			if( pt > PtBinEdges[i] && pt < PtBinEdges[i+1] )
			{
				ptbin = i;
				break;
			}
		}
	}

	return ptbin;
}

Int_t DYAnalyzer::Find_muon_EtaBin_ID(Double_t eta)
{
	const Int_t nEtaBins = 4;
	Double_t EtaBinEdges[nEtaBins+1] = {0, 0.9, 1.2, 2.1, 2.4};

	Int_t etabin = 9999;

	for(Int_t i=0; i<nEtaBins; i++)
	{
		if( fabs(eta) > EtaBinEdges[i] && fabs(eta) < EtaBinEdges[i+1] )
		{
			etabin = i;
			break;
		}
	}

	return etabin;
}
















