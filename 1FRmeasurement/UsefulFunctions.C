#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TVector2.h>
#include "/u/user/joonblee/Headers/Object.h"
#include "/u/user/joonblee/Headers/NtupleHandle.h"
#include "/u/user/joonblee/Headers/RoCorr/RoccoR.cc"

#include <iostream>

#define M_Z 91.2

#define trig "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
#define trig_ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"

double second_isolation(vector<Muon> MuonCollection, Muon mu) {
  double isolation = mu.trkiso;
  double dR = 9999;
  for(unsigned i=0; i<MuonCollection.size(); i++) {
    Muon smu = MuonCollection[i];
    if( !(smu.testMuonIDv2() && smu.acceptance(5,2.5)) ) continue;
    dR = TMath::Sqrt( (mu.Inner_eta - smu.Inner_eta)*(mu.Inner_eta - smu.Inner_eta) + TVector2::Phi_mpi_pi(mu.Inner_phi - smu.Inner_phi)*TVector2::Phi_mpi_pi(mu.Inner_phi - smu.Inner_phi) );
    if( 0.01 < dR && dR < 0.3 ) {
      isolation = mu.trkiso - smu.Inner_pT/mu.Pt;
      if(isolation < 0) isolation = 0;
      break;
    }
  }
  return isolation;
}


bool tracker_muon_test(vector<Muon> MuonCollection, Muon mu) {
  bool test = false;
  if( mu.isGlobalMuon ) test = true;
  else if( mu.isTrackerMuon ) {
    for(unsigned i=0; i<MuonCollection.size(); i++) {
      Muon mu_ = MuonCollection[i];
      if( mu.Momentum.DeltaR(mu_.Momentum) < 0.3 // dR(trk,glb) < 0.3
          && mu_.acceptance(5,2.7)               // min cut of glb
          && mu_.isGlobalMuon
          && mu_.testMuonIDv2() ) {
        test = true;
        break;
      }
    }
  }
  return test;
}

vector<GenLepton> DYspliter(NtupleHandle* ntuple){
  vector<GenLepton> GenLeptonCollection;
  for(int j=0; j<ntuple->gnpair; j++) {
    GenLepton genlep; genlep.FillFromNtuple( ntuple , j );
    if( genlep.isMuon()
        && genlep.isPromptFinalState ) {// acceptance + fromHardProcessFinalState
     GenLeptonCollection.push_back( genlep );
    }
  }
  return GenLeptonCollection;
}

vector<GenLepton> DiMuFromMergedJet(NtupleHandle* ntuple){
  vector<GenLepton> GenLeptonCollection;
  for(int j=0; j<ntuple->gnpair; j++) {
    GenLepton genlep; genlep.FillFromNtuple( ntuple , j );
    if( genlep.isMuon()
        && genlep.isLastCopy ) {// acceptance + fromHardProcessFinalState
     GenLeptonCollection.push_back( genlep );
    }
  }
  return GenLeptonCollection;
}



bool QCD_suppression(vector<Muon> pBox) {
  bool test = true;
  for(unsigned i=0; i<pBox.size(); i++) {
    for(unsigned j=i+1; j<pBox.size(); j++) {
      if( pBox[i].charge == pBox[j].charge ) continue;
      if( (pBox[i].Momentum + pBox[j].Momentum).M() < 4 ) {
        test = false;
        break;
      }
    }
    if( !test ) break;
  }
  return test;
}


bool JPsi_suppression(vector<Muon> pBox) {
  bool test = true;
  for(unsigned i=0; i<pBox.size(); i++) {
    for(unsigned j=i+1; j<pBox.size(); j++) {
      if( pBox[i].charge == pBox[j].charge ) continue;
			double M = (pBox[i].Momentum + pBox[j].Momentum).M();
      if( (3.0 < M && M < 3.2) /*|| M < 0.4*/ ) {
        test = false;
        break;
      }
    }
    if( !test ) break;
  }
  return test;
}


void time() {
  time_t t = time(0);
  struct tm * now = localtime( & t );
  cout<<(now->tm_year + 1900)<<"-"
      <<(now->tm_mon + 1)<<"-"
      <<now->tm_mday<<" : "
      <<now->tm_hour<<"-"
      <<now->tm_min<<"-"
      <<now->tm_sec
      <<endl;
}


bool boost_selection(vector<Muon> muons, vector<Muon> muons_) {
	muons.insert(muons.end(), muons_.begin(), muons_.end());

	for(int k=0; k<muons.size(); k++) {
		Muon m = muons[k];
		for(int l=k+1; l<muons.size(); l++) {
			Muon m_ = muons[l];

			if( (m.Momentum + m_.Momentum).M() < 4
					&& m.Momentum.DeltaR( m_.Momentum ) < 0.3 )
				return true;
		}
	}
	return false;
}


bool DYVirtualGammaSelector(NtupleHandle* ntuple) {
  bool TT = false;
  GenLepton mu, mu_;
  double pt = 0, pt_ = 0;
  for(int i = ntuple->gnpair - 1; i >= 0; i--) {
    GenLepton m; m.FillFromNtuple( ntuple , i );
    if( m.ID == 22 && m.fromHardProcessFinalState ) {
      return false;
    }
    if( abs( m.ID ) == 13 && m.isPromptFinalState && !m.fromHardProcessFinalState ) {
      TT = true;
    }
  }

  if( !TT ) return false;

  for(int i = ntuple->gnpair - 1; i >= 0; i--) {
    GenLepton m; m.FillFromNtuple( ntuple , i );
    if( !( abs(m.ID) == 13 && m.isPrompt ) ) continue;
    if( m.isPromptFinalState && !m.fromHardProcessFinalState ) {
      if( m.ID == 13 ) {
        mu = m;
        pt = m.Mother_Pt;
      }
      else if( m.ID == -13 ) {
        mu_ = m;
        pt_ = m.Mother_Pt;
      }
    }
    if( pt == m.Pt ) {
      mu = m; pt = m.Mother_Pt;
    }
    else if( pt_ == m.Pt ) {
      mu_ = m; pt_ = m.Mother_Pt;
    }

    if( mu.Mother == 22 || mu_.Mother == 22 )
      return false;
  }
  return true;
}


void RochesterCorrection(NtupleHandle* ntuple, TString dataSet, TRandom* r1, vector<Muon> *add_MuonCollection, vector<Muon> *add_CorrMuonCollection) {
  add_MuonCollection->clear(); add_CorrMuonCollection->clear();

	RoccoR rc("../Headers/RoCorr/rcdata.2016.v3");

  for(int j=0; j<ntuple->nMuon; j++) {
    Muon mu; mu.FillFromNtuple( ntuple , j );

		if( mu.Pt < 4 || fabs(mu.eta) > 2.7 || !(mu.isGlobalMuon || mu.isTrackerMuon) ) continue;

    add_MuonCollection->push_back(mu);

		Double_t rndm[2]; r1->RndmArray(2, rndm);
    double SF = 0;
    Int_t s, m;

    if(dataSet.Index("data") >= 0)
      SF = rc.kScaleDT(mu.charge, mu.Pt, mu.eta, mu.phi, s=0, m=0);
    else 
      SF = rc.kScaleAndSmearMC(mu.charge, mu.Pt, mu.eta, mu.phi, mu.trackerLayers, rndm[0], rndm[1], s=0, m=0);

    mu.Pt = SF*mu.Pt;

    TLorentzVector v;
    v.SetPtEtaPhiM( mu.Pt , mu.eta , mu.phi , M_Mu );

    mu.Momentum = v;

    add_CorrMuonCollection->push_back(mu);
  }
}








