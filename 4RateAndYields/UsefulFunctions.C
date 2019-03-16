#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TVector2.h>
#include "/u/user/joonblee/Headers/Object.h"
#include "/u/user/joonblee/Headers/NtupleHandle.h"

#include <iostream>

#define M_Z 91.2

#define trig "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
#define trig_ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"

double second_isolation(NtupleHandle* ntuple, Muon mu) {
  double isolation = mu.trkiso;
  double dR = 9999;
  for(unsigned i=0; i<ntuple->nMuon; i++) {
    Muon smu; smu.FillFromNtuple(ntuple,i);
    if( !(smu.testMuonIDv2() && smu.acceptance(5,2.5)) ) continue;
    /*dR = smu.Momentum.DeltaR(mu.Momentum);*/
    dR = TMath::Sqrt( (mu.Inner_eta - smu.Inner_eta)*(mu.Inner_eta - smu.Inner_eta) + TVector2::Phi_mpi_pi(mu.Inner_phi - smu.Inner_phi)*TVector2::Phi_mpi_pi(mu.Inner_phi - smu.Inner_phi) );
    if( 0.01 < dR && dR < 0.3 ) {
      isolation = mu.trkiso - smu.Inner_pT/mu.Pt;
      if(isolation < 0) isolation = 0;
      break;
    }
  }
  return isolation;
}

bool tracker_muon_test(NtupleHandle* ntuple, Muon mu) {
  bool test = false;
  if( mu.isGlobalMuon ) test = true;
  else if( mu.isTrackerMuon ) {
    for(unsigned i=0; i<ntuple->nMuon; i++) {
      Muon mu_; mu_.FillFromNtuple(ntuple,i);
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
      if( (2.9 < M && M < 3.3) /*|| M < 0.4*/ ) {
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

/*
double FR_template(Muon mu) {
  double pt = mu.Pt; double eta = mu.eta;
  double fakerate = -999;

  const int ptbinnum = 15;
  double ptbin[ptbinnum+1] = {5,6,7,8,9,10,12,15,20,25,32,40,52,65,100,500};
  const int ptbinnum_endcap = 12;
  double ptbin_endcap[ptbinnum_endcap+1] = {5,6,7,8,9,10,12,15,17,25,43,80,500};

  if( fabs(eta) < 1.2 ) { // Barrel
    double FR[] = {0.0626778,0.0681467,0.0679842,0.0693651,0.0636564,0.0604347,0.0537839,0.0451948,0.0280992,0.0303409,0.0456232,0.0439918,0.0570109,0.0554287,0.0363951};

    for(int i=0; i < ptbinnum; i++) {
      if( ptbin[i] < pt && pt < ptbin[i+1] ) fakerate = FR[i];
    }
    if( 100 < pt ) fakerate = 0.0554;
  }
  else { // Endcap
    double FR[] = {0.0753232,0.0758815,0.0796902,0.0768849,0.0759609,0.068662,0.0691455,0.0561826,0.0549706,0.0318213,0.0670364,0.0358052};

    for(int i=0; i < ptbinnum_endcap; i++) {
      if( ptbin_endcap[i] < pt && pt < ptbin_endcap[i+1] ) fakerate = FR[i];
    }
    if( 100 < pt ) fakerate = 0.067;
  }

  return fakerate;
}
*/

double FR_template(vector<Muon> muons) {
	double fakerate = 1;

  const int ptbinnum = 15;
  double ptbin[ptbinnum+1] = {5,6,7,8,9,10,12,15,20,25,32,40,52,65,100,500};
  const int ptbinnum_endcap = 12;
  double ptbin_endcap[ptbinnum_endcap+1] = {5,6,7,8,9,10,12,15,17,25,43,80,500};

	double FR_barrel[] = {0.0626778,0.0681467,0.0679842,0.0693651,0.0636564,0.0604347,0.0537839,0.0451948,0.0280992,0.0303409,0.0456232,0.0439918,0.0570109,0.0554287,0.0554};
	double FR_endcap[] = {0.0753232,0.0758815,0.0796902,0.0768849,0.0759609,0.068662,0.0691455,0.0561826,0.0549706,0.0318213,0.0670364,0.067};

	for(int i=0; i<muons.size(); i++) {
  	double pt = muons[i].Pt; double eta = muons[i].eta;
	  if( fabs(eta) < 1.2 ) { // Barrel
	    for(int i=0; i < ptbinnum; i++) {
	      if( ptbin[i] < pt && pt < ptbin[i+1] ) {
					fakerate = fakerate*FR_barrel[i]/(1-FR_barrel[i]);
					break;
				}
	    }
	    if( 500 < pt /*&& fakerate == -999*/ ) fakerate = fakerate*FR_barrel[ptbinnum-1]/(1-FR_barrel[ptbinnum-1]);
	  }
	  else { // Endcap
	    for(int i=0; i < ptbinnum_endcap; i++) {
	      if( ptbin_endcap[i] < pt && pt < ptbin_endcap[i+1] ) {
					fakerate = fakerate*FR_endcap[i]/(1-FR_endcap[i]);
					break;
	    	}
			}
	    if( 500 < pt /*&& fakerate == -999*/ ) fakerate = fakerate*FR_endcap[ptbinnum_endcap-1]/(1-FR_endcap[ptbinnum_endcap-1]);
	  }
	}

  return fakerate;
}

double Boosted_FR_template(vector<Muon> muons) {
  double fakerate = 1;

  //const int ptbinnum = 11;
  //double ptbin[ptbinnum+1] = {5,7,10,15,20,25,32,40,52,65,100,500};
  //const int ptbinnum_endcap = 8;
  //double ptbin_endcap[ptbinnum_endcap+1] = {5,7,10,15,20,25,43,100,500};

  //double FR_barrel[] = {0.0926386,0.117703,0.110768,0.109505,0.112307,0.0780903,0.0648539,0.0436538,0.0555524,0.00693524,0.0117617};
  //double FR_endcap[] = {0.050443,0.0653745,0.0730131,0.0575,0.0575,0.0297531,0.0142844,0.00989454};


  const int ptbinnum = 7;
  double ptbin[ptbinnum+1] = {5,7,10,15,22,40,100,500};
  const int ptbinnum_endcap = 7;
  double ptbin_endcap[ptbinnum_endcap+1] = {5,7,10,15,22,40,100,500};


  double FR_barrel[] = {0.0926386,0.117703,0.110768,0.111212,0.0815894,0.0370332,0.0117617};
  double FR_endcap[] = {0.050443,0.0653745,0.0730131,0.0547601,0.0377863,0.0153047,0.00989454};

  for(int i=0; i<muons.size(); i++) {
    double pt = muons[i].Pt; double eta = muons[i].eta;
    if( fabs(eta) < 1.2 ) { // Barrel
      for(int i=0; i < ptbinnum; i++) {
        if( ptbin[i] < pt && pt < ptbin[i+1] ) {
          fakerate = fakerate*FR_barrel[i]/(1-FR_barrel[i]);
          break;
        }
      }
      if( 500 < pt ) fakerate = fakerate*FR_barrel[ptbinnum-1]/(1-FR_barrel[ptbinnum-1]);
    }
    else { // Endcap
      for(int i=0; i < ptbinnum_endcap; i++) {
        if( ptbin_endcap[i] < pt && pt < ptbin_endcap[i+1] ) {
          fakerate = fakerate*FR_endcap[i]/(1-FR_endcap[i]);
          break;
        }
      }
      if( 500 < pt ) fakerate = fakerate*FR_endcap[ptbinnum_endcap-1]/(1-FR_endcap[ptbinnum_endcap-1]);
    }
  }

  return fakerate;
}


bool boost_selection(vector<Muon> muons, vector<Muon> muons_) {
	muons.insert(muons.end(), muons_.begin(), muons_.end());

/*
	for(int i=0; i<muons.size(); i++) {
		Muon mu = muons[i];
		for(int j=i+1; j<muons.size(); j++) {
			Muon mu_ = muons[j];
			if( mu.charge != mu_.charge ) {
				double mass = (mu.Momentum + mu_.Momentum).M();
				if( 40 < mass && mass < 120 ) {
*/
					for(int k=0; k<muons.size(); k++) {
						//if( k==i || k==j ) continue;
						Muon m = muons[k];
						for(int l=k+1; l<muons.size(); l++) {
							//if( l==i || l==j ) continue;
							Muon m_ = muons[l];

							if( (m.Momentum + m_.Momentum).M() < 4
									&& m.Momentum.DeltaR( m_.Momentum ) < 0.3 )
								return true;
						}
					}
//				}
//			}
//		}
//	}

	return false;
}


bool DYSoftGammaSelector(NtupleHandle* ntuple) {
	bool NoPromptFSGamma = true;
	for(int i = 0; i < ntuple->nGenOthers; i++) {
		GenOthers gn; gn.FillFromNtuple( ntuple , i );
		if( gn.ID == 22 && gn.isPromptFinalState ) {
			NoPromptFSGamma = false;
			break;
		}
	}
	if( NoPromptFSGamma ) {
		cout<<endl<<"No promptFS Gamma"<<endl;
		return true;
	}

	GenLepton PFSmu, PFSmu_;
	GenLepton MotherMu, MotherMu_;
	int nPFSmu = 0; bool TMom1 = false; bool TMom2 = false;

	for(int i = ntuple->gnpair - 1; i >= 0; i--) {
		GenLepton m; m.FillFromNtuple( ntuple , i );
		if( m.ID == 13 && m.isPromptFinalState && !m.fromHardProcessFinalState) {
			PFSmu = m;
			nPFSmu++;
			if( abs(m.Mother) == 13 ) TMom1 = true;
			else MotherMu = m;
		}
		else if( m.ID == -13 && m.isPromptFinalState && !m.fromHardProcessFinalState) {
			PFSmu_ = m;
			nPFSmu++;
			if( abs(m.Mother) == 13 ) TMom2 = true;
			else MotherMu_ = m;
		}
		if( nPFSmu == 2 ) break;
	}

	if( nPFSmu == 2 ) {
		if( TMom1 ) {
			bool TMom3 = false;
			for(int i = ntuple->gnpair - 1; i >= 0; i--) {
				GenLepton m; m.FillFromNtuple( ntuple , i );
				if( m.Pt == PFSmu.Mother_Pt ) {
					MotherMu = m;
					if( abs(m.Mother) == 13 ) TMom3 = true;
					break;
				}
			}
			if( TMom3 ) {
				for(int i = ntuple->gnpair - 1; i >= 0; i--) {
					GenLepton m; m.FillFromNtuple( ntuple , i );
					if( m.Pt == MotherMu.Mother_Pt ) {
						MotherMu = m;
						if( abs(m.Mother) != 13 ) break;
					}
				}
			}
		}
		if( TMom2 ) {
			bool TMom4 = false;
      for(int i = 0; i < ntuple->gnpair; i++) {
        GenLepton m; m.FillFromNtuple( ntuple , i );
        if( m.Pt == PFSmu_.Mother_Pt ) {
          MotherMu_ = m;
          break;
        }
      }
      if( TMom4 ) {
        for(int i = ntuple->gnpair - 1; i >= 0; i--) {
          GenLepton m; m.FillFromNtuple( ntuple , i );
          if( m.Pt == MotherMu_.Mother_Pt ) {
            MotherMu_ = m;
            if( abs(m.Mother) != 13 ) break;
          }
        }
      }
		}
	}
	else {
		cout<<endl<<"nPFSmu != 2"<<endl;
		return false;
	}

	for(int i = 0; i < ntuple->nGenOthers; i++) {
		GenOthers gm; gm.FillFromNtuple( ntuple , i );
		if( gm.ID == 22 && gm.Pt < 10
				&& gm.Pt == MotherMu.Pt && gm.Pt == MotherMu_.Pt )
			return true;
	}
	return false;
}


