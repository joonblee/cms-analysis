#ifndef analysis_H
#define analysis_H

#include <TLorentzVector.h>
#include "NtupleMaker.h"

const double muon_mass     = 0.1056583715;
const double electron_mass = 0.000510998;

class PhysicsEvent : public NtupleEvent {
public:
  bool TriggerSelection( std::string trigger_ = "HLT_IsoMu24_v") {
    for( std::vector<NtupleTrigger>::const_iterator trigger = triggers.begin(); trigger != triggers.end(); ++trigger ) {
      if( trigger->name.find(trigger_) == 0 ) {
        return trigger->isFired;
      }
    }
    return false;
  }
};

/////////////////////////////////////////////////////////////////////////
// genparticles' mass //

class PhysicsGen : public NtupleGenParticle {
public:
  bool mu_acceptance(double pt_){
    if( fromHardProcessFinalState && fabs(id)==13 && pt > pt_ && fabs(eta) < 2.4 ) return true;
    else return false;
  }
  bool el_acceptance(double pt_) {
    if( fromHardProcessFinalState && fabs(id)==11 && pt > pt_ && (fabs(eta) < 1.4442 || (1.556 < fabs(eta) && fabs(eta) < 2.5)) ) return true;
    else return false;
  }
  bool el_candidate() {
    if( fromHardProcessFinalState && fabs(id)==11 ) return true;
    else return false;
  }
  TLorentzVector momentum(){
    TLorentzVector momentum_;
    momentum_.SetPtEtaPhiM(pt,eta,phi,muon_mass);
    return momentum_;
  }
  TVector3 xyz() {
    TVector3 xyz_;
    xyz_.SetXYZ(px,py,pz);
    return xyz_;
  }
};

/////////////////////////////////////////////////////////////////////////
class PhysicsElectron : public NtupleElectron {
public:
  bool acceptance(double pt_) { // accpetance //
    if( pt > pt_ && (fabs(eta)<1.4442 || (fabs(eta)>1.556 && fabs(eta)<2.5))) return true;
    else return false;
  } // acceptance //
  TLorentzVector momentum(){
    TLorentzVector momentum_;
    momentum_.SetPtEtaPhiM(pt,eta,phi,electron_mass);
    return momentum_;
  }

  bool accWithMediumId(double pt_) { // passMediumId //
    if( pt>pt_ ) {
      if(fabs(eta)<1.4442) {
        if(fabs(dEtaIn)<0.00311 && fabs(dPhiIn)<0.103 && sigmaIetaIeta<0.00998 && hOverE<0.253 && isoRho<0.0695 && fabs(ooEmooP)<0.134 && expectedMissingInnerHits<=1 && passConversionVeto) return true;
        else return false;
      }
      else if(fabs(eta)>1.556 && fabs(eta)<2.5) {
        if(fabs(dEtaIn)<0.00609 && fabs(dPhiIn)<0.045 && sigmaIetaIeta<0.0298 && hOverE<0.0878 && isoRho<0.0821 && fabs(ooEmooP)<0.13 && expectedMissingInnerHits<=1 && passConversionVeto) return true;
        else return false;
      }
      else return false;
    } 
    else return false;
  }

  bool accWMediumIdNm1Cut(double pt_ ,int ncut) {
    bool tResult = true;
    vector<bool> result;
    if( pt>pt_ ) {
      if(fabs(eta)<1.4442) {
        if( ncut!=0 ) {
          if(hOverE<0.253) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=1 ) {
          if(fabs(ooEmooP)<0.134) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=2 ) {
          if(expectedMissingInnerHits<=1) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=3 ) {
          if(passConversionVeto) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=4 ) {
          if(fabs(dPhiIn)<0.103) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=5 ) {
          if(fabs(dEtaIn)<0.00311) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=6 ) {
          if(sigmaIetaIeta<0.00998) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=7 ) {
          if(isoRho<0.0695) result.push_back(true); else result.push_back(false);
        }
      }
      else if(fabs(eta)>1.556 && fabs(eta)<2.5) {
        if( ncut!=0 ) {
          if(hOverE<0.0878) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=1 ) {
          if(fabs(ooEmooP)<0.13) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=2 ) {
          if(expectedMissingInnerHits<=1) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=3 ) {
          if(passConversionVeto) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=4 ) {
          if(fabs(dPhiIn)<0.045) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=5 ) {
          if(fabs(dEtaIn)<0.00609) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=6 ) {
          if(sigmaIetaIeta<0.0298) result.push_back(true); else result.push_back(false);
        }
        if( ncut!=7 ) {
          if(isoRho<0.0821) result.push_back(true); else result.push_back(false);
        }
      }
      else result.push_back(false);
    }
    else result.push_back(false);
    for(int i = 0; i != result.size(); i++) {
      tResult = tResult*result[i];
    }
    return tResult;
  }

  bool passMediumIdFlow(int flow) { // cut flow efficiency //
    if( fabs(eta)<1.4442 ) {
      if( flow == 0 ) return true;
      else if( !(hOverE<0.253) ) return false; else if( flow == 1 ) return true;
      else if( !(fabs(ooEmooP)<0.134) ) return false; else if( flow == 2 ) return true;
      else if( !(expectedMissingInnerHits<=1) ) return false; else if( flow == 3 ) return true;
      else if( !(passConversionVeto) ) return false; else if( flow == 4 ) return true;
      else if( !(fabs(dPhiIn)<0.103) ) return false; else if( flow == 5 ) return true;
      else if( !(fabs(dEtaIn)<0.00311) ) return false; else if( flow == 6 ) return true;
      else if( !(sigmaIetaIeta<0.00998) ) return false; else if( flow == 7 ) return true;
      else if( !(isoRho<0.0695) ) return false; else if( flow == 8 ) return true;
      else return false;
    }
    else if( (fabs(eta)>1.556 && fabs(eta)<2.5) ) {
      if( flow == 0 ) return true;
      else if( !(hOverE<0.0878) ) return false; else if( flow == 1 ) return true;
      else if( !(fabs(ooEmooP)<0.13) ) return false; else if( flow == 2 ) return true;
      else if( !(expectedMissingInnerHits<=1) ) return false; else if( flow == 3 ) return true;
      else if( !(passConversionVeto) ) return false; else if( flow == 4 ) return true;
      else if( !(fabs(dPhiIn)<0.045) ) return false; else if( flow == 5 ) return true;
      else if( !(fabs(dEtaIn)<0.00609) ) return false; else if( flow == 6 ) return true;
      else if( !(sigmaIetaIeta<0.0298) ) return false; else if( flow == 7 ) return true;
      else if( !(isoRho<0.0821) ) return false; else if( flow == 8 ) return true;
      else return false;
    }
    else return false;
  } // pass Medium Id //
}; // PhsicsElectron //

class PhysicsMuon : public NtupleMuon {
public:
/*
  bool acceptance(double pt_, double eta_){
    if( pt > pt_ && fabs(eta) < eta_ ) return true;
    else return false;
  }
*/
  bool acceptance(double pt_){
    if( pt > pt_ && fabs(eta) < 2.4 ) return true;
    else return false;
  }
  TLorentzVector momentum(){
    TLorentzVector momentum_;
    momentum_.SetPtEtaPhiM(pt,eta,phi,muon_mass);
    return momentum_;
  }
  TVector3 xyz() {
    TVector3 xyz_;
    xyz_.SetXYZ(px,py,pz);
    return xyz_;
  }

  bool highPtMuonId(int ncut = 0) {
//    if( isGlobalMuon && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX)<0.2 && muonBestTrack_ptError/muonBestTrack_pt<0.3 )
    if(ncut==10) return true;
    else if(fabs(dxyVTX)>0.2||fabs(dzVTX)>0.5) return false; else if(ncut==1) return true;
    else if(nValidPixelHits<1) return false; else if(ncut==2) return true;
    else if(nTrackerLayers<6) return false; else if(ncut==3) return true;
    else if(muonBestTrack_ptError/muonBestTrack_pt>0.3) return false; else if(ncut==4) return true;
    else if(!isGlobalMuon) return false; else if(ncut==5) return true;
    else if(nValidMuonHits<1) return false; else if(ncut==6) return true;
    else if(nMatchedStations<2) return false; else if(ncut==7) return true;
    else if( isGlobalMuon && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX)<0.2 && fabs(dzVTX)<0.5 && muonBestTrack_ptError/muonBestTrack_pt<0.3 ) return true;
    else 
    return false;
  }

  bool looseMuonId() {
    if( isPFMuon && (isGlobalMuon||isTrackerMuon) && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5 )
      return true;
    else
      return false;
  }

  bool tightMuonId(int ncut = 0) {
/*                if( isPFMuon && isGlobalMuon && normalizedChi2 < 10 && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5 )
                        return true;
*/
    if(ncut==10) return true;
    else if(fabs(dxyVTX)>0.2) return false; else if(ncut==1) return true;
    else if(fabs(dzVTX)>0.5) return false; else if(ncut==2) return true;
    else if(nValidPixelHits<1) return false; else if(ncut==3) return true;
    else if(nTrackerLayers<6) return false; else if(ncut==4) return true;
    else if(!isPFMuon) return false; else if(ncut==5) return true;    
    else if(nMatchedStations<2) return false; else if(ncut==6) return true;
    else if(nValidMuonHits<1) return false; else if(ncut==7) return true;
    else if(normalizedChi2>10) return false; else if(ncut==8) return true;
    else if(!isGlobalMuon) return false; else if(ncut==9) return true;
    else if(isPFMuon && isGlobalMuon && normalizedChi2 < 10 && nValidMuonHits>0 && nMatchedStations>1 && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5) return true;
    else
    return false;
  }
  bool trackerTightMuonId() {
   if( isPFMuon && nValidPixelHits>0 && nTrackerLayers>5 && fabs(dxyVTX) < 0.2 && fabs(dzVTX) < 0.5 )
     return true;
   else
     return false;
  }
  bool trackerTightMuonIdFlow(int ncut) {
    if(ncut==0) return true;
    else if(!(isGlobalMuon||isTrackerMuon)) return false; else if(ncut==1) return true;
    else if(fabs(dxyVTX)>0.2) return false; else if(ncut==2) return true;
    else if(fabs(dzVTX)>0.5) return false; else if(ncut==3) return true;
    else if(nValidPixelHits<1) return false; else if(ncut==4) return true;
    else if(nTrackerLayers<6) return false; else if(ncut==5) return true;
    else if(!isPFMuon) return false; else if(ncut==6) return true;
    else return true;
  }
  bool trackerTightMuonIdNm1Cut(int ncut) {
    bool tResult = true;
    vector<bool> result;
    if( isGlobalMuon||isTrackerMuon ) result.push_back(true); else result.push_back(false);
    if( ncut!=0 ) {
      if(fabs(dxyVTX)<0.2) result.push_back(true); else result.push_back(false);
    }
    if( ncut!=1 ) {
      if(fabs(dzVTX)<0.5) result.push_back(true); else result.push_back(false);
    }
    if( ncut!=2 ) {
      if(nValidPixelHits>0) result.push_back(true); else result.push_back(false);
    }
    if( ncut!=3 ) {
      if(nTrackerLayers>5) result.push_back(true); else result.push_back(false);
    }
    if( ncut!=4 ) {
      if(isPFMuon) result.push_back(true); else result.push_back(false);
    }
    for(int i = 0; i != result.size(); i++) {
      tResult = tResult*result[i];
    }
    return tResult;
  }

  bool isolation(double iso) {
    if( isolationR03_sumpt/pt < iso) return true;
    else return false;
  }
};

class Tester {
public:
  int dR_tester(double deltaR) {
    if     (             deltaR<0.01) return 0;
    else if(0.01<deltaR&&deltaR<0.02) return 1;
    else if(0.02<deltaR&&deltaR<0.05) return 2;
    else if(0.05<deltaR&&deltaR<0.1 ) return 3;
    else if(0.1 <deltaR             ) return 4;
    else return 5;
  }
};

#endif

