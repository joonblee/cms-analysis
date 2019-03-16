import ROOT

from Validation.RecoTrack.plotting.ntuple import *

def main():
  ntuple = TrackingNtuple("trackingNtuple.root")

  iEvent = 0
  totTP = 0
  totTrk = 0

  print "Start iteration ..."
  for event in ntuple:
    vertices = event.vertices()

    ### --------- 1. event selector ----------- ###
    ### Only choose the events
    ### with 4 leptons   
    ### --------------------------------------- ###
    tps = event.trackingParticles()

    nLeptons = 0
    for tp in tps:
      if (tp.pdgId() == 13 or tp.pdgId() == -13 or tp.pdgId() == 11 or tp.pdgId() == -11) and tp.pt() > 10:
        nLeptons += 1
    if not nLeptons == 4: continue

    iEvent += 1
    if iEvent > 1: break

    print "\n\nevent", iEvent
    # -------------------------------------- 1. event selector END

    ### ---------- 2. sim track ---------- ###
    ### links from tracking particles      ###
    ### to tracks                          ###
    ### ---------------------------------- ###
    print "\n --- sim track --- "
    totTP += len(tps)

    nMatchedTP = 0
    nDupTP = 0
    for tp in tps:
      if (tp.pdgId() == 13 or tp.pdgId() == -13 or tp.pdgId() == 11 or tp.pdgId() == -11) and tp.pt() > 10:        
        print "TP ID", tp.pdgId(), "-> (",tp.pt(),tp.eta(),tp.phi(),")"
        if tp.nMatchedTracks() == 0: continue
        nMatchedTP += 1
        if tp.nMatchedTracks() >= 2:
          nDupTP += 1
          print "  - There are 2 tracks matched w/ 1 TP"

        #for hit in tp.hits():
        #  print hit.layerStr(), hit.detId(), " x", hit.x(), " y", hit.y(), " tof", hit.tof()
    print "# of TP : ", len(tps)
    print "# of matched TP (w/ reco track) : ", nMatchedTP
    print "# of duplicated TP              : ", nDupTP
    # -------------------- 2. sim track END

    ### ---------- 3. reco track ---------- ###
    ### links from reco tracks to TP        ###
    ### ----------------------------------- ###
    print "\n --- reco track --- "
    tracks = event.tracks()
    totTrk += len(tracks)

    for track in tracks:
      
      if track.nMatchedTrackingParticles == 0: continue
      for tpInfo in track.matchedTrackingParticleInfos():
        tp = tpInfo.trackingParticle()
        if not ((tp.pdgId() == 13 or tp.pdgId() == -13 or tp.pdgId() == 11 or tp.pdgId() == -11) and tp.pt() > 10): continue
        print "trk    -> (",track.pt(),track.eta(),track.phi(),")"
        print "TP ID", tp.pdgId(), "-> (",tp.pt(),tp.eta(),tp.phi(),")"
        if tp.nMatchedTracks() > 1:
          print "  - track is duplicated"
        if tp.parentVertex().nSourceTrackingParticles() > 0:
          print "  - track is secondary"

        nRecHit = 0
        nSimHit = 0
        if ntuple.hasHits():
          pix_simTrkIdx = set()
          for hit in track.pixelHits():
            nRecHit+=1
            
            if hit.nSimHits() >= 1:
              nSimHit+=1
              for simHit in hit.simHits():
                simHitMatTP = simHit.trackingParticle()
                print "----> sim hit matched tracking particle pt :",simHitMatTP.pt()
        print "  - reco hits in tracks :",nRecHit
        print "  - sim  hits matched   :",nSimHit

    print "# of trk : ", len(tracks)
    # -------------------- 3. reco track END

    ### hits
    if ntuple.hasHits():
      for hit in event.pixelHits():
        for track in hit.tracks():
          print "hit trk pt :",track.pt()

    """
    tracks = event.tracks()
    ntracks = len(tracks)
    tot_ntracks += ntracks
    tot_pv_ntracks += vertices[0].nTracks()
    """

  print "total # of evts      : ", iEvent
  print "total # of tp        : ", totTP
  print "ave # of tp in 1 evt : ", float(totTP)/iEvent
  print "total # of trk       : ", totTrk
  print "ave # of trk ~       : ", float(totTrk)/iEvent

if __name__ == "__main__":
  main()


