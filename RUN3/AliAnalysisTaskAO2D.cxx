#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisManager.h"

//#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVEventHandler.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliESDInputHandler.h"
#include "AliAO2DInputHandler.h"

#include "AliAnalysisTaskAO2D.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

ClassImp(AliAnalysisTaskAO2D)

//________________________________________________________________________
AliAnalysisTaskAO2D::AliAnalysisTaskAO2D(const char *name) 
  : AliAnalysisTaskSE(name)
{
  // Constructor
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskAO2D::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fListOut = new TList();
  fListOut->SetOwner();
  fListOut->SetName("listHistos");

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);

  fHistQ = new TH1F("fHistQ", "TPC clusters: Q distribution", 1000, 0, 10000);
  fHistQ->GetXaxis()->SetTitle("Q");
  fHistQ->GetYaxis()->SetTitle("dN/dQ");
  fHistQ->SetMarkerStyle(kFullCircle);

  fHistNTPCCl = new TH1F("fHistNTPCCl", "Number of TPC clusters", 160, -0.5, 159.5);
  fHistNTPCCl->GetXaxis()->SetTitle("n. TPC Cl.");
  fHistNTPCCl->GetYaxis()->SetTitle("dN/d(n. TPC Cl)");
  fHistNTPCCl->SetMarkerStyle(kFullCircle);

  fHistNESDtracks = new TH1F("fHistNESDtracks", "Number of ESD tracks", 1000, -0.5, 999.5);
  fHistNESDtracks->GetXaxis()->SetTitle("n. ESD tracks");
  fHistNESDtracks->GetYaxis()->SetTitle("dN/d(n. ESD tracks)");
  fHistNESDtracks->SetMarkerStyle(kFullCircle);

  fListOut->Add(fHistPt);
  fListOut->Add(fHistQ);
  fListOut->Add(fHistNTPCCl);
  fListOut->Add(fHistNESDtracks);

  PostData(1, fListOut);
}

//________________________________________________________________________
void AliAnalysisTaskAO2D::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliAO2DInputHandler *ih = dynamic_cast<AliAO2DInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (ih) {
    UserExecAO2D();
    return;
  }

  Int_t nESDtracks = fInputEvent->GetNumberOfTracks();
  Printf("There are %d tracks in this event", nESDtracks);

  fHistNESDtracks->Fill(nESDtracks);

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < nESDtracks; iTracks++) {
    const AliVTrack* track = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(iTracks));
    fHistPt->Fill(track->Pt());
    fHistNTPCCl->Fill(track->GetTPCNcls());
  } //track loop 

  // Post output data.
  PostData(1, fListOut);
  fEv++;
}  

void AliAnalysisTaskAO2D::UserExecAO2D()
{
   AliAO2DInputHandler *ih = dynamic_cast<AliAO2DInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   
   //using Event_t = AliAO2DTypes::Event_t;
   //using Track_t = AliAO2DTypes::Track_t;

   auto event = ih->CurrentEvent();
   auto tracks = ih->GetTracks();
   int id = 0;
   printf("Current event %lu: %lld tracks\n", ih->GetCurrentEntry(), tracks.size());

   //printf("  %d: %d %g\n", id++, tracks[0].fCollisionsID, tracks[0].fX);
   for (const auto &track : tracks) {
      double pt = TMath::Abs(1./track.fSigned1Pt);
      fHistPt->Fill(pt);
   }
}

//________________________________________________________________________
void AliAnalysisTaskAO2D::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fListOut = dynamic_cast<TList*> (GetOutputData(1)); 

  if (fListOut) {
    fHistPt = dynamic_cast<TH1F*>(fListOut->FindObject("fHistPt")); 
    if (!fHistPt) {
      Printf("ERROR: fHistPt not available");
      return;
    }
   
    TCanvas *c1 = new TCanvas("AliAnalysisTaskAO2D","Pt",10,10,510,510);
    c1->cd(1)->SetLogy();
    fHistPt->DrawCopy("E");
  }
  else {
    Printf("In Terminate: no TList found");
  }

}
