/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <TTree.h>
#include <TFile.h>

#include "AliAO2DInputHandler.h"
#include <AliLog.h>

ClassImp(AliAO2DInputHandler)

static Option_t *gAO2DDataType = "AO2D";

//______________________________________________________________________________
AliAO2DInputHandler::AliAO2DInputHandler() :
    AliInputEventHandler()
{
  /// Default constructor
  for (int i = 0; i < AliAO2DTypes::kTrees; ++i) {
    fTrees[i] = nullptr;
    fReadTree[i] = false;
  }
}

//______________________________________________________________________________
AliAO2DInputHandler::AliAO2DInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title)
{
    /// Constructor
  for (int i = 0; i < AliAO2DTypes::kTrees; ++i) {
    fTrees[i] = nullptr;
    fReadTree[i] = false;
  }
}

//______________________________________________________________________________
AliAO2DInputHandler::~AliAO2DInputHandler() 
{
/// Destructor

}

//______________________________________________________________________________
Bool_t AliAO2DInputHandler::Init(TTree* tree, Option_t* opt)
{
    /// Initialisation necessary for each new tree

  fTree = tree;
  if (!fTree) return kFALSE;
  fTree->GetEntries();
  const char *treeName[AliAO2DTypes::kTrees] = { "O2events", "O2tracks", "O2calo",  "O2caloTrigger", "O2muon", "O2muoncls", "O2zdc", "O2vzero", "O2v0s", "O2cascades", "O2tof", "O2kine" };
  TTree *ttree = fTree->GetTree();
  if (!ttree) ttree = fTree; // is this ok?
  TString filename(ttree->GetCurrentFile()->GetName());
  auto dir = ttree->GetDirectory();
  fTrees[0] = ttree;
  fFileEntry = -1;
  for (int i = 1; i < AliAO2DTypes::kTrees; ++i)
    fTrees[i] = (TTree*)dir->Get(treeName[i]);

  ConnectEventTree(fTrees[0]);
  ConnectTracksTree(fTrees[1]);
    
  // Connect all trees with local objects
    

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAO2DInputHandler::BeginEvent(Long64_t entry)
{
  /// Begin event
  fEntry = entry;
  fFileEntry++;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAO2DInputHandler::Notify(const char* path)
{
  /// Notifiction of directory change    
  TTree *ttree = fTree->GetTree();
  if (!ttree) ttree = fTree; // is this ok?
  TString filename(ttree->GetCurrentFile()->GetName());
  AliInfo(Form("Moving to file %s", filename.Data()));
  
  return kTRUE;
}

//______________________________________________________________________________
void AliAO2DInputHandler::ConnectEventTree(TTree *tree) const
{
  TBranch *br;
  // Connect trees with local objects
  br = tree->GetBranch("fStart");
  br->SetAddress((void*)fEvent.fStart);
  br = tree->GetBranch("fNentries");
  br->SetAddress((void*)fEvent.fNentries);
  br = tree->GetBranch("fRunNumber");
  br->SetAddress((void*)&fEvent.fRunNumber);
  br = tree->GetBranch("fEventId");
  br->SetAddress((void*)&fEvent.fEventId);
  br = tree->GetBranch("fX");
  br->SetAddress((void*)&fEvent.fX);
  br = tree->GetBranch("fY");
  br->SetAddress((void*)&fEvent.fY);
  br = tree->GetBranch("fZ");
  br->SetAddress((void*)&fEvent.fZ);
  br = tree->GetBranch("fCovXX");
  br->SetAddress((void*)&fEvent.fCovXX);
  br = tree->GetBranch("fCovXY");
  br->SetAddress((void*)&fEvent.fCovXY);
  br = tree->GetBranch("fCovXZ");
  br->SetAddress((void*)&fEvent.fCovXZ);
  br = tree->GetBranch("fCovYY");
  br->SetAddress((void*)&fEvent.fCovYY);
  br = tree->GetBranch("fCovYZ");
  br->SetAddress((void*)&fEvent.fCovYZ);
  br = tree->GetBranch("fCovZZ");
  br->SetAddress((void*)&fEvent.fCovZZ);
  br = tree->GetBranch("fChi2");
  br->SetAddress((void*)&fEvent.fChi2);
  br = tree->GetBranch("fN");
  br->SetAddress((void*)&fEvent.fN);
  br = tree->GetBranch("fEventTime");
  br->SetAddress((void*)&fEvent.fEventTime);
  br = tree->GetBranch("fEventTimeRes");
  br->SetAddress((void*)&fEvent.fEventTimeRes);
  br = tree->GetBranch("fEventTimeMask");
  br->SetAddress((void*)&fEvent.fEventTimeMask);
}

//______________________________________________________________________________
void AliAO2DInputHandler::ConnectTracksTree(TTree *tree) const
{
  TBranch *br;
  // Connect trees with local objects
  br = tree->GetBranch("fCollisionsID");
  br->SetAddress((void*)&fReadTrack.fCollisionsID);
  br = tree->GetBranch("fX");
  br->SetAddress((void*)&fReadTrack.fX);
  br = tree->GetBranch("fAlpha");
  br->SetAddress((void*)&fReadTrack.fAlpha);
  br = tree->GetBranch("fY");
  br->SetAddress((void*)&fReadTrack.fY);
  br = tree->GetBranch("fZ");
  br->SetAddress((void*)&fReadTrack.fZ);
  br = tree->GetBranch("fSnp");
  br->SetAddress((void*)&fReadTrack.fSnp);
  br = tree->GetBranch("fTgl");
  br->SetAddress((void*)&fReadTrack.fTgl);
  br = tree->GetBranch("fSigned1Pt");
  br->SetAddress((void*)&fReadTrack.fSigned1Pt);
  br = tree->GetBranch("fCYY");
  br->SetAddress((void*)&fReadTrack.fCYY);
  br = tree->GetBranch("fCZY");
  br->SetAddress((void*)&fReadTrack.fCZY);
  br = tree->GetBranch("fCZZ");
  br->SetAddress((void*)&fReadTrack.fCZZ);
  br = tree->GetBranch("fCSnpY");
  br->SetAddress((void*)&fReadTrack.fCSnpY);
  br = tree->GetBranch("fCSnpZ");
  br->SetAddress((void*)&fReadTrack.fCSnpZ);
  br = tree->GetBranch("fCSnpSnp");
  br->SetAddress((void*)&fReadTrack.fCSnpSnp);
  br = tree->GetBranch("fCTglY");
  br->SetAddress((void*)&fReadTrack.fCTglY);
  br = tree->GetBranch("fCTglZ");
  br->SetAddress((void*)&fReadTrack.fCTglZ);
  br = tree->GetBranch("fCTglSnp");
  br->SetAddress((void*)&fReadTrack.fCTglSnp);
  br = tree->GetBranch("fCTglTgl");
  br->SetAddress((void*)&fReadTrack.fCTglTgl);
  br = tree->GetBranch("fC1PtY");
  br->SetAddress((void*)&fReadTrack.fC1PtY);
  br = tree->GetBranch("fC1PtZ");
  br->SetAddress((void*)&fReadTrack.fC1PtZ);
  br = tree->GetBranch("fC1PtSnp");
  br->SetAddress((void*)&fReadTrack.fC1PtSnp);
  br = tree->GetBranch("fC1PtTgl");
  br->SetAddress((void*)&fReadTrack.fC1PtTgl);
  br = tree->GetBranch("fC1Pt21Pt2");
  br->SetAddress((void*)&fReadTrack.fC1Pt21Pt2);
  br = tree->GetBranch("fTPCinnerP");
  br->SetAddress((void*)&fReadTrack.fTPCinnerP);
  br = tree->GetBranch("fFlags");
  br->SetAddress((void*)&fReadTrack.fFlags);
  br = tree->GetBranch("fITSClusterMap");
  br->SetAddress((void*)&fReadTrack.fITSClusterMap);
  br = tree->GetBranch("fTPCncls");
  br->SetAddress((void*)&fReadTrack.fTPCncls);
  br = tree->GetBranch("fTRDntracklets");
  br->SetAddress((void*)&fReadTrack.fTRDntracklets);
  br = tree->GetBranch("fITSchi2Ncl");
  br->SetAddress((void*)&fReadTrack.fITSchi2Ncl);
  br = tree->GetBranch("fTPCchi2Ncl");
  br->SetAddress((void*)&fReadTrack.fTPCchi2Ncl);
  br = tree->GetBranch("fTRDchi2");
  br->SetAddress((void*)&fReadTrack.fTRDchi2);
  br = tree->GetBranch("fTOFchi2");
  br->SetAddress((void*)&fReadTrack.fTOFchi2);
  br = tree->GetBranch("fTPCsignal");
  br->SetAddress((void*)&fReadTrack.fTPCsignal);
  br = tree->GetBranch("fTRDsignal");
  br->SetAddress((void*)&fReadTrack.fTRDsignal);
  br = tree->GetBranch("fTOFsignal");
  br->SetAddress((void*)&fReadTrack.fTOFsignal);
  br = tree->GetBranch("fLength");
  br->SetAddress((void*)&fReadTrack.fLength);
}

//______________________________________________________________________________
void AliAO2DInputHandler::ReadTracks()
{
  Int_t ntracks = fEvent.fNentries[AliAO2DTypes::kTracks];
  Int_t start = fEvent.fStart[AliAO2DTypes::kTracks];
  fTracks.resize(ntracks);

  for (Int_t entry = start; entry < start + ntracks; ++entry) {
    fTrees[AliAO2DTypes::kTracks]->GetEntry(entry);
    fTracks[entry - start] = fReadTrack;
  }
  fReadTree[AliAO2DTypes::kTracks] = true;
}

void AliAO2DInputHandler::ReadMCtracks() {}
void AliAO2DInputHandler::ReadTOFclusters() {}
void AliAO2DInputHandler::ReadCaloData() {}
void AliAO2DInputHandler::ReadCaloTrigger() {}
void AliAO2DInputHandler::ReadMUONTracks() {}
void AliAO2DInputHandler::ReadMUONclusters() {}
void AliAO2DInputHandler::ReadV0s() {}
void AliAO2DInputHandler::ReadCascades() {}

//______________________________________________________________________________
Bool_t AliAO2DInputHandler::FinishEvent()
{
  /// Finish event
  for (int i = 0; i < AliAO2DTypes::kTrees; ++i)
    fReadTree[i] = kFALSE;
  // if (fEvent) fEvent->Reset();
  return kTRUE;
}

//______________________________________________________________________________
Option_t *AliAO2DInputHandler::GetDataType() const
{
/// Returns handled data type.

   return gAO2DDataType;
}
