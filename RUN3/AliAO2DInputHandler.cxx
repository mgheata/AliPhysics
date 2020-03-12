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
  TTree *ttree = fTree->GetTree();
  if (!ttree) ttree = fTree; // is this ok?
  TString filename(ttree->GetCurrentFile()->GetName());
  auto dir = ttree->GetDirectory();
  fTrees[0] = ttree;
  fFileEntry = -1;
  for (int i = 1; i < AliAO2DTypes::kTrees; ++i)
    fTrees[i] = (TTree*)dir->Get(AliAO2DTypes::kTreeName[i]);

  ConnectEventTree(fTrees[AliAO2DTypes::kEvents]);
  assert(fTrees[AliAO2DTypes::kTracks]);
  ConnectTracksTree(fTrees[AliAO2DTypes::kTracks]);
  // Check if centrality tree exists
  if (fTrees[AliAO2DTypes::kCentrality]) {
    fHasCentrality = true;
    ConnectCentralityTree(fTrees[AliAO2DTypes::kCentrality]);
  }
  // Check if multiplicity tree exists
  if (fTrees[AliAO2DTypes::kMultiplicity]) {
    fHasMultiplicity = true;
    ConnectMultiplicityTree(fTrees[AliAO2DTypes::kMultiplicity]);
  }
    
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
void AliAO2DInputHandler::ConnectCentralityTree(TTree *tree) const
{
  TBranch *br;
  // Connect trees with local objects
  br = tree->GetBranch("fQuality");
  br->SetAddress((void*)&fCentrality.fQuality);
  br = tree->GetBranch("fCentralityV0M");
  br->SetAddress((void*)&fCentrality.fCentralityV0M);
  br = tree->GetBranch("fCentralityV0A");
  br->SetAddress((void*)&fCentrality.fCentralityV0A);
  br = tree->GetBranch("fCentralityV0A0");
  br->SetAddress((void*)&fCentrality.fCentralityV0A0);
  br = tree->GetBranch("fCentralityV0A123");
  br->SetAddress((void*)&fCentrality.fCentralityV0A123);
  br = tree->GetBranch("fCentralityV0C");
  br->SetAddress((void*)&fCentrality.fCentralityV0C);
  br = tree->GetBranch("fCentralityV0A23");
  br->SetAddress((void*)&fCentrality.fCentralityV0A23);
  br = tree->GetBranch("fCentralityV0C01");
  br->SetAddress((void*)&fCentrality.fCentralityV0C01);
  br = tree->GetBranch("fCentralityV0S");
  br->SetAddress((void*)&fCentrality.fCentralityV0S);
  br = tree->GetBranch("fCentralityV0MEq");
  br->SetAddress((void*)&fCentrality.fCentralityV0MEq);
  br = tree->GetBranch("fCentralityV0AEq");
  br->SetAddress((void*)&fCentrality.fCentralityV0AEq);
  br = tree->GetBranch("fCentralityV0CEq");
  br->SetAddress((void*)&fCentrality.fCentralityV0CEq);
  br = tree->GetBranch("fCentralityFMD");
  br->SetAddress((void*)&fCentrality.fCentralityFMD);
  br = tree->GetBranch("fCentralityTRK");
  br->SetAddress((void*)&fCentrality.fCentralityTRK);
  br = tree->GetBranch("fCentralityTKL");
  br->SetAddress((void*)&fCentrality.fCentralityTKL);
  br = tree->GetBranch("fCentralityCL0");
  br->SetAddress((void*)&fCentrality.fCentralityCL0);
  br = tree->GetBranch("fCentralityCL1");
  br->SetAddress((void*)&fCentrality.fCentralityCL1);
  br = tree->GetBranch("fCentralityCND");
  br->SetAddress((void*)&fCentrality.fCentralityCND);
  br = tree->GetBranch("fCentralityZNA");
  br->SetAddress((void*)&fCentrality.fCentralityZNA);
  br = tree->GetBranch("fCentralityZNC");
  br->SetAddress((void*)&fCentrality.fCentralityZNC);
  br = tree->GetBranch("fCentralityZPA");
  br->SetAddress((void*)&fCentrality.fCentralityZPA);
  br = tree->GetBranch("fCentralityZPC");
  br->SetAddress((void*)&fCentrality.fCentralityZPC);
  br = tree->GetBranch("fCentralityNPA");
  br->SetAddress((void*)&fCentrality.fCentralityNPA);
  br = tree->GetBranch("fCentralityV0MvsFMD");
  br->SetAddress((void*)&fCentrality.fCentralityV0MvsFMD);
  br = tree->GetBranch("fCentralityTKLvsV0M");
  br->SetAddress((void*)&fCentrality.fCentralityTKLvsV0M);
  br = tree->GetBranch("fCentralityZEMvsZDC");
  br->SetAddress((void*)&fCentrality.fCentralityZEMvsZDC);
  br = tree->GetBranch("fCentralityV0Mtrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0Mtrue);
  br = tree->GetBranch("fCentralityV0Atrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0Atrue);
  br = tree->GetBranch("fCentralityV0Ctrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0Ctrue);
  br = tree->GetBranch("fCentralityV0MEqtrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0MEqtrue);
  br = tree->GetBranch("fCentralityV0AEqtrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0AEqtrue);
  br = tree->GetBranch("fCentralityV0CEqtrue");
  br->SetAddress((void*)&fCentrality.fCentralityV0CEqtrue);
  br = tree->GetBranch("fCentralityFMDtrue");
  br->SetAddress((void*)&fCentrality.fCentralityFMDtrue);
  br = tree->GetBranch("fCentralityTRKtrue");
  br->SetAddress((void*)&fCentrality.fCentralityTRKtrue);
  br = tree->GetBranch("fCentralityTKLtrue");
  br->SetAddress((void*)&fCentrality.fCentralityTKLtrue);
  br = tree->GetBranch("fCentralityCL0true");
  br->SetAddress((void*)&fCentrality.fCentralityCL0true);
  br = tree->GetBranch("fCentralityCL1true");
  br->SetAddress((void*)&fCentrality.fCentralityCL1true);
  br = tree->GetBranch("fCentralityCNDtrue");
  br->SetAddress((void*)&fCentrality.fCentralityCNDtrue);
  br = tree->GetBranch("fCentralityZNAtrue");
  br->SetAddress((void*)&fCentrality.fCentralityZNAtrue);
  br = tree->GetBranch("fCentralityZNCtrue");
  br->SetAddress((void*)&fCentrality.fCentralityZNCtrue);
  br = tree->GetBranch("fCentralityZPAtrue");
  br->SetAddress((void*)&fCentrality.fCentralityZPAtrue);
  br = tree->GetBranch("fCentralityZPCtrue");
  br->SetAddress((void*)&fCentrality.fCentralityZPCtrue);
}

//______________________________________________________________________________
void AliAO2DInputHandler::ConnectMultiplicityTree(TTree *tree) const
{
  TBranch *br;
  // Connect trees with local objects
  br = tree->GetBranch("fNtracks");
  br->SetAddress((void*)&fMultiplicity.fNtracks);
  br = tree->GetBranch("fNsingle");
  br->SetAddress((void*)&fMultiplicity.fNsingle);
  br = tree->GetBranch("fDPhiWindow2");
  br->SetAddress((void*)&fMultiplicity.fDPhiWindow2);
  br = tree->GetBranch("fDThetaWindow2");
  br->SetAddress((void*)&fMultiplicity.fDThetaWindow2);
  br = tree->GetBranch("fDPhiShift");
  br->SetAddress((void*)&fMultiplicity.fDPhiShift);
  br = tree->GetBranch("fNStdDev");
  br->SetAddress((void*)&fMultiplicity.fNStdDev);
  br = tree->GetBranch("fFiredChips");
  br->SetAddress((void*)fMultiplicity.fFiredChips);
  br = tree->GetBranch("fITSClusters");
  br->SetAddress((void*)fMultiplicity.fITSClusters);
  br = tree->GetBranch("fCentroidXY");
  br->SetAddress((void*)fMultiplicity.fCentroidXY);
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

//______________________________________________________________________________
void AliAO2DInputHandler::ReadCentrality()
{
  Int_t start = fEvent.fStart[AliAO2DTypes::kCentrality];
  fTrees[AliAO2DTypes::kCentrality]->GetEntry(start);
}

//______________________________________________________________________________
void AliAO2DInputHandler::ReadMultiplicity()
{
  Int_t start = fEvent.fStart[AliAO2DTypes::kMultiplicity];
  fTrees[AliAO2DTypes::kMultiplicity]->GetEntry(start);
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
