#include "AliExFlowTask.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TRandom.h>
#include <TVector3.h>
#include <THnSparse.h>
#include <TProfile2D.h>


// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "AliAODVZERO.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliAODv0.h"
//#include "AliMultSelection.h"
#include "AliAODcascade.h"

// STL includes
#include <iostream>
#include <ctime>
#include <sys/time.h>
using std::cout;
using std::endl;

ClassImp(AliExFlowTask)
//_____________________________________________________________________________
AliExFlowTask::AliExFlowTask():
  AliAnalysisTaskSE(),
  fAOD(0),
  fVtxCut(10.0),
  fFilterbit(128),
  fEtaCut(0.8),
  fNoClus(70),
  fMinPt(0.2),
  fMaxPt(20.0),
  fNHarm(2.),
  fEtaGap(0),
  fListOfObjects(0),
  fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0),
  fResvn(0)
{
    
    for (Int_t i = 0; i < 10; i++){

        fVnAllA[i] = 0;
        fVnAllC[i] = 0;
      
    }

}

//______________________________________________________________________________
AliExFlowTask::AliExFlowTask(const char *name):
    AliAnalysisTaskSE(name),
fAOD(0),
fVtxCut(10.0),
fFilterbit(128),
fEtaCut(0.8),
fNoClus(70),
fMinPt(0.2),
fMaxPt(20.0),
fNHarm(2.),
fEtaGap(0),
fListOfObjects(0),
fVtx(0), fVtxBeforeCuts(0), fVtxAfterCuts(0),
fResvn(0)
{
    
    for (Int_t i = 0; i < 10; i++){
        
        fVnAllA[i] = 0;
        fVnAllC[i] = 0;
        
    }

  // Output slot #1 writes into a TTree
  DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliExFlowTask::~AliExFlowTask()
{
  // Destructor
  if (fListOfObjects) 
    delete fListOfObjects;
  
}

//______________________________________________________________________________
void AliExFlowTask::UserCreateOutputObjects()
{ 
 
    OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner();


    const Int_t nCenB = 10;
    Float_t cenBins[nCenB+1] = {0, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    
    const Int_t nPtB = 36;
    Double_t ptBins[nPtB+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0, 30.0, 40.0, 50.0};

    
    fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
    fListOfObjects->Add(fVtx);
        
    fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
    fListOfObjects->Add(fVtxBeforeCuts);
        
    fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
    fListOfObjects->Add(fVtxAfterCuts);

    
    fResvn = new TProfile("fResvn", "; centrality percentile; resolution", nCenB, cenBins);
    fListOfObjects->Add(fResvn);

    for (Int_t ic = 0; ic < nCenB; ic++){
            
        fVnAllA[ic] = new TProfile(Form("fVnAllA_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
        fListOfObjects->Add(fVnAllA[ic]);
            
        fVnAllC[ic] = new TProfile(Form("fVnAllC_%d", ic), "; p_{T} (GeV/c); v_{n}", nPtB, ptBins);
        fListOfObjects->Add(fVnAllC[ic]);
        
    }
    
  fAO2Dhandler = dynamic_cast<AliAO2DInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  assert(fAO2Dhandler != nullptr);

  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliExFlowTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  auto vtx = fAO2Dhandler->CurrentEvent();
  auto tracks = fAO2Dhandler->GetTracks();

  Float_t zvtx = GetVertex(&vtx);
   
  if(zvtx< -990){

    fVtx->Fill(0);
      
  } else {

      fVtx->Fill(1);
      fVtxBeforeCuts->Fill(zvtx);
      
      if (TMath::Abs(zvtx) < fVtxCut)
          Analyze(&vtx, zvtx);
          
  }
    
  // Post output data.
  PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliExFlowTask::Analyze(AliAO2DTypes::Vertex_t *vtx, Float_t vtxZ)
{  

    //Centrality
    Float_t v0Centr = -100.;

    /* // Multiplicity selection object not available, so use the V0M centrality percentile
    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) aod->FindListObject("MultSelection");
    if( !MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return;
    } else {
        v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
    }
    */
    auto centralityObj = fAO2Dhandler->GetCentrality();
    if (!centralityObj) return;
    v0Centr = centralityObj->GetCentralityPercentile("V0M");


    Short_t centrCode = -10;
    if ((v0Centr >= 0) && (v0Centr < 5.))
        centrCode = 0;
    else if ((v0Centr >= 5.) && (v0Centr < 10.))
        centrCode = 1;
    else if ((v0Centr >= 10.) && (v0Centr < 20.))
        centrCode = 2;
    else if ((v0Centr >= 20.) && (v0Centr < 30.))
        centrCode = 3;
    else if ((v0Centr >= 30.) && (v0Centr < 40.))
        centrCode = 4;
    else if ((v0Centr >= 40.) && (v0Centr < 50.))
        centrCode = 5;
    else if ((v0Centr >= 50.) && (v0Centr < 60.))
        centrCode = 6;
    else if ((v0Centr >= 60.) && (v0Centr < 70.))
        centrCode = 7;
    else if ((v0Centr >= 70.) && (v0Centr < 80.))
        centrCode = 8;
    else if ((v0Centr >= 80.) && (v0Centr < 90.))
        centrCode = 9;
    
    if (centrCode < 0)
        return;
    
    
    fVtxAfterCuts->Fill(vtxZ);

    auto tracks = fAO2Dhandler->GetTracks();    
    
    const Int_t nTracks = tracks.size();
    

    Double_t Qxan = 0, Qyan = 0;
    Double_t Qxcn = 0, Qycn = 0;
    Int_t sumMta = 0, sumMtc = 0;
    
    for (Int_t loop1 = 0; loop1 < 2; loop1++){

        for (Int_t it1 = 0; it1 < nTracks; it1++) {
            
            auto const &trk1 = tracks[it1];

            // no filter bit
            /*
            if (!(aodTrk1->TestFilterBit(fFilterbit)))
                continue;
            */
            if ((TMath::Abs(trk1.Eta()) >= fEtaCut) || (trk1.fTPCncls < fNoClus) || (trk1.Pt() < fMinPt) || (trk1.Pt() >= fMaxPt))
                continue;
            
            if (loop1 == 0){
                
                if (trk1.Eta() > fEtaGap){
                    
                    Qxan += TMath::Cos(fNHarm*trk1.Phi());
                    Qyan += TMath::Sin(fNHarm*trk1.Phi());
                    
                    sumMta++;
                    
                }
                
                if (trk1.Eta() < -fEtaGap){
                    
                    Qxcn += TMath::Cos(fNHarm*trk1.Phi());
                    Qycn += TMath::Sin(fNHarm*trk1.Phi());
                    
                    sumMtc++;
                }
               
                
            } else {
            
  
                
                if (trk1.Eta() < -fEtaGap && sumMta > 0){
                //if (aodTrk1->Eta() < -fEtaGap){
        
                    Double_t vnSPA = (TMath::Cos(fNHarm*trk1.Phi())*Qxan + TMath::Sin(fNHarm*trk1.Phi())*Qyan)/(Double_t)sumMta;
                    //Double_t vnSPA = (TMath::Cos(fNHarm*aodTrk1->Phi())*Qxan + TMath::Sin(fNHarm*aodTrk1->Phi())*Qyan);

                        
                    fVnAllA[centrCode]->Fill(trk1.Pt(), vnSPA);
                    
                }
                
                
                if (trk1.Eta() > fEtaGap && sumMtc > 0){
                //if (aodTrk1->Eta() > fEtaGap){
                    
                    Double_t vnSPC = (TMath::Cos(fNHarm*trk1.Phi())*Qxcn + TMath::Sin(fNHarm*trk1.Phi())*Qycn)/(Double_t)sumMtc;
                    //Double_t vnSPC = (TMath::Cos(fNHarm*aodTrk1->Phi())*Qxcn + TMath::Sin(fNHarm*aodTrk1->Phi())*Qycn);

                    
                    fVnAllC[centrCode]->Fill(trk1.Pt(), vnSPC);
                    
                }
                
                
            }
            
        }
        
    }
    

    if (sumMta > 0 && sumMtc > 0){
        
        Double_t resvn = (Qxan*Qxcn + Qyan*Qycn)/(Double_t)sumMta/(Double_t)sumMtc;
        fResvn->Fill(v0Centr, resvn);
        
    }
    
    //Double_t resvn = (Qxan*Qxcn + Qyan*Qycn);
    //fResvn->Fill(v0Centr, resvn);
    
}

//_____________________________________________________________________________
Float_t AliExFlowTask::GetVertex(AliAO2DTypes::Vertex_t *vtx) const
{

  Float_t vtxz = -999;
  if (vtx->fN < 2)
    return vtxz;

  // The SPD vertex information is missing in AO2D, so we can use only a simplified version
  return vtx->fZ;
    
    
  /*
    //new vertex selection
    const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
    const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();
    
    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)
        return vtxz; // one of vertices is missing
    
    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    
    double dz = vtTrc->GetZ() - vtSPD->GetZ();
    
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;
    
    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
        return vtxz; // bad vertexing
    
    
    TString vtxTyp = vtSPD->GetTitle();
    Double_t zRes = TMath::Sqrt(covSPD[5]);
    if ((vtxTyp.Contains("vertexer:Z")) && (zRes>0.25) && (vtSPD->GetNContributors() < 20))
        return vtxz; // bad vertexing
    
    
  vtxz = vtTrc->GetZ();

  return vtxz;
  */
}

//_____________________________________________________________________________
void AliExFlowTask::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
