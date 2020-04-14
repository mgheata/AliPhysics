#ifndef ALIEXFLOWTASK_H
#define ALIEXFLOWTASK_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TString.h>
#include <THnSparse.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliVHeader.h>
#include <AliVVertex.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliPIDResponse.h>
#include <AliAO2DInputHandler.h>


class AliExFlowTask : public AliAnalysisTaskSE {
 public:
  AliExFlowTask();
  AliExFlowTask(const char *name);

  virtual ~AliExFlowTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  void   UserExecESD(Option_t *option);
  void   UserExecAOD(Option_t *option);
  void   UserExecAO2D(Option_t *option);
  virtual void   Terminate(Option_t *); 

  
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetMinPt() { return fMinPt; }   


    virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
    virtual void  SetFilterbit(UInt_t filterbit){fFilterbit = filterbit;}
    virtual void  SetNoClus(Int_t noclus){fNoClus = noclus;}
    virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
    virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
    virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
    virtual void  SetNHarmonic(Double_t nHarm){fNHarm = nHarm;}
    virtual void  SetEtaGap(Float_t etaGap){fEtaGap = etaGap;}

 private:
    Float_t GetVertex(AliESDEvent* esd) const;
    Float_t GetVertex(AliAODEvent* aod) const;
    Float_t GetVertex(AliAO2DTypes::Vertex_t *vtx) const;

    void Analyze(AliESDEvent* esd, Float_t vtxZ);
    void Analyze(AliAODEvent* aod, Float_t vtxZ);
    void Analyze(AliAO2DTypes::Vertex_t *vtx, Float_t vtxZ);


    Int_t        fIev = 0;       // event number
    AliAO2DInputHandler *fAO2Dhandler = nullptr; // AO2D input handler
    AliAODEvent* fAOD = nullptr;                //! AOD object
    AliESDEvent* fESD = nullptr;
    
    // Cuts and options
    Double_t     fVtxCut;             // Vtx cut on z position in cm
    UInt_t       fFilterbit;          // filter bit
    Double_t     fEtaCut;             // Eta cut used to select particles
    Int_t        fNoClus;	          // No of TPC clusters
    Double_t     fMinPt;              // Min pt - for histogram limits
    Double_t     fMaxPt;              // Max pt - for histogram limits
    Double_t     fNHarm;              // harmonic number
    Float_t      fEtaGap;             // eta gap for gap 1 -> fEtaGap = 0.5


    //
    // Output objects
    //
    TList*        fListOfObjects;     //! Output list of objects
    TH1I*         fVtx;               //! Event vertex info
    TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
    TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts

    TProfile*     fResvn;         //! correlation V0A-V0C for resolution

    
    TProfile*     fVnAllA[10];             //! vn V0A all
    TProfile*     fVnAllC[10];             //! vn V0C all
    
  ClassDef(AliExFlowTask, 1);    //Example flow analysis
};

#endif
