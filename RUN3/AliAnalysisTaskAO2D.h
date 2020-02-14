#ifndef AliAnalysisTaskAO2D_cxx
#define AliAnalysisTaskAO2D_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliVEvent;
class TList;

#include "AliAnalysisTaskSE.h"
#include "AliAO2Dtypes.h"

class AliAnalysisTaskAO2D : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskAO2D() {}
  AliAnalysisTaskAO2D(const char *name);
  virtual ~AliAnalysisTaskAO2D() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  void UserExecAO2D();
  virtual void   Terminate(Option_t *);
  
private:
  TH1F*             fHistPt = nullptr;         // Pt spectrum
  Int_t             fEv = -1;
  TH1F*             fHistQ = nullptr;          // TPC clusters Q spectrum
  TList*            fListOut = nullptr;        // output list
  TH1F*             fHistNTPCCl = nullptr;     // histo with the number of TPC clusters
  TH1F*             fHistNESDtracks = nullptr; // histo with number of ESD tracks
 
  AliAnalysisTaskAO2D(const AliAnalysisTaskAO2D&) = delete; // not implemented
  AliAnalysisTaskAO2D& operator=(const AliAnalysisTaskAO2D&) = delete; // not implemented
  
  ClassDef(AliAnalysisTaskAO2D, 1); // example of analysis
};

#endif
