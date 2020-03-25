#include "AliExFlowTask.h"
AliExFlowTask* AddTaskExFlow(UInt_t filterbit = 768,
                             Double_t etaCut = 0.8,
                             Double_t vtxCut = 10.,
                             Double_t minPt = 0.2,
                             Double_t maxPt = 10.0,
                             Int_t noclus = 70,
                             Double_t nHarm = 2.,
                             Float_t etaGap = 0.1,
                             TString uniqueID = "")
{
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddExFlowTask", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskData", "This task requires an input event handler");
    return NULL;
  }  

    

  // Create the task and configure it 
  //========================================================================
  AliExFlowTask* taskD = new AliExFlowTask(Form("taskFlow_%s", uniqueID.Data()));
    
    taskD->SetDebugLevel(3);
    taskD->SetFilterbit(filterbit);
    taskD->SetNoClus(noclus);
    taskD->SetEtaCut(etaCut);
    taskD->SetVtxCut(vtxCut);
    taskD->SetMinPt(minPt);
    taskD->SetMaxPt(maxPt);
    taskD->SetNHarmonic(nHarm);
    taskD->SetEtaGap(etaGap);

    mgr->AddTask(taskD);
  
  
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer* cout = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("vn_SP_TPC.root:%s", uniqueID.Data()));
    mgr->ConnectInput (taskD, 0, cinput);
    mgr->ConnectOutput(taskD, 1, cout);
  
    // Return task pointer at the end
    return taskD;
    
}
