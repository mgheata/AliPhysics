R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <AliAO2DInputHandler.h>
#include <RUN3/AddTaskExFlow.C>


TChain* CreateChain(const char *xmlfile, const char *type="O2events");
TChain *CreateLocalChain(const char *txtfile, const char *type, int nfiles);

void testExFlowTask(int type = 1)
{
   if (type < 0 || type > 2) {
     printf("Error: supported types are: 0=ESD 1=AOD 2=AO2D\n");
     return;
   }
   const char *anatype[3] = {"esdEvent", "aodEvent", "O2events"};

 //  TGrid::Connect("alien:");

   AliAnalysisManager *mgr = new AliAnalysisManager("test_ExFlow");
   AliInputEventHandler *handler = nullptr;
   TChain *chain = nullptr;
   switch (type) {
     case 0:
       handler = new AliESDInputHandler();
       chain = CreateLocalChain("wn.txt", "esdTree", 10);
       break;
     case 1:
       handler = new AliAODInputHandler();
       chain = CreateLocalChain("wnaod.txt", "aodTree", 10);
       break;
     case 2:
       handler = new AliAO2DInputHandler();
       chain = CreateLocalChain("wnao2d.txt", "O2events", 10);
       break;
   };
   if (!chain) return;
   
   mgr->SetInputEventHandler(handler);
      
   auto taskFlow1 = AddTaskExFlow(768, 0.8, 10.0, 0.2, 50.0, 70, 2., 0.5, "v2_def_gap1");
   auto taskFlow2 = AddTaskExFlow(768, 0.8, 10.0, 0.2, 50.0, 70, 3., 0.5, "v3_def_gap1");
   auto taskFlow3 = AddTaskExFlow(768, 0.8, 10.0, 0.2, 50.0, 70, 4., 0.5, "v4_def_gap1");

   
   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();

   mgr->SetDebugLevel(1);
   mgr->StartAnalysis("localfile", chain, 123456789, 0);
}

TChain *CreateLocalChain(const char *txtfile, const char *type, int nfiles)
{
   TString treename = type;
   // treename.ToLower();
   // treename += "Tree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   // Open the file
   ifstream in;
   in.open(txtfile);
   Int_t count = 0;
    // Read the input list of files and add them to the chain
   TString line;
   TChain *chain = new TChain(treename);
   while (in.good())
   {
      in >> line;
      if (line.IsNull() || line.BeginsWith("#")) continue;
      if (count++ == nfiles) break;
      TString esdFile(line);
      TFile *file = TFile::Open(esdFile);
      if (file && !file->IsZombie()) {
         chain->Add(esdFile);
         file->Close();
      } else {
         Error("GetChainforTestMode", "Skipping un-openable file: %s", esdFile.Data());
      }   
   }
   in.close();
   if (!chain->GetListOfFiles()->GetEntries()) {
       Error("CreateLocalChain", "No file from %s could be opened", txtfile);
       delete chain;
       return nullptr;
   }
   return chain;
}

//________________________________________________________________________________
TChain* CreateChain(const char *xmlfile, const char *type)
{
// Create a chain using url's from xml file
   TString filename;
   Int_t run = 0;
   TString treename = type;
   //treename.ToLower();
   //treename += "Tree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   TGridCollection *coll = gGrid->OpenCollection(xmlfile);
   if (!coll) {
      ::Error("CreateChain", "Cannot create an AliEn collection from %s", xmlfile);
      return NULL;
   }
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TChain *chain = new TChain(treename);
   coll->Reset();
   while (coll->Next()) {
      filename = coll->GetTURL();
      if (mgr) {
         Int_t nrun = AliAnalysisManager::GetRunFromAlienPath(filename);
         if (nrun && nrun != run) {
            printf("### Run number detected from chain: %d\n", nrun);
            mgr->SetRunFromPath(nrun);
            run = nrun;
         }
      }
      chain->Add(filename);
   }
   if (!chain->GetNtrees()) {
      ::Error("CreateChain", "No tree found from collection %s", xmlfile);
      return NULL;
   }
   return chain;
}
