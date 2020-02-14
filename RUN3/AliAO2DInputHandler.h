#ifndef ALIAO2DINPUTHANDLER_H
#define ALIAO2DINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \class AliAO2DInputHandler
/// \brief AO2D Input Handler realisation of the AliInputEventHandler interface
///
/// \author Mihaela Gheata, CERN

#include <vector>
#include <array>

#include <AliInputEventHandler.h>

#include "AliAO2Dtypes.h"

class AliAO2DInputHandler : public AliInputEventHandler {

public:
    AliAO2DInputHandler();
    AliAO2DInputHandler(const char* name, const char* title);
    virtual ~AliAO2DInputHandler();
    virtual Bool_t       Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliVEventHandler::Notify();}
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    Option_t            *GetDataType() const;
    Long64_t             GetCurrentEntry() const { return fEntry; }
    Long64_t             GetFileEntry() const { return fFileEntry; }
    // Get the statistics object (currently TH2). Option can be BIN0.
    virtual TObject     *GetStatistics(Option_t *option="") const {return nullptr;}

    AliAO2DTypes::Event_t const &CurrentEvent() const {return fEvent;}
    AliAO2DTypes::MCvtx_t const &CurrentMCEvent()  const {return fMCEvent;}

    std::vector<AliAO2DTypes::Track_t> const &GetTracks()
    {
      if (!fReadTree[AliAO2DTypes::kTracks]) ReadTracks();
      return fTracks;
    }

    std::vector<AliAO2DTypes::MCparticle_t> const &GetMCtracks()
    {
      if (!fReadTree[AliAO2DTypes::kKinematics]) ReadMCtracks();
      return fMCtracks;
    }  

   private:
    AliAO2DInputHandler(const AliAO2DInputHandler& handler) = delete;
    AliAO2DInputHandler& operator=(const AliAO2DInputHandler& handler) = delete;

    void ConnectEventTree(TTree *tree) const;
    void ConnectTracksTree(TTree *tree) const;
    /// \brief Functions reading on demand different arrays
    void ReadTracks();
    void ReadMCtracks();
    void ReadTOFclusters();
    void ReadCaloData();
    void ReadCaloTrigger();
    void ReadMUONTracks();
    void ReadMUONclusters();
    void ReadV0s();
    void ReadCascades();

   private:
    std::array<TTree*, AliAO2DTypes::kTrees> fTrees;        //!<! Array of trees in the AO2D file
    bool fReadTree[AliAO2DTypes::kTrees] = {false};         //!<! Flags for read trees

    AliAO2DTypes::Event_t     fEvent;                       //!<! Pointer to the event
    AliAO2DTypes::MCvtx_t     fMCEvent;                     //!<! Pointer to the MCEvent

    AliAO2DTypes::Track_t     fReadTrack;                   //!<! Track connected to the tree
    std::vector<AliAO2DTypes::Track_t> fTracks;             //!<! Vector of tracks for the current event

    AliAO2DTypes::MCparticle_t fReadMCparticle;             //!<! MCparticle connected to the tree
    std::vector<AliAO2DTypes::MCparticle_t> fMCtracks;      //!<! Vactor of MC tracks for the current event
    
    AliAO2DTypes::TOFcluster_t fReadTOFcluster;             //!<! TOF cluster connected to the tree
    std::vector<AliAO2DTypes::TOFcluster_t> fTOFclusters;   //!<! Vector of TOF clusters

    AliAO2DTypes::CaloData_t fReadCaloData;                 //!<! Calo data connected to the tree
    std::vector<AliAO2DTypes::CaloData_t> fCaloData;        //!<! Vector of calo data
    
    AliAO2DTypes::CaloTrigger_t  fReadCaloTrigger;          //!<! CALO trigger connected to the tree
    std::vector<AliAO2DTypes::CaloTrigger_t> fCaloTrigger;  //!<! Vector of calo trigger data

    AliAO2DTypes::MuonTrack_t fReadMuonTrack;               //!<! MUON track connected to the tree
    std::vector<AliAO2DTypes::MuonTrack_t> fMuonTracks;     //!<! Vector of muon tracks for the current event
    
    AliAO2DTypes::MUONcluster_t fReadMUONcluster;           //!<! MUON cluster connected to the tree
    std::vector<AliAO2DTypes::MUONcluster_t> fMUONclusters; //!<! Vector of MUON clusters
    
    AliAO2DTypes::ZDCdata_t fZDCdata;                       //!<! ZDC data
    AliAO2DTypes::VZEROdata_t fVZEROdata;                   //!<! VZERO data
    
    AliAO2DTypes::V0_t fReadV0;                             //!<! V0 connected to the tree
    std::vector<AliAO2DTypes::V0_t> fV0s;                   //!<! V0s
    
    AliAO2DTypes::Cascade_t fCascade;                       //!<! Cascade connected to the tree
    std::vector<AliAO2DTypes::Cascade_t> fCascades;         //!<! Cascades

    Long64_t fEntry = -1;                                   //!<! Current read entry
    Long64_t fFileEntry = -1;                               //!<! Current read entry in the current file
    bool fUseMC = false;                                    //!<! Use MC

    ClassDef(AliAO2DInputHandler, 1);
};

#endif
