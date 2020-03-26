/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAO2Dtypes_H
#define AliAO2Dtypes_H

namespace AliAO2DTypes {

enum TreeIndex { // Index of the output trees
  kEvents = 0,  // O2events
  kTracks,      // O2tracks
  kCalo,        // O2calo
  kCaloTrigger, // O2caloTrigger
  kMuon,        // O2muon
  kMuonCls,     // O2muoncls
  kZdc,         // O2zdc
  kVzero,       // O2vzero
  kV0s,         // O2v0s
  kCascades,    // O2cascades
  kTOF,         // O2tof
  kKinematics,  // O2kine
  kTrees
};

struct Vertex_t {
    // Start indices and numbers of elements for data in the other trees matching this vertex.
    // Needed for random access of collision-related data, allowing skipping data discarded by the user
    Int_t     fStart[kTrees]    = {0}; /// Start entry indices for data in the other trees matching this vertex
    Int_t     fNentries[kTrees] = {0}; /// Numbers of entries for data in the other trees matching this vertex
    // Event data
    Int_t     fRunNumber;       /// Run Number (added in case of multirun skimming)
    ULong64_t fEventId = 0u;    /// Event (collision) unique id. Contains period, orbit and bunch crossing numbers
    // Primary vertex position
    Float_t  fX = -999.f;       /// Primary vertex x coordinate
    Float_t  fY = -999.f;       /// Primary vertex y coordinate
    Float_t  fZ = -999.f;       /// Primary vertex z coordinate
    // Primary vertex covariance matrix
    Float_t  fCovXX = 999.f;    /// cov[0]
    Float_t  fCovXY = 0.f;      /// cov[1]
    Float_t  fCovXZ = 0.f;      /// cov[2]
    Float_t  fCovYY = 999.f;    /// cov[3]
    Float_t  fCovYZ = 0.f;      /// cov[4]
    Float_t  fCovZZ = 999.f;    /// cov[5]
    // Quality parameters
    Float_t  fChi2;             /// Chi2 of the vertex
    UInt_t   fN;                /// Number of contributors

    // The calculation of event time certainly will be modified in Run3
    // The prototype below can be switched on request
    Float_t fEventTime = -999.f;    /// Event time (t0) obtained with different methods (best, T0, T0-TOF, ...)
    Float_t fEventTimeRes = -999.f; /// Resolution on the event time (t0) obtained with different methods (best, T0, T0-TOF, ...)
    UChar_t fEventTimeMask = 0u;    /// Mask with the method used to compute the event time (0x1=T0-TOF,0x2=T0A,0x3=TOC) for each momentum bins  
};

using Event_t = Vertex_t;

struct Track_t {
  // Track data

  Int_t   fCollisionsID;    /// The index of the collision vertex in the TF, to which the track is attached
  // In case we need connection to TOF clusters, activate next lines
  // Int_t   fTOFclsIndex;     /// The index of the associated TOF cluster
  // Int_t   fNTOFcls;         /// The number of TOF clusters

  // Coordinate system parameters
  Float_t fX = -999.f;          /// X coordinate for the point of parametrisation
  Float_t fAlpha = -999.f;      /// Local <--> global coor.system rotation angle

  // Track parameters
  Float_t fY = -999.f;          /// fP[0] local Y-coordinate of a track (cm)
  Float_t fZ = -999.f;          /// fP[1] local Z-coordinate of a track (cm)
  Float_t fSnp = -999.f;        /// fP[2] local sine of the track momentum azimuthal angle
  Float_t fTgl = -999.f;        /// fP[3] tangent of the track momentum dip angle
  Float_t fSigned1Pt = -999.f;  /// fP[4] 1/pt (1/(GeV/c))

  // Covariance matrix
  Float_t fCYY = -999.f;        /// fC[0]
  Float_t fCZY = -999.f;        /// fC[1]
  Float_t fCZZ = -999.f;        /// fC[2]
  Float_t fCSnpY = -999.f;      /// fC[3]
  Float_t fCSnpZ = -999.f;      /// fC[4]
  Float_t fCSnpSnp = -999.f;    /// fC[5]
  Float_t fCTglY = -999.f;      /// fC[6]
  Float_t fCTglZ = -999.f;      /// fC[7]
  Float_t fCTglSnp = -999.f;    /// fC[8]
  Float_t fCTglTgl = -999.f;    /// fC[9]
  Float_t fC1PtY = -999.f;      /// fC[10]
  Float_t fC1PtZ = -999.f;      /// fC[11]
  Float_t fC1PtSnp = -999.f;    /// fC[12]
  Float_t fC1PtTgl = -999.f;    /// fC[13]
  Float_t fC1Pt21Pt2 = -999.f;  /// fC[14]

  // Additional track parameters
  Float_t fTPCinnerP = -999.f;  /// Full momentum at the inner wall of TPC for dE/dx PID

  // Track quality parameters
  ULong64_t fFlags = 0u;        /// Reconstruction status flags

  // Clusters
  UChar_t fITSClusterMap = 0u;  /// ITS map of clusters, one bit per a layer
  UShort_t fTPCncls = 0u;       /// number of clusters assigned in the TPC
  UChar_t fTRDntracklets = 0u;  /// number of TRD tracklets used for tracking/PID (TRD/TOF pattern)

  // Chi2
  Float_t fITSchi2Ncl = -999.f; /// chi2/Ncl ITS
  Float_t fTPCchi2Ncl = -999.f; /// chi2/Ncl TPC
  Float_t fTRDchi2 = -999.f;    /// chi2 TRD match (?)
  Float_t fTOFchi2 = -999.f;    /// chi2 TOF match (?)

  // PID
  Float_t fTPCsignal = -999.f;  /// dE/dX TPC
  Float_t fTRDsignal = -999.f;  /// dE/dX TRD
  Float_t fTOFsignal = -999.f;  /// TOFsignal
  Float_t fLength = -999.f;     /// Int.Lenght @ TOF
};                              //  structure to keep track information

struct MCvtx_t {
  Int_t     fRunNumber;         /// Run Number (added in case of multirun skimming)
  // MC information on the event
  Short_t fGeneratorsID = 0u;   /// Generator ID used for the MC
  Float_t fX = -999.f;          /// Primary vertex x coordinate from MC
  Float_t fY = -999.f;          /// Primary vertex y coordinate from MC
  Float_t fZ = -999.f;          /// Primary vertex z coordinate from MC
  Float_t fT = -999.f;          /// Time of the collision from MC
};                              // MC vertices

// Track labels
struct MCparticle_t {
  // Int_t fLabel = -1;         /// Track label
  // Int_t fTOFLabel[3] = {-1}; /// Label of the track matched to TOF

  Int_t   fCollisionsID;        /// The index of the MC collision vertex

  // MC information (modified version of TParticle
  Int_t fPdgCode    = -99999;   /// PDG code of the particle
  Int_t fStatusCode = -99999;   /// generation status code
  Int_t fMother[2]   = { 0 };   /// Indices of the mother particles
  Int_t fDaughter[2] = { 0 };   /// Indices of the daughter particles
  Float_t fWeight    = 1;       /// particle weight from the generator or ML

  Float_t fPx = -999.f;         /// x component of momentum
  Float_t fPy = -999.f;         /// y component of momentum
  Float_t fPz = -999.f;         /// z component of momentum
  Float_t fE  = -999.f;         /// Energy (covers the case of resonances, no need for calculated mass)

  Float_t fVx = -999.f;         /// x of production vertex
  Float_t fVy = -999.f;         /// y of production vertex
  Float_t fVz = -999.f;         /// z of production vertex
  Float_t fVt = -999.f;         /// t of production vertex
  // We do not use the polarisation so far
};                              // MC particles from the kinematics tree

struct TOFcluster_t {
  // TOF clusters
  // PH: Do we store the TOF information per track?
  Int_t fTOFChannel = -1;       /// Index of the matched channel
  Short_t fTOFncls = -1;        /// Number of matchable clusters of the track
  Float_t fDx = -999.f;         /// Residual along x
  Float_t fDz = -999.f;         /// Residual along z
  Float_t fToT = -999.f;        /// ToT
  Float_t fLengthRatio = -999.f; /// Ratio of the integrated track length @ TOF to the cluster with respect to the matched cluster
};                              // structure to keep TOF clusters

struct CaloData_t {
  // Calorimeter data (EMCAL & PHOS)

  Int_t   fCollisionsID;        /// The index of the collision vertex in the TF, to which the track is attached

  Short_t fCellNumber = -1;     /// Cell absolute Id. number
  Float_t fAmplitude = -999.f;  /// Cell amplitude (= energy!)
  Float_t fTime = -999.f;       /// Cell time
  Char_t fCellType = -1;        /// EMCAL: High Gain: 0 / Low Gain: 1 / TRU: 2 / LEDmon 3 (see DataFromatsEMCAL/Constants.h)
  Char_t fType = -1;            /// Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
};                              //  structure to keep EMCAL info

struct CaloTrigger_t {
  // Calorimeter trigger data (EMCAL & PHOS)
  Int_t   fCollisionsID;        /// The index of the collision vertex in the TF, to which the track is attached
  Short_t fFastorAbsID = - 1;   /// FastOR absolute ID
  Float_t fL0Amplitude = -1.f;  /// L0 amplitude (ADC) := Peak Amplitude
  Float_t fL0Time = -1.f;       /// L0 time
  Int_t fL1TimeSum = -1;        /// L1 amplitude (ADC) := Integral over L0 time samples
  Char_t fNL0Times = -1;        /// Number of L0 times
  Int_t fTriggerBits = 0;       /// Online trigger bits
  Char_t fType = -1;            /// Calorimeter type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
};                              // structure to keep calo trigger info

struct MuonTrack_t {
  // MUON track data

  Int_t   fCollisionsID;        /// The index of the collision vertex, to which the muon is attached
  // In case we need connection to muon clusters, activate next lines
  // Int_t   fClusterIndex;     /// The index of the associated MUON clusters
  // Int_t   fNclusters;        /// The number of MUON clusters

  /// Parameters at vertex
  Float_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge 
  Float_t fThetaX;              ///< Angle of track at vertex in X direction (rad)
  Float_t fThetaY;              ///< Angle of track at vertex in Y direction (rad)
  Float_t fZ;                   ///< Z coordinate (cm)
  Float_t fBendingCoor;         ///< bending coordinate (cm)
  Float_t fNonBendingCoor;      ///< non bending coordinate (cm)

  /// Reduced covariance matrix of UNCORRECTED track parameters, ordered as follow:      <pre>
  /// [0] =  <X,X>
  /// [1] =<X,ThetaX>  [2] =<ThetaX,ThetaX>
  /// [3] =  <X,Y>     [4] =  <Y,ThetaX>     [5] =  <Y,Y>
  /// [6] =<X,ThetaY>  [7] =<ThetaX,ThetaY>  [8] =<Y,ThetaY>  [9] =<ThetaY,ThetaY>
  /// [10]=<X,InvP_yz> [11]=<ThetaX,InvP_yz> [12]=<Y,InvP_yz> [13]=<ThetaY,InvP_yz> [14]=<InvP_yz,InvP_yz>  </pre>
  Float_t fCovariances[15];     ///< \brief reduced covariance matrix of parameters AT FIRST CHAMBER

  /// Global tracking info
  Float_t fChi2;                ///< chi2 in the MUON track fit
  Float_t fChi2MatchTrigger;    ///< chi2 of trigger/track matching
};                              // structure to keep muons information

struct MUONcluster_t {
  // Muon cluster data
    
  Int_t   fMuonsID;             /// The index of the muon track to which the clusters are attached
  Float_t fX;                   ///< cluster X position
  Float_t fY;                   ///< cluster Y position
  Float_t fZ;                   ///< cluster Z position
  Float_t fErrX;                ///< transverse position errors
  Float_t fErrY;                ///< transverse position errors
  Float_t fCharge;              ///< cluster charge
  Float_t fChi2;                ///< cluster chi2
};                              // structure to keep muon clusters information

struct ZDCdata_t {
  // ZDC: it is not clear what is the minimal set of information (PH)

  Int_t     fCollisionsID;      /// The index of the collision vertex

  Float_t   fZEM1Energy;      ///< E in ZEM1
  Float_t   fZEM2Energy;   ///< E in ZEM2

  Float_t   fZNCTowerEnergy[5]; ///< E in 5 ZNC sectors - high gain chain
  Float_t   fZNATowerEnergy[5]; ///< E in 5 ZNA sectors - high gain chain
  Float_t   fZPCTowerEnergy[5]; ///< E in 5 ZPC sectors - high gain chain
  Float_t   fZPATowerEnergy[5]; ///< E in 5 ZPA sectors - high gain chain

  Float_t   fZNCTowerEnergyLR[5]; ///< E in 5 ZNC sectors - low gain chain
  Float_t   fZNATowerEnergyLR[5]; ///< E in 5 ZNA sectors - low gain chain
  Float_t   fZPCTowerEnergyLR[5]; ///< E in 5 ZPC sectors - low gain chain
  Float_t   fZPATowerEnergyLR[5]; ///< E in 5 ZPA sectors - low gain chain

  Float_t   fZDCTDCCorrected[32][4]; /// ZDC TDC data in ns corrected 4 phase shift

  UChar_t   fFired;             /// Bits: 0 - ZNA, 1 - ZNC, 2 - ZPA, 3 - ZPC, 4 - ZEM1, 5 - ZEM2
};                              // structure to keep ZDC information

struct VZEROdata_t {
  /// VZERO as proxy for FIT

  Int_t   fCollisionsID;        /// The index of the collision vertex

  Float_t fAdc[64];             ///  adc for each channel
  Float_t fTime[64];            ///  time for each channel
  Float_t fWidth[64];           ///  time width for each channel
  ULong64_t fBBFlag;            ///  BB Flags from Online V0 Electronics
  ULong64_t fBGFlag;            ///  BG Flags from Online V0 Electronics
};                              // structure to keep VZERO information

struct V0_t {
 /// V0s (Ks, Lambda)

 Int_t fPosTrackID;             // Positive track ID
 Int_t fNegTrackID;             // Negative track ID
};                              // structure to keep v0sinformation

struct Cascade_t {
  /// Cascades

  Int_t fV0sID;                 // V0 ID
  Int_t fTracksID;              // Bachelor track ID
};                              // structure to keep cascades information

} // namespace AliAO2DTypes
#endif

