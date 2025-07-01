// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file candidateSelectorXic0XicpToHadronic.cxx
/// \brief Selection of Xic0 and Xicp candidates 
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

#include <vector>
#include <string>

// Mandatory includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// For selection
#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h" // -> findBin function

// Related to ML 

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorXic0XicpToHadronic {

  // cursors to fill the tables
  struct : ProducesGroup {
    Produces<aod::HfSelXic0ToXiPi> hfSelXic0ToXiPiCandidate;
    Produces<aod::HfMlXic0ToXiPi> hfMlXic0ToXiPiCandidate;
    Produces<aod::HfSelXicpToXiPiPi> hfSelXicpToXiPiPiCandidate;
    Produces<aod::HfMlXicpToXiPiPi> hfMlXicpToXiPiPiCandidate;
  } cursors;

  struct : ConfigurableGroup {

    Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
    Configurable<float> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
    // Topological cuts Xic0
    Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic0_xicp_to_hadronic::Cuts[0], 
	  hf_cuts_xic0_xicp_to_hadronic::NBinsPt, 
	  hf_cuts_xic0_xicp_to_hadronic::NCutVars,
	  hf_cuts_xic0_xicp_to_hadronic::labelsPt,
	  hf_cuts_xic0_xicp_to_hadronic::labelsCutVar}, "Xic0 candidate selection per pT bin"};
    Configurable<std::vector<double>> binsPt{"binsPt", std::vector{hf_cuts_xic0_xicp_to_hadronic::vecBinsPt}, "pT bin limits"};
     // Topological cuts Xicp
    Configurable<LabeledArray<double>> cutsXicp{"Xicp cuts", {hf_cuts_xic0_xicp_to_hadronic::CutsXicp[0], 
	  hf_cuts_xic0_xicp_to_hadronic::NBinsPtXicp, 
	  hf_cuts_xic0_xicp_to_hadronic::NCutVarsXicp,
	  hf_cuts_xic0_xicp_to_hadronic::labelsPtXicp,
	  hf_cuts_xic0_xicp_to_hadronic::labelsCutVarXicp}, "Xicp candidate selection per pT bin"};
    Configurable<std::vector<double>> binsPtXicp{"binsPtXicp", std::vector{hf_cuts_xic0_xicp_to_hadronic::vecBinsPtXicp}, "pT bin limits"};
    // QA switch
    Configurable<bool> activateQA{"activateQA", true, "Flag to enable QA histograms"};
    // Enable PID
    Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
    Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept status::NotApplicable [(NotApplicable for one detector) and NotApplicable or Conditional for the other)] in PID selection"};
    // TPC PID
    Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
    Configurable<float> ptPidTpcMax{"ptPidTpcMax", 20., "Lower bound of track pT for TPC PID"};
    Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma on TPC only"};
    Configurable<float> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
    // TOF PID
    Configurable<float> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
    Configurable<float> ptPidTofMax{"ptPidTofMax", 20., "Lower bound of track pT for TOF PID"};
    Configurable<float> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma on TOF only"};
    Configurable<float> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
    // ML interface
    // [empty for now]
    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "https://alice-ccdb.cern.ch", "url of the ccdb repository"};

  } configs;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;

  std::vector<float> outputMlXicToXiPiPi = {};
  
  using TracksPidWithSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    // Initialize TrackSelector
    if (configs.usePid) {
      // pion
      selectorPion.setRangePtTpc(configs.ptPidTpcMin, configs.ptPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-configs.nSigmaTpcMax, configs.nSigmaTpcMax);
      selectorPion.setRangeNSigmaTpcCondTof(-configs.nSigmaTpcCombinedMax, configs.nSigmaTpcCombinedMax);
      selectorPion.setRangePtTof(configs.ptPidTofMin, configs.ptPidTofMax);
      selectorPion.setRangeNSigmaTof(-configs.nSigmaTofMax, configs.nSigmaTofMax);
      selectorPion.setRangeNSigmaTofCondTpc(-configs.nSigmaTofCombinedMax, configs.nSigmaTofCombinedMax);
      // proton 
      selectorProton.setRangePtTpc(configs.ptPidTpcMin, configs.ptPidTpcMax);
      selectorProton.setRangeNSigmaTpc(-configs.nSigmaTpcMax, configs.nSigmaTpcMax);
      selectorProton.setRangeNSigmaTpcCondTof(-configs.nSigmaTpcCombinedMax, configs.nSigmaTpcCombinedMax);
      selectorProton.setRangePtTof(configs.ptPidTofMin, configs.ptPidTofMax);
      selectorProton.setRangeNSigmaTof(-configs.nSigmaTofMax, configs.nSigmaTofMax);
      selectorProton.setRangeNSigmaTofCondTpc(-configs.nSigmaTofCombinedMax, configs.nSigmaTofCombinedMax);
    }

    // Initializing QA histogram
    if (configs.activateQA) {
      constexpr int kNBinsSelections = 1 + SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + SelectionStep::RecoMl] = "Skims & Topological & PID & ML selections";

      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections+0.5, ""};

      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {kTH2F, {axisSelections, {(std::vector<double>)configs.binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
			
      for (int iBin=0; iBin<kNBinsSelections; ++iBin) {
	registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin+1, labels[iBin].data());
      }
    }
		
    // Initializing ML
    // [ empty for now ]
  }

  // Conjugate independent topological cuts
  // \param candidate is Xic0 candidate
  // \return true if candidate passes all cuts
  template <typename T1>
  bool selectionTopol(const T1& hfCandXic0)
  {
    auto candPt = hfCandXic0.pt();
    int ptBin = findBin(configs.binsPt, candPt);
    if (ptBin==-1) {
      return false;
    }
	
    // check if the candidate pt is within analysis range
    if (candPt < configs.ptCandMin || candPt >= configs.ptCandMax) {
      return false;
    }

    // check if candidate Xic0  mass is within a defined mass window
    if (std::abs(hfCandXic0.invMassXic0() - o2::constants::physics::MassXiC0) > configs.cuts->get(ptBin, "m")) {
      return false;
    }

    // cosine of pointing angle
    if (hfCandXic0.cpa() <= configs.cuts->get(ptBin, "cosine pointing angle")) {
      return false;
    }

    // cosine of pointing angle xy
    if (hfCandXic0.cpaXY() <= configs.cuts->get(ptBin, "cosine pointing angle XY")) {
      return false;	
    }
		
#if 1 
    // maximum decay length
    if (hfCandXic0.decayLength() > configs.cuts->get(ptBin, "max decay length")) {
      return false;
    }

    // maximum decay length XY
    if (hfCandXic0.decayLengthXY() > configs.cuts->get(ptBin, "max decay length XY")) {
      return false;
    }

    // candidate chi2PCA
    if (hfCandXic0.chi2PCA() > configs.cuts->get(ptBin, "chi2PCA")) {
      return false;
    }

    // maximum DCA of daughters of xic0
    if (std::abs(hfCandXic0.impactParameter0()) > configs.cuts->get(ptBin, "max impParXY Xi") ||
	std::abs(hfCandXic0.impactParameter1()) > configs.cuts->get(ptBin, "max impParXY Pi")) {
      return false;
    }

    // cuts on daughter pT
    if (hfCandXic0.ptProng0() < configs.cuts->get(ptBin, "pT Xi") ||
	hfCandXic0.ptProng1() < configs.cuts->get(ptBin, "pT Pi")) {
      return false;
    }
#endif	
    // candidate reached by here passed all the topoplodical cuts
    return true;
    
  } // end selection toplogy Xic0

  // Apply Xic0 PID selection
  // \param pidTrackPi0 PID status of trackPi(Prong1 of Xic0 candidate)
  // \param pidTrackPr PID status of trackPr(positive daughter from V0 candidate
  // \param pidTrackPiLam PID status of trackPiLam (negative daughter from V0 candidate)
  // \param pidTrackPiXi PID status of trackPiXi (bachelor candidate of cascade candidate)
  // \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  // \return true if prongs of Xic0 candidate passes all selections	
  bool selectionPID(TrackSelectorPID::Status const pidTrackPi0,
		    TrackSelectorPID::Status const pidTrackPr,
		    TrackSelectorPID::Status const pidTrackPiLam,
		    TrackSelectorPID::Status const pidTrackPiXi,
		    bool const acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && (pidTrackPi0 != TrackSelectorPID::Accepted ||
				    pidTrackPr != TrackSelectorPID::Accepted ||
				    pidTrackPiLam != TrackSelectorPID::Accepted ||
				    pidTrackPiXi != TrackSelectorPID::Accepted)) {
      return false;
    }

    if (acceptPIDNotApplicable && (pidTrackPi0 == TrackSelectorPID::Rejected ||
				   pidTrackPr == TrackSelectorPID::Rejected ||
				   pidTrackPiLam == TrackSelectorPID::Rejected ||
				   pidTrackPiXi == TrackSelectorPID::Rejected)) {
      return false;
    }
	
    return true;
  }

  // process Xic0 function
  void process(aod::HfCandXic0 const& hfCandsXic0,
	       TracksPidWithSel const&)
  {
    for (auto const& hfCandXic0 : hfCandsXic0) {

      int statusXic0ToXiPi=0;
      auto ptCandXic0 = hfCandXic0.pt();

      if (configs.activateQA) {
	registry.fill(HIST("hSelections"), 1, ptCandXic0);	
      }

      // No hfflag -> By default skim selected
      SETBIT(statusXic0ToXiPi, SelectionStep::RecoSkims);
      if (configs.activateQA) {
	registry.fill(HIST("hSelections"), 2+SelectionStep::RecoSkims, ptCandXic0);
      }
	
#if 1 
      // Topological cuts
      if (!selectionTopol(hfCandXic0)) {
	cursors.hfSelXic0ToXiPiCandidate(statusXic0ToXiPi);
	continue;
      }
      SETBIT(statusXic0ToXiPi, SelectionStep::RecoTopol);
      if (configs.activateQA) {
	registry.fill(HIST("hSelections"), 2+SelectionStep::RecoTopol, ptCandXic0);
      }
#endif	

      // Track-level PID selection
      if (configs.usePid) {
			
	// Get tracks
	auto trackPi = hfCandXic0.pi_as<TracksPidWithSel>();
	auto trackV0Pos = hfCandXic0.posTrack_as<TracksPidWithSel>();
	auto trackV0Neg = hfCandXic0.negTrack_as<TracksPidWithSel>();
	auto trackPiFromXi = hfCandXic0.bachelor_as<TracksPidWithSel>();

	// assign proton and pion hypothesis by cascade sign. By default, Xi sign is considered negative
	auto trackPrFromLambda = trackV0Pos;
	auto trackPiFromLambda = trackV0Neg;
	if (hfCandXic0.cascSign() > 0) {
					
	  // Brief sign check
	  if (trackPiFromXi.sign() < 0) {
	    LOG(info) << ">>>>>>>>>>>>>>> Cascade sign and Pi from cascade sign doesn't match";
	  }

	  trackPrFromLambda = trackV0Neg;
	  trackPiFromLambda = trackV0Pos;
	}
				
	// PID info
	TrackSelectorPID::Status pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
	TrackSelectorPID::Status pidTrackPiFromXi = selectorPion.statusTpcAndTof(trackPiFromXi);
	TrackSelectorPID::Status pidTrackPrFromLambda = selectorProton.statusTpcAndTof(trackPrFromLambda);
	TrackSelectorPID::Status pidTrackPiFromLambda = selectorPion.statusTpcAndTof(trackPiFromLambda);

	if (!selectionPID(pidTrackPi, pidTrackPiFromXi, pidTrackPrFromLambda, pidTrackPiFromLambda, configs.acceptPIDNotApplicable.value)) {
	  cursors.hfSelXic0ToXiPiCandidate(statusXic0ToXiPi);
	  continue;
	}

	SETBIT(statusXic0ToXiPi, SelectionStep::RecoPID);
	if (configs.activateQA) {
	  registry.fill(HIST("hSelections"), 2+SelectionStep::RecoPID, ptCandXic0);
	}

	cursors.hfSelXic0ToXiPiCandidate(statusXic0ToXiPi); // -> This part will be moved into if (configs.applyMl) part later
      }

    }// candidate loop
  } // end process Xic0

  ////////////////////////////////////////////////////////////////////////////////////// Xicp begins
   template <typename T2>
  bool selectionTopolXicp(const T2& hfCandXic)
  {
    auto candpT = hfCandXic.pt();
    int pTBin = findBin(configs.binsPtXicp, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < configs.ptCandMin || candpT >= configs.ptCandMax) {
      return false;
    }
      
    // check candidate mass is within a defined mass window
    if (std::abs(hfCandXic.invMassXicPlus() - o2::constants::physics::MassXiCPlus) > configs.cutsXicp->get(pTBin, "m")) {
      return false;
    }
  
    // cosine of pointing angle
    if (hfCandXic.cpa() <= configs.cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // cosine of pointing angle XY
    if (hfCandXic.cpaXY() <= configs.cutsXicp->get(pTBin, "cos pointing angle XY")) {
      return false;
    }

    // candidate maximum decay length
    if (hfCandXic.decayLength() > configs.cutsXicp->get(pTBin, "max decay length")) {
      return false;
    }

    // candidate maximum decay length XY
    if (hfCandXic.decayLengthXY() > configs.cutsXicp->get(pTBin, "max decay length XY")) {
      return false;
    }

    // candidate chi2PC
    if (hfCandXic.chi2PCA() > configs.cutsXicp->get(pTBin, "chi2PCA")) {
      return false;
    }

    // maximum DCA of daughters
   
    if ((std::abs(hfCandXic.impactParameter0()) > configs.cutsXicp->get(pTBin, "max impParXY Xi")) ||
	(std::abs(hfCandXic.impactParameter1()) > configs.cutsXicp->get(pTBin, "max impParXY Pi0")) ||
	(std::abs(hfCandXic.impactParameter2()) > configs.cutsXicp->get(pTBin, "max impParXY Pi1"))) {
      return false;
    }

    // cut on daughter pT
        
    if (hfCandXic.ptProng0() < configs.cutsXicp->get(pTBin, "pT Xi") ||
	hfCandXic.ptProng1() < configs.cutsXicp->get(pTBin, "pT Pi0") ||
	hfCandXic.ptProng2() < configs.cutsXicp->get(pTBin, "pT Pi1")) {
      return false;
    }
      
    return true;
  } // end selection toplogy Xicp
 
  bool selectionPidXicp(TrackSelectorPID::Status const pidTrackPi0,
			TrackSelectorPID::Status const pidTrackPi1,
			TrackSelectorPID::Status const pidTrackPr,
			TrackSelectorPID::Status const pidTrackPiLam,
			TrackSelectorPID::Status const pidTrackPiXi,
			bool const acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && (pidTrackPi0 != TrackSelectorPID::Accepted || pidTrackPi1 != TrackSelectorPID::Accepted || pidTrackPr != TrackSelectorPID::Accepted || pidTrackPiLam != TrackSelectorPID::Accepted || pidTrackPiXi != TrackSelectorPID::Accepted)) {
      return false;
    }
    if (acceptPIDNotApplicable && (pidTrackPi0 == TrackSelectorPID::Rejected || pidTrackPi1 == TrackSelectorPID::Rejected || pidTrackPr == TrackSelectorPID::Rejected || pidTrackPiLam == TrackSelectorPID::Rejected || pidTrackPiXi == TrackSelectorPID::Rejected)) {
      return false;
    }
    return true;
    
  }
  
  void processXicp(aod::HfCandXic const& hfCandsXic,
		   TracksPidWithSel const&)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicpToXiPiPi = 0;
      
      outputMlXicToXiPiPi.clear(); 

      auto ptCandXic = hfCandXic.pt();

      if (configs.activateQA) {
        registry.fill(HIST("hSelections"), 1, ptCandXic);
      }

      // No hfflag -> by default skim selected
      SETBIT(statusXicpToXiPiPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusXicToXiPiPi = 1
      if (configs.activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandXic);
      }

      // topological cuts
      if (!selectionTopolXicp(hfCandXic)) {
        cursors.hfSelXicpToXiPiPiCandidate(statusXicpToXiPiPi);
	continue;
      }
      SETBIT(statusXicpToXiPiPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusXicToXiPiPi = 3
      if (configs.activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandXic);
      }

      // track-level PID selection
      if (configs.usePid) {
        auto trackPi0 = hfCandXic.pi0_as<TracksPidWithSel>();
        auto trackPi1 = hfCandXic.pi1_as<TracksPidWithSel>();
        auto trackV0PosDau = hfCandXic.posTrack_as<TracksPidWithSel>();
        auto trackV0NegDau = hfCandXic.negTrack_as<TracksPidWithSel>();
        auto trackPiFromXi = hfCandXic.bachelor_as<TracksPidWithSel>();
        // assign proton and pion hypothesis to V0 daughters
        auto trackPr = trackV0PosDau;
        auto trackPiFromLam = trackV0NegDau;
        if (hfCandXic.sign() < 0) {
          trackPr = trackV0NegDau;
          trackPiFromLam = trackV0PosDau;
        }
        // PID info
        TrackSelectorPID::Status pidTrackPi0 = selectorPion.statusTpcAndTof(trackPi0);
        TrackSelectorPID::Status pidTrackPi1 = selectorPion.statusTpcAndTof(trackPi1);
        TrackSelectorPID::Status pidTrackPr = selectorProton.statusTpcAndTof(trackPr);
        TrackSelectorPID::Status pidTrackPiLam = selectorPion.statusTpcAndTof(trackPiFromLam);
        TrackSelectorPID::Status pidTrackPiXi = selectorPion.statusTpcAndTof(trackPiFromXi);

        if (!selectionPidXicp(pidTrackPi0, pidTrackPi1, pidTrackPr, pidTrackPiLam, pidTrackPiXi, configs.acceptPIDNotApplicable.value)) {
          cursors.hfSelXicpToXiPiPiCandidate(statusXicpToXiPiPi);
	  continue;
        }
        SETBIT(statusXicpToXiPiPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusXicToXiPiPi = 7
        if (configs.activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandXic);
        }
      }

      // ML selection
      // [ empty for now ]

      cursors.hfSelXicpToXiPiPiCandidate(statusXicpToXiPiPi);
    }
  }

}; // end struct candidate selector

    WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
    {
      return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXic0XicpToHadronic>(cfgc)};
    }
