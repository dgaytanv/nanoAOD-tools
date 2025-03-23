from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib

class jetvetomapProducer(Module):
    def __init__(self, isMC, era, corrName, veto_map_name = "jetvetomap"):
        """Module to apply veto maps that veto out events with important jets in the "hot" or 
        "cold" zones. These maps should be applied similarly both on Data and MC, to keep the
        phase-spaces equal. 
        Parameters:
        """
        self.isMC = isMC
        self.corrName = corrName
        pogdir = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/"
        fname = os.path.join(pogdir, f"POG/JME/{era}/jetvetomaps.json.gz")
        self.veto_map_name = veto_map_name
        print('jetvetomap file', fname)
        self.evaluator_JVMAP= correctionlib.CorrectionSet.from_file(fname)
        print('jetvetomap name', corrName)
        self.evaluator_VETO = self.evaluator_JVMAP[corrName]
        self.lenVar = "nJet"
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if '16' in self.corrName or '17' in self.corrName or '18' in self.corrName:
            self.out.branch("Jet_veto_flag", "I", lenVar=self.lenVar)
        else:
            self.out.branch("Flag_JetVetoed", "I", title="Event veto flag from Jet Veto Map")
        
    
    def fixPhi(self, phi):
        epsilon = 1e-6  # Small offset to avoid boundary issues
        if phi > np.pi:
            #print(f"phi {phi} is greater than pi. Setting phi to pi - epsilon.")
            phi = np.pi - epsilon
        elif phi < -np.pi:
            #print(f"phi {phi} is less than -pi. Setting phi to -pi + epsilon.")
            phi = -np.pi + epsilon
        return phi

    def analyze(self, event):
        '''nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        - jets that don't overlap with PF muon (dR < 0.2)
        '''
        jets = Collection(event, "Jet")
        
        if '16' in self.corrName or '17' in self.corrName or '18' in self.corrName:
            jets_veto_flag = []
            for i, jet in enumerate(jets):
                veto_flag = 0
                if (jet.pt> 15 and (jet.jetId ==2 or jet.jetId ==6) and (jet.chEmEF + jet.neEmEF)<0.9 and jet.muonIdx1 == -1 and jet.muonIdx2 == -1):
                    phi = self.fixPhi(jet.phi)
                    veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi)
                    if veto_map_value > 0:
                        veto_flag=1
                jets_veto_flag.append(veto_flag)
            self.out.fillBranch("Jet_veto_flag", jets_veto_flag)
            return True #add veto flag to jets but don't veto events
        else:
            veto_flag = 0
            for i, jet in enumerate(jets):
                if (jet.pt> 15 and (jet.jetId ==2 or jet.jetId ==6) and (jet.chEmEF + jet.neEmEF)<0.9 and jet.muonIdx1 == -1 and jet.muonIdx2 == -1):
                    # Correct phi and evaluate veto map
                    phi = self.fixPhi(jet.phi)
                    veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi)

                    # Check if the jet is vetoed
                    if veto_map_value > 0:
                        veto_flag = 1  # Set flag if a vetoed jet is found
                        break  # Break out of the loop since we only need one veto to trigger
            # Fill the branch with the veto result
            self.out.fillBranch("Flag_JetVetoed", veto_flag)
            if veto_flag == 0:
                return True
            else:
                return False #veto event if one jet is found for run 3
            
            
jetvetomapUL2017 = lambda: jetvetomapProducer(
    True, "2017_UL", "Summer19UL17_V2")
jetvetomapUL2018 = lambda: jetvetomapProducer(
    True, "2018_UL", "Summer19UL18_V1")
jetvetomap2022 = lambda: jetvetomapProducer(
    True, "2022_Summer22", "Summer22_23Sep2023_RunCD_V1")
#for testing at nanoaod script locally only