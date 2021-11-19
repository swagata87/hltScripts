#!/usr/bin/env python

#RUN it like this
#python3 filename.py  inputFile.root   -o=outFile.root 

from array import array
import re
import argparse
import sys
import math
from DataFormats.FWLite import Events, Handle
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
import json
import FWCore.ParameterSet.Config as cms
from CondCore.CondDB.CondDB_cfi import *
from PhysicsTools.PythonAnalysis import *
import time

def getListFilterPassedObj(filterName,hltsevt):
        eg_trig_objs = []
        #get a list of trigger objects that passed a given filter
        filterIndex = getFilterIndex(hltsevt,filterName)
        if (filterIndex < hltsevt.sizeFilters() ):
                for filterKey in hltsevt.filterKeys(filterIndex):
                        obj = hltsevt.getObjects()[filterKey]
                        eg_trig_objs.append(obj)
        return eg_trig_objs


#from https://github.com/cms-egamma/EgammaDAS2020/blob/solutions/test/egFWLiteExercise2a.py
def match_trig_objs(eta,phi,trig_objs,max_dr=0.1):    
    max_dr2 = max_dr*max_dr
    matched_objs = [obj for obj in trig_objs if ROOT.reco.deltaR2(eta,phi,obj.eta(),obj.phi()) < max_dr2]
    return matched_objs


#from https://github.com/Sam-Harper/usercode/blob/100XNtup/SHNtupliser/test/checkTrigsAOD.py
def getFilterIndex(trigEvt,filterName):
    for index in range(0,trigEvt.sizeFilters()):
        if filterName==trigEvt.filterLabel(index):
                #print 'filtername found' 
                return index
    return trigEvt.sizeFilters()


def get_genparts(genparts,pid=11,antipart=True,status=1):
    """                                                                                                                                                                               returns a list of the gen particles matching the given criteria from hard process                                                                                                 might not work for all generators as depends on isHardProcess()                                                                                                                   """
    selected = []
    if genparts==None:
        return selected

    for part in genparts:
        pdg_id = part.pdgId()
        if pdg_id == pid or (antipart and abs(pdg_id) == abs(pid)):
            if part.isHardProcess():
                if status == 1:
                    selected.append(part)
    return selected


def match_to_gen(eta,phi,genparts,pid=11,antipart=True,max_dr=0.1,status=1):
    """                                                                                                                                                                               Matches an eta,phi to gen level particle from the hard process                                                                                                                    might not work for all generaters as depends on isHardProcess()                                                                                                                   """
    best_match = None
    best_dr2 = max_dr*max_dr
    selected_parts = get_genparts(genparts,pid,antipart,status)
    for part in selected_parts:
        dr2 = ROOT.reco.deltaR2(eta,phi,part.eta(),part.phi())
        if dr2 < best_dr2:
            best_match = part
            best_dr2 = dr2
    return best_match,best_dr2


if __name__ == "__main__":

    oldargv = sys.argv[:]
    sys.argv = [ '-b-' ]
    sys.argv = oldargv
    ROOT.gSystem.Load("libFWCoreFWLite.so");
    ROOT.gSystem.Load("libDataFormatsFWLite.so");
    ROOT.FWLiteEnabler.enable()

    parser = argparse.ArgumentParser(description='example e/gamma HLT analyser')
    parser.add_argument('in_filename',nargs="+",help='input filename')
    parser.add_argument("-o", "--output", help="Directs the output to a name of your choice")

    args = parser.parse_args()

    ele_handle, ele_label = Handle("std::vector<trigger::EgammaObject>"), "hltEgammaHLTExtra"
    #ele_handle, ele_label = Handle("std::vector<reco::EgTrigSumObj>"), "hltEgammaHLTExtra"
    hlt_handle, hlt_label = Handle("edm::TriggerResults"), "TriggerResults::HLTX"
    hltevt_handle, hltevt_label = Handle("trigger::TriggerEvent"), "hltTriggerSummaryAOD::HLTX"
    gen_handle, gen_label = Handle("std::vector<reco::GenParticle>"), "genParticles"

    to_remove=[]
    for file in args.in_filename:
        file_temp = ROOT.TFile(file)
        if ( file_temp.IsZombie() ) :
            to_remove.append(file)

    new_list = [x for x in args.in_filename if x not in to_remove]

    events = Events(new_list)

    den_ele_eta_ele32 = ROOT.TH1D("den_ele_eta",";#eta;nEntries",40,-4,4)
    den_ele_pt_EB = ROOT.TH1D("den_ele_pt_EB",";pT;nEntries",100,0,500)
    den_ele_pt_EE = ROOT.TH1D("den_ele_pt_EE",";pT;nEntries",100,0,500)

    num_ele_eta_hltEle32WPTightGsfTrackIsoFilter = ROOT.TH1D("num_ele_eta_hltEle32WPTightGsfTrackIsoFilter",";#eta;nEntries",40,-4,4)
    num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EB = ROOT.TH1D("num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EB",";pT;nEntries",100,0,500)
    num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EE = ROOT.TH1D("num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EE",";pT;nEntries",100,0,500)
    
    for event_nr,event in enumerate(events):

        event.getByLabel(ele_label,ele_handle)
        event.getByLabel(hlt_label,hlt_handle)
        event.getByLabel(hltevt_label,hltevt_handle)
        event.getByLabel(gen_label,gen_handle)

        eles = ele_handle.product()
        hlts = hlt_handle.product()
        hltsevt = hltevt_handle.product()
        genobjs=gen_handle.product()

        trigdict=event.object().triggerNames(hlts).triggerNames()

        eg_trig_objs_hltEle32WPTightGsfTrackIsoFilter = getListFilterPassedObj("hltEle32WPTightGsfTrackIsoFilter",hltsevt)

        for eg in eles:

          #check if the electron matches with any trigger object that passed a given filter 
          matched_objs_hltEle32WPTightGsfTrackIsoFilter = match_trig_objs(eg.eta(),eg.phi(),eg_trig_objs_hltEle32WPTightGsfTrackIsoFilter)
          nmatch_hltEle32WPTightGsfTrackIsoFilter = len(matched_objs_hltEle32WPTightGsfTrackIsoFilter)
          gen_match_ele = match_to_gen(eg.eta(),eg.phi(),gen_handle.product(),pid=11)[0]

          if (gen_match_ele):

                  ## Fill denominators
                  if (abs(eg.eta()) < 1.44 ): den_ele_pt_EB.Fill(eg.pt())
                  if (abs(eg.eta()) > 1.56 ): den_ele_pt_EE.Fill(eg.pt())
                  if (eg.pt()>=40.0): 
                          den_ele_eta_ele32.Fill(eg.eta())

                  ## Fill numerators
                  if (nmatch_hltEle32WPTightGsfTrackIsoFilter>0) :
                      if (abs(eg.eta()) < 1.44 ): num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EB.Fill(eg.pt())
                      if (abs(eg.eta()) > 1.56 ): num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EE.Fill(eg.pt())
                      if (eg.pt()>=40.0): 
                        num_ele_eta_hltEle32WPTightGsfTrackIsoFilter.Fill(eg.eta())


    output_file = TFile( args.output, 'recreate' )

    den_ele_eta_ele32.Write()
    den_ele_pt_EB.Write()
    den_ele_pt_EE.Write()

    num_ele_eta_hltEle32WPTightGsfTrackIsoFilter.Write()
    num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EB.Write()
    num_ele_pt_hltEle32WPTightGsfTrackIsoFilter_EE.Write()

    output_file.Close()
