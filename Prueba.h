#ifndef analysis_Prueba_h
#define analysis_Prueba_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"
#include <TH1F.h>
#include <TFile.h>
#include <vector>

namespace MA5
{
class Prueba : public AnalyzerBase
{
  INIT_ANALYSIS(Prueba,"Prueba")

  TH1F* plot_deltaR_b1b2;
  TH1F* plot_PT_leptons;
  TH1F* plot_PT_b1;
  TH1F* plot_MET;
  TH1F* plot_sdETA_b1b2;

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private: 
  std::vector<const MCParticleFormat*> lepton_final_state_array;
  std::vector<const MCParticleFormat*> neutrino_final_state_array;
  std::vector<const MCParticleFormat*> j_final_state_array;
  std::vector<const MCParticleFormat*> b_final_state_array;
  std::vector<const MCParticleFormat*> j1_final_state_array;
  std::vector<const MCParticleFormat*> j2_final_state_array;
  std::vector<const MCParticleFormat*> b1_final_state_array;
  std::vector<const MCParticleFormat*> b2_final_state_array;
  std::vector<const MCParticleFormat*> b3_final_state_array;
  std::vector<const MCParticleFormat*> b4_final_state_array;

  MAbool is_lepton_final_state(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-15)&&(part->pdgid()!=-13)&&(part->pdgid()!=-11)&&(part->pdgid()!=11)&&(part->pdgid()!=13)&&(part->pdgid()!=15) ) return false;
     return true; }
  
  MAbool is_neutrino_final_state(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-16)&&(part->pdgid()!=-14)&&(part->pdgid()!=-12)&&(part->pdgid()!=12)&&(part->pdgid()!=14)&&(part->pdgid()!=16) ) return false;
     return true; }

  MAbool is_j_final_state(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-4)&&(part->pdgid()!=-3)&&(part->pdgid()!=-2)&&(part->pdgid()!=-1)&&(part->pdgid()!=1)&&(part->pdgid()!=2)&&(part->pdgid()!=3)&&(part->pdgid()!=4)&&(part->pdgid()!=21) ) return false;
     return true; }
  
  MAbool is_b_final_state(const MCParticleFormat* part) const {
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-5)&&(part->pdgid()!=5) ) return false;
     return true; }


};
}

#endif
