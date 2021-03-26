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
  bool notMerged;
  bool partiallyMerged;
  bool fullyMerged;
  int index_b;
  MA5::MALorentzVector dijet;
  MA5::MALorentzVector b_dijet;
  MA5::MALorentzVector reconstructed_W;

  MAbool is_lepton_final_state(const MCParticleFormat* part){
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-15)&&(part->pdgid()!=-13)&&(part->pdgid()!=-11)&&(part->pdgid()!=11)&&(part->pdgid()!=13)&&(part->pdgid()!=15) ) return false;
     return true; }
  
  MAbool is_neutrino_final_state(const MCParticleFormat* part){
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-16)&&(part->pdgid()!=-14)&&(part->pdgid()!=-12)&&(part->pdgid()!=12)&&(part->pdgid()!=14)&&(part->pdgid()!=16) ) return false;
     return true; }

  MAbool is_j_final_state(const MCParticleFormat* part){
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-4)&&(part->pdgid()!=-3)&&(part->pdgid()!=-2)&&(part->pdgid()!=-1)&&(part->pdgid()!=1)&&(part->pdgid()!=2)&&(part->pdgid()!=3)&&(part->pdgid()!=4)&&(part->pdgid()!=21) ) return false;
     return true; }
  
  MAbool is_b_final_state(const MCParticleFormat* part){
     if ( part==0 ) return false;
     if ( !PHYSICS->Id->IsFinalState(part) ) return false;
     if ( (part->pdgid()!=-5)&&(part->pdgid()!=5) ) return false;
     return true; }

  void jets_Merge(std::vector<const MCParticleFormat*> j1, std::vector<const MCParticleFormat*> j2, std::vector<const MCParticleFormat*> b, MA5::MALorentzVector &dijet, MA5::MALorentzVector &b_dijet, MA5::MALorentzVector &reconstructed_W){

     const MALorentzVector j1_p = j1[0]->momentum();
     const MALorentzVector j2_p = j2[0]->momentum();
     dijet = j1_p + j2_p;
     Double_t dr_dijet = j1_p.DeltaR(j2_p);
     //std::cout << dijet.Mag() << std::endl;
     if (dr_dijet > 0.8)
     {
         notMerged = true;
         partiallyMerged = false ;
         fullyMerged = false;
     }
  
     else 
     {
         notMerged = false;
         Double_t smallest_dr_b_dijet = 999999.0;
         for (MAuint32 k=0;k<b.size();k++)
         {
              const MALorentzVector b_p = b[k]->momentum();
              Double_t dr_b_dijet = b_p.DeltaR(dijet);
              if (dr_b_dijet < smallest_dr_b_dijet)
              {
                  smallest_dr_b_dijet = dr_b_dijet;
                  index_b = k;
              }
         }

         if (smallest_dr_b_dijet < 1.0)
         {  
            partiallyMerged = false ;
            fullyMerged = true;
            b_dijet = b[index_b]->momentum() + dijet;
         }

         else
         {
            partiallyMerged = true ;
            fullyMerged = false;
            reconstructed_W = dijet;
            b_dijet = b[index_b]->momentum() + dijet;
         }
   
     }



     }


};
}

#endif
