// $Id: TupleToolMuonIso.h 
#ifndef TUPLETOOLMUONISOt_H 
#define TUPLETOOLMUONISOt_H 1

// Include files
#include "Kernel/IParticleTupleTool.h" 
#include "DecayTreeTupleBase/TupleToolBase.h"

#include "Kernel/IDVAlgorithm.h"
#include "Kernel/IRelatedPVFinder.h"
#include "Event/RecVertex.h"
#include "LoKi/ParticleCuts.h"
#include "LoKi/AParticleCuts.h"
#include "LoKi/ParticleContextCuts.h"
#include "TH1D.h"
/** @class TupleToolMuonIso TupleToolMuonIso.h

 * 
 *  @author Fatima Soomro
 *  @date   2011-08-08
 */

class IDistanceCalculator;
class IDVAlgorithm;
class IPVReFitter;

class TupleToolMuonIso_LC: public TupleToolBase, virtual public IParticleTupleTool {
public:
  /// Standard constructor
  TupleToolMuonIso_LC( const std::string& type, 
              const std::string& name,
              const IInterface* parent );
  /// Loop over differnt conesizes and fill the variables into the tuple
  virtual StatusCode fill( const LHCb::Particle*
                           , const LHCb::Particle*
                           , const std::string&
                           , Tuples::Tuple& );
  
  
  StatusCode fillIsolation(const LHCb::Particle *part, std::string, Tuples::Tuple&);

  StatusCode CDFIsolation(const LHCb::Particle* B, 
                          const LHCb::Particle::ConstVector  decprod,
                          std::string, Tuples::Tuple&);
  
  
  StatusCode fillMCTruth(const LHCb::Particle *part, std::string, Tuples::Tuple&);
  const LHCb::MCParticle* originof( const LHCb::MCParticle* ) ;
  
  virtual ~TupleToolMuonIso_LC( ); ///< Destructor
  
  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode finalize ();


  
protected:

private:

  //std::vector<const LHCb::Particle*> m_decayParticles;

  int is_tau;
  int is_muon;
  

  double m_pvdis_h;
  double m_pvdis;
  double m_svdis_h;
  double m_svdis;
  int m_tracktype;
  double m_ips;
  double m_fc;
  double m_angle;



  double m_count_mu1;
  double m_count_mu2;
  double m_count_mu3;



  double m_count_mup_Giampi; 
  double m_count_mum_Giampi; 
  double m_count_mup_Fatima;
  double m_count_mum_Fatima;
  double m_count_mum;
  double m_count_mup;
  double m_count_p;

  double m_count_mum_f;
  double m_count_mup_f;
  double m_isoCDFyuri;
  double m_CDFIso;
  double m_cosnk;
  double m_doca_iso;
  double m_doca_tr;
  double m_doca;

  bool m_isMC;
  double m_mcancestorpid1;
  double m_mcancestorpid2;
  double m_mcancestorpid;
  double m_endVertices1;
  double m_endVertices2;
  double m_endVertices;
  double m_typeDecay1 ;   double m_zDecay1 ;
  double m_typeDecay2 ;   double m_zDecay2 ;
  double m_typeOrigin1;   double m_zOrigin1;
  double m_typeOrigin2;   double m_zOrigin2;
  double m_typeOrigin; double m_zOrigin;
  //Particle2MCLinker* m_pLinker;
  std::string m_particlePaths;
  double m_MCI ;
  double m_MCPV ;
  IParticle2MCAssociator* m_p2mcAssoc;
  std::string m_p2mcAssocType;
  
  std::string m_TracksPath;
  std::string m_ParticlePath;
  std::string m_PVInputLocation;		\



  const IRelatedPVFinder* pFinder; 
  IDVAlgorithm* m_dva;
  IPVReFitter* m_pvReFitter;
  const IDistanceCalculator* m_dist;

  std::vector<int> getIso(const LHCb::Particle*);
  void IsHltGood(Gaudi::XYZPoint o, Gaudi::XYZVector p, Gaudi::XYZPoint o_mu,Gaudi::XYZVector p_mu, Gaudi::XYZPoint PV, bool& hltgood, double& fc);
  double pointer (Gaudi::XYZVector vtx, Gaudi::XYZVector P_tr,  Gaudi::XYZVector P_mu);
  double gettheta(const LHCb::Particle* vdau1, const LHCb::Particle* vdau2);
  double getphi(const LHCb::Particle* vdau1, const LHCb::Particle* vdau2);
  double ratio( double p1, double p2);
  double IsClose(const LHCb::Particle* p1, const LHCb::Particle* p2);
  void InCone(Gaudi::XYZPoint o1, Gaudi::XYZVector p1,Gaudi::XYZPoint o2,
              Gaudi::XYZVector p2, Gaudi::XYZPoint& vtx, double& doca, double& angle);
  void closest_point(Gaudi::XYZPoint o,Gaudi::XYZVector p, Gaudi::XYZPoint o_mu,Gaudi::XYZVector p_mu, 
                     Gaudi::XYZPoint& close1, Gaudi::XYZPoint& close2, Gaudi::XYZPoint& vertex, bool& fail);
  double arcosine (Gaudi::XYZVector, Gaudi::XYZVector);
  
};



#endif // TUPLETOOLMUONVARIABLES_H

