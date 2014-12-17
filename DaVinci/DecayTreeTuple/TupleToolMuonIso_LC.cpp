// Include files
#include "GaudiKernel/ToolFactory.h" 
#include "Event/Particle.h"  
// kernel 
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Kernel/IParticle2MCAssociator.h"
#include <Kernel/IDistanceCalculator.h>
#include "Kernel/IPVReFitter.h"
//#include "Kernel/IOnOffline.h"
#include "Kernel/IDVAlgorithm.h"
#include <Kernel/GetIDVAlgorithm.h>

// MC stuff
#include "Event/GenHeader.h" 
#include "Event/MCHeader.h" 
#include "TH1F.h"
#include "TH1D.h"
#include "TupleToolMuonIso_LC.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include <TROOT.h>
#include <TObject.h>
#include "TH1D.h"
//#include "TMVA/Reader.h"
#include "Kernel/IRelatedPVFinder.h"
#include "Math/Boost.h"

// STL                                                                                                                                
#include <string>
#include <vector>
// GSL                                                                                                                                
#include "gsl/gsl_math.h"
//from Gaudi                                                                                                                          
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IAlgTool.h"

//from Kernel                                                                                                                         
#include <Kernel/GetIDVAlgorithm.h>
#include "Kernel/DaVinciAlgorithm.h"
// from Event                                                                                                                         
#include "Event/Track.h"
#include "Event/RecVertex.h"
// Local                                                                                                                              

#include "Kernel/DVAlgorithm.h"
#include "Kernel/DaVinciAlgorithm.h"
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
// from Gaudi                                                                                                                         
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
// from Event                                                                                                                         
#include "Event/RecHeader.h"
#include "Event/Track.h"
#include "Event/RecVertex.h"
#include "Event/MCParticle.h"
#include "Kernel/IJetMaker.h"
#include <iostream>     // std::cout                                                                                                  
#include <algorithm>    // std::sort                                                                                                  
#include <vector>       // std::vector        



 
//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolMuonIso for Kstarmumu
//
// 2012-10-28 : Justine Serrano 
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_TOOL_FACTORY( TupleToolMuonIso_LC );


const LHCb::RecVertex* closestPV(const LHCb::Particle* part);


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolMuonIso_LC::TupleToolMuonIso_LC( const std::string& type,
                                                  const std::string& name,
                                                  const IInterface* parent) 
  : TupleToolBase ( type, name , parent )
  , m_dva(0)
  , m_dist(0)


{
  declareInterface<IParticleTupleTool>(this);
  declareProperty( "tracktype", m_tracktype = 3,
                   "Set the type of tracks which are considered inside the cone (default = 3)");
  
  declareProperty("ParticlePath",   m_ParticlePath="/Event/Phys/StdAllNoPIDsPions/Particles");
  declareProperty("PVInputLocation",m_PVInputLocation = LHCb::RecVertexLocation::Primary);
  declareProperty ( "TrackContainer",   m_TracksPath = LHCb::TrackLocation::Default); // default is "Rec/Track/Best "
  declareProperty("angle"     , m_angle  = 0.27     ); // 
  declareProperty("fc"        , m_fc  = 0.60     ); // 
  declareProperty("doca_iso"  , m_doca_iso  = 0.13     ); // 
  declareProperty("ips"       , m_ips  = 3.0     ); // 
  declareProperty("svdis"     , m_svdis  = -0.15     ); // 
  declareProperty("svdis_h"   , m_svdis_h  = 30.     ); // 
  declareProperty("pvdis"     , m_pvdis  = 0.5     ); // 
  declareProperty("pvdis_h"   , m_pvdis_h  = 40.    ); // 
  declareProperty("isMC", m_isMC = true); 
  declareProperty( "IP2MCPAssociatorType", m_p2mcAssocType =  "DaVinciSmartAssociator");


}

//=============================================================================
// Destructor
//=============================================================================
TupleToolMuonIso_LC::~TupleToolMuonIso_LC() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode TupleToolMuonIso_LC::initialize() {
  StatusCode sc = TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;
  
  
  if(m_isMC){
    m_p2mcAssoc = tool<IParticle2MCAssociator>(m_p2mcAssocType, this);
    if (0==m_p2mcAssoc) return Error("Couldn't get MC associator!! ",  StatusCode::FAILURE);  
  }
  
  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc() ) ;
  if (0==m_dva) return Error("Couldn't get parent DVAlgorithm", 
                             StatusCode::FAILURE);  

  Assert( m_p2mcAssoc
          , "The DaVinci smart associator hasn't been initialized!");

  //if ( "" == m_PVInputLocation){
  //  const IOnOffline* oo = tool<IOnOffline>("OnOfflineTool",this);
  //  m_PVInputLocation = oo->primaryVertexLocation();
  //  info() << "Fatima, pay attention. Will be looking for PVs at " << m_PVInputLocation << endmsg ;
    //debug() << " Will be looking for PVs at " << m_PVInputLocation << endmsg ;
  // } 
  
  m_dist = m_dva->distanceCalculator ();
  if( !m_dist ){
    Error("Unable to retrieve the IDistanceCalculator tool");
    return StatusCode::FAILURE;
  }

  m_pvReFitter = tool<IPVReFitter>("AdaptivePVReFitter", this );
  if(! m_pvReFitter) {
    fatal() << "Unable to retrieve AdaptivePVReFitter" << endreq;
    return StatusCode::FAILURE;
  }
  
  
  //pFinder = tool<IRelatedPVFinder>("GenericParticle2PVRelator__p2PVWithIPChi2_OfflineDistanceCalculatorName_",this);
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  return StatusCode::SUCCESS;

}


//=============================================================================
// Fill the tuple
//=============================================================================
StatusCode TupleToolMuonIso_LC::fill( const LHCb::Particle *top, const LHCb::Particle *part,
                                          const std::string &  	head, Tuples::Tuple &  	tuple	) {

  is_tau=0;
  is_muon=0;
  const LHCb::MCParticle* mclink = 0;
  mclink=m_p2mcAssoc->relatedMCP(part);
  if(mclink) {is_tau=is_tau=abs(mclink->particleID().pid()); }
  else is_tau=0;
  // if(abs(part->particleID().pid())==15) is_tau=1;;



  const std::string prefix=fullName(head);
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Fill" << endmsg;
  

  
  LHCb::Particle::ConstVector decprod;
  

  const SmartRefVector< LHCb::Particle > Bdaughters = part->daughters();
  int countDaughters =0;
  int count_protons=0;
  int count_muon=0;
  
    
  for( SmartRefVector< LHCb::Particle >::const_iterator idau = Bdaughters.begin() ; idau != Bdaughters.end() ; ++idau){
    if (  abs((*idau)->particleID().pid()) == 13 || abs((*idau)->particleID().pid())==2212 ) countDaughters++;
    if (  abs((*idau)->particleID().pid()) == 13 ) count_muon++;
    if (  abs((*idau)->particleID().pid()) == 2212 ) count_protons++;
    
    const LHCb::MCParticle* mclink1 = 0;
    const LHCb::Particle* part = *(idau);
    mclink1=m_p2mcAssoc->relatedMCP(part);
    if(mclink1)
      {	
	std::cout<<"Muon ID "<< mclink1->particleID().pid()<<std::endl;
	if( abs(mclink1->particleID().pid()) ==13)  is_muon+=1; 

      }
  }
 
 if  (countDaughters ==3 && abs(part->particleID().pid()) ==4122 && count_protons==1 && count_muon==2 )
   //  if(1)
   {
     StatusCode fillIsolation_mu = fillIsolation(part, prefix, tuple);
     if (!fillIsolation_mu) return fillIsolation_mu;
     
     StatusCode scfillIso= CDFIsolation( part, decprod, prefix, tuple); 
     if (!scfillIso) return scfillIso;
     return StatusCode::SUCCESS;
   }
 
  else return StatusCode::SUCCESS;
  
}


//-------------------------------------------------------------------------
// isolation variables
//-------------------------------------------------------------------------
StatusCode TupleToolMuonIso_LC::fillIsolation(const LHCb::Particle *part,
						 const std::string prefix,
					   Tuples::Tuple& tuple ){ // part here is B
  
  m_count_mum = 0;
  m_count_mup = 0;
  m_count_p=0;


  const LHCb::VertexBase* goodPV =  m_dva->bestVertex(part); //getRelatedPV(part);
  Gaudi::XYZVector PV_pos;
    debug() << " will fill isolation info  "  << endmsg ;

  if (NULL==goodPV) {
    warning() << "No best PV for " << part->particleID().pid() << endmsg ;
    return StatusCode::SUCCESS;
  }
  
  PV_pos = goodPV->position();
  
  const LHCb::Particle::ConstVector pv = part->daughtersVector();
  
  int idx = 0;
  for (LHCb::Particle::ConstVector::const_iterator ipart=pv.begin();ipart!=pv.end();++ipart) 
    {
      debug() << "daughter ID " << (*ipart)->particleID().pid() << endmsg ;
      //      if ( NULL==(*ipart)->proto() ) continue;
      idx++;
    }
  
  if(idx != 3 )  {
    return StatusCode::SUCCESS;
  }
  
  bool test =true; 
  std::string prefix1, prefix2, prefix3; 

  //-------------------  Giampi's isolation: -------------------//
  
  
  debug()<<"  just before calling getIso "<<endreq;
  

  std::vector<int> iso5 = getIso( part);  // part =B
  debug()<<"  just after calling getIso "
	 <<"GIAMPI  iso "<<iso5[0]<<"  "<<iso5[1]<<endreq;
  
  if (iso5[0]==-999) return StatusCode::FAILURE;
  
  const LHCb::Particle* dau1 = part->daughtersVector().at(0);
  const LHCb::Particle* dau2 = part->daughtersVector().at(1);
  const LHCb::Particle* dau3 = part->daughtersVector().at(2);

  
  //  const LHCb::Particle::ConstVector parts = jpsi->daughtersVector();  
  prefix1 = "proton";
  prefix2 = "mu_1";
  prefix3 = "mu_2";
  ///////////////////////////////////////////////////////////////////////
  m_count_mu1 = static_cast<float>( iso5[0] );                                                                               
  m_count_mu2 = static_cast<float>( iso5[1] );           
  m_count_mu3 = static_cast<float>( iso5[2] );
  

  

  /*
  if (jpsi->daughtersVector().at(0)->charge() < 0) {
    
    m_count_mum_Giampi = static_cast<float>( iso5[0] );
    m_count_mup_Giampi = static_cast<float>( iso5[1] );

  }
  else{
   
    m_count_mum_Giampi = static_cast<float>( iso5[1] );
    m_count_mup_Giampi = static_cast<float>( iso5[0] );

  }
  */
  debug()<<"filling tuple with Giampi's isolation vars "<<  m_count_mup_Giampi<<"  "<< m_count_mum_Giampi <<endmsg;
  
  test &= tuple->column( prefix1+"_isolation", m_count_mu1);
  test &= tuple->column( prefix2+"_isolation", m_count_mu2);
  test &= tuple->column( prefix3+"_isolation", m_count_mu3);

  //-------------------  Giampi's isolation ends -------------------//
  
  if(test) return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
  
}//fillIsolation(part)

//=============================================================================
// method to calculate Giampi's isolation
//=============================================================================
std::vector<int>  TupleToolMuonIso_LC::getIso(const LHCb::Particle* tau){
  
  const LHCb::VertexBase *PV = m_dva->bestVertex(tau);  
  const LHCb::VertexBase *SV = tau->endVertex();
  //  const LHCb::Particle::ConstVector parts = B->daughtersVector();  
  
  const LHCb::MCParticle* mclink1 = 0;
  mclink1=m_p2mcAssoc->relatedMCP(tau);
  //  std::cout<<mclink1<<std::endl;  
  std::vector<int> iso(3, -9999);

  const LHCb::Particle* dau1 = tau->daughtersVector().at(0);
  const LHCb::Particle* dau2 = tau->daughtersVector().at(1);
  const LHCb::Particle* dau3 = tau->daughtersVector().at(2);
  /////////////////////////////////////////////////////////////////////////////////
    
  
  LHCb::Particle::ConstVector::const_iterator i3mu_part;
  ROOT::Math::SMatrix<double, 3, 3> p_3mu;
  ROOT::Math::SMatrix<double, 3, 3> o_3mu;
  ROOT::Math::SVector<double, 3> pt3mu;

  const LHCb::Particle::ConstVector parts = tau->daughtersVector();
  
  
    
  //Loop over kstar particles, get their simple kinematics
  int itrack=0;

  int j=0;
  bool wasmuon=false;
  
  for ( i3mu_part = parts.begin(); i3mu_part != parts.end(); i3mu_part++) {
    const LHCb::ProtoParticle * proto =  (*i3mu_part)->proto();
    const LHCb::Track* track = proto->track();
    if( abs((*i3mu_part)->particleID().pid() ) ==2212 ) j==0; // proton
    if( (abs((*i3mu_part)->particleID().pid()) ==13)&&(!wasmuon) ) {j==1; wasmuon=true;  } // mu-
    if( (abs((*i3mu_part)->particleID().pid()) ==13)&&(wasmuon)) j==2; // mu+  
    
    p_3mu(0,j) = track->momentum().x();
    p_3mu(1,j) = track->momentum().y();
    p_3mu(2,j) = track->momentum().z();
    o_3mu(0,j) = track->position().x();
    o_3mu(1,j) = track->position().y();
    o_3mu(2,j) = track->position().z();

    pt3mu[j] = sqrt(pow(p_3mu(0,j),2)+pow(p_3mu(1,j),2));
    j++;
  }
  
  if(parts.size() !=3) {
    return iso;
  }
  
  StatusCode sc = StatusCode::SUCCESS;
  
  LHCb::Particles* allparts = get<LHCb::Particles>(m_ParticlePath);
  //LHCb::Particles* allparts = get<LHCb::Particles>(m_TracksPath);
  if (!allparts) {
    error() << " Failed to get particles container "
            <<  m_ParticlePath << endmsg;
    
    return iso;
  }
  
  
  Gaudi::XYZPoint PosPV = PV->position();
  Gaudi::XYZPoint PosSV = SV->position();
  
  int i = 0;
    ROOT::Math::SVector<int, 3> iso5;
    iso5[0]=0; // proton
    iso5[1]=0; // mu- 
    iso5[2]=0; // mu+
  
  
  



  double fc_min=10000000000000;
  double angle_min=10000000000000;
  double pvdis_min=10000000000000;
  double svdis_min=10000000000000;
  double doca_min=10000000000000;
  double ips_min=10000000000000;
  double pt_min=10000000000000;
  double p_min=10000000000000;
  ///////////////////////////////////////////////////////
  double fc_max=-10000000000000;
  double angle_max=-10000000000000;
  double pvdis_max=-10000000000000;
  double svdis_max=-10000000000000;
  double doca_max=-10000000000000;
  double ips_max=-10000000000000;
  double pt_max=-10000000000000;
  double p_max=-10000000000000;




for ( int i2 = 0; i2 < 3; i2++) {
    bool hltgood = false;
    double fc = 0.;
    Gaudi::XYZPoint omu(o_3mu(0,i2),o_3mu(1,i2),o_3mu(2,i2));
    Gaudi::XYZVector pmu(p_3mu(0,i2),p_3mu(1,i2),p_3mu(2,i2));
    
    
    //Loop over all particles  
    LHCb::Particles::const_iterator ip;
    int j=0;  // for counting ntracks
    for ( ip = allparts->begin(); ip != allparts->end() ; ++ip) {
      
      
      const LHCb::ProtoParticle * proto =  (*ip)->proto();
      const LHCb::Track* track = proto->track();
      const LHCb::Particle*  cand = (*ip);
      Gaudi::XYZPoint o(track->position());
      Gaudi::XYZVector p(track->momentum());
      bool isInList = 0;
      double pt = p.Rho();
      double ptot = p.R();
      int charge= proto->charge();
      
      float clone = proto->track()->info(LHCb::Track::CloneDist,100000);
      float ghost = proto->track()->ghostProbability();
      
      
      if (track->type()!=m_tracktype)   continue; 
      j++; 
      
      // skip if other particle is in input list 
      if (ratio(pt, pt3mu[0]) < 0.0001 || ratio(pt,pt3mu[1]) <0.0001 || ratio(pt,pt3mu[2]) <0.0001 ) {  continue; }
      
      IsHltGood(o, p, omu, pmu ,PosPV, hltgood, fc); 
      
      // find doca and angle between input and other tracks
      Gaudi::XYZPoint vtx(0.,0.,0.);
      double doca(-1.);          
      double angle(-1.);
      InCone(omu,pmu, o, p, vtx, doca, angle);
      
      // find impact parameters, distances from secondary and primary vertex
      double imp = 0.;
      double impchi2 = 0.;
      double ips,pvdis,svdis, IP;
      ips=100000.;
      LHCb::RecVertex::Container::const_iterator iv;
      LHCb::RecVertex::Container* verts = NULL;
      
      if(exist<LHCb::RecVertex::Container>(m_PVInputLocation))
	{
	  verts = get<LHCb::RecVertex::Container>(m_PVInputLocation);
	}
      for ( iv = verts->begin(); iv != verts->end(); iv++) {
	m_dist->distance(&(*cand),(*iv),imp,impchi2);
	if (impchi2<ips) ips = impchi2;
      }      
      ips=sqrt(ips);
      IP=imp;
      //      double trkchi2 = 
      double deltaZvtx = (vtx.z()-PosPV.z());
      double trkchi2 = cand->proto()->track()->chi2PerDoF();
      pvdis = (vtx.z()-PosPV.z())/fabs(vtx.z()-PosPV.z())*(vtx-PosPV).R();
      svdis = (vtx.z()-PosSV.z())/fabs(vtx.z()-PosSV.z())*(vtx-PosSV).R();
      double distPV = sqrt(   (vtx.z()-PosPV.z())*(vtx.z()-PosPV.z()) + (vtx.y()-PosPV.y())*(vtx.y()-PosPV.y()) + (vtx.x()-PosPV.x())*(vtx.x()-PosPV.x()) );
      double distSV = sqrt(   (vtx.z()-PosSV.z())*(vtx.z()-PosSV.z()) + (vtx.y()-PosSV.y())*(vtx.y()-PosSV.y()) +(vtx.x()-PosSV.x())*(vtx.x()-PosSV.x()) );
      
      
      

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //double ipChi2 = track->
     

      // non-isolating criteria #5
      if (angle <m_angle && fc<m_fc && doca<m_doca_iso && ips>m_ips && 
          svdis>m_svdis && svdis<m_svdis_h && pvdis>m_pvdis && pvdis<m_pvdis_h 
          && track->type()==m_tracktype) {
	iso5[i2] += 1;
      }// if 
      if(hltgood)
	{
	  /*
	  fc_h.Fill(fc);
	  angle_h.Fill(angle);
	  pvdis_h.Fill(pvdis);
	  svdis_h.Fill(svdis);
	  doca_h.Fill(doca);
	  ips_h.Fill(ips);
	  pt_h.Fill(pt);
	  p_h.Fill(ptot);
	  if(fc<fc_min) fc_min=fc;
	  if(angle<angle_min) angle_min=angle;
	  if(pvdis< pvdis_min)  pvdis_min= pvdis;
          if(svdis< svdis_min)  svdis_min= svdis;
	  if(doca< doca_min)  doca_min= doca;
	  if(ips< ips_min)  ips_min= ips;
	  if(pt< pt_min)  pt_min= pt;
	  if(ptot< p_min)  p_min= ptot;
	  ////////////////////////////////////////////////
	  if(fc>fc_max) fc_max=fc;
          if(angle>angle_max) angle_max=angle;
          if(pvdis> pvdis_max)  pvdis_max= pvdis;
          if(svdis> svdis_max)  svdis_max= svdis;
          if(doca> doca_max)  doca_max= doca;
          if(ips> ips_max)  ips_max= ips;
          if(pt> pt_max)  pt_max= pt;
          if(ptot> p_max)  p_max= ptot;
	  */
	}



      /*


      Tuple tuple = nTuple( "IsolationTwoBody") ;
      tuple->column("is_tau", is_tau);
      tuple->column("is_muon", is_muon);
      tuple -> column ( "muon_number" , i2) ;
      tuple -> column ( "track" , itrack) ;
      tuple -> column ( "angle" , angle) ;
      tuple -> column ( "fc" , fc) ;
      tuple -> column ( "trckchi2" , trkchi2) ;
      tuple -> column ( "doca" , doca) ;
      tuple -> column ( "ip" , IP) ;
      tuple -> column ( "ips" , ips) ;
      tuple -> column ( "svdis" , svdis) ;
      tuple -> column ( "pvdis" , pvdis) ;
      tuple -> column ( "ishltgood" , hltgood ) ;
      tuple -> column ( "clone" , clone ) ;
      tuple -> column ( "ghost" , ghost ) ;
      tuple -> column ( "pttot" , pt3mu[0]+pt3mu[1]+pt3mu[2] ) ;
      tuple -> column ( "ptot" , ptot ) ;
      tuple -> column ( "pt" , pt ) ;
      tuple -> column ( "distPV" , distPV ) ;
      tuple -> column ( "distSV" , distSV ) ;


      tuple->write();
      */
      
      itrack++;
    }// allparts
    
  } //i2
/*
  Tuple tuple_mum = nTuple( "Mu_minus") ;
  


  tuple_mum->column("is_tau", is_tau);
  tuple_mum->column("is_muon", is_muon);



  tuple_mum->column("fc_mean", fc_h.GetMean());
  tuple_mum->column("fc_RMS", fc_h.GetRMS());
  tuple_mum->column("fc_min", fc_min);
  tuple_mum->column("fc_max", fc_max);

  tuple_mum->column("angle_mean", angle_h.GetMean());
  tuple_mum->column("angle_RMS", angle_h.GetRMS());
  tuple_mum->column("angle_min", angle_min);
  tuple_mum->column("angle_max", angle_max);

  tuple_mum->column("pvdis_mean", pvdis_h.GetMean());
  tuple_mum->column("pvdis_RMS", pvdis_h.GetRMS());
  tuple_mum->column("pvdis_min", pvdis_min);
  tuple_mum->column("pvdis_max", pvdis_max);

  tuple_mum->column("svdis_mean", svdis_h.GetMean());
  tuple_mum->column("svdis_RMS", svdis_h.GetRMS());
  tuple_mum->column("svdis_min", svdis_min);
  tuple_mum->column("svdis_max", svdis_max);

  tuple_mum->column("doca_mean", doca_h.GetMean());
  tuple_mum->column("doca_RMS", doca_h.GetRMS());
  tuple_mum->column("doca_min", doca_min);
  tuple_mum->column("doca_max", doca_max);

  tuple_mum->column("ips_mean", ips_h.GetMean());
  tuple_mum->column("ips_RMS", ips_h.GetRMS());
  tuple_mum->column("ips_min", ips_min);
  tuple_mum->column("ips_max", ips_max);

  tuple_mum->column("pt_mean", pt_h.GetMean());
  tuple_mum->column("pt_RMS", pt_h.GetRMS());
  tuple_mum->column("pt_min", pt_min);
  tuple_mum->column("pt_max", pt_max);

  tuple_mum->write();
*/  
//////////////////////////////////////////////////////////////


 

  
  

  
    iso[0] = iso5[0] ;
    iso[1] = iso5[1] ;  
    iso[2] = iso5[2] ;

  
  return iso;

}


//=============================================================================
// IsHLTGood method,used by isolation calculation
//=============================================================================
void  TupleToolMuonIso_LC::IsHltGood(Gaudi::XYZPoint o,Gaudi::XYZVector p,
					 Gaudi::XYZPoint o_mu,Gaudi::XYZVector 
					 p_mu, Gaudi::XYZPoint PV, bool& hltgood,
					 double& fc) {
  
  Gaudi::XYZVector rv;
  Gaudi::XYZPoint vtx;
  Gaudi::XYZPoint close;
  Gaudi::XYZPoint close_mu;
  bool fail(false);

  closest_point(o,p,o_mu,p_mu,close,close_mu,vtx,fail);

  if (fail) {
    fc = -1.;
    hltgood = -1;
  }
  else {
    double pete = p.Rho();
    rv = vtx - PV;
    double DOCA_b = (close-close_mu).R();
    double ZZ = rv.z();
    fc = pointer(rv,p,p_mu);
    hltgood=( (DOCA_b<0.2) && (ZZ>0.) && (ZZ<30.) && (fc<0.4) && (pete>2.) );
  }
}


double TupleToolMuonIso_LC::pointer (Gaudi::XYZVector vertex,
                                 Gaudi::XYZVector p, Gaudi::XYZVector
                                 p_mu)  {
  double pt=p.Rho()+p_mu.Rho();
  Gaudi::XYZVector ptot(p+p_mu);
  double temp = arcosine(vertex,ptot);
  double  num=ptot.R()*sin(temp);
  double  den=(num+pt);
  double fc = num/den;
  return fc;
}






//=============================================================================
// Other functions needed by isolation 
//=============================================================================

double TupleToolMuonIso_LC::getphi(const LHCb::Particle* vdau1, const LHCb::Particle* vdau2){
  double dphi = vdau1->momentum().Phi() - vdau2->momentum().Phi();
  return dphi;
}


//=============================================================================
double TupleToolMuonIso_LC::gettheta(const LHCb::Particle* vdau1, const LHCb::Particle* vdau2){
  
  double dtheta = vdau1->momentum().Eta() -  vdau2->momentum().Eta();
  return dtheta;
}


//=============================================================================
double TupleToolMuonIso_LC::ratio( double p1, double p2){  
  return fabs(p1 -p2)*(1./fabs(p1+p2)); 
}



//=============================================================================
double TupleToolMuonIso_LC::IsClose(const LHCb::Particle* p1,const LHCb::Particle* p2) {
  
  double deta = gettheta(p1,p2);
  double dphi = getphi(p1,p2);
  return sqrt(deta*deta+dphi*dphi);
}


//=============================================================================
void TupleToolMuonIso_LC::closest_point(Gaudi::XYZPoint o,Gaudi::XYZVector p,
					    Gaudi::XYZPoint o_mu,Gaudi::XYZVector p_mu, 
					    Gaudi::XYZPoint& close1, 
					    Gaudi::XYZPoint& close2, 
					    Gaudi::XYZPoint& vertex, bool& fail) {
  
  
  Gaudi::XYZVector v0(o - o_mu);
  Gaudi::XYZVector v1(p.unit());
  Gaudi::XYZVector v2(p_mu.unit());
  Gaudi::XYZPoint temp1(0.,0.,0.);
  Gaudi::XYZPoint temp2(0.,0.,0.);
  fail = false;
  double  d02 = v0.Dot(v2); 
  double  d21 = v2.Dot(v1); 
  double  d01 = v0.Dot(v1);
  double  d22 = v2.Dot(v2);
  double  d11 = v1.Dot(v1); 
  double  denom = d11 * d22 - d21 * d21; 
  if (fabs(denom) <= 0.) {
    close1 = temp1;
    close2 = temp2;
    fail = true;
  }
  else {
    double numer = d02 * d21 - d01 * d22; 
    double mu1 = numer / denom;            
    double mu2 = (d02 + d21 * mu1) / d22; 
    close1 = o+v1*mu1;
    close2 = o_mu+v2*mu2;
  }
  vertex = (close1+(close2-close1)*0.5);
}


double TupleToolMuonIso_LC::arcosine(Gaudi::XYZVector p1,Gaudi::XYZVector p2) {
  
  double num    = (p1.Cross(p2)).R();
  double den    = p1.R()*p2.R();
  double seno   = num/den;
  double coseno = p1.Dot(p2)/den;
  double alpha  = asin(fabs(seno));
  if (coseno < 0 ) {
    alpha = ROOT::Math::Pi() - alpha;
  }
  return alpha;
}



//=============================================================================
void TupleToolMuonIso_LC::InCone(Gaudi::XYZPoint o1, 
				     Gaudi::XYZVector p1,Gaudi::XYZPoint o2, 
				     Gaudi::XYZVector p2,
				     Gaudi::XYZPoint& vtx, double& 
				     doca, double& angle){
  
  Gaudi::XYZPoint rv;
  Gaudi::XYZPoint close;
  Gaudi::XYZPoint close_mu;
  bool fail(false);
  closest_point(o1,p1,o2,p2,close,close_mu,vtx, fail);
  if (fail) {
    doca =-1.;
    angle=-1.;
  } 
  else {  
    doca = (close-close_mu).R();
    angle = arcosine(p1,p2);
  }
}

//=============================================================================
StatusCode TupleToolMuonIso_LC::CDFIsolation(const LHCb::Particle* B,
						 const LHCb::Particle::ConstVector  pions,
						 const std::string prefix, 
						 Tuples::Tuple& tuple ){
  
  bool test = true;
  
  LHCb::Particles*  parts = get<LHCb::Particles>(m_ParticlePath);
  if (!parts) {
    error() << " Failed to get particles container " << endmsg;
    return StatusCode::SUCCESS;
  }

  
  double pt_bs = B->momentum().rho(); //pt();
  double iso_giampi = 0.0;
     
  for(LHCb::Particles::const_iterator ipp=parts->begin();ipp!=parts->end();ipp++){
    const LHCb::ProtoParticle *proto = (*ipp)->proto();
    if(proto) {
      const LHCb::Track* atrack = proto->track();
      if(atrack) {
        
        Gaudi::XYZVector ptrack =  ((*ipp)->proto()->track()->momentum()); 
        //double ptrack = atrack->p(); //
        double pttrack = ptrack.rho(); //   (*ipp)->pt();  //  
        int overlap=0;
        for ( LHCb::Particle::ConstVector::const_iterator ipion = pions.begin(); ipion != pions.end() ; ++ipion) {
          double ptmu =  (*ipion)->momentum().rho();
          if  (ratio(pttrack, ptmu) < 0.0001) {overlap=1; break;}
        }
        if (overlap==1) continue;
        
        
        double deta      =  B->momentum().Eta() - (*ipp)->momentum().Eta() ; 
        double delta_phi =  B->momentum().Phi() - (*ipp)->momentum().Phi(); 
        delta_phi = TMath::Abs(delta_phi);
        
        if (delta_phi > TMath::Pi() )  delta_phi = 2*TMath::Pi() - delta_phi;   
        
        double rad_cone = TMath::Sqrt(  (delta_phi*delta_phi + deta*deta) );
        if (  rad_cone <=1.0)
          {
            iso_giampi += pttrack;
          }
      } // atrack
    } //proto
  } // ipp
  
  iso_giampi = pt_bs/(iso_giampi+pt_bs);
  test &= tuple->column( prefix+"_CDFiso", iso_giampi);
  
  return test;
  
}




//==============================
// MC Ancestor Finding Routine:
//==============================
const LHCb::MCParticle* TupleToolMuonIso_LC::originof( const LHCb::MCParticle* product ) {
  const LHCb::MCParticle* mother = product->mother();
  if ( (!mother) || product->particleID().hasBottom() ) return product;
  else return originof( mother );
}

  


//=============================================================================
// Finalize
//=============================================================================
StatusCode TupleToolMuonIso_LC::finalize(){
  
  return TupleToolBase::finalize();
  

}
////////////////////////////////////////////////////////////////////////////////////
/*
const LHCb::RecVertex* TupleToolMuonIso_LC::closestPV(const LHCb::Particle* part)
{
  const LHCb::RecVertex* bestPV = 0;
  double smallest = 999999.0;
  std::vector<double> chi2_map;

  LHCb::RecVertices* PVs =  get<LHCb::RecVertices>(LHCb::RecVertexLocation::Primary);
  for(LHCb::RecVertices::const_iterator ipv = PVs->begin(); ipv != PVs->end(); ipv++ ) {
    const LHCb::RecVertex* tmppv = *ipv;
    double tmpip = 9999.0;
    double tmpipchi2 = 1000.0;
    StatusCode sc = m_dist()->distance(part,tmppv,tmpip,tmpipchi2);
    if(sc.isSuccess()) {
      chi2_map.push_back(tmpipchi2);
      if(tmpipchi2 < smallest) {
        bestPV = tmppv;
        smallest = tmpipchi2;
      }
    }

  }
  // sort vector                                                                                                                    
  std::sort (chi2_map.begin(), chi2_map.end(), myobject);
  for(int i=0;i<(int)chi2_map.size();++i)
    {
      if(i==0) continue;
      if(chi2_map[i] < m_PV_separation) return NULL;

    }

  return bestPV;
}
*/
