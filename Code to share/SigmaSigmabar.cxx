// This is v0.4.
// All functional. Improvements of efficiency sought.
// 
// Benjamin Verbeek, Hefei, 2022-11-29

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"



#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "SigmaSigmabarAlg/SigmaSigmabar.h"

//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/SecondVertexFit.h"	// added
#include "ParticleID/ParticleID.h"
#include "McTruth/McParticle.h"

#include <vector>

const double mpi = 0.139570;  // pion mass [GeV]
const double mproton  = 0.938272;  // Proton mass, from PDG 2022 [GeV]
const double mlambda  = 1.115683; // Lambda mass, PDG 2022
const double mSigma0  = 1.192642; // Sigma0 mass, PDG 2022

const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};

const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTP;

int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut6,Ncut30,Ncut31;  // Initializing cut-counters. OK. /BV
int Npi, Nproton; // Initializing counters for pions and protons

/////////////////////////////////////////////////////////////////////////////

SigmaSigmabar::SigmaSigmabar(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  
  //Declare the properties  
  declareProperty("Vr0cut", m_vr0cut=10.0);  // Cylinder cut condition
  declareProperty("Vz0cut", m_vz0cut=20.0); // From memo is 20.0 /BV 221129.

  declareProperty("EnergyThreshold", m_energyThreshold=0.04); // Fake gamma cut
  declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
  declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
  declareProperty("GammaAngleCut", m_gammaAngleCut=20.0);

  // Flag: use or not (1/0)
  declareProperty("Test4C", m_test4C = 1);
  declareProperty("CheckDedx", m_checkDedx = 1);
  declareProperty("CheckTof",  m_checkTof = 1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode SigmaSigmabar::initialize(){
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  
  //
  //  BOOK ALL OUTPUT NTUPLES
  //
  StatusCode status;

  if(m_test4C==1) {
    NTuplePtr nt4(ntupleSvc(), "FILE1/fit4c");
    if ( nt4 ) m_tuple4 = nt4;
    else {
      m_tuple4 = ntupleSvc()->book ("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple4 )    {
        //OLD
	status = m_tuple4->addItem ("chi2",         m_chi1);
  // status = m_tuple4->addItem ("mLambda",      m_mLambda);
  // status = m_tuple4->addItem ("mLambdabar",   m_mLambdabar);
  // status = m_tuple4->addItem ("pLambda",      m_pLambda);
  // status = m_tuple4->addItem ("pLambdabar",   m_pLambdabar);
  // status = m_tuple4->addItem ("pg1",          m_pg1);
  // status = m_tuple4->addItem ("pg2",          m_pg2);
  // status = m_tuple4->addItem ("m_mSigma0",    m_mSigma0);
  // status = m_tuple4->addItem ("m_mSigmabar0", m_mSigmabar0);
  // status = m_tuple4->addItem ("m_pSigma0",    m_pSigma0);
  // status = m_tuple4->addItem ("m_pSigmabar0", m_pSigmabar0);

////
  // Masses -- use .m() on momentum
  //status = m_tuple4->addItem ("m_m_p",         m_m_p);
  //status = m_tuple4->addItem ("m_m_pbar",      m_m_pbar);
  // status = m_tuple4->addItem ("m_m_pip",       m_m_pip);
  // status = m_tuple4->addItem ("m_m_pim",       m_m_pim);
  status = m_tuple4->addItem ("m_m_Lambda",    m_m_Lambda);
  status = m_tuple4->addItem ("m_m_Lambdabar", m_m_Lambdabar);
  status = m_tuple4->addItem ("m_m_Sigma0",    m_m_Sigma0);
  status = m_tuple4->addItem ("m_m_Sigmabar0", m_m_Sigmabar0);
  status = m_tuple4->addItem ("m_m_gg",        m_m_gg);

  // Full momentum -- use .rho()
  status = m_tuple4->addItem ("m_p_p",         m_p_p);
  status = m_tuple4->addItem ("m_p_pbar",      m_p_pbar);
  status = m_tuple4->addItem ("m_p_pip",       m_p_pip);
  status = m_tuple4->addItem ("m_p_pim",       m_p_pim);
  status = m_tuple4->addItem ("m_p_Lambda",    m_p_Lambda);
  status = m_tuple4->addItem ("m_p_Lambdabar", m_p_Lambdabar);
  status = m_tuple4->addItem ("m_p_Sigma0",    m_p_Sigma0);
  status = m_tuple4->addItem ("m_p_Sigmabar0", m_p_Sigmabar0);
  status = m_tuple4->addItem ("m_p_g1",        m_p_g1);
  status = m_tuple4->addItem ("m_p_g2",        m_p_g2);
  status = m_tuple4->addItem ("m_p_gg",        m_p_gg);

  // momentum x-direction -- use .px()
  status = m_tuple4->addItem ("m_px_p",         m_px_p);
  status = m_tuple4->addItem ("m_px_pbar",      m_px_pbar);
  status = m_tuple4->addItem ("m_px_pip",       m_px_pip);
  status = m_tuple4->addItem ("m_px_pim",       m_px_pim);
  status = m_tuple4->addItem ("m_px_Lambda",    m_px_Lambda);
  status = m_tuple4->addItem ("m_px_Lambdabar", m_px_Lambdabar);
  status = m_tuple4->addItem ("m_px_Sigma0",    m_px_Sigma0);
  status = m_tuple4->addItem ("m_px_Sigmabar0", m_px_Sigmabar0);
  status = m_tuple4->addItem ("m_px_g1",        m_px_g1);
  status = m_tuple4->addItem ("m_px_g2",        m_px_g2);

  // momentum y-direction -- use .py()
  status = m_tuple4->addItem ("m_py_p",         m_py_p);
  status = m_tuple4->addItem ("m_py_pbar",      m_py_pbar);
  status = m_tuple4->addItem ("m_py_pip",       m_py_pip);
  status = m_tuple4->addItem ("m_py_pim",       m_py_pim);
  status = m_tuple4->addItem ("m_py_Lambda",    m_py_Lambda);
  status = m_tuple4->addItem ("m_py_Lambdabar", m_py_Lambdabar);
  status = m_tuple4->addItem ("m_py_Sigma0",    m_py_Sigma0);
  status = m_tuple4->addItem ("m_py_Sigmabar0", m_py_Sigmabar0);
  status = m_tuple4->addItem ("m_py_g1",        m_py_g1);
  status = m_tuple4->addItem ("m_py_g2",        m_py_g2);
  
  // momentum z-direction -- use .pz()
  status = m_tuple4->addItem ("m_pz_p",         m_pz_p);
  status = m_tuple4->addItem ("m_pz_pbar",      m_pz_pbar);
  status = m_tuple4->addItem ("m_pz_pip",       m_pz_pip);
  status = m_tuple4->addItem ("m_pz_pim",       m_pz_pim);
  status = m_tuple4->addItem ("m_pz_Lambda",    m_pz_Lambda);
  status = m_tuple4->addItem ("m_pz_Lambdabar", m_pz_Lambdabar);
  status = m_tuple4->addItem ("m_pz_Sigma0",    m_pz_Sigma0);
  status = m_tuple4->addItem ("m_pz_Sigmabar0", m_pz_Sigmabar0);
  status = m_tuple4->addItem ("m_pz_g1",        m_pz_g1);
  status = m_tuple4->addItem ("m_pz_g2",        m_pz_g2);

  // Energy -- use .e() ?
  status = m_tuple4->addItem ("m_E_p",         m_E_p);
  status = m_tuple4->addItem ("m_E_pbar",      m_E_pbar);
  status = m_tuple4->addItem ("m_E_pip",       m_E_pip);
  status = m_tuple4->addItem ("m_E_pim",       m_E_pim);
  status = m_tuple4->addItem ("m_E_Lambda",    m_E_Lambda);
  status = m_tuple4->addItem ("m_E_Lambdabar", m_E_Lambdabar);
  status = m_tuple4->addItem ("m_E_Sigma0",    m_E_Sigma0);
  status = m_tuple4->addItem ("m_E_Sigmabar0", m_E_Sigmabar0);
  status = m_tuple4->addItem ("m_E_g1",        m_E_g1);
  status = m_tuple4->addItem ("m_E_g2",        m_E_g2);
  status = m_tuple4->addItem ("m_E_gg",        m_E_gg);


  // Also save all after 6C.

  // Save after 2nd vx-fit
////

  // for topoana
  status = m_tuple4->addItem("indexmc",       m_idxmc, 0, 100);
  status = m_tuple4->addItem("pdgid",         m_idxmc, m_pdgid);
  status = m_tuple4->addItem("motheridx",     m_idxmc, m_motheridx);

      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // test 4C


// Comment in head-file and all assignments. Add new output to tuple! 

  //
  //--------end of book--------
  //

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * *  END INITIALIZE * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * BEGIN EXECUTE * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

StatusCode SigmaSigmabar::execute() {
  
  // Messages
  std::cout << "execute()" << std::endl;


  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  int event=eventHeader->eventNumber();
  log << MSG::DEBUG <<"run, evtnum = "
      << runNo << " , "
      << event <<endreq;
  cout<<"event "<<event<<endl;
  Ncut0++;  // Counts total number of events

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  //  log << MSG::INFO << "get event tag OK" << endreq;
    log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
      << evtRecEvent->totalCharged() << " , "
      << evtRecEvent->totalNeutral() << " , "
      << evtRecEvent->totalTracks() <<endreq;


// NOT WORKING... What are the particles?
 
  if( runNo < 0){ 
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
    if( mcParticleCol ){
      int m_numParticle = 0;
      bool psipDecay = false;
      int rootIndex = -1;
      int nmc_lambda=0, nmc_lambdabar=0, nmc_pp=0, nmc_nbar=0, nmc_pi01=0, nmc_pi02 = 0, nmc_pm = 0, nmc_n = 0;
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++) {
        if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        if ((*iter_mc)->particleProperty()==443 ) {
                        psipDecay = true;
                        rootIndex = (*iter_mc)->trackIndex();
        }
        if (!psipDecay) continue;
        int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
        int pdgid = (*iter_mc)->particleProperty();
        m_pdgid[m_numParticle] = pdgid;
        m_motheridx[m_numParticle] = mcidx;
        m_numParticle +=1;
      }
      m_idxmc = m_numParticle;
    }
  }
  

  // Ncut_mc++; // what is this? Initialize or ignore?
  //
  //  ######### ########## ########## ########## ########## ########## #########
  // ######## START FIRST CUT: NUMBER OF CHARGED TRACKS AND TOTAL CHARGE + REGION OF ORGIN ########
  // ########## Here, we look for charge conservation (sum 0), 4 charged tracks and only tracks
  // ########## from near the interaction point
  //  ######### ########## ########## ########## ########## ########## #########
  //

  // Some initializations:
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  //
  // check x0, y0, z0, r0
  // suggest cut: |z0|<5 && r0<1
  //
  Vint iGood, ipip, ipim, ip, ipbar; // CHNGD: add also for protons: ip, ipbar
  iGood.clear();
  ipip.clear();
  ipim.clear();
  ip.clear();
  ipbar.clear();
  Vp4 ppip, ppim, pp, ppbar; // CHNGD: add also for protons. pp, ppbar
  ppip.clear();
  ppim.clear();
  pp.clear();
  ppbar.clear();
  
  // Start cut 1: nCharge = 4, totCharge = 0
  int nCharge = 0;

  Hep3Vector xorigin(0,0,0);
  HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
  //if (m_reader.isRunNumberValid(runNo)) {
   IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double* dbv = vtxsvc->PrimaryVertex();  // Get primary vertex.
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
  //    HepVector dbv = m_reader.PrimaryVertex(runNo);
  //    HepVector vv = m_reader.SigmaPrimaryVertex(runNo);
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
      vx.setX(dbv[0]);
      vx.setY(dbv[1]);
      vx.setZ(dbv[2]);
      Evx[0][0]=vv[0]*vv[0];
      Evx[1][1]=vv[1]*vv[1];
      Evx[2][2]=vv[2]*vv[2];
  }
  VertexParameter vx_db;
	vx_db.setVx(vx);
	vx_db.setEvx(Evx);

  // cout << "Hello SigmaSigmabar world! v0.3.2. output Sigmamass - from VSCode" << endl; // edited!
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){ // loop over charged tracks
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;  // get track
    if(!(*itTrk)->isMdcTrackValid()) continue;  // check if MDC track is valid
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack(); // get MDC track
    double pch=mdcTrk->p();
    double x0=mdcTrk->x();
    double y0=mdcTrk->y();
    double z0=mdcTrk->z();
    double phi0=mdcTrk->helix(1);
    double xv=xorigin.x();
    double yv=xorigin.y();
    double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
    double cost = cos(mdcTrk->theta());


    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
    VFHelix helixip(point0,a,Ea); 
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
    double  Rvphi0=vecipa[1];       


  // THROW AWAY TRACKS NOT ORIGINATING FROM INTERACTION REGION  
    if(fabs(Rvz0) >= m_vz0cut) continue;  // CHNGD: from 10.0 to m_vz0cut
    if(fabs(Rvxy0) >= 10.0) continue;  // TODO: Here? Angle? // NOTE: Xiaorong claims angle cut is superfluous anyway. /BV 11-30
    if(fabs(cost) >= 0.93 ) continue; 
    
    iGood.push_back(i); // OK, add to list of good charged tracks
    nCharge += mdcTrk->charge();  // add the charge (+1 or -1)
  }
  
  
  //
  // Finish Good Charged Track Selection (we expect pi+, pi-, p, pbar)
  //
  int nGood = iGood.size();
  log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
  if((nGood != 4)||(nCharge!=0)){ // CHNGD: Here, we expect 4 charged tracks with total charge = 0
    return StatusCode::SUCCESS;   // Failed test? Return ok and move to next event
  }

  Ncut1++;  // Programme didn't return? Cool! Add to counter of items that came this far.
  //
  //  ######### ########## ########## ########## ########## ########## #########
  // ######### END FIRST CUT: NUMBER OF CHARGED TRACKS AND TOTAL CHARGE #########
  //  ######### ########## ########## ########## ########## ########## #########
  // ########## START SECOND CUT: GAMMA
  // ########## Here, we look for two or more gammas.
  // ######### ########## ########## ########## ########## ########## #########
  
  // Start gamma finder, nGam >= 2
  Vint iGam;
  iGam.clear();
  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) { // Iterate over unchaarged tracks
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;  // Uncharged are appended after charged.
    if(!(*itTrk)->isEmcShowerValid()) continue; // Black-box. (Keep) If not valid, continue.
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) { // loop over charged tracks
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      //      double ctht = extpos.cosTheta(emcpos);
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if(angd < dang){ // Update if closer
        dang = angd;
        dthe = thed;
        dphi = phid;
      }
    }
    if(dang>=200) continue; // CUT: Angle between charged track and gamma must be less than 200
    double eraw = emcTrk->energy(); // ^ TODO: ok limit?
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);

    if(eraw < m_energyThreshold) continue; // CUT: Energy of gamma must be below threshold
    if(fabs(dang) < m_gammaAngleCut) continue; // CUT: Angle between charged track and gamma must be less than m_gammaAngleCut
    //                                                ^ TODO ok limit?
    // good photon cut will be set here
    //
    iGam.push_back(i);  // Save gamma to list of good gammas
  }
  
  //
  // Finish Good Photon Selection
  //
  int nGam = iGam.size();

  log << MSG::DEBUG << "num Good Photon " << nGam  << " , " <<evtRecEvent->totalNeutral()<<endreq;
  if(nGam<2){
    return StatusCode::SUCCESS;
  }
  Ncut2++;
  //  ######### ########## ########## ########## ########## ########## #########
  // ######### END SCOND CUT: NUMBER OF OK GAMMAS >= 2 #########
  //  ######### ########## ########## ########## ########## ########## #########
 
  //
  // Assign 4-momentum to each photon
  // 

  Vp4 pGam;
  pGam.clear();
  for(int i = 0; i < nGam; i++) { // loop all passed gamma
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i]; 
    RecEmcShower* emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector ptrk;
    ptrk.setPx(eraw*sin(the)*cos(phi));
    ptrk.setPy(eraw*sin(the)*sin(phi));
    ptrk.setPz(eraw*cos(the));
    ptrk.setE(eraw);

//    ptrk = ptrk.boost(-0.011,0,0);// boost to cms

    pGam.push_back(ptrk); // 4-momentum of each photon added to vector
  }
  cout<<"before pid"<<endl;

  //  ######### ########## ########## ########## ########## ########## #########
  // ######### Particle Identification (PID) and momentum assignments #########
  //  ######### ########## ########## ########## ########## ########## #########
//  ######### ########## ########## ########## ########## ########## #########
  // ######### START THIRD CUT: Exact check of decay products (PID) #########
  //  ######### ########## ########## ########## ########## ########## #########

  //
  // Assign 4-momentum to each charged track
  //
  ParticleID *pid = ParticleID::instance(); // Black-box to me but okay.
  for(int i = 0; i < nGood; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];   // iterate good tracks
    //    if(pid) delete pid;
    
    RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
    bool isPion = ((mdcTrk->p()) < 0.45) ? true : false; // Value taken from plot/memo // Seems to work fine. One switches.

    // If change, keep prob-checks here? 
    //if((pid->probPion() < 0.001) && (pid->probProton() < 0.001)) continue; // Edited. If too unlikely, throw it away.
    // removed
    isPion ? Npi++ : Nproton++; // Count pions and protons

    RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
    RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion
// NOTE: Wrong

    // moved out from below if-else
    HepLorentzVector ptrk;
    ptrk.setPx(mdcKalTrk->px());
    ptrk.setPy(mdcKalTrk->py());
    ptrk.setPz(mdcKalTrk->pz());
    double p3 = ptrk.mag();
    // Add this momentum to output... curious is two separate peaks:


    // Characterize!
    if((mdcKalTrk->charge() > 0) && isPion ) { // if positive & pion, its pi+
      ipip.push_back(iGood[i]); // dep
      ptrk.setE(sqrt(p3*p3+mpi*mpi)); // This is why proton mass is needed. TODO: define in intialize
      ppip.push_back(ptrk); // 4-momentum of each pion added to vector

    } else if ((mdcKalTrk->charge() < 0) && isPion) {  // if not positive and is pion, its pi-
      ipim.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mpi*mpi));
      ppim.push_back(ptrk); // 4-momentum of each pion added to vector

      // Added: now check for proton in the same way:
    } else if ((mdcKalTrk->charge() > 0) && !isPion) { // if positive & proton, its p+
      ip.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mproton*mproton));
      pp.push_back(ptrk); // 4-momentum of each proton added to vector

    // and for pbar
    } else if ((mdcKalTrk->charge() < 0) && !isPion) {  // if not positive and is proton, its pbar-
      ipbar.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mproton*mproton)); 
      ppbar.push_back(ptrk); // 4-momentum of each proton added to vector
    }
//      ptrk = ptrk.boost(-0.011,0,0);//boost to cms
  } // characterized tracks and momenta now stored to respective lists.


  int npip = ipip.size();
  int npim = ipim.size(); // Add for proton/-bar. TODO.
  int np = ip.size();
  int npbar = ipbar.size();
  // cout << "Here come npip, npim, np, npbar: " << npip << npim << np << npbar << endl; // edited!

  if(npip*npim*np*npbar != 1) return StatusCode::SUCCESS;  // Here we check if there is only one pion and one anti-pion. TODO: Add proton/-bar.
          // Benefit of writing product is unclear to me. Harder to read and more computationally 
          // expensive, is my guess. TODO.  We could probably just use momentum split.

  Ncut3++;

//  ######### ########## ########## ########## ########## ########## #########
  // ######### END THIRD CUT: Exact check of decay products #########
  //  ######### ########## ########## ########## ########## ########## #########
  // ########## START KINEMATIC FIT. #########
  //######### ########## ########## ########## ########## ########## #########
  // TODO: What now? Secondary vx fit to get pion + proton to virtual lambda?


  // Construct virtual lambda. Check if pi+ and pbar origin from same vertex
  // Check if pi- and p origin from same vertex

  //
  // Loop each gamma pair, check ppi0 and pTot 
  //

  HepLorentzVector pTot;
  for(int i = 0; i < nGam - 1; i++){
    for(int j = i+1; j < nGam; j++) {
      HepLorentzVector p2g = pGam[i] + pGam[j];
      pTot = ppip[0] + ppim[0];
      pTot += p2g;
      // m_m2gg = p2g.m();
      // m_etot = pTot.e();
      // m_tuple3 -> write();
    }
  }
  
  Vp4 p4Lambdavtx, p4pvtx, p4pbarvtx, p4pimvtx, p4pipvtx;
  p4Lambdavtx.clear();
  p4pbarvtx.clear();
  p4pvtx.clear();
  p4pimvtx.clear();
  p4pipvtx.clear();
  //Vdouble  decayL_lambdabar, decayLerr_lambdabar, chisq_lambdabar;
  //decayLerr_lambdabar.clear();
  //decayL_lambdabar.clear();
  //chisq_lambdabar.clear();
  //VWTP wlambdabar_vertex;
  //wlambdabar_vertex.clear();
  Vdouble  decayL_Lambda, decayLerr_Lambda, chisq_Lambda;
  decayLerr_Lambda.clear();
  decayL_Lambda.clear();
  chisq_Lambda.clear();
  VWTP wLambda_vertex;
  wLambda_vertex.clear();
  HepPoint3D cPoint;

  // point to start of pions. TODO: also for protons?
  RecMdcKalTrack *pipTrk  = (*(evtRecTrkCol->begin()+ipip[0]))->mdcKalTrack();
  RecMdcKalTrack *pimTrk  = (*(evtRecTrkCol->begin()+ipim[0]))->mdcKalTrack();
  RecMdcKalTrack *pTrk    = (*(evtRecTrkCol->begin()+ip[0]))->mdcKalTrack();
  RecMdcKalTrack *pbarTrk = (*(evtRecTrkCol->begin()+ipbar[0]))->mdcKalTrack();

  WTrackParameter wvpipTrk, wvpimTrk, wvpTrk, wvpbarTrk; // TODO: also for protons? see below
  wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
  wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());
  wvpTrk   = WTrackParameter(mproton,  pTrk->getZHelixP(),  pTrk->getZErrorP());
  wvpbarTrk= WTrackParameter(mproton, pbarTrk->getZHelixP(), pbarTrk->getZErrorP());


  // Initialize (origin + large error)
  // HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx2(3, 0);
  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx2[0][0] = bx*bx;
  Evx2[1][1] = by*by;
  Evx2[2][2] = bz*bz;

  // Set these initial values to the vertex object
  VertexParameter vxpar;
  vxpar.setVx(vx);

  vxpar.setEvx(Evx2);  // Previously wrong error
  // fit for Lambda from pi- and p
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpTrk);
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);
  //if(!vtxfit->Fit(0)) return SUCCESS; // If the fit fails, I guess? Is it assuming from origin?
  // Primary VERTEX FIT
  if(!vtxfit->Fit(0)) return StatusCode::SUCCESS;	// else, ... (to keep in scope)
    vtxfit->Swim(0);
    vtxfit->BuildVirtualParticle(0);
    WTrackParameter wLambda = vtxfit->wVirtualTrack(0);
    VertexParameter vtxLambda = vtxfit->vpar(0);
    WTrackParameter wtrkproton = vtxfit->wtrk(0);
    p4pvtx.push_back(wtrkproton.p());
    WTrackParameter wtrkpion = vtxfit->wtrk(1);
    p4pimvtx.push_back(wtrkpion.p());

    // Secondary VERTEX FIT
    SecondVertexFit *vtxfit2 = SecondVertexFit::instance();
    vtxfit2->init();
    vtxfit2->setPrimaryVertex(vx_db); // Average interaction point from database
    vtxfit2->AddTrack(0, wLambda);
    vtxfit2->setVpar(vtxLambda);

    if(vtxfit2->Fit()) {
      HepLorentzVector p4Lambda = vtxfit2->p4par();

      p4Lambdavtx.push_back(p4Lambda);
      decayL_Lambda.push_back(vtxfit2->decayLength());
      decayLerr_Lambda.push_back(vtxfit2->decayLengthError());
      chisq_Lambda.push_back(vtxfit2->chisq());
      wLambda_vertex.push_back(vtxfit2->wpar());
      cPoint = vtxfit2->crossPoint();
    } else {
			return StatusCode::SUCCESS;
		}
  
  Ncut30++;
  
  Vp4 p4Lambdabarvtx;
	p4Lambdabarvtx.clear();
  Vdouble  decayL_Lambdabar, decayLerr_Lambdabar, chisq_Lambdabar;
  decayLerr_Lambdabar.clear();
  decayL_Lambdabar.clear();
  chisq_Lambdabar.clear();
  VWTP wLambdabar_vertex;
  wLambdabar_vertex.clear();

  // fit for Lambdabar from pi+ and pbar
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpbarTrk);
  vtxfit->AddTrack(1,  wvpipTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);

  // Primary vertex fit for Lambdabar
  if(!vtxfit->Fit(0)) return StatusCode::SUCCESS;
    vtxfit->Swim(0);
    vtxfit->BuildVirtualParticle(0);
    WTrackParameter wLambdabar = vtxfit->wVirtualTrack(0);
    VertexParameter vtxLambdabar = vtxfit->vpar(0);
    WTrackParameter wtrkprotonbar = vtxfit->wtrk(0);
    p4pbarvtx.push_back(wtrkprotonbar.p());
    WTrackParameter wtrkpionbar = vtxfit->wtrk(1);
    p4pipvtx.push_back(wtrkpionbar.p());

    // Secondary vtx fit
    //SecondVertexFit *vtxfit2 = SecondVertexFit::instance(); // Still in same scope.
    vtxfit2->init();
    vtxfit2->setPrimaryVertex(vx_db); // What is this?
    vtxfit2->AddTrack(0, wLambdabar);
    vtxfit2->setVpar(vtxLambdabar);
    if(vtxfit2->Fit()) {
      HepLorentzVector p4Lambdabar = vtxfit2->p4par();

      p4Lambdabarvtx.push_back(p4Lambdabar);
      decayL_Lambdabar.push_back(vtxfit2->decayLength());
      decayLerr_Lambdabar.push_back(vtxfit2->decayLengthError());
      chisq_Lambdabar.push_back(vtxfit2->chisq());
      wLambdabar_vertex.push_back(vtxfit2->wpar()); // push_back not needed
      //cPoint = vtxfit2->crossPoint(); // only done first time. Why?
    } else {
			return StatusCode::SUCCESS;
		}

  Ncut31++; // Counts passed Lambdabar vertex fit

  // Now I have wLambdabar and wLambda (virtual tracks)
  // Kinematic fit with photons...

  //KinematicFit * kmfit = KinematicFit::instance();
  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();

  //  
  //  Apply Kinematic 4C fit
  // 
  cout<<"before 4c"<<endl;
  if(m_test4C==1) {
//    double ecms = 3.097;
    HepLorentzVector ecms(0.034,0,0,3.097); // TODO: Check

    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    for(int i = 0; i < nGam-1; i++) { // iterate over pairs of gammas
      RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
      for(int j = i+1; j < nGam; j++) {
        RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
        kmfit->init();
        kmfit->AddTrack(0, wLambda_vertex[0]);
        kmfit->AddTrack(1, wLambdabar_vertex[0]);
        kmfit->AddTrack(2, 0.0, g1Trk);
        kmfit->AddTrack(3, 0.0, g2Trk); // also protons here then...
        kmfit->AddFourMomentum(0, ecms);
        bool oksq = kmfit->Fit();
        if(oksq) {
          double chi2 = kmfit->chisq();
          if(chi2 < chisq) {  // better fit? Update gamma choice.
            chisq = chi2;
            ig1 = iGam[i];
            ig2 = iGam[j];
          }
        }
            }
    }
    
    // Get the best values
    if(chisq < 200) {   // Memo suggests less than 100...

      RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
      RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
      kmfit->init();
      kmfit->AddTrack(0, wLambda);
      kmfit->AddTrack(1, wLambdabar);
      kmfit->AddTrack(2, 0.0, g1Trk); // <--- Här läggs gamma1-spåret in
      kmfit->AddTrack(3, 0.0, g2Trk); // <--- Här läggs gamma2-spåret in
      kmfit->AddFourMomentum(0, ecms);
      bool oksq = kmfit->Fit(); // Unnecessary re-check...? REMOVE
      if(oksq) {

  // NEW:
  HepLorentzVector pLambda = kmfit->pfit(0);
  HepLorentzVector pLambdabar = kmfit->pfit(1);
  HepLorentzVector pg1 = kmfit->pfit(2);
  HepLorentzVector pg2 = kmfit->pfit(3);

  // NEW
  m_m_Lambda = pLambda.m();
  m_m_Lambdabar = pLambdabar.m();
  m_p_Lambda = pLambda.mag();
  m_p_Lambdabar = pLambdabar.mag();
  m_p_g1 = pg1.mag();
  m_p_g2 = pg2.mag();

  m_m_gg = (pg1 + pg2).m();

  // Also find the best sigma candidate:
  // Try combining the lambdas with each gamma and check which produces lowest chi2
  HepLorentzVector pLambda_g1 = (kmfit->pfit(0)) + (kmfit->pfit(2));
  HepLorentzVector pLambda_g2 = (kmfit->pfit(0)) + (kmfit->pfit(3));
  HepLorentzVector pLambdabar_g1 = (kmfit->pfit(1)) + (kmfit->pfit(2));
  HepLorentzVector pLambdabar_g2 = (kmfit->pfit(1)) + (kmfit->pfit(3));

  double mLambda_g1 = pLambda_g1.m();
  double mLambda_g2 = pLambda_g2.m();
  double mLambdabar_g1 = pLambdabar_g1.m();
  double mLambdabar_g2 = pLambdabar_g2.m();

  // Find the best combination by taking the sum of the absolute difference from the nominal mass squared
  double chi2_Lambda_g1_Lambdabar_g2 = pow(mLambda_g1 - mSigma0, 2) + pow(mLambdabar_g2 - mSigma0, 2);
  double chi2_Lambda_g2_Lambdabar_g1 = pow(mLambda_g2 - mSigma0, 2) + pow(mLambdabar_g1 - mSigma0, 2);

  // Pick the best combination (smallest chi2) and store it
  if(chi2_Lambda_g1_Lambdabar_g2 < chi2_Lambda_g2_Lambdabar_g1) {
    m_m_Sigma0 = mLambda_g1;
    m_m_Sigmabar0 = mLambdabar_g2;
    m_p_Sigma0 = pLambda_g1.mag();
    m_p_Sigmabar0 = pLambdabar_g2.mag();
  } else {
    m_m_Sigma0 = mLambda_g2;
    m_m_Sigmabar0 = mLambdabar_g1;
    m_p_Sigma0 = pLambda_g2.mag();
    m_p_Sigmabar0 = pLambdabar_g1.mag();
  }
  // save also momentum

	m_chi1 = kmfit->chisq();

	m_tuple4->write();
        Ncut4++;    // What survived the reconstruction...
      } else {
        return StatusCode::SUCCESS;
      }
    } else {
      return StatusCode::SUCCESS;
    } // Just a safeguard.
  }


Ncut6++;
} // end execute()



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode SigmaSigmabar::finalize() {
  cout<<"total number:         "<<Ncut0<<endl;
  cout<<"nGood==4, nCharge==0: "<<Ncut1<<endl;
  cout<<"nGam>=2:              "<<Ncut2<<endl;
  cout<<"Pass Pid:             "<<Ncut3<<endl;
  cout<<"Pass lambda tag:      "<<Ncut30<<endl;
  cout<<"Pass lambdabar tag:   "<<Ncut31<<endl;
  cout<<"Pass 4C:              "<<Ncut4<<endl;
  cout<<"J/psi->Sig0 Sigbar0:  "<<Ncut6<<endl;
  cout<<"Npi, Nproton, mom.:   "<< Npi << " " << Nproton << endl;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  return StatusCode::SUCCESS;
}

