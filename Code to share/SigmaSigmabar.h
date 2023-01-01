#ifndef Physics_Analysis_SigmaSigmabar_H
#define Physics_Analysis_SigmaSigmabar_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "TH1.h"
//#include "VertexFit/ReadBeamParFromDb.h"


class SigmaSigmabar : public Algorithm {

public:
  SigmaSigmabar(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:

  //ReadBeamParFromDb m_reader;
  // Declare r0, z0 cut for charged tracks
  double m_vr0cut;
  double m_vz0cut;

  //Declare energy, dphi, dthe cuts for fake gamma's
  double m_energyThreshold;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaAngleCut;

  // Flag: use or not (1/0)
  int m_test4C;
  int m_test5C;

  // Flag: use or not (1/0)
  int m_checkDedx;
  int m_checkTof;

  // define Ntuples here

  // NTuple::Tuple*  m_tuple1;      // charged track vertex
  // NTuple::Item<double>  m_vx0;
  // NTuple::Item<double>  m_vy0;
  // NTuple::Item<double>  m_vz0;
  // NTuple::Item<double>  m_vr0;
  // NTuple::Item<double>  m_rvxy0;
  // NTuple::Item<double>  m_rvz0;
  // NTuple::Item<double>  m_rvphi0;

  // NTuple::Tuple*  m_tuple2;      // fake photon
  // NTuple::Item<double>  m_dthe;
  // NTuple::Item<double>  m_dphi;
  // NTuple::Item<double>  m_dang;
  // NTuple::Item<double>  m_eraw;

  // NTuple::Tuple*  m_tuple3;     // rhopi: raw mgg, etot
  // NTuple::Item<double>  m_m2gg;
  // NTuple::Item<double>  m_etot;

  NTuple::Tuple*  m_tuple4;     // rhopi 4C
  // Chi2 from kinematic fit
  NTuple::Item<double>  m_chi1; 
  // Masses
  NTuple::Item<double>  m_m_p;
  NTuple::Item<double>  m_m_pbar;
  NTuple::Item<double>  m_m_pip;
  NTuple::Item<double>  m_m_pim;
  NTuple::Item<double>  m_m_Lambda;
  NTuple::Item<double>  m_m_Lambdabar;
  NTuple::Item<double>  m_m_Sigma0; 
  NTuple::Item<double>  m_m_Sigmabar0;
  NTuple::Item<double>  m_m_gg;  // g1 + g2

  // Full momentum (absolute value) -- use .rho()
  NTuple::Item<double>  m_p_p;
  NTuple::Item<double>  m_p_pbar;
  NTuple::Item<double>  m_p_pip;
  NTuple::Item<double>  m_p_pim;
  NTuple::Item<double>  m_p_Lambda;
  NTuple::Item<double>  m_p_Lambdabar;
  NTuple::Item<double>  m_p_Sigma0; 
  NTuple::Item<double>  m_p_Sigmabar0;
  NTuple::Item<double>  m_p_g1;
  NTuple::Item<double>  m_p_g2;
  NTuple::Item<double>  m_p_gg;  // g1 + g2

  // momentum x-direction -- use .px()
  NTuple::Item<double>  m_px_p;
  NTuple::Item<double>  m_px_pbar;
  NTuple::Item<double>  m_px_pip;
  NTuple::Item<double>  m_px_pim;
  NTuple::Item<double>  m_px_Lambda;
  NTuple::Item<double>  m_px_Lambdabar;
  NTuple::Item<double>  m_px_Sigma0; 
  NTuple::Item<double>  m_px_Sigmabar0;
  NTuple::Item<double>  m_px_g1;
  NTuple::Item<double>  m_px_g2;
  NTuple::Item<double>  m_px_gg;  // g1 + g2

  // momentum y-direction -- use .py()
  NTuple::Item<double>  m_py_p;
  NTuple::Item<double>  m_py_pbar;
  NTuple::Item<double>  m_py_pip;
  NTuple::Item<double>  m_py_pim;
  NTuple::Item<double>  m_py_Lambda;
  NTuple::Item<double>  m_py_Lambdabar;
  NTuple::Item<double>  m_py_Sigma0; 
  NTuple::Item<double>  m_py_Sigmabar0;
  NTuple::Item<double>  m_py_g1;
  NTuple::Item<double>  m_py_g2;
  NTuple::Item<double>  m_py_gg;  // g1 + g2

  // momentum z-direction -- use .pz()
  NTuple::Item<double>  m_pz_p;
  NTuple::Item<double>  m_pz_pbar;
  NTuple::Item<double>  m_pz_pip;
  NTuple::Item<double>  m_pz_pim;
  NTuple::Item<double>  m_pz_Lambda;
  NTuple::Item<double>  m_pz_Lambdabar;
  NTuple::Item<double>  m_pz_Sigma0; 
  NTuple::Item<double>  m_pz_Sigmabar0;
  NTuple::Item<double>  m_pz_g1;
  NTuple::Item<double>  m_pz_g2;
  NTuple::Item<double>  m_pz_gg;  // g1 + g2  // NOTE: not necessary?

  // Energy
  NTuple::Item<double>  m_E_p;
  NTuple::Item<double>  m_E_pbar;
  NTuple::Item<double>  m_E_pip;
  NTuple::Item<double>  m_E_pim;
  NTuple::Item<double>  m_E_Lambda;
  NTuple::Item<double>  m_E_Lambdabar;
  NTuple::Item<double>  m_E_Sigma0; 
  NTuple::Item<double>  m_E_Sigmabar0;
  NTuple::Item<double>  m_E_g1;
  NTuple::Item<double>  m_E_g2;
  NTuple::Item<double>  m_E_gg;  // g1 + g2

  // That's 68 outputs... large files & time consuming?

  // Xiaorong suggests: 
  // p_x, p_y, p_z of all particles + Energy + invariant mass (for photons: combined)
  // for p, pbar, pi+, pi-, lambda, lambdabar, sigma, sigmabar, g1, g2, gg

  // For topoana
  NTuple::Item<int>     m_idxmc;  // (m_indexmc) ?? Check required input from other email
  NTuple::Array<int>    m_pdgid;   
  NTuple::Array<int>    m_motheridx;


  // NTuple::Tuple*  m_tuple5;     // rhopi 5C
  // NTuple::Item<double>  m_chi2;
  // NTuple::Item<double>  m_mrh0;
  // NTuple::Item<double>  m_mrhp;
  // NTuple::Item<double>  m_mrhm;

  // NTuple::Tuple*  m_tuple6;    // photons
  // NTuple::Item<double>  m_fcos;
  // NTuple::Item<double>  m_elow;

  // NTuple::Tuple* m_tuple7;    // dE/dx
  // NTuple::Item<double> m_ptrk;
  // NTuple::Item<double> m_chie;
  // NTuple::Item<double> m_chimu;
  // NTuple::Item<double> m_chipi;
  // NTuple::Item<double> m_chik;
  // NTuple::Item<double> m_chip;
  // NTuple::Item<double> m_probPH;
  // NTuple::Item<double> m_normPH;
  // NTuple::Item<double> m_ghit;
  // NTuple::Item<double> m_thit;

  // NTuple::Tuple* m_tuple8;   // endcap tof
  // NTuple::Item<double> m_ptot_etof;
  // NTuple::Item<double> m_cntr_etof;
  // NTuple::Item<double> m_te_etof;
  // NTuple::Item<double> m_tmu_etof;
  // NTuple::Item<double> m_tpi_etof;
  // NTuple::Item<double> m_tk_etof;
  // NTuple::Item<double> m_tp_etof;
  // NTuple::Item<double> m_ph_etof;
  // NTuple::Item<double> m_rhit_etof;
  // NTuple::Item<double> m_qual_etof;

  // NTuple::Tuple* m_tuple9;  // barrel inner tof
  // NTuple::Item<double> m_ptot_btof1;
  // NTuple::Item<double> m_cntr_btof1;
  // NTuple::Item<double> m_te_btof1;
  // NTuple::Item<double> m_tmu_btof1;
  // NTuple::Item<double> m_tpi_btof1;
  // NTuple::Item<double> m_tk_btof1;
  // NTuple::Item<double> m_tp_btof1;
  // NTuple::Item<double> m_ph_btof1;
  // NTuple::Item<double> m_zhit_btof1;
  // NTuple::Item<double> m_qual_btof1;

  // NTuple::Tuple* m_tuple10;  // barrel outer tof
  // NTuple::Item<double> m_ptot_btof2;
  // NTuple::Item<double> m_cntr_btof2;
  // NTuple::Item<double> m_te_btof2;
  // NTuple::Item<double> m_tmu_btof2;
  // NTuple::Item<double> m_tpi_btof2;
  // NTuple::Item<double> m_tk_btof2;
  // NTuple::Item<double> m_tp_btof2;
  // NTuple::Item<double> m_ph_btof2;
  // NTuple::Item<double> m_zhit_btof2;
  // NTuple::Item<double> m_qual_btof2;
  
  // NTuple::Tuple* m_tuple11;  // Particle ID info.
  // NTuple::Item<double> m_ptrk_pid;
  // NTuple::Item<double> m_cost_pid;
  // NTuple::Item<double> m_dedx_pid;
  // NTuple::Item<double> m_tof1_pid;
  // NTuple::Item<double> m_tof2_pid;
  // NTuple::Item<double> m_prob_pid;
  // NTuple::Item<double> m_pproton_pid;
  // NTuple::Item<double> m_ppion_pid;

  // TH1D *hist;
};

#endif 
