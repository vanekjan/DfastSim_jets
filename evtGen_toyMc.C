/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for phi --> K+K-
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TParticlePDG.h"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "StarEvtGenDecayer.h"

//for dca from helix
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"

using namespace std;

//----------------------FROM PYTHIA FAST-SIM-------------------------------------------------------------------------------

TTree* outTree;

float SV_x, SV_y, SV_z;
float MCdecayLength;

float pi_E, pi_px, pi_py, pi_pz, pi_pt, pi_ch;
float pi_eta, pi_phi;

float K_E, K_px, K_py, K_pz, K_pt, K_ch;
float K_eta, K_phi;

TFile* result; //output file in root format

ifstream inFile; //input text files

ofstream outFile; //output text file

string outFileName = "D0.toyMc.root"; //default output file name (do not change - see submit XML)

//-------------------------------------------------------------------------------------------------------------------------

void initEvtGen();
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TVector3 v00);
void decayAndFill(int PDG_id, TLorentzVector* b, double const weight, TClonesArray& daughters); //Decay particle - PDG ID
void bookObjects();
void write();

StarEvtGenDecayer* starEvtGenDecayer = NULL;

//============== main  program ==================

//centrality range in %, default range is 0% to 80%, pTspectrum = 0 - flat pT, 1 - Levy (weight function)
//specifically set in run_evtGent_toyMc.C
void evtGen_toyMc() 
{
  cout<<"Starting EvtGen"<<endl;
  initEvtGen();
  cout<<"initEvtGen() done..."<<endl;

  gRandom->SetSeed();
  bookObjects();

  TClonesArray ptl("TParticle", 10); //array of 10 TParticles for daughters

  string Input_kine_line; //to load input kinematics from a text file

  //to store input kinematics in array: 
  //#N(0)	pid(1)	stat(2) 	E(3) 	Px(4)	Py(5)	Pz(6)	Eta(7)	Phi(8)
  double InKinematics[9]; 

  while(getline(inFile, Input_kine_line))
  {
    stringstream Input_kine_line_stream(Input_kine_line);

    //loop over kinemaric variables on one line of the input file (kinematics of one mother particle)
    for(unsigned int i = 0; i < 9; i++)
    {
      Input_kine_line_stream>>InKinematics[i];
    }    

    //TLorentzVector(px, py, pz, E);
    TLorentzVector* b_d = new TLorentzVector(InKinematics[4], InKinematics[5], InKinematics[6], InKinematics[3]);

    cout<<InKinematics[1]<<endl;

    decayAndFill(InKinematics[1], b_d, 1., ptl);

  }


  cout<<"Write!"<<endl;
  write();
  cout<<"Written!"<<endl;



  return;

}


void decayAndFill(int PDG_id, TLorentzVector* b, double const weight, TClonesArray& daughters)//decay for D+
{

  starEvtGenDecayer->Decay(PDG_id, b); //Decay particle

  starEvtGenDecayer->ImportParticles(&daughters); //get daughters from decay

  //daugter four-momenta
  TLorentzVector kMom;
  TLorentzVector p1Mom;

  TVector3 v00; //decay vertex

  int nTrk = daughters.GetEntriesFast();

  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (abs(ptl0->GetPdgCode()))
    {
      case 321:
        ptl0->Momentum(kMom);
        v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
        break;
      case 211:
        ptl0->Momentum(p1Mom);
        break;
      default:
        break;
    }
  }

  //clear daughters for next decay
  daughters.Clear();

  //save decay
  fill(PDG_id, b, weight, kMom, p1Mom, v00);

}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TVector3 v00)
{
  
  

   
  // reconstruct

  //float const kDca = dca(kMom.Vect(), v00, vertex);
  //float const pDca = dca(p1Mom.Vect(), v00, vertex);

  //secondary vertex
  SV_x = v00.x();
  SV_y = v00.y();
  SV_z = v00.z();

  MCdecayLength = v00.Mag();  


  //pion
  pi_E = p1Mom.E();
  pi_px = p1Mom.Px();
  pi_py = p1Mom.Py();
  pi_pz = p1Mom.Pz();
  pi_pt = p1Mom.Pt();

  pi_eta = p1Mom.Eta();
  pi_phi = p1Mom.Phi();

  //kaon
  K_E = kMom.E();
  K_px = kMom.Px();
  K_py = kMom.Py();
  K_pz = kMom.Pz();
  K_pt = kMom.Pt();

  K_eta = kMom.Eta();
  K_phi = kMom.Phi();


  outTree->Fill();



}


//___________
void bookObjects()
{

  
  cout << "Loading input files ..." << endl;
  cout<<endl;

  inFile.open("./input/d0_jetscape.txt");

  if(!inFile.is_open())
  {
    cout<<"Failed to open input file!"<<endl;
    return;
  } 

  
  result = new TFile(outFileName.c_str(), "recreate");
  //result->SetCompressionLevel(1);
  result->cd();


  outTree = new TTree("ntp", "ntp");

  //event
  outTree->Branch("SV_x", &SV_x, "SV_x/F"); 
  outTree->Branch("SV_y", &SV_y, "SV_y/F"); 
  outTree->Branch("SV_z", &SV_z, "SV_z/F"); 

  outTree->Branch("MCdecayLength", &MCdecayLength, "MCdecayLength/F");

  //pion
  outTree->Branch("pi_ch", &pi_ch, "pi_ch/F");   
  outTree->Branch("pi_E", &pi_E, "pi_E/F");  
  outTree->Branch("pi_px", &pi_px, "pi_px/F");
  outTree->Branch("pi_py", &pi_py, "pi_py/F");
  outTree->Branch("pi_pz", &pi_pz, "pi_pz/F");  
  outTree->Branch("pi_pt", &pi_pt, "pi_pt/F");

  outTree->Branch("pi_phi", &pi_phi, "pi_phi/F");
  outTree->Branch("pi_eta", &pi_eta, "pi_eta/F");

  //kaon
  outTree->Branch("K_ch", &K_ch, "K_ch/F");
  outTree->Branch("K_E", &K_E, "K_E/F");  
  outTree->Branch("K_px", &K_px, "K_px/F");
  outTree->Branch("K_py", &K_py, "K_py/F");
  outTree->Branch("K_pz", &K_pz, "K_pz/F");  
  outTree->Branch("K_pt", &K_pt, "K_pt/F");

  outTree->Branch("K_phi", &K_phi, "K_phi/F");
  outTree->Branch("K_eta", &K_eta, "K_eta/F");

  

  cout << "Done with loading all files ..." << endl;
}
//_______________________________________________________________________________________________________________________


float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaHelix(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex, float charge) //DCA to PV - calcualted same way as in StPicoTrack
{
  //all position vectors used for calcualtion with helix transformed form mum to cm

  StThreeVectorF p_work(p.X(), p.Y(), p.Z()); //create StThreeVectorF of momentum
  StThreeVectorF pos_work(pos.X()/1e4, pos.Y()/1e4, pos.Z()/1e4); //create StThreeVectorF of position

  StThreeVectorF vertex_work(vertex.X()/1e4, vertex.Y()/1e4, vertex.Z()/1e4); //create StThreeVectorF of PV


  StPhysicalHelixD helix_work(p_work, pos_work, -5*kilogauss , charge); //true helix with origin near SV

  helix_work.moveOrigin(helix_work.pathLength(vertex_work)); // move origin to DCA to PV

  StThreeVectorF origin = helix_work.origin(); //new helix origin at DCA


  return (origin - vertex_work).mag()*1e4;

}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

  return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 newPos(pos);
  newPos.SetZ(0);

  TVector3 newP(p);
  newP.SetZ(0);

  TVector3 newVertex(vertex);
  newVertex.SetZ(0);

  TVector3 posDiff = newPos - newVertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
  return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  if (sin(p.Theta()) == 0) return 0;
  else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}


void write()
{
  result->cd();
  outTree->Write();
  result->Close();
}

void initEvtGen()
{
  cout<<"initEvtGen start..."<<endl;
  EvtRandomEngine* eng = 0;
  eng = new EvtSimpleRandomEngine();
  cout<<"setting random engine..."<<endl;
  EvtRandom::setRandomEngine((EvtRandomEngine*)eng);
  cout<<"done"<<endl;
  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

  EvtExternalGenList genList;
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();

  TString Decay_DEC="StRoot/StarGenerator/EvtGen1_06_00/DECAY.DEC";
  TString Evt_pdl="StRoot/StarGenerator/EvtGen1_06_00/evt.pdl";
  EvtGen *myGenerator=new EvtGen(Decay_DEC,Evt_pdl,(EvtRandomEngine*)eng,radCorrEngine, &extraModels);
  starEvtGenDecayer=new StarEvtGenDecayer(myGenerator);
  starEvtGenDecayer->SetDecayTable("D0.PHSP.DEC");
  starEvtGenDecayer->SetDebug(0);
}
