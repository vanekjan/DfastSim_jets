void run_evtGent_toyMc(int nEvents = 1e3)
{

  gSystem->Load("libStTableUtilities.so");
  gSystem->Load("libTable");
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");        // new addition 22jul99

  gSystem->Load( "libVMC.so" );
  gSystem->Load( "libSt_g2t.so" );
  gSystem->Load( "libSt_geant_Maker.so" );
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load("libHepMC2_06_09.so");
  gSystem->Load("libPythia8_1_86.so");
  gSystem->Load("libPhotos3_61.so");
  gSystem->Load("libTauola1_1_5.so");
  gSystem->Load(".sl73_gcc485/lib/libEvtGen1_06_00.so");

  gInterpreter->AddIncludePath("StRoot/StarGenerator/EvtGen1_06_00");
  gInterpreter->AddIncludePath("$STAR/StRoot/StarGenerator/HepMC2_06_09");

  StChain chain("myChain");

  gROOT->ProcessLine(Form(".x evtGen_toyMc.C+()"));

  chain.Finish();

  return;
}
