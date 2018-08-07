//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script 
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//
// assume that the environment variable MU2E_TEST_RELEASE points to the 
//  Jul 08 2014 P.Murat
//------------------------------------------------------------------------------

#include <time.h>
#include <TStyle.h>

void rootlogon(){
  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
                                // the line below tells rootcling where to look for 
				// the include files

   // from ROOT plain style
   myStyle->SetCanvasBorderMode(0);
   myStyle->SetPadBorderMode(0);
   myStyle->SetPadColor(0);
   myStyle->SetCanvasColor(0);
   myStyle->SetTitleColor(1);
   myStyle->SetStatColor(0);

   myStyle->SetLabelSize(0.03,"xyz"); // size of axis values

   myStyle->SetHistLineWidth(2);
   myStyle->SetHistLineColor(kBlue+1);

   myStyle->SetLineWidth(1);

   //myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
   //myStyle->SetErrorX(0.001);

   //myStyle->SetPadTickX(0);
   //myStyle->SetPadTickY(0);


   myStyle->SetFuncColor(kRed+1);
   myStyle->SetFuncWidth(3);
   //myStyle->SetLineColor(kRed+1);

   myStyle->SetTitleFont(62, "X");
   myStyle->SetTitleFont(62, "Y");

   myStyle->SetLabelSize(0.05, "X");
   myStyle->SetTitleSize(0.05, "X");

   myStyle->SetLabelSize(0.05, "Y");
   myStyle->SetTitleSize(0.05, "Y");

   myStyle->SetTitleOffset(0.95, "X");

   // default canvas positioning
   myStyle->SetCanvasDefX(900);
   myStyle->SetCanvasDefY(20);
   //myStyle->SetCanvasDefH(550);
   //myStyle->SetCanvasDefW(540);

   myStyle->SetPadBottomMargin(0.1);
   myStyle->SetPadTopMargin(0.1);
   myStyle->SetPadLeftMargin(0.1);
   myStyle->SetPadRightMargin(0.1);
   myStyle->SetPadTickX(1);
   myStyle->SetPadTickY(1);
   myStyle->SetFrameBorderMode(0);

   myStyle->SetTitleBorderSize(0);
   myStyle->SetOptTitle(0);

   // Din letter
   myStyle->SetPaperSize(21, 28);

   myStyle->SetStatBorderSize(0);
   myStyle->SetStatX(0.85);
   myStyle->SetStatY(0.85);
   myStyle->SetStatFont(42);
   myStyle->SetOptStat(111110);// Show overflow and underflow as well
   myStyle->SetOptFit(1111);
   myStyle->SetPalette(1);

   // apply the new style
   gROOT->SetStyle("MyStyle"); //uncomment to set this style
   gROOT->ForceStyle(); // use this style, not the one saved in root files

  gInterpreter->AddIncludePath("./include");
  gInterpreter->AddIncludePath(gSystem->Getenv("CLHEP_INC"));
  gInterpreter->AddIncludePath(Form("%s/include",gSystem->Getenv("ROOTSYS")));


  gInterpreter->AddIncludePath(Form("%s/tex/cdfnotes",
				    gSystem->Getenv("HOME")));

  //  gSystem->SetMakeSharedLib("cd $BuildDir ; g++ -c -g $Opt -pipe -m32 -Wall -W -Woverloaded-virtual -fPIC -pthread $IncludePath $SourceFiles ;  g++ -g $ObjectFiles -shared -Wl,-soname,$LibName.so -m32 $LinkedLibs -o $SharedLib");
//-----------------------------------------------------------------------------
// load in ROOT physics vectors and event generator libraries
//-----------------------------------------------------------------------------
  gSystem->Load("$ROOTSYS/lib/libEG.so");
  //  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("$ROOTSYS/lib/libMinuit.so");
  gSystem->Load("$ROOTSYS/lib/libFumili.so");
  //  gSystem->Load("$ROOTSYS/lib/libTree.so");
  //  gSystem->Load("$ROOTSYS/lib/libRuby.so");
//-----------------------------------------------------------------------------
//  check batch mode
//-----------------------------------------------------------------------------
  const char* opt ;
  int batch_mode = 0;

  int nargs = gApplication->Argc();

  for (int i=1; i<nargs; i++) {
    opt  = gApplication->Argv(i);
    if (strcmp(opt,"-b") == 0) {
      batch_mode = 1;
      break;
    }
  }

  printf("   batch_mode = %i\n",batch_mode);
//-----------------------------------------------------------------------------
// always need libStntuple_loop, but the other 2 libs should be loaded in 
// only if we're running bare root
//-----------------------------------------------------------------------------
  const char* exec_name = gApplication->Argv(0);
 
  printf(" nargs = %2i exec_name = %s\n",nargs, exec_name);

  if (exec_name) {
    if (strstr(exec_name,"root.exe") != 0) {
//-----------------------------------------------------------------------------
// assume STNTUPLE  analysis job
//-----------------------------------------------------------------------------
      if (batch_mode == 1) gSystem->Load("$ROOTSYS/lib/libGui.so");
//-----------------------------------------------------------------------------
// Mu2e Offline libraries
//-----------------------------------------------------------------------------
//     //     gSystem->Load("$MU2E_BASE_RELEASE/lib/libmu2e_Mu2eInterfaces.so");
//     //     gSystem->Load("$MU2E_BASE_RELEASE/lib/libmu2e_CalorimeterGeom.so");
// 
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_base.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_obj.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_loop.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_alg.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_ana.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_val.so");
      
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libmuminus_eplus_ana.so");
      
					// print overflows/underflows in the stat box
      gStyle->SetOptStat(11111111);
					// print fit results in the stat box
      gStyle->SetOptFit(1110);
      TArrow::SetDefaultArrowSize(0.015);
    }
    else if (strstr(exec_name,"mu2e.NNN") != 0) {
      //      gSystem->Load("$MU2E_BASE_RELEASE/lib/libmurat_obj.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_val.so");
      //      gSystem->Load("$MU2E_BASE_RELEASE/lib/libmurat_plot.so");
    }
//-----------------------------------------------------------------------------
//  databases
//-----------------------------------------------------------------------------
//   gSystem->Load("libStntuple_oracle.so");

    if (gSystem->Exec("ls $HOME/root_macros/set_style.C &> /dev/null") == 0) {
      gInterpreter->ExecuteMacro("$HOME/root_macros/set_style.C");
    }
  }
//-----------------------------------------------------------------------------
// report the process ID which simplifies debugging
//-----------------------------------------------------------------------------
  printf(" process ID: %i\n",gSystem->GetPid());
  TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
  gInterpreter->ProcessLine(".! ps | grep root");
}
