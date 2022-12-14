///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "muminus_eplus/ana/scripts/modules.hh"
#include "Stntuple/scripts/global_vars.h"

def_name mutoeplus_001("mutoeplus_ana001");
//def_name track_ce ("track_ana001");
//def_name track_dio("track_ana001");
///////////////////////////////////////////////////////////////////////////////



//-----------------------------------------------------------------------------
void  mutoeplus_ana001(int GenCode=28, double NEvents=-1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_ptrk = (TMuminuseplusAna001Module*) g.x->AddModule("TMuminuseplusAna001Module",0);  
  //  m_ptrk->SetDebugBit(43, 1);
  m_ptrk->SetNormalization(NEvents);
  m_ptrk->SetGeneratorCode(GenCode);

}
