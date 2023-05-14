#!/usr/bin/env python
#
import os, re, string, sys, subprocess
Import('env')
sys.path.append(os.getenv("MUSE_WORK_DIR")+'/site_scons')
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
x = subprocess.call(os.getenv("MUSE_WORK_DIR")+'/muminus_eplus/scripts/build_config_muse muminus_eplus',shell=True)
muminus_eplus_env = env.Clone()

muminus_eplus_env['CPPPATH' ].append('-I'+os.environ['MUSE_WORK_DIR']+'/include');
muminus_eplus_env['CXXFLAGS'].append('-I'+os.environ['MUSE_WORK_DIR']+'/include');
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
exec(open(os.environ['MUSE_WORK_DIR']+"/site_scons/stntuple_site_init.py").read())

from stntuple_helper import *

muminus_eplus_env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen})
muminus_eplus_env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('muminus_eplus_env')
Export('stntuple_helper')
