project('shyfem',
        ['fortran','c','cpp'],
        default_options : [
                          'buildtype=release', # {plain, debug, debugoptimized, release, minsize, custom}
                          'warning_level=0',  # From 0 = none to 3 = highesti
                          ]
        )

_compiler = meson.get_compiler('fortran')
_Ccompiler = meson.get_compiler('c')
_CPPcompiler = meson.get_compiler('cpp')
#######################################################################
#----------------- BEGINNING OF THE CUSTOMISABLE AREA -----------------
#######################################################################

#------------------- eventually setup optimization arguments -------------
_bt = get_option('buildtype')
if _bt == 'debugoptimized' or _bt == 'release'
   if _compiler.get_id() == 'intel'
    add_global_arguments([  '-march=native',
                            '-funroll-loops',
                            '-xhost',
                            '-xcore-avx2',
                            '-axCORE-AVX2,AVX'
                          ],language : ['fortran'] )
  elif _compiler.get_id() == 'gcc'
    add_global_arguments([  '-mcpu=native',
                            '-ffast-math',
                            '-ftree-vectorize',
                            '-ftree-vectorizer-verbose=0',
                            '-funroll-loops',
                          ],language : ['fortran'] )
  endif
endif
if _bt == 'debug'
  if _compiler.get_id() == 'intel'
    add_global_arguments([ 
                          '-trace',
                          '-qopt-report=5',
                          '-qopt-report-phase=all',
                          '-diag-disable 10397', #disables info on qoptr files creation
                          '-traceback',
                          '-fp-stack-check',
                          '-fpe0',
                          '-ftrapuv',
                          '-CB',#'-fp-model source','-debug all',
                          ],language : ['fortran'] )
  else
    add_global_arguments([
                          '-fbacktrace','-fbounds-check',
                          '-Wall','-Wextra',
                          '-fimplicit-none',
                          '-fcheck=all',
                          '-Wtabs',
                          '-ffpe-trap=zero,invalid,overflow',
                          ],language : ['fortran'] )
  endif

endif

#######################################################################
#-------------------- END OF THE CUSTOMISABLE AREA --------------------
#######################################################################

#--------------- check eventual incompatibilities ---------------------------
if get_option('use_AmgX')==true and get_option('use_PETSc')==false 
   error('In order to use AmgX you must activate the PETSc solver')
endif
count_solvers=0
if get_option('use_SPK')==true
  count_solvers+=1
endif
if get_option('use_PETSc')==true
  count_solvers+=1
endif
if get_option('use_paradiso')==true
  count_solvers+=1
endif
if count_solvers != 1
   error('You must activate no more and no less than 1 solver among PETSc, SparseKit and Paradiso')
endif

###############################################################################
message('------------------- preparing global settings --------------------')
###############################################################################

if get_option('Verbose')
  warning('Verbose option is activated, it induces verbose shyfem output slowing down A LOT the shyfem program')
endif
if get_option('test_zeta')
warning('test_zeta option is activated, it creates an output file of the solution vector at every timestep which slows down the shyfem program')
endif

#--------------- setting up pragma precompilation directives ----------------
_config_file = files('pragma_directives.h.meson_input')
_config_data = configuration_data()

_config_data.set('Verbose',get_option('Verbose'))
_config_data.set('test_zeta',get_option('test_zeta'))

configure_file(output:'pragma_directives.h',
               configuration:_config_data,
               input:_config_file)

if _compiler.get_id() == 'intel'
    add_global_arguments('-fpp',language : ['fortran'] )
elif _compiler.get_id() == 'gcc'
  add_global_arguments('-cpp',language : ['fortran'] )
endif

#----------------setting up mpi dependency -----------------------------------
if get_option('use_PETSc') == true or get_option('use_MPI') == true
  _with_mpi=true
  message('Compiling with mpi dependency')
  _mpifortran_dep = dependency('mpi',language:'fortran')

  _code = '''#include "mpi.h" '''
  _compiles = _compiler.compiles(_code, name:'can include mpi' )
  if not _compiles
      error('\n Error Cannot find header mpi.h')
  endif
  _mpi_dependencies = [_mpifortran_dep]
  _mpi_src_files = ['fem3d/mod_shympi_node.f','fem3d/mod_shympi_node_internal.f']
else
  _with_mpi=false
  _mpi_src_files = ['fem3d/mod_shympi_dummy.f']
  _mpi_dependencies = []
endif  

###############################################################################
message('------------------- preparing shyfem build target ----------------')
###############################################################################
_shyfem_dependencies= []
_shyfem_include_directories=[]
_shyfem_source_files=[]
_shyfem_fortran_args = []

_shyfem_source_files += _mpi_src_files
_shyfem_dependencies += _mpi_dependencies

#------------------- petsc solver ----------------------------------
if get_option('use_PETSc') == true
   message('Selecting the PETSc solver')
   # create list of dependencies, containing for each one the array : [library_name, installation_path, include_code, compiler]
   _lib_properties=[['PETSc','path_petsc','''#include "petsc/finclude/petsc.h" ''',_compiler] ] #----------------PETSc------------------
   if get_option('use_AmgX') == true
     _lib_properties +=[['amgxsh','path_amgx','''# include <amgx_c.h> ''',_Ccompiler]]          #----------------AmgX------------------
     _lib_properties +=[['AmgXWrapper','path_amgxwrapper','''# include <AmgXSolver.hpp>''',_CPPcompiler]] #---------AmgXWrapper----------
     _shyfem_source_files += ['amgx-c-wrapper/amgx-c-wrapper/src/amgx_c_wrapper.cpp']   #---------------------- c-wrapper to c++ AmgXWrapper ----
   else
     _shyfem_source_files +=['fem3d/amgx_c_wrapper_dummy.f']
   endif
   foreach _libproperty : _lib_properties
       _lib=_libproperty[0] 
       _opt=_libproperty[1]
       _dir=get_option(_opt)
       _code=_libproperty[2]
       _comp=_libproperty[3] 
       if _dir != 'none'
          message('Using user-provided '+_lib+' directory '+_dir)
          _shyfem_dependencies += [ _compiler.find_library(_lib,dirs:[ _dir+'/lib/', _dir+'/lib64/' ] )  ]
          _shyfem_include_directories += [ include_directories(_dir+'/include/') ]
       else
          _shyfem_dependencies += [dependency(_lib,
                                    not_found_message:_lib+' dependency couldnt be found using pkg-config,'
                                         +' meson cannot proceed because '
                                         +' both library and include paths are required to use this libray\.' 
                                         +' Please provide installation path through command line '
                                         +' -D'+_opt+'=PATH_TO_INSTALLATION_DIRECTORY'
                                         +' or add file .pc path to environmental variable PKG_CONFIG_PATH'
                                         +' (path might be $<LIBNAME>_HOME/lib/pkgconfig)',
                                    )]
       endif
       _compiles = _comp.compiles(_code, name:'can include '+_lib,include_directories:_shyfem_include_directories,dependencies:_shyfem_dependencies) 
       
       if not _compiles
           error('Failed to compile code >> '+_code)
       endif
   endforeach
   _shyfem_source_files+=[
                  'fem3d/simsys_petsc.f', 
                  'fem3d/mod_petsc_system.f', 
                  'fem3d/mod_petsc_global.f', 
                  'fem3d/subcoo.f', # still needed for function loccoo3d used in newnohydro 
                  'fem3d/spk_blas.f',     # needed by spk_linpack
                  'fem3d/spk_linpack.f',  # still needed for function dsvdc used in tide_utils.f 
                 ] 
   #-------------------- cuda ------------------------------------------
   if get_option('num_GPU')>0
     _cudalibs=['cublas','cudart','cusparse','cusolver']
     _cuda_dep=[]
     if get_option('path_cuda') != 'none'
       _cuda_dir=get_option('path_cuda')
       foreach _culib : _cudalibs
          _cuda_dep += [ _compiler.find_library(_culib,dirs:[ _cuda_dir+'/lib/', _cuda_dir+'/lib64/' ] )  ]
          _shyfem_include_directories+=include_directories(_cuda_dir+'/include/')
       endforeach
     else
       foreach _culib : _cudalibs
          _dep = dependency(_culib,
                           not_found_message: _culib+' dependency couldnt be found using pkg-config,'
                                              +' now trying with find_library as we only need the library path,' 
                                              +'  if it does not work please provide installation path in var "cuda_dir"'
                                              +' or add cuda library path to environmental variable "LIBRARY_PATH"',
                           required:false)
          if not _dep.found()
            _dep=_compiler.find_library(_culib)
          endif
          _cuda_dep += [ _dep ]
       endforeach
     endif
     _shyfem_dependencies+=_cuda_dep
   endif
#------------------- sparsekit solver ----------------------------------
elif get_option('use_SPK') == true
   message('Selecting the sparsekit solver')
   _shyfem_source_files+=[
                  'fem3d/mod_system.f', 
                  'fem3d/simsys_spk.f',
                  'fem3d/subspk.f', 
                  'fem3d/subcoo.f', 
                  'fem3d/subcooutil.f', 
                  'fem3d/spk_ilut.f', 
                  'fem3d/spk_iters.f', 
                  'fem3d/spk_itaux.f', 
                  'fem3d/spk_coo.f', 
                  'fem3d/spk_itpack.f', 
                  'fem3d/spk_matvec.f', 
                  'fem3d/spk_blas.f', 
                  'fem3d/spk_linpack.f', 
                 ] 
#------------------- pardiso solver ----------------------------------
elif get_option('use_paradiso') == true
   message('Selecting the paradiso solver')
   _shyfem_source_files += [ 'simsys_pard.f', 'subpard.f', 'subcoo.o']
   _shyfem_fortran_args += ['./bin/find_pardlib.sh' ]
endif


#------------------- gotm dependency -----------------------------------
subdir('femgotm') # defines gotm_dep
_shyfem_dependencies+=[_gotm_dep]

#------------------- add common source files ---------------------------
_shyfem_source_files += [
             'fem3d/shyfem.f', 
             'fem3d/new3di.f', 
             'fem3d/newrog.f', 
             'fem3d/newbcl.f', 
             'fem3d/newbcl0.f', 
             'fem3d/newini.f', 
             'fem3d/newtra.f', 
             'fem3d/newchao.f', 
             'fem3d/newcon.f', 
             'fem3d/newconz.f', 
             'fem3d/newcon_omp.f', 
             'fem3d/newchk.f', 
             'fem3d/newfix.f', 
             'fem3d/newtvd.f', 
             'fem3d/newstab.f', 
             'fem3d/newsig.f', 
             'fem3d/newnudge.f', 
             'fem3d/newpoi.f', 
             'fem3d/newbfm.f', 
             'fem3d/newlevels.f', 
             'fem3d/mod_hydro.f', 
             'fem3d/mod_hydro_vel.f', 
             'fem3d/mod_hydro_print.f', 
             'fem3d/mod_hydro_baro.f', 
             'fem3d/mod_diff_visc_fric.f', 
             'fem3d/mod_roughness.f', 
             'fem3d/mod_ts.f', 
             'fem3d/mod_area.f', 
             'fem3d/mod_bound_dynamic.f', 
             'fem3d/mod_diff_aux.f', 
             'fem3d/mod_gotm_aux.f', 
             'fem3d/mod_bnd_aux.f', 
             'fem3d/mod_depth.f', 
             'fem3d/mod_geom_dynamic.f', 
             'fem3d/mod_internal.f', 
             'fem3d/mod_nohyd.f', 
             'fem3d/mod_nudging.f', 
             'fem3d/mod_bclfix.f', 
             'fem3d/mod_keps.f', 
             'fem3d/mod_fluidmud.f', 
             'fem3d/mod_sinking.f', 
             'fem3d/mod_waves.f', 
             'fem3d/mod_meteo.f', 
             'fem3d/mod_geom.f', 
             'fem3d/mod_conz.f', 
             'fem3d/mod_bound_geom.f', 
             'fem3d/mod_bnd.f', 
             'fem3d/mod_tvd.f', 
             'fem3d/mod_bndo.f', 
             'fem3d/mod_system.f', 
             'fem3d/mod_layer_thickness.f', 
             'fem3d/mod_subset.f', 
             'fem3d/mod_lagrange.f', 
             'fem3d/submpi_ghost.f', 
             'fem3d/submpi_basin.f', 
             'fem3d/submpi_custom.f', 
             'fem3d/submpi_buffer.f', 
             'fem3d/submpi_time.f', 
             'fem3d/shympi_general.f', 
             'fem3d/shympi_debug.f', 
             'fem3d/submpi_assert.f', 
             'fem3d/subver.f', 
             'fem3d/subcst.f', 
             'fem3d/subsys.f', 
             'fem3d/subnsh.f', 
             'fem3d/subtime.f', 
             'fem3d/suboutput.f', 
             'fem3d/suboutputd.f', 
             'fem3d/subbas.f', 
             'fem3d/subnev.f', 
             'fem3d/subnsu.f', 
             'fem3d/subn35.f', 
             'fem3d/subdep.f', 
             'fem3d/subdef.f', 
             'fem3d/subfil.f', 
             'fem3d/subnls.f', 
             'fem3d/subpar3.f', 
             'fem3d/subintp.f', 
             'fem3d/subsss.f', 
             'fem3d/subround.f', 
             'fem3d/subscn.f', 
             'fem3d/subssv.f', 
             'fem3d/subssm.f', 
             'fem3d/subext.f', 
             'fem3d/subexta.f', 
             'fem3d/subous.f', 
             'fem3d/subousa.f', 
             'fem3d/subbnd.f', 
             'fem3d/subnos.f', 
             'fem3d/subwat.f', 
             'fem3d/subn11.f', 
             'fem3d/subdry.f', 
             'fem3d/subcus.f', 
             'fem3d/subcon1.f', 
             'fem3d/subtide.f', 
             'fem3d/tide_utils.f', 
             'fem3d/subres.f', 
             'fem3d/subuti.f', 
             'fem3d/subutil.f', 
             'fem3d/sublin.f', 
             'fem3d/subvola.f', 
             'fem3d/subgeo.f', 
             'fem3d/debug.f', 
             'fem3d/submod.f', 
             'fem3d/subele.f', 
             'fem3d/subqfxt.f', 
             'fem3d/subqfxf.f', 
             'fem3d/subqfxm1.f', 
             'fem3d/subqfxm2.f', 
             'fem3d/subqfxm3.f', 
             'fem3d/subqfxm4.f', 
             'fem3d/subqfxm5.f', 
             'fem3d/subqfxm6.f', 
             'fem3d/subqfxu1.f', 
             'fem3d/subqfxu2.f', 
             'fem3d/subqfxu3.f', 
             'fem3d/subqfxu4.f', 
             'fem3d/submud_dummy.f', 
             'fem3d/subgotm.f', 
             'fem3d/subbndo.f', 
             'fem3d/subbnds.f', 
             'fem3d/intp.f', 
             'fem3d/subssf.f', 
             'fem3d/bdist.f', 
             'fem3d/links.f', 
             'fem3d/sublnka.f', 
             'fem3d/sublnkd.f', 
             'fem3d/sublnks.f', 
             'fem3d/sublnku.f', 
             'fem3d/sedi3d.f', 
             'fem3d/sedtrans05.f', 
             'fem3d/subwaves.f', 
             'fem3d/substress.f', 
             'fem3d/subhisto.f', 
             'fem3d/lagrange_main.f', 
             'fem3d/lagrange_back.f', 
             'fem3d/lagrange_dif.f', 
             'fem3d/lagrange_cont.f', 
             'fem3d/lagrange_init.f', 
             'fem3d/lagrange_inout.f', 
             'fem3d/lagrange_util.f', 
             'fem3d/lagrange_utils.f', 
             'fem3d/lagrange_larve.f', 
             'fem3d/lagrange_flux.f', 
             'fem3d/lagrange_track.f', 
             'fem3d/lagrange_util_tr.f', 
             'fem3d/lagrange_decay.f', 
             'fem3d/lagrange_sedim.f', 
             'fem3d/lagrange_connect.f', 
             'fem3d/lagrange_vertical.f', 
             'fem3d/lagrange_stk.f', 
             'fem3d/lagrange_read.f', 
             'fem3d/line_util.f', 
             'fem3d/subdts.f', 
             'fem3d/loading.f', 
             'fem3d/subrst.f', 
             'fem3d/subnosa.f', 
             'fem3d/subdif.f', 
             'fem3d/subreg.f', 
             'fem3d/subfind.f', 
             'fem3d/atoxi3d.f', 
             'fem3d/atoxi.f', 
             'fem3d/kepsd.f', 
             'fem3d/keps_util.f', 
             'fem3d/newexpl.f', 
             'fem3d/subrnd.f', 
             'fem3d/submet.f', 
             'fem3d/subrgf.f', 
             'fem3d/subcoord.f', 
             'fem3d/subsrt.f', 
             'fem3d/subfvl.f', 
             'fem3d/subtrace.f', 
             'fem3d/subproj.f', 
             'fem3d/subeos.f', 
             'fem3d/sigmautil.f', 
             'fem3d/subwrt.f', 
             'fem3d/subvintp.f', 
             'fem3d/subfemfile.f', 
             'fem3d/subtsuvfile.f', 
             'fem3d/subflxa.f', 
             'fem3d/subflx3d.f', 
             'fem3d/subflx.f', 
             'fem3d/subflxu.f', 
             'fem3d/subboxa.f', 
             'fem3d/newnohydro.f', 
             'fem3d/mod_offline.f', 
             'fem3d/suboff.f', 
             'fem3d/subfemintp.f', 
             'fem3d/subtsfile.f', 
             'fem3d/submeteo2.f', 
             'fem3d/subiso8601.f', 
             'fem3d/ecological_dummy.f', 
             'fem3d/newbfm_dummy.f', 
             'fem3d/subomp_dummy.f', 
             'fem3d/subgrd.f', 
             'fem3d/grdutil.f', 
             'fem3d/subxi.f', 
             'fem3d/subshy.f', 
             'fem3d/subshyutil.f', 
             'fem3d/substrings.f', 
             'fem3d/subbit10.f', 
             'fem3d/shyutil.f', 
             'fem3d/subfloodfill.f', 
             'fem3d/subssed.f', 
             'fem3d/subclo.f', 
             'fem3d/sub2linear.f', 
             'fem3d/subdates.f', 
             'fem3d/subpenta.f', 
             'fem3d/submir.f', 
             'fem3d/mod_arrays.f', 
             'fem3d/subfifo.f', 
             'fem3d/subww3.f', 
             'fem3d/newclose.f', 
             'fem3d/mercury_dummy.f', 
             'fem3d/subtimeutil.f',
             'fem3d/subice.f',
             'fem3d/subicemodel.f90',
             'fem3d/mod_test_zeta.f'
           ] 



_shyfem_exe=executable('shyfem',
           _shyfem_source_files ,
           fortran_args : _shyfem_fortran_args,
           include_directories : _shyfem_include_directories,
           dependencies : _shyfem_dependencies,
           )

subdir('tests_meson')
###############################################################################
message('----------------- preparing shyparts build target ----------------')
###############################################################################

#------------------- metis partitioning -----------------------------
_shyparts_dependencies = []
if _with_mpi
  if get_option('path_metis') != 'none'
     _metis_dir=get_option('path_metis') 
     _shyparts_dependencies += [ _compiler.find_library('metis',dirs:[ _metis_dir+'/lib/' ] )  ]
  else
     _metis_dep = dependency('metis',
                             not_found_message:'metis dependency couldnt be found using pkg-config,' 
                                                +' now trying with find_library as we only need the library path,' 
                                                +'  if it does not work please provide installation path in var "_metis_dir"'
                                                +' or add metis library path to environmental variable "LIBRARY_PATH"',
                             required:false)
     if not _metis_dep.found()
       _metis_dep= _compiler.find_library('metis')
     endif
     _shyparts_dependencies += [ _metis_dep ]
  endif
  _code = '''#include "metis.h" '''
  _compiles = _compiler.compiles(_code, name:'can include metis' )
  if not _compiles
      error('\n Error Cannot find header metis.h')
  endif
  #--------------------------------------------------------------------

  executable('shyparts',
             [
              'fem3d/shyparts_metis.f', 
              'fem3d/basutil.f', 
              'fem3d/basutil_resol.f', 
              'fem3d/basutil_hsigma.f', 
              'fem3d/basutil_bathy.f', 
              'fem3d/basutil_interp.f', 
              'fem3d/basutil_smooth.f', 
              'fem3d/basutil_box.f', 
              'fem3d/basutil_part.f', 
              'fem3d/basutil_custom.f', 
              'fem3d/subnos.f', 
              'fem3d/sub2linear.f', 
              'fem3d/sublnka.f', 
              'fem3d/sublnku.f', 
              'fem3d/links.f', 
              'fem3d/subpar3.f', 
              'fem3d/subnls.f', 
              'fem3d/subsss.f', 
              'fem3d/subround.f', 
              'fem3d/subscn.f', 
              'fem3d/subfil.f', 
              'fem3d/subdef.f', 
              'fem3d/subnsu.f', 
              'fem3d/subssv.f', 
              'fem3d/subbas.f', 
              'fem3d/subsys.f', 
              'fem3d/subapn.f', 
              'fem3d/subnsa.f', 
              'fem3d/subdum.f', 
              'fem3d/subdts.f', 
              'fem3d/subdep.f', 
              'fem3d/debug.f', 
              'fem3d/subreg.f', 
              'fem3d/subnev.f', 
              'fem3d/subfemfile.f', 
              'fem3d/subfemutil.f', 
              'fem3d/subtsfile.f', 
              'fem3d/subvintp.f', 
              'fem3d/subver.f', 
              'fem3d/sigmautil.f', 
              'fem3d/suboutput.f', 
              'fem3d/suboutputd.f', 
              'fem3d/subclo.f', 
              'fem3d/mod_depth.f', 
              'fem3d/mod_hydro.f', 
              'fem3d/mod_hydro_print.f', 
              'fem3d/mod_geom.f', 
              'fem3d/mod_geom_dynamic.f', 
              'fem3d/mod_bound_geom.f', 
              'fem3d/mod_bnd.f', 
              'fem3d/mod_hydro_baro.f', 
              'fem3d/newlevels.f', 
              'fem3d/subgrd.f', 
              'fem3d/grdutil.f', 
              'fem3d/subsrt.f', 
              'fem3d/subgeo.f', 
              'fem3d/substrings.f', 
              'fem3d/subiso8601.f', 
              'fem3d/submpi_ghost.f', 
              'fem3d/submpi_basin.f', 
              'fem3d/submpi_custom.f', 
              'fem3d/submpi_buffer.f', 
              'fem3d/submpi_time.f', 
              'fem3d/shympi_general.f', 
              'fem3d/shympi_debug.f', 
              'fem3d/subfloodfill.f',
             'fem3d/mod_shympi_dummy.f', 
             ],
             dependencies : _shyparts_dependencies 
            )
endif #with mpi
###############################################################################
message('------------------- preparing shypart build target ---------------')
###############################################################################
  
executable('shypart',
           [
            'fem3d/shypart.f', 
            'fem3d/basutil.f', 
            'fem3d/basutil_resol.f', 
            'fem3d/basutil_hsigma.f', 
            'fem3d/basutil_bathy.f', 
            'fem3d/basutil_interp.f', 
            'fem3d/basutil_smooth.f', 
            'fem3d/basutil_box.f', 
            'fem3d/basutil_part.f', 
            'fem3d/basutil_custom.f', 
            'fem3d/subnos.f', 
            'fem3d/sub2linear.f', 
            'fem3d/sublnka.f', 
            'fem3d/sublnku.f', 
            'fem3d/links.f', 
            'fem3d/subpar3.f', 
            'fem3d/subnls.f', 
            'fem3d/subsss.f', 
            'fem3d/subround.f', 
            'fem3d/subscn.f', 
            'fem3d/subfil.f', 
            'fem3d/subdef.f', 
            'fem3d/subnsu.f', 
            'fem3d/subssv.f', 
            'fem3d/subbas.f', 
            'fem3d/subsys.f', 
            'fem3d/subapn.f', 
            'fem3d/subnsa.f', 
            'fem3d/subdum.f', 
            'fem3d/subdts.f', 
            'fem3d/subdep.f', 
            'fem3d/debug.f', 
            'fem3d/subreg.f', 
            'fem3d/subnev.f', 
            'fem3d/subfemfile.f', 
            'fem3d/subfemutil.f', 
            'fem3d/subtsfile.f', 
            'fem3d/subvintp.f', 
            'fem3d/subver.f', 
            'fem3d/sigmautil.f', 
            'fem3d/suboutput.f', 
            'fem3d/suboutputd.f', 
            'fem3d/subclo.f', 
            'fem3d/mod_depth.f', 
            'fem3d/mod_hydro.f', 
            'fem3d/mod_hydro_print.f', 
            'fem3d/mod_geom.f', 
            'fem3d/mod_geom_dynamic.f', 
            'fem3d/mod_bound_geom.f', 
            'fem3d/mod_bnd.f', 
            'fem3d/mod_hydro_baro.f', 
            'fem3d/newlevels.f', 
            'fem3d/subgrd.f', 
            'fem3d/grdutil.f', 
            'fem3d/subsrt.f', 
            'fem3d/subgeo.f', 
            'fem3d/substrings.f', 
            'fem3d/subiso8601.f', 
            'fem3d/submpi_ghost.f', 
            'fem3d/submpi_basin.f', 
            'fem3d/submpi_custom.f', 
            'fem3d/submpi_buffer.f', 
            'fem3d/submpi_time.f', 
            'fem3d/shympi_general.f', 
            'fem3d/shympi_debug.f', 
            'fem3d/subfloodfill.f', 
             'fem3d/mod_shympi_dummy.f', ]
          )

###############################################################################
message('------------------- preparing shybas build target ----------------')
###############################################################################

executable('shybas',
           [
            'fem3d/shybas.f', 
            'fem3d/basutil.f', 
            'fem3d/basutil_resol.f', 
            'fem3d/basutil_hsigma.f', 
            'fem3d/basutil_bathy.f', 
            'fem3d/basutil_interp.f', 
            'fem3d/basutil_smooth.f', 
            'fem3d/basutil_box.f', 
            'fem3d/basutil_part.f', 
            'fem3d/basutil_custom.f', 
            'fem3d/subnos.f', 
            'fem3d/sub2linear.f', 
            'fem3d/sublnka.f', 
            'fem3d/sublnku.f', 
            'fem3d/links.f', 
            'fem3d/subpar3.f', 
            'fem3d/subnls.f', 
            'fem3d/subsss.f', 
            'fem3d/subround.f', 
            'fem3d/subscn.f', 
            'fem3d/subfil.f', 
            'fem3d/subdef.f', 
            'fem3d/subnsu.f', 
            'fem3d/subssv.f', 
            'fem3d/subbas.f', 
            'fem3d/subsys.f', 
            'fem3d/subapn.f', 
            'fem3d/subnsa.f', 
            'fem3d/subdum.f', 
            'fem3d/subdts.f', 
            'fem3d/subdep.f', 
            'fem3d/debug.f', 
            'fem3d/subreg.f', 
            'fem3d/subnev.f', 
            'fem3d/subfemfile.f', 
            'fem3d/subfemutil.f', 
            'fem3d/subtsfile.f', 
            'fem3d/subvintp.f', 
            'fem3d/subver.f', 
            'fem3d/sigmautil.f', 
            'fem3d/suboutput.f', 
            'fem3d/suboutputd.f', 
            'fem3d/subclo.f', 
            'fem3d/mod_depth.f', 
            'fem3d/mod_hydro.f', 
            'fem3d/mod_hydro_print.f', 
            'fem3d/mod_geom.f', 
            'fem3d/mod_geom_dynamic.f', 
            'fem3d/mod_bound_geom.f', 
            'fem3d/mod_bnd.f', 
            'fem3d/mod_hydro_baro.f', 
            'fem3d/newlevels.f', 
            'fem3d/subgrd.f', 
            'fem3d/grdutil.f', 
            'fem3d/subsrt.f', 
            'fem3d/subgeo.f', 
            'fem3d/substrings.f', 
            'fem3d/subiso8601.f', 
            'fem3d/submpi_ghost.f', 
            'fem3d/submpi_basin.f', 
            'fem3d/submpi_custom.f', 
            'fem3d/submpi_buffer.f', 
            'fem3d/submpi_time.f', 
            'fem3d/shympi_general.f', 
            'fem3d/shympi_debug.f',
            'fem3d/mod_shympi_dummy.f',
           ]
          )


###############################################################################
message('------------------- preparing shypre build target ----------------')
###############################################################################
executable('shypre',
           [
            'fem3d/shypre.f', 
            'fem3d/subver.f', 
            'fem3d/subros.f', 
            'fem3d/subcmv.f', 
            'fem3d/subdef.f', 
            'fem3d/subbas.f', 
            'fem3d/subfil.f', 
            'fem3d/subpar3.f', 
            'fem3d/subsrt.f', 
            'fem3d/subsss.f', 
            'fem3d/subround.f', 
            'fem3d/subscn.f', 
            'fem3d/subnsu.f', 
            'fem3d/subapn.f', 
            'fem3d/subnsa.f', 
            'fem3d/subnls.f', 
            'fem3d/subsys.f', 
            'fem3d/subdum.f', 
            'fem3d/subgrd.f', 
            'fem3d/grdutil.f', 
            'fem3d/subclo.f', 
            'fem3d/substrings.f', 
            'fem3d/mod_geom_dynamic.f', 
           ]
          )

