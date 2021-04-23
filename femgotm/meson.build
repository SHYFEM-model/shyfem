###############################################################################
message('- - - - - [ preparing dependency libgotm build target ] - - - - -')
###############################################################################
_i = include_directories('.') # list of the object returned by the function include directories

_l = library('gotm',
             [
                'tridiagonal.F90',
                'eqstate.F90',
                'util.F90',
                'turbulence.F90',
                'kpp.F90',
                'diff_face.F90',
                'fem_gotm_interface.F90',
                'gotm_lib_version.F90',
                'production.F90',
                'alpha_mnb.F90',
                'r_ratio.F90',
                'tkeeq.F90',
                'q2over2eq.F90',
                'tkealgebraic.F90',
                'kbeq.F90',
                'kbalgebraic.F90',
                'dissipationeq.F90',
                'genericeq.F90',
                'lengthscaleeq.F90',
                'algebraiclength.F90',
                'epsbalgebraic.F90',
                'variances.F90',
                'ispralength.F90',
                'potentialml.F90',
                'cmue_a.F90',
                'cmue_b.F90',
                'cmue_c.F90',
                'cmue_d.F90',
                'cmue_ma.F90',
                'cmue_rf.F90',
                'cmue_sg.F90',
                'compute_cpsi3.F90',
                'compute_rist.F90',
                'fk_craig.F90',
                'internal_wave.F90',
             ],
             include_directories: _i)

_gotm_dep = declare_dependency(include_directories : _i,
                              link_with: _l)

