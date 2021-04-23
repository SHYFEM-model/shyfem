###############################################################################
message('- - - - - [ preparing dependency libamgx_c_wrapper ] - - - - -')
###############################################################################
_i = include_directories(['src/',_amgx_inc,_amgxwrapper_inc])

_l = library('amgx_c_wrapper',
             [
             ],
             dependencies:[_petsc_dep,_amgx_dep,_amgxwrapper_dep],
             include_directories: _i)

_amgx_c_wrapper_dep = declare_dependency(include_directories : _i,
                              link_with: _l)

