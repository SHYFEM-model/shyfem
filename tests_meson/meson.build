_prefix = meson.current_source_dir()
message('- - - - - [ meson configuring running tests ] - - - - - - - - - - ')

_types=['EXPLICIT']
if get_option('use_SPK')
    _types+=['SPK_CPU']
endif

_with_amgx='false'
_with_petsc='false'
if get_option('use_PETSc')
    _types+=['PETSc_CPU']
    _with_petsc='cpu'
    if get_option('use_AmgX')
       _with_amgx='true'
    endif
endif

_num_GPU=get_option('num_GPU')

_run_and_check_success = find_program('run_and_check_success.sh')
_create_basin = find_program('create_basin.sh')
_synch = find_program('synch.sh')

if get_option('use_PETSc') == true or get_option('use_MPI') == true
  _mpi_procs=['0001','0002']
else
  _mpi_procs=['0001']
endif

if _num_GPU > 0
  _types+=['PETSc_GPU']
  if get_option('use_AmgX')==true
    _types+=['PETSc_AMGX']
  endif
endif   

foreach _np : _mpi_procs
  test('create .bas file for '+_np+' mpi procs',
          _create_basin,
          args:[
                _np,
                'mesh.grd',
                _prefix+'/',
               ],
          workdir:meson.current_build_dir(),
          priority:100, 
          is_parallel : true)
endforeach             

test('--------------- barrier insuring .bas files are created --------------',
    _synch,
    args:[],
    workdir:meson.current_build_dir(),
    priority:99, 
    is_parallel : false)
_priority=98

foreach _type : _types
  foreach _np : _mpi_procs
    if _type=='EXPLICIT'
        _impexplicit='EXPLICIT'
    else
        _impexplicit='IMPLICIT'
    endif
    if _type=='EXPLICIT' or _type=='SPK_CPU'
        _case=''
    else
        _case='_'+_type
    endif
    if get_option('test_zeta') == true
      _threshold='1E-08'
      _test_name='run shyfem zeta tests for '+_type+' solver using '+_np+' mpi procs'
    else
      _threshold='none'
      _test_name='run shyfem success tests for '+_type+' solver using '+_np+' mpi procs'
    endif
    if _type=='PETSc_GPU' and  _np != '0001'
      message('only one thread can be used with PETSc_GPU')
    else
      test(_test_name,
          _run_and_check_success,
          args:[
                _np,
                _shyfem_exe.full_path(),
                _impexplicit+_case+'.'+_np+'.str',
                'test_zeta.'+ _np + '.'+_impexplicit,
                'test_zeta.'+ _np + '.'+_impexplicit,
                _prefix+'/',
                _threshold,
               ],
          workdir:meson.current_build_dir(),
          priority:_priority,
          is_parallel : true)
    endif
  endforeach      
  _priority=_priority-1
  if get_option('test_zeta')==true          
     test('--------------- barrier insuring zeta output files wont overwrite --------------',
       _synch,
       args:[],
       workdir:meson.current_build_dir(),
       priority:_priority, 
       is_parallel : false)
      _priority=_priority-1
  endif
endforeach             
