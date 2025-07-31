su2_err_src = ['SU2_ERR.cpp']
if get_option('enable-normal')

  su2_err = executable('SU2_ERR',
                      su2_err_src,
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep, su2_itp_dep],
                      cpp_args :[default_warning_flags, su2_cpp_args])

endif
