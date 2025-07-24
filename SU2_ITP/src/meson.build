su2_itp_src = ['SU2_ITP.cpp']
if get_option('enable-normal')

  su2_itp = executable('SU2_ITP',
                      su2_itp_src,
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep],
                      cpp_args :[default_warning_flags, su2_cpp_args])

endif
