su2_itp_src = ['SU2_ITP.cpp']
itp_lib_src = ['interpolation.cpp']
if get_option('enable-normal')

  su2_itp_lib = static_library('su2_itp', itp_lib_src,
                              dependencies: [su2_deps, common_dep, su2_cfd_dep],
                              cpp_args: [default_warning_flags, su2_cpp_args])
  su2_itp_dep = declare_dependency(link_with: su2_itp_lib,
                                  include_directories: include_directories('..', '../include'))

  su2_itp = executable('SU2_ITP',
                      su2_itp_src + itp_lib_src,
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep],
                      cpp_args :[default_warning_flags, su2_cpp_args])

endif
