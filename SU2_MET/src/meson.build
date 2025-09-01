su2_met_src = ['SU2_MET.cpp']
if get_option('enable-normal')

  su2_met = executable('SU2_MET',
                      su2_met_src,
                      install: true,
		                  dependencies: [su2_deps, common_dep, su2_cfd_dep],
		                  cpp_args : [default_warning_flags, su2_cpp_args])

endif
