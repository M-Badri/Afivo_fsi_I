!< FiNeR, Fortran INI ParseR and generator.

module m_config
!< FiNeR, Fortran INI ParseR and generator.
use finer_backend
use finer_file_ini_t

implicit none
private
public :: err_option_name
public :: err_option_vals
public :: err_option
public :: err_section_name
public :: err_section_options
public :: err_section
public :: err_source_missing
public :: t_config
public :: t_config_autotest
endmodule m_config
