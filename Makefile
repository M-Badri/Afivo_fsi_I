EXE = Afivo_fsi_I
FC = gfortran
LD = gfortran
EXT_LIB_DIR = external_libs

IDIR = -Isrc/cpp_macros -I$(EXT_LIB_DIR)/hypre/include -I$(EXT_LIB_DIR)/silo/include 
CFLAGS = -Wall -cpp -fopenmp -w -ffast-math -flto  -J$(OBJS_DIR) $(IDIR)
LFLAGS = -O2 -fopenmp 
LIBS = -L$(EXT_LIB_DIR)/hypre/lib -L$(EXT_LIB_DIR)/silo/lib -lHYPRE -lsilo


SRCS_f90d1 = \
m_af_advance.f90 \
m_af_all.f90 \
m_af_core.f90 \
m_af_flux_schemes.f90 \
m_af_ghostcell.f90 \
m_af_interp.f90 \
m_af_multigrid.f90 \
m_af_output.f90 \
m_af_particles.f90 \
m_af_prolong.f90 \
m_af_restrict.f90 \
m_af_types.f90 \
m_af_utils.f90 \
m_coarse_solver.f90 \
m_dielectric.f90 \
m_mg_types.f90 \
m_vtk.f90 

SRCS_f90d2 = \
m_boundary_condition.f90 

SRCS_f90d3 = \
m_helmholtz.f90 

SRCS_f90d4 = \
immersed_boundary.f90 \
immersed_boundary_external_force.f90 \
immersed_boundary_object.f90 

SRCS_f90d5 = \
incomp_flow_field.f90 \
incomp_flow_field_interfaces.f90 

SRCS_f90d6 = \
initial_condition.f90 

SRCS_f90d7 = \
m_math_lib.f90 

SRCS_f90d8 = \
mls.f90 

SRCS_f90d9 = \
navier_stoeks_solver.f90 

SRCS_f90d10 = \
projection.f90 

SRCS_f90d11 = \
rbf_interp_2d.f90 

SRCS_f90d12 = \
m_solver_parameters.f90 

SRCS_f90d13 = \
predefined_move.f90 

SRCS_f90d14 = \
write.f90 

SRCS_f90d15 = \
main.f90 

SRCS_F90d16 = \
befor64.F90 \
befor64_pack_data_m.F90 

SRCS_F90d17 = \
face.F90 

SRCS_F90d18 = \
penf.F90 \
penf_b_size.F90 \
penf_global_parameters_variables.F90 \
penf_stringify.F90 

SRCS_f90d19 = \
m_config.f90 

SRCS_f90d20 = \
finer_backend.f90 \
finer_file_ini_t.f90 \
finer_section_t.f90 

SRCS_F90d20 = \
finer_option_t.F90 

SRCS_F90d21 = \
stringifor.F90 \
stringifor_string_t.F90 

OBJS_f90d1 = \
m_af_advance.o \
m_af_all.o \
m_af_core.o \
m_af_flux_schemes.o \
m_af_ghostcell.o \
m_af_interp.o \
m_af_multigrid.o \
m_af_output.o \
m_af_particles.o \
m_af_prolong.o \
m_af_restrict.o \
m_af_types.o \
m_af_utils.o \
m_coarse_solver.o \
m_dielectric.o \
m_mg_types.o \
m_vtk.o 

OBJS_f90d2 = \
m_boundary_condition.o 

OBJS_f90d3 = \
m_helmholtz.o 

OBJS_f90d4 = \
immersed_boundary.o \
immersed_boundary_external_force.o \
immersed_boundary_object.o 

OBJS_f90d5 = \
incomp_flow_field.o \
incomp_flow_field_interfaces.o 

OBJS_f90d6 = \
initial_condition.o 

OBJS_f90d7 = \
m_math_lib.o 

OBJS_f90d8 = \
mls.o 

OBJS_f90d9 = \
navier_stoeks_solver.o 

OBJS_f90d10 = \
projection.o 

OBJS_f90d11 = \
rbf_interp_2d.o 

OBJS_f90d12 = \
m_solver_parameters.o 

OBJS_f90d13 = \
predefined_move.o 

OBJS_f90d14 = \
write.o 

OBJS_f90d15 = \
main.o 

OBJS_F90d16 = \
befor64.o \
befor64_pack_data_m.o 

OBJS_F90d17 = \
face.o 

OBJS_F90d18 = \
penf.o \
penf_b_size.o \
penf_global_parameters_variables.o \
penf_stringify.o 

OBJS_f90d19 = \
m_config.o 

OBJS_f90d20 = \
finer_backend.o \
finer_file_ini_t.o \
finer_section_t.o 

OBJS_F90d20 = \
finer_option_t.o 

OBJS_F90d21 = \
stringifor.o \
stringifor_string_t.o 

SRC_DIR_f90d1 = src/afivo_fsi_lib/afivo/

SRC_DIR_f90d2 = src/afivo_fsi_lib/boundary_conditon/

SRC_DIR_f90d3 = src/afivo_fsi_lib/helmholtz/

SRC_DIR_f90d4 = src/afivo_fsi_lib/immersed_boundary/

SRC_DIR_f90d5 = src/afivo_fsi_lib/incomp_field/

SRC_DIR_f90d6 = src/afivo_fsi_lib/initial_condition/

SRC_DIR_f90d7 = src/afivo_fsi_lib/math_lib/

SRC_DIR_f90d8 = src/afivo_fsi_lib/mls/

SRC_DIR_f90d9 = src/afivo_fsi_lib/navier_stokes_solver/

SRC_DIR_f90d10 = src/afivo_fsi_lib/projection/

SRC_DIR_f90d11 = src/afivo_fsi_lib/radial_basis/

SRC_DIR_f90d12 = src/afivo_fsi_lib/solver_parameters/

SRC_DIR_f90d13 = src/afivo_fsi_lib/user_defined_procedures/

SRC_DIR_f90d14 = src/afivo_fsi_lib/write/

SRC_DIR_f90d15 = src/main/

SRC_DIR_F90d16 = src/third_party_lib/BeFor64/

SRC_DIR_F90d17 = src/third_party_lib/FACE/

SRC_DIR_F90d18 = src/third_party_lib/PENF/

SRC_DIR_f90d19 = src/third_party_lib/config/

SRC_DIR_f90d20 = src/third_party_lib/finer/

SRC_DIR_F90d20 = src/third_party_lib/finer/

SRC_DIR_F90d21 = src/third_party_lib/stringfor/
OBJS_DIR = obj/Release/src/afivo_fsi_lib/afivo/
EXE_DIR = bin/Release/


VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR):$(SRC_DIR_f90d2):$(OBJS_DIR):$(SRC_DIR_f90d3):$(OBJS_DIR):$(SRC_DIR_f90d4):$(OBJS_DIR):$(SRC_DIR_f90d5):$(OBJS_DIR):$(SRC_DIR_f90d6):$(OBJS_DIR):$(SRC_DIR_f90d7):$(OBJS_DIR):$(SRC_DIR_f90d8):$(OBJS_DIR):$(SRC_DIR_f90d9):$(OBJS_DIR):$(SRC_DIR_f90d10):$(OBJS_DIR):$(SRC_DIR_f90d11):$(OBJS_DIR):$(SRC_DIR_f90d12):$(OBJS_DIR):$(SRC_DIR_f90d13):$(OBJS_DIR):$(SRC_DIR_f90d14):$(OBJS_DIR):$(SRC_DIR_f90d15):$(OBJS_DIR):$(SRC_DIR_F90d16):$(OBJS_DIR):$(SRC_DIR_F90d17):$(OBJS_DIR):$(SRC_DIR_F90d18):$(OBJS_DIR):$(SRC_DIR_f90d19):$(OBJS_DIR):$(SRC_DIR_f90d20):$(OBJS_DIR):$(SRC_DIR_F90d20):$(OBJS_DIR):$(SRC_DIR_F90d21):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1) $(OBJS_f90d2) $(OBJS_f90d3) $(OBJS_f90d4) $(OBJS_f90d5) $(OBJS_f90d6) $(OBJS_f90d7) $(OBJS_f90d8) $(OBJS_f90d9) $(OBJS_f90d10) $(OBJS_f90d11) $(OBJS_f90d12) $(OBJS_f90d13) $(OBJS_f90d14) $(OBJS_f90d15) $(OBJS_F90d16) $(OBJS_F90d17) $(OBJS_F90d18) $(OBJS_f90d19) $(OBJS_f90d20) $(OBJS_F90d20) $(OBJS_F90d21))

all : $(EXE)

$(EXE) : $(OBJS_f90d1) $(OBJS_f90d2) $(OBJS_f90d3) $(OBJS_f90d4) $(OBJS_f90d5) $(OBJS_f90d6) $(OBJS_f90d7) $(OBJS_f90d8) $(OBJS_f90d9) $(OBJS_f90d10) $(OBJS_f90d11) $(OBJS_f90d12) $(OBJS_f90d13) $(OBJS_f90d14) $(OBJS_f90d15) $(OBJS_F90d16) $(OBJS_F90d17) $(OBJS_F90d18) $(OBJS_f90d19) $(OBJS_f90d20) $(OBJS_F90d20) $(OBJS_F90d21)
	@mkdir -p $(EXE_DIR)
	$(LD) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d2):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d2)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d3):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d3)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d4):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d4)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d5):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d5)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d6):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d6)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d7):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d7)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d8):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d8)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d9):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d9)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d10):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d10)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d11):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d11)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d12):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d12)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d13):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d13)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d14):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d14)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d15):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d15)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_F90d16):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90d16)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_F90d17):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90d17)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_F90d18):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90d18)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_f90d19):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d19)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d20):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d20)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_F90d20):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90d20)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_F90d21):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90d21)$(@:.o=.F90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# File dependencies
m_af_advance.o: \
    m_af_advance.f90 \
    m_af_types.o
m_af_all.o: \
    m_af_all.f90 \
    m_af_advance.o \
    m_af_core.o \
    m_af_flux_schemes.o \
    m_af_ghostcell.o \
    m_af_interp.o \
    m_af_multigrid.o \
    m_af_output.o \
    m_af_particles.o \
    m_af_prolong.o \
    m_af_restrict.o \
    m_af_types.o \
    m_af_utils.o \
    m_coarse_solver.o \
    m_dielectric.o \
    m_mg_types.o
m_af_core.o: \
    m_af_core.f90 \
    m_af_ghostcell.o \
    m_af_prolong.o \
    m_af_restrict.o \
    m_af_types.o \
    m_af_utils.o
m_af_flux_schemes.o: \
    m_af_flux_schemes.f90 \
    m_af_core.o \
    m_af_ghostcell.o \
    m_af_restrict.o \
    m_af_types.o
m_af_ghostcell.o: \
    m_af_ghostcell.f90 \
    m_af_prolong.o \
    m_af_types.o
m_af_interp.o: \
    m_af_interp.f90 \
    m_af_types.o \
    m_af_utils.o
m_af_multigrid.o: \
    m_af_multigrid.f90 \
    m_af_core.o \
    m_af_ghostcell.o \
    m_af_prolong.o \
    m_af_restrict.o \
    m_af_types.o \
    m_af_utils.o \
    m_coarse_solver.o \
    m_mg_types.o
m_af_output.o: \
    m_af_output.f90 \
    m_af_core.o \
    m_af_interp.o \
    m_af_types.o \
    m_vtk.o
m_af_particles.o: \
    m_af_particles.f90 \
    m_af_ghostcell.o \
    m_af_restrict.o \
    m_af_types.o \
    m_af_utils.o
m_af_prolong.o: \
    m_af_prolong.f90 \
    m_af_types.o
m_af_restrict.o: \
    m_af_restrict.f90 \
    m_af_types.o
m_af_types.o: \
    m_af_types.f90
m_af_utils.o: \
    m_af_utils.f90 \
    m_af_types.o
m_coarse_solver.o: \
    m_coarse_solver.f90 \
    m_af_ghostcell.o \
    m_af_types.o \
    m_mg_types.o
m_dielectric.o: \
    m_dielectric.f90 \
    m_af_ghostcell.o \
    m_af_types.o \
    m_af_utils.o
m_mg_types.o: \
    m_mg_types.f90 \
    m_af_types.o
m_vtk.o: \
    m_vtk.f90
m_boundary_condition.o: \
    m_boundary_condition.f90 \
    m_af_all.o \
    m_solver_parameters.o
m_helmholtz.o: \
    m_helmholtz.f90 \
    m_af_all.o \
    m_boundary_condition.o \
    m_config.o \
    incomp_flow_field.o \
    initial_condition.o \
    m_solver_parameters.o
immersed_boundary.o: \
    immersed_boundary.f90 \
    m_af_all.o \
    m_config.o \
    incomp_flow_field.o \
    mls.o \
    predefined_move.o \
    rbf_interp_2d.o \
    m_solver_parameters.o \
    write.o
immersed_boundary_external_force.o: \
    immersed_boundary_external_force.f90 \
    immersed_boundary.o
immersed_boundary_object.o: \
    immersed_boundary_object.f90 \
    immersed_boundary.o
incomp_flow_field.o: \
    incomp_flow_field.f90 \
    m_af_prolong.o \
    m_boundary_condition.o \
    m_config.o \
    initial_condition.o \
    m_solver_parameters.o
incomp_flow_field_interfaces.o: \
    incomp_flow_field_interfaces.f90 \
    m_af_all.o \
    m_config.o
initial_condition.o: \
    initial_condition.f90 \
    m_af_all.o \
    m_solver_parameters.o
m_math_lib.o: \
    m_math_lib.f90
mls.o: \
    mls.f90 \
    m_math_lib.o
navier_stoeks_solver.o: \
    navier_stoeks_solver.f90 \
    face.o \
    m_af_all.o \
    m_config.o \
    m_helmholtz.o \
    immersed_boundary.o \
    incomp_flow_field.o \
    projection.o \
    m_solver_parameters.o
projection.o: \
    projection.f90 \
    m_af_all.o \
    m_boundary_condition.o \
    m_config.o \
    incomp_flow_field.o \
    initial_condition.o \
    m_solver_parameters.o
rbf_interp_2d.o: \
    rbf_interp_2d.f90 \
    m_math_lib.o
m_solver_parameters.o: \
    m_solver_parameters.f90 \
    m_af_all.o
predefined_move.o: \
    predefined_move.f90 \
    m_af_all.o \
    m_solver_parameters.o
write.o: \
    write.f90 \
    face.o \
    m_af_all.o \
    m_solver_parameters.o
main.o: \
    main.f90 \
    m_af_all.o \
    m_config.o \
    incomp_flow_field.o \
    navier_stoeks_solver.o
befor64.o: \
    befor64.F90 \
    befor64_pack_data_m.o \
    penf.o
befor64_pack_data_m.o: \
    befor64_pack_data_m.F90 \
    penf.o
face.o: \
    face.F90
penf.o: \
    penf.F90 \
    penf_b_size.o \
    penf_global_parameters_variables.o \
    penf_stringify.o
penf_b_size.o: \
    penf_b_size.F90 \
    penf_global_parameters_variables.o
penf_global_parameters_variables.o: \
    penf_global_parameters_variables.F90
penf_stringify.o: \
    penf_stringify.F90 \
    penf_b_size.o \
    penf_global_parameters_variables.o
m_config.o: \
    m_config.f90 \
    finer_backend.o \
    finer_file_ini_t.o
finer_backend.o: \
    finer_backend.f90 \
    penf.o
finer_file_ini_t.o: \
    finer_file_ini_t.f90 \
    finer_backend.o \
    finer_option_t.o \
    finer_section_t.o \
    penf.o \
    stringifor.o
finer_option_t.o: \
    finer_option_t.F90 \
    finer_backend.o \
    penf.o \
    stringifor.o
finer_section_t.o: \
    finer_section_t.f90 \
    finer_backend.o \
    finer_option_t.o \
    penf.o \
    stringifor.o
stringifor.o: \
    stringifor.F90 \
    penf.o \
    stringifor_string_t.o
stringifor_string_t.o: \
    stringifor_string_t.F90 \
    befor64.o \
    penf.o

