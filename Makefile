EXEC := 2D_NS_Solver
EXEC_DIRECTORY := ./obj/
EXECUTABLES := 2D_NS_Solver

BASE_PATH := ./
SOURCE_PATH := ./src/
#IPATH  := ./
#LPATH  := ./
#LIBPATH := /usr/local/lib
#LIBPATH := 
#LIBS := -lm
IFLAGS := -I $(IPATH)
LFLAGS :=  $(LIBS) -L $(LIBPATH) #-lgsl -lgslcblas -lm -l fnccheck -L$(CFDRC)/lib -lDTF
DFLAGS := #-ggdb3 # -g
OFLAGS := #-O3 #-fno-strict-aliasing -funroll-loops -Wstrict-aliasing #-O5 #-O3  -O2 -O
WFLAGS := -Wall #-Wno-long-long
EFLAGS := #-shared-intel #-pg
AFLAGS := # -m64
CFLAGS := $(EFLAGS) -c $(AFLAGS) $(DFLAGS) $(OFLAGS) $(WFLAGS)
CC     := g++ # i586-mingw32msvc-g++
STRIP  := #strip # i586-mingw32msvc-strip

OBJS := $(BASE_PATH)$(SOURCE_PATH)avg_conv_flux1.o\
	$(BASE_PATH)$(SOURCE_PATH)avg_conv_flux2.o\
	$(BASE_PATH)$(SOURCE_PATH)AvgFlux1.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_cut.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_eulerwall.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_farfiled.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_inflow.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_prescribed_inflow.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_symmetry.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_transmitive.o\
	$(BASE_PATH)$(SOURCE_PATH)bc_wall.o\
	$(BASE_PATH)$(SOURCE_PATH)boundary_conditions.o\
	$(BASE_PATH)$(SOURCE_PATH)conv_var_diffs.o\
	$(BASE_PATH)$(SOURCE_PATH)dependent_variables.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_llf1.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_llf2_prim.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers2_prim.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_h1.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_h2_conv.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_h2_prim.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_le1.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_le2_conv.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_le2_prim.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_movers_nwsc.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_roe1.o\
	$(BASE_PATH)$(SOURCE_PATH)diss_roe2_prim.o\
	$(BASE_PATH)$(SOURCE_PATH)ECCS.o\
	$(BASE_PATH)$(SOURCE_PATH)flux_boundary.o\
	$(BASE_PATH)$(SOURCE_PATH)flux_conv_defns.o\
	$(BASE_PATH)$(SOURCE_PATH)flux_viscous.o\
	$(BASE_PATH)$(SOURCE_PATH)flux_wall.o\
	$(BASE_PATH)$(SOURCE_PATH)forces.o\
	$(BASE_PATH)$(SOURCE_PATH)grad_face.o\
	$(BASE_PATH)$(SOURCE_PATH)grid_computations.o\
	$(BASE_PATH)$(SOURCE_PATH)ic_steady_flow.o\
	$(BASE_PATH)$(SOURCE_PATH)ic_unsteady_flow.o\
	$(BASE_PATH)$(SOURCE_PATH)imp_res_smoo.o\
	$(BASE_PATH)$(SOURCE_PATH)ini_flow_dimensional.o\
	$(BASE_PATH)$(SOURCE_PATH)iniflow_nondimen.o\
	$(BASE_PATH)$(SOURCE_PATH)initialize_grads.o\
	$(BASE_PATH)$(SOURCE_PATH)initialize_variables.o\
	$(BASE_PATH)$(SOURCE_PATH)kfds_1st_order.o\
	$(BASE_PATH)$(SOURCE_PATH)kfds_2nd_order.o\
	$(BASE_PATH)$(SOURCE_PATH)limiter_reference.o\
	$(BASE_PATH)$(SOURCE_PATH)main.o\
	$(BASE_PATH)$(SOURCE_PATH)read_grid.o\
	$(BASE_PATH)$(SOURCE_PATH)read_grid_top.o\
	$(BASE_PATH)$(SOURCE_PATH)read_restart.o\
	$(BASE_PATH)$(SOURCE_PATH)read_solver_input.o\
	$(BASE_PATH)$(SOURCE_PATH)residue_cal.o\
	$(BASE_PATH)$(SOURCE_PATH)solver.o\
	$(BASE_PATH)$(SOURCE_PATH)time_step.o\
	$(BASE_PATH)$(SOURCE_PATH)var_diffs.o\
	$(BASE_PATH)$(SOURCE_PATH)write_solution.o\
	

$(EXEC): $(OBJS)
	$(CC) $(EFLAGS) $(AFLAGS) $(OBJS)  $(LFLAGS) $(IFLAGS) -o $(EXEC_DIRECTORY)$@

%.o: %.cpp
	$(CC) $(CFLAGS) `pwd`/$< -o $@
#               @-makedepend -I$(IPATH) -I$(LPATH) -a -fMake_dependencies $< > /dev/null 2>&1

clean:
	rm -f $(OBJS)
	rm -f $(EXEC_DIRECTORY)$(EXEC)
#               @-touch Make_dependencies
#               @-makedepend -fMake_dependencies __arbitrary_file > /dev/null 2>&1
archive:
	ar rcs $(EXEC_DIRECTORY)$(EXEC).a $(OBJS)

#-include Make_dependencies

