#	File name		: Makefile
#	Date			: December 2022
#	Version			: 3.0.0
#	Author			: Adhika
#	Co-Author		: Angga

DEST			= .

COMPILER		= g++

OPT_STD 		= -std=c++17

OPT_O 			= -O2 # -O3: newest optimizer but can be slower

OPT_WARNING		= -Wall

OPT_OMP			= -fopenmp

OPT_MODEL		= #-mcmodel=medium # For a large dataCurrently no model

OPTFLAG			= $(OPT_STD) $(OPT_O) $(OPT_WARNING) $(OPT_OMP) $(OPT_MODEL)
 
MAKEFILE		= Makefile

PROGRAM			= main

COREs			= 20	# Number of core use for parallel programming

# List of .cpp source file path
SRCS			= main.cpp\
				global.cpp\
				Utils.cpp\
				src/save_data/summary_log.cpp\
				src/save_data/state_data.cpp\
				src/geometry/2D_objects_generator.cpp\
				src/geometry/3D_objects_generator.cpp\
				src/geometry/geometry_utils.cpp\
				src/geometry/def_vel_var.cpp\
				src/geometry/moving_body.cpp\
				src/geometry/geometry.cpp\
				src/grid_block/gridNode.cpp\
				src/grid_block/generateGrid.cpp\
				src/grid_block/gridNodeNgh.cpp\
				src/grid_block/gridNodeAdapt.cpp\
				src/grid_block/gridNodeAdaptUtils.cpp\
				src/initialization/initialization.cpp\
				src/initialization/init2D.cpp\
				src/initialization/init3D.cpp\
				src/neighbor/neighbor_utilities/base_grid.cpp\
				src/neighbor/direct/direct_find.cpp\
				src/neighbor/link_list/link_list.cpp\
				src/neighbor/spatial_hash/spatial_hash.cpp\
				src/neighbor/inter_search/inter_search.cpp\
				src/neighbor/cell_list/cell_list_utils.cpp\
				src/neighbor/cell_list/cell_list_init.cpp\
				src/neighbor/cell_list/cell_list_ngh.cpp\
				src/neighbor/cell_list/cell_list_adaptive.cpp\
				src/neighbor/neighbor.cpp\
				src/adaptation/adaptation.cpp\
				src/remeshing/redistribute_particles.cpp\
				src/remeshing/remeshing.cpp\
				src/LSMPS/LSMPSa.cpp\
				src/LSMPS/LSMPSb.cpp\
				src/FMM/treeCell.cpp\
				src/FMM/fast3DFMM.cpp\
				src/FMM/fmm2D.cpp\
				src/FMM/fmm3D.cpp\
				src/velocity_calculation/velocity_calc.cpp\
				src/velocity_calculation/regularization_function_2d.cpp\
				src/velocity_calculation/biotsavart_direct_2d.cpp\
				src/velocity_calculation/biotsavart_fmm_2d.cpp\
				src/velocity_calculation/subroutine_all_fmm_2d.cpp\
				src/pressure_poisson/pressure_poisson.cpp\
				src/penalization/penalization_def.cpp\
				src/penalization/penalization_calc.cpp\
				src/penalization/penalization.cpp\
				src/force_calculation/force_calc.cpp\
				src/force_calculation/force_output.cpp\
				src/force_calculation/force_linear_impulse.cpp\
				src/advection/main_advection.cpp\
				src/diffusion/main_diffusion.cpp\
				src/stretching/main_stretching.cpp\
				# src/solver_poisson/solver.cpp\
				# src/solver_poisson/solver_poisson.cpp\
				# src/diffusion/pse.cpp\
				# src/FMM/treeCellNew.cpp\
				# src/DC_operator/vandermonde.cpp\
				# src/DC_operator/dc_pse.cpp\
				# src/DC_operator/dc_gradient.cpp\
				# src/remeshing/sign_particles.cpp\
				# src/remeshing/basic_remeshing.cpp\
				# src/remeshing/particle_adaption/multiblock_adaption.cpp\
				# src/remeshing/particle_adaption/particle_adaption.cpp\
				# src/remeshing/particle_adaption/imr_removal_insertion.cpp\
				# src/remeshing/particle_adaption/imr_resolution_field.cpp\
				# src/remeshing/particle_adaption/imr_resolution_interp.cpp\
				# src/remeshing/particle_adaption/imr_searching.cpp\
				# src/remeshing/particle_adaption/imr_set_properties.cpp\
				# src/remeshing/particle_adaption/imr_steepest_descent.cpp\

# List of target .o object file path
OBJS			= $(SRCS:.cpp=.o)

# List of target .o object file path
.cpp.o:
			$(COMPILER) $(OPTFLAG) -c $*.cpp -o $*.o

$(PROGRAM):	$(OBJS)
			$(COMPILER) $(OPTFLAG) -o $(PROGRAM) $(OBJS)
			@echo "\n--------------- DONE COMPILE --------------"

# List of all method
all:		$(PROGRAM)

set_core:
			export OMP_NUM_THREADS=$(COREs)

run:
			@echo "--------------- RUN PROGRAM ---------------"
			./$(PROGRAM)

clean:
			@rm -f $(OBJS)
			@rm -f $(PROGRAM)
			@echo "----------------- CLEANED -----------------"

data_clear:
			@while [ -z "$$CONTINUE" ]; do \
			read -r -p "Do you really want to delete data files [Y/N] ? " CONTINUE; \
			done ; \
			if [ $$CONTINUE != "y" ] && [ $$CONTINUE != "Y" ]; then \
			echo "Exiting." ; exit 1 ; \
			fi
			@rm -f output/*.dat output/*.csv
			@rm -f *.dat *.csv
			@echo "-------------- DATA DELETED ---------------"
