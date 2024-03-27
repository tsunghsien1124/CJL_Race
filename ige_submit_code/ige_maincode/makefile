#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
FC = ifort               						# fortran compiler
FC_MPI = mpif90               					# fortran compiler MPI
FFLAGS = -O3 -check bounds -check pointers -check uninit -u -zero -unroll		# flags: more flags: -mcmodel medium -shared-intel
FFLAGS_DEBUG = -g -check bounds					# ./debug.exe

#------------------------------------------------------------------------------
# directories
#------------------------------------------------------------------------------
FORLIB = /home/sylee/now/fortranlib
LIBDIR = lib
OBJDIR = obj

#------------------------------------------------------------------------------
# program name and files: only need to change this part.
#------------------------------------------------------------------------------
PROGRAM = ige

# fortranlib files
LIB =	nrutil.f90 \
	spcfcn.f90 \
	random.f90 \
	matinv.f90 \
	bfgs_pqn.f90 \
	simplex.f90 \
	brent.f90 \
	rouwenhorst.f90 \
	ols.f90 \
	qsort_module.f90 \

#	newuoa.f90 \
#	sort.f90 \
#	bobyqa_tim.f90 \


# main files
SRC =	global.f90 \
	functions.f90 \
	printresults.f90 \
	grids.f90 \
	old_parent.f90 \
	mid_parent.f90 \
	sec_parent.f90 \
	pri_parent.f90 \
	new_parent.f90 \
	not_parent.f90 \
	valuefuncs.f90 \
	simulation.f90 \
	compute_iter.f90 \
	distance_func.f90 \
	lifetime_stats.f90 \
	bt_module.f90 \
	polex.f90 \
	main.f90

#------------------------------------------------------------------------------
# lib.o and obj files
#------------------------------------------------------------------------------
LIB_O =	$(patsubst $(FORLIB)/%.f90,$(LIBDIR)/%.o,$(addprefix /home/sylee/now/fortranlib/,$(LIB)))
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.f90=.o)))

#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------
# 1. separate library.
library: buildlib $(LIB_O)
$(LIBDIR)/%.o: $(FORLIB)/%.f90
	$(FC_MPI) -c $< $(FFLAGS) -o $@ -module $(LIBDIR)

# 2. main code
$(PROGRAM): buildobj $(OBJ)
	$(FC_MPI) $(LIB_O) $(OBJ) $(FFLAGS) -o $@

$(OBJDIR)/%.o: %.f90
	$(FC_MPI) -c $< $(FFLAGS) -o $@ -module $(OBJDIR) -I$(LIBDIR)

# 3. can do all
all: library $(PROGRAM)

# 4. clean clutter, if any
clean:
	rm -rf *.mod *.o *.exe *.ipo *.il *~ \#* PI* $(PROGRAM) $(OBJDIR)

# 5. clean everything
cleanall: clean
	rm -rf $(LIBDIR) $(OBJDIR)

#------------------------------------------------------------------------------
# lib and obj directory
#------------------------------------------------------------------------------
buildlib:
	@$(call make-lib)
buildobj:
	@$(call make-obj)

define make-lib
	rm -rf $(LIBDIR)
	mkdir -p $(LIBDIR)
	cd $(LIBDIR)
endef
define make-obj
	rm -rf $(OBJDIR)
	mkdir -p $(OBJDIR)
	cd $(OBJDIR)
endef
