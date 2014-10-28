#!/usr/bin/make
#----------------------------------------------------------------------------------------------------------------------------------
# make init

# shell
SHELL = /bin/bash
# no verbose
$(VERBOSE).SILENT:
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# User options
COMPILER = intel
DEBUG    = no
F03STD   = no
OPTIMIZE = no
OPENMP   = no
BIGEIN   = no
R16P     = no
DEXE     = ~/bin/

.PHONY : DEFAULTRULE
DEFAULTRULE: $(DEXE)XnPatches

.PHONY : help
help:
	@echo
	@echo -e '\033[1;31m Make options of XnPatches code\033[0m'
	@echo
	@echo -e '\033[1;31m Compiler choice: COMPILER=$(COMPILER)\033[0m\033[1m => default\033[0m'
	@echo -e '\033[1;31m  COMPILER=gnu  \033[0m\033[1m => GNU gfortran           \033[0m'
	@echo -e '\033[1;31m  COMPILER=intel\033[0m\033[1m => Intel Fortran         \033[0m'
	@echo
	@echo -e '\033[1;31m Compiling options\033[0m'
	@echo -e '\033[1;31m  DEBUG=yes(no)   \033[0m\033[1m => on(off) debug                  (default $(DEBUG))\033[0m'
	@echo -e '\033[1;31m  F03STD=yes(no)  \033[0m\033[1m => on(off) check standard fortran (default $(F03STD))\033[0m'
	@echo -e '\033[1;31m  OPTIMIZE=yes(no)\033[0m\033[1m => on(off) optimization           (default $(OPTIMIZE))\033[0m'
	@echo -e '\033[1;31m  OPENMP=yes(no)  \033[0m\033[1m => on(off) OpenMP directives      (default $(OPENMP))\033[0m'
	@echo -e '\033[1;31m  BIGEIN=yes(no)  \033[0m\033[1m => on(off) Big Endian input files (default $(BIGEIN))\033[0m'
	@echo
	@echo -e '\033[1;31m Preprocessing options\033[0m'
	@echo -e '\033[1;31m  R16P=yes(no)\033[0m\033[1m => on(off) definition of real with "128 bit" (default $(R16P))\033[0m'
	@echo
	@echo -e '\033[1;31m Executable directory\033[0m'
	@echo -e '\033[1;31m  DEXE="your_path"\033[0m\033[1m => directory where exe is placed (default $(DEXE))\033[0m'
	@echo
	@echo -e '\033[1;31m External libraries\033[0m'
	@echo -e '\033[1;31m  TECIO=yes(no)\033[0m\033[1m => on(off) Tecplot IO library linking (default $(TECIO))\033[0m'
	@echo
	@echo -e '\033[1;31m Provided Rules\033[0m'
	@echo -e '\033[1;31m  Defualt rule    =>\033[0m\033[1m $(DEXE)XnPatches\033[0m'
	@echo -e '\033[1;31m  help            =>\033[0m\033[1m printing this help message\033[0m'
	@echo -e '\033[1;31m  $(DEXE)XnPatches =>\033[0m\033[1m building OFF code\033[0m'
	@echo -e '\033[1;31m  cleanobj        =>\033[0m\033[1m cleaning compiled object\033[0m'
	@echo -e '\033[1;31m  cleanmod        =>\033[0m\033[1m cleaning .mod files\033[0m'
	@echo -e '\033[1;31m  cleanmsg        =>\033[0m\033[1m cleaning make-log massage files\033[0m'
	@echo -e '\033[1;31m  cleanexe        =>\033[0m\033[1m cleaning executable files\033[0m'
	@echo -e '\033[1;31m  clean           =>\033[0m\033[1m running cleanobj, cleanmod and cleanmsg\033[0m'
	@echo -e '\033[1;31m  cleanall        =>\033[0m\033[1m running clean and cleanexe\033[0m'
	@echo -e '\033[1;31m  tar             =>\033[0m\033[1m creating a tar archive of the project\033[0m'
	@echo
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# directory & file
DSRC  = ./src/
DOBJ  = ./obj/
DMOD  = ./mod/
DLIB  = ./lib/
VPATH = $(DSRC) $(DOBJ) $(DMOD) $(DLIB)

MKDIRS = $(DOBJ) $(DMOD) $(DEXE)
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
  LIBS = $(DLIB)64bit/tecio64.a $(DLIB)64bit/libstdc++64.5.0.7.so
else
  LIBS = $(DLIB)32bit/tecio.a $(DLIB)32bit/libstdc++.5.0.7.so
endif
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# compiling and linking options
# compiler specific rules
# GNU
WRN_GNU = -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized
CHK_GNU = -fcheck=all
DEB_GNU = -fmodule-private -ffree-line-length-132 -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan #-fno-range-check  ,precision,denormal,underflow
STD_GNU = -std=f2003 -fall-intrinsics
OMP_GNU = -fopenmp
OPT_GNU = -O3
BGE_GNU =
PRF_GNU =
# Intel
WRN_INT = -warn all
CHK_INT = -check all
DEB_INT = -debug all -extend-source 132 -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -traceback -gen-interfaces
STD_INT = -std03
OMP_INT = -openmp
OPT_INT = -O3 -ipo -ipo-jobs4 -vec-report1
BGE_INT = -convert big_endian
PRF_INT = #-p
# setting rules according user options
ifeq "$(COMPILER)" "gnu"
  FC = gfortran
  OPTSC = -cpp -c -J$(DMOD) -static -fprotect-parens -fno-realloc-lhs
  OPTSL =
  WRN = $(WRN_GNU)
  CHK = $(CHK_GNU)
  DEB = $(DEB_GNU)
  STD = $(STD_GNU)
  OMP = $(OMP_GNU)
  OPT = $(OPT_GNU)
  BGE = $(BGE_GNU)
  PRF = $(PRF_GNU)
endif
ifeq "$(COMPILER)" "intel"
  FC = ifort
  OPTSC = -cpp -c -module $(DMOD) -static#-assume protect_parens -assume norealloc_lhs -fp-model source
  OPTSL =
  WRN = $(WRN_INT)
  CHK = $(CHK_INT)
  DEB = $(DEB_INT)
  STD = $(STD_INT)
  OMP = $(OMP_INT)
  OPT = $(OPT_INT)
  BGE = $(BGE_INT)
  PRF = $(PRF_INT)
endif
ifeq "$(DEBUG)" "yes"
  PREPROC := $(PREPROC) -DDEBUG
  OPTSC := $(OPTSC) -O0 -C -g $(WRN) $(CHK) $(DEB)
  OPTSL := $(OPTSL) -O0 -C -g $(WRN) $(CHK) $(DEB)
endif
ifeq "$(F03STD)" "yes"
  OPTSC := $(OPTSC) $(STD)
  OPTSL := $(OPTSL) $(STD)
endif
ifeq "$(OPTIMIZE)" "yes"
  OPTSC := $(OPTSC) $(OPT)
  OPTSL := $(OPTSL) $(OPT)
endif
ifeq "$(OPENMP)" "yes"
  PREPROC := $(PREPROC) -DOPENMP
  OPTSC := $(OPTSC) $(OMP)
  OPTSL := $(OPTSL) $(OMP)
endif
ifeq "$(BIGEIN)" "yes"
  OPTSC := $(OPTSC) $(BGE)
  OPTSL := $(OPTSL) $(BGE)
endif
# pre-processing options
# R16 precision
R16PCHK = (Unknown R16P switch) Used default R16P=no
ifeq "$(R16P)" "no"
  R16PCHK = (Known R16P switch) Used R16P=$(R16P)
endif
ifeq "$(R16P)" "yes"
  R16PCHK = (Known R16P switch) Used R16P=$(R16P)
  PREPROC := $(PREPROC) -Dr16p
endif
OPTSC := $(OPTSC) $(PREPROC)
OPTSL := $(OPTSL) $(PREPROC)

WHICHFC = $(shell which $(FC))

PRINTCHK = "\\033[1;31m Compiler used \\033[0m\\033[1m $(COMPILER) => $(WHICHFC)\\033[0m \n\
            \\033[1;31mSource dir    \\033[0m\\033[1m $(DSRC)\\033[0m \n\
            \\033[1;31mLibraries     \\033[0m\\033[1m $(LIBS)\\033[0m \n \n\
            \\033[1;31m Debug         \\033[0m\\033[1m $(DEBUG)\\033[0m \n\
            \\033[1;31m F-standard    \\033[0m\\033[1m $(F03STD)\\033[0m \n\
            \\033[1;31m Optimize      \\033[0m\\033[1m $(OPTIMIZE)\\033[0m \n\
            \\033[1;31m OpenMP        \\033[0m\\033[1m $(OPENMP)\\033[0m \n\
            \\033[1;31m R16P          \\033[0m\\033[1m $(R16PCHK)\\033[0m\n\
	    \\033[1;31m DEXE          \\033[0m\\033[1m $(DEXE)\\033[0m"
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
.PHONY : PRINTINFO
.NOTPARALLEL : PRINTINFO
PRINTINFO:
	@echo | tee make.log
	@echo -e $(PRINTCHK) | tee -a make.log
	@echo | tee -a make.log
	@echo -e '\033[1;31m Compiling options\033[0m' | tee -a make.log
	@echo -e '\033[1m [$(OPTSC)]\033[0m' | tee -a make.log
	@echo | tee -a make.log
	@echo -e '\033[1;31m Linking options \033[0m' | tee -a make.log
	@echo -e '\033[1m [$(OPTSL)]\033[0m' | tee -a make.log
	@echo | tee -a make.log

.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@

.PHONY : cleanobj
cleanobj:
	@echo -e "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(DOBJ)

.PHONY : cleanmod
cleanmod:
	@echo -e "\033[1;31m deleting mods \033[0m" | tee -a make.log
	@rm -fr $(DMOD)

.PHONY : cleanexe
cleanexe:
	@echo -e "\033[1;31m deleting exes \033[0m" | tee -a make.log
	@rm -f $(addprefix $(DEXE),$(EXES))

.PHONY : cleanmsg
cleanmsg:
	@rm -f diagnostic_messages
	@rm -f error_messages

.PHONY : clean
clean: cleanobj cleanmod cleanmsg

.PHONY : cleanall
cleanall: clean cleanexe

.PHONY : tar
tar: cleanall
	@echo -e "\033[1;31m Creating tar archive of the code \033[0m" | tee make.log
	@rm -rf XnPatches
	@mkdir -p XnPatches
	@cp -rL lib util src makefile XnPatches/
	@tar czf XnPatches.tgz XnPatches
	@rm -rf XnPatches
#----------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------
# rules of linking and compiling
COTEXT  = -e "\033[1;31m Compiling\033[0m\033[1m $(<F)\033[0m"
LITEXT  = -e "\033[1;31m Assembling\033[0m\033[1m $@\033[0m"
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

$(DEXE)XnPatches : PRINTINFO $(MKDIRS) $(DOBJ)xnpatches.o
	@rm -f $(filter-out $(DOBJ)xnpatches.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) XnPatches

$(DOBJ)block_variables.o: src/Block_Variables.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_command_line_interface.o: src/Data_Type_Command_Line_Interface.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_os.o : Data_Type_OS.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_postprocess.o: src/Data_Type_PostProcess.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_command_line_interface.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_vector.o : Data_Type_Vector.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ir_precision.o : IR_Precision.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_base64.o : Lib_Base64.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_io_misc.o : Lib_IO_Misc.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_os.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_tec.o: src/Lib_TEC.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)block_variables.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_vtk_io.o : Lib_VTK_IO.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_base64.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)xnpatches.o : XnPatches.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)block_variables.o \
	$(DOBJ)data_type_command_line_interface.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)lib_tec.o \
	$(DOBJ)lib_vtk_io.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#-----------------------------------------------------------------------------------------------------------------------------------
