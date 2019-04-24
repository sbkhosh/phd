###############################################################################
###                                                                         ###
###                         Master makefile                                 ###
###                                                                         ###
###############################################################################

#--- Debug --------------------------------------------------------------------

DEBUG = no

##############################################################################

ifeq ($(MACHINE),)
 ifeq ($(findstring eslogin,$(HOSTNAME)),eslogin)
  MACHINE = archer
 endif
 ifeq ($(findstring polaris,$(HOSTNAME)),polaris)
  MACHINE = polaris
 endif
 ifeq ($(MACHINE),)
  MACHINE = new_machine
 endif
endif

##############################################################################

#check if user wants to compile inc3d for the present prob directory only
pathsearch = $(firstword $(wildcard $(addsuffix /$(1),$(subst :, ,$(PWD)))))
MK_LOCAL:= $(wildcard make_local)
ifeq ($(MK_LOCAL),make_local)
   #executable should be local
   BASEDIR = $(PWD)/..
   EXECDIR  = $(PWD)
else
   pathsearch = $(firstword $(wildcard $(addsuffix /$(1),$(subst :, ,$(PWD)))))
   CHECKBASE := $(call pathsearch,src)
   ifeq ($(CHECKBASE),)
      #call from prob directory
      BASEDIR = $(PWD)/..
   else
      #call from inc3d base directory
      BASEDIR = $(PWD)
   endif
   EXECDIR  = $(BASEDIR)/bin
endif

#set dependend paths
MKFILE  = $(EXECDIR)/$(MACHINE).mk

RUNDIR = $(PWD)
SRCDIR  = $(BASEDIR)/src
EXTDIR  = $(BASEDIR)/external
OBJDIR  = $(EXECDIR)/obj_$(MACHINE)
FFILES  = $(SRCDIR)/files.mk

EXEC = inc3d_$(MACHINE)
LIBN = libinc3d

export DEBUG EXEC LIBN FFILES BASEDIR SRCDIR OBJDIR\
       EXECDIR RUNDIR EXTDIR NP_OS SVNDEF MAKE MKFILE MACHINE

##############################################################################

all: $(EXEC)

checkpath:
	@test -d $(EXECDIR) || mkdir -p $(EXECDIR)
	@test -f $(MKFILE) || cp $(BASEDIR)/makefiles/$(MACHINE)/*.mk $(EXECDIR)
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
ifneq ($(RUNDIR),$(BASEDIR))
#copy submit scripts for calls from prob directory
	@diff $(BASEDIR)/makefiles/$(MACHINE) $(RUNDIR) | grep -e 'makefiles/$(MACHINE)' | \
	  grep -v -e svn -e '.mk' | awk -F ': ' '{print $$2}' > list
	@for file in `cat list`; do cp $(BASEDIR)/makefiles/$(MACHINE)/$$file $(RUNDIR) ; done
	@rm list
endif

ifeq ($(MK_LOCAL),make_local)
	@echo compiling a INCOMPACT3D executable in $(RUNDIR)
	@test ! -L $(RUNDIR)/$(EXEC) || (rm -f $(RUNDIR)/$(EXEC))
else
	@echo using INCOMPACT3D executable in bin $(RUNDIR)
ifneq ($(RUNDIR),$(BASEDIR))
	@test ! -d $(RUNDIR)/$(MACHINE) || (echo removing $(RUNDIR)/$(MACHINE) directory && \
	   rm -r $(RUNDIR)/$(MACHINE))
	@test ! -d $(RUNDIR)/obj_$(MACHINE) || (echo removing $(RUNDIR)/obj_$(MACHINE) directory && \
	   rm -r $(RUNDIR)/obj_$(MACHINE))
	@rm -f $(RUNDIR)/$(EXEC)
	@ln -s $(EXECDIR)/$(EXEC) $(RUNDIR)
endif
endif

#------- Create executable ----- #
$(EXEC): checkpath
	@echo Calling gmake in subdir
	$(MAKE) -f $(MKFILE) $(EXECDIR)/$(EXEC)

#------- Create INCOMPACT3D library ---- #
lib:	checkpath
	@echo Calling gmake in subdir $(ADDSTR)
	$(MAKE) -f $(MKFILE) $(EXECDIR)/$(LIBN).a

run: $(EXEC)
	$(MAKE) -f $(MKFILE) run

#-----------------------------------------------------------------------------

.PHONY: clean distclean mach_wrapper install run lib all checkpath

clean:
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.a
	-rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.flc
	-rm -f $(EXECDIR)/$(EXEC) $(EXECDIR)/*.a
	-rm -f opari.rc opari.tab.c *.opari.inc opari.tab.o
	-rm -f $(SRCDIR)/*.mod.f90
	-rm -f $(RUNDIR)/*.mod

distclean: clean
	-rm -f $(OBJDIR)/*~

mach_wrapper:
	@echo $(MACHINE)

$(OBJDIR)/%.o: checkpath
	$(MAKE) -f $(MKFILE) $@

forcheck:
	$(MAKE) -f $(MKFILE) forcheck

show_srcnames:
	$(MAKE) -f $(MKFILE) show_srcnames
