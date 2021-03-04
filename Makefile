###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
ALL_MODULES := no
HEAT_CONDUCTION := yes
NAVIER_STOKES := yes
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

ifeq ($(NAVIER_STOKES),yes)
        FLUID_PROPERTIES            := yes
				RAY_TRACING                 := yes
        RDG                         := yes
        HEAT_CONDUCTION             := yes
endif

# Use the SQUIRREL submodule if it exists and SQUIRREL_DIR is not set
SQUIRREL_SUBMODULE    := $(CURDIR)/squirrel
ifneq ($(wildcard $(SQUIRREL_SUBMODULE)/Makefile),)
  SQUIRREL_DIR        ?= $(SQUIRREL_SUBMODULE)
else
  SQUIRREL_DIR        ?= $(shell dirname `pwd`)/squirrel
endif

# squirrel
APPLICATION_DIR    := $(SQUIRREL_DIR)
APPLICATION_NAME   := squirrel
include            $(FRAMEWORK_DIR)/app.mk

# Use the CRANE submodule if it exists and CRANE_DIR is not set
CRANE_SUBMODULE := $(CURDIR)/crane
ifneq ($(wildcard $(CRANE_SUBMODULE)/Makefile),)
  CRANE_DIR        ?= $(CRANE_SUBMODULE)
else
  CRANE_DIR        ?= $(shell dirname `pwd`)/crane
endif

# crane
APPLICATION_DIR    := $(CRANE_DIR)
APPLICATION_NAME   := crane
include            $(FRAMEWORK_DIR)/app.mk

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := zapdos
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
