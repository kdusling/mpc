# ------------------------------------------------
# Makefile
#
# Author: Kevin Dusling
# Date  : 2016-26-04
#
# ------------------------------------------------

# get gsl library flags
GSL_LIBS = $(shell gsl-config --libs)
GSL_CFLAGS = $(shell gsl-config --cflags)

# change these to set the proper directories where each files should be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

EXCLUDES := $(SRCDIR)/tabulate.c  

SOURCES  := $(filter-out $(EXCLUDES), $(wildcard $(SRCDIR)/*.c))
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

TARGET   = main.exe
CC       = gcc 
CFLAGS   = -std=c99 -O3 -Wall -I. $(GSL_CFLAGS)
LINKER   = gcc 
LFLAGS   = -Wall -I. -lm -I. $(GSL_LIBS)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $(LFLAGS) $(OBJECTS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONEY: clean
clean:
	@rm -f $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	@rm -f $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

