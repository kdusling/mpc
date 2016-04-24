# ------------------------------------------------
# Makefile
#
# Author: Kevin Dusling
# Date  : 2016-26-04
#
# ------------------------------------------------

# project name (generate executable with this name)
TARGET   = main.exe

# get gsl library flags
GSL_LIBS = $(shell gsl-config --libs)
GSL_CFLAGS = $(shell gsl-config --cflags)

CC       = gcc
# compiling flags here
CFLAGS   = -std=c99 -O3 -Wall -I. $(GSL_CFLAGS)

LINKER   = gcc -o
# linking flags here
LFLAGS   = -Wall -I. -lm $(GSL_LIBS)

# change these to set the proper directories where each files shoould be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm       = rm -f


$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONEY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

