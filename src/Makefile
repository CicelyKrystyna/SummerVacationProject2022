# define compiler options
CPP = g++
CPPFLAGS = -Wall -g -O3 -DEBUG  -fpermissive 

# macros
RM = rm -f

# libraries
LDFLAGS = 

# paths to include files (out of current directory)
INCLUDES =

# output
EXEC=cell++

# input
INPUT=main.cpp \
	CoupledModel.cpp Fibre.cpp Cell.cpp Vessel.cpp \
	PDE.cpp  Param.cpp Mesh.cpp


DIRS := ./
FILES:= $(INPUT) $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
SRCS := $(notdir $(INPUT))
# rule to create object files
OBJS := $(SRCS:.cpp=.o)


.PHONY: all clean

all: $(EXEC) 

$(EXEC): $(OBJS) 
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(OBJS) $(LDFLAGS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@


clean: 	
	rm -f *.o cell++ 
