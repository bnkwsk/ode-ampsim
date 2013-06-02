program_NAME := simulator
program_SRC_C_DIR := src/c
program_SRC_CXX_DIR := src/cpp
program_GEN_CXX_DIR := $(program_SRC_CXX_DIR)/generated
program_INT_DIR := intermediate
program_BIN_DIR := bin
program_C_SRCS := $(wildcard $(program_SRC_C_DIR)/*.c)
program_GEN_CXX_SRCS := $(program_GEN_CXX_DIR)/FirstBlock.cpp.part $(program_GEN_CXX_DIR)/SecondBlock.cpp.part $(program_GEN_CXX_DIR)/ThirdBlock.cpp.part
program_CXX_SRCS := $(wildcard $(program_SRC_CXX_DIR)/*.cpp)
program_HXX_SRCS := $(wildcard $(program_SRC_CXX_DIR)/*.hpp)
program_C_OBJS := $(program_C_SRCS:$(program_SRC_C_DIR)/%.c=$(program_INT_DIR)/%.o)
program_CXX_OBJS := $(program_CXX_SRCS:$(program_SRC_CXX_DIR)/%.cpp=$(program_INT_DIR)/%.o)
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := lib/flens/ lib/ .
program_LIBRARY_DIRS := 
program_LIBRARIES :=

CXX := g++-4.8

CXXFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) -pipe -std=c++11 -O3 -pthread -lpthread `pkg-config gtkmm-3.0 --cflags` -DNDEBUG -flto -march=native -ggdb3
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
LFLAGS += -lsndfile `pkg-config gtkmm-3.0 --libs` -lOpenCL -ltcmalloc -flto

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	mkdir -p $(program_BIN_DIR)
	$(LINK.cc) $(program_OBJS) -o $(program_BIN_DIR)/$(program_NAME) $(LFLAGS)

$(program_GEN_CXX_SRCS):
	mkdir -p $(program_GEN_CXX_DIR)
	maxima < src/maxima/dae.max > /dev/null

$(program_OBJS): $(program_INT_DIR)/%.o : $(program_SRC_CXX_DIR)/%.cpp $(program_GEN_CXX_SRCS) $(program_HXX_SRCS)
	mkdir -p $(program_INT_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	@- $(RM) -r $(program_GEN_CXX_DIR)

distclean: clean
