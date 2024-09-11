SHELL:=/bin/bash

GCOM_DIR := $(shell git rev-parse --show-toplevel)
CXX=g++
# stdc++fs is for <filesystem>
# -lz -lgzstream is for gzip compression 
# -lboost_serialization for reading cache stat binary file
LIBS=-lstdc++fs -lgzstream -lboost_serialization -lz -lsassparser 
LIB_DIRS=-L$(GCOM_DIR)/third-party/gzstream -L$(GCOM_DIR)/third-party/zlib-1.2.12/lib -L$(GCOM_DIR)/third-party/boost_1_86_0/stage/lib -Wl,-rpath=$(GCOM_DIR)/third-party/boost_1_86_0/stage/lib
TARGET=GCoM
# TODO : differentiate build and bin directory per make target? OBJS should be redefined
BIN_DIR=$(GCOM_DIR)/bin
BIULD_DIR=$(GCOM_DIR)/build

SRCS=$(shell find src third-party/accel-sim/addrdec -type f -name *.cc)
INCLUDES=-I$(GCOM_DIR)/third-party/gzstream -I$(GCOM_DIR)/third-party/zlib-1.2.12/include -I$(GCOM_DIR)/third-party/boost_1_86_0/ -I$(GCOM_DIR)/src/include -I$(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/trace-parser -I$(GCOM_DIR)/third-party/accel-sim/addrdec
OBJS=$(SRCS:%.cc=$(BIULD_DIR)/%.o)
HEADS=$(shell find src -type f -name *.h)

.PHONY=all default makedirs sassparser clean

all: default

default: CXXFLAG=-std=c++17 -O3 -g -Werror
default: DEBUG=release
default: makedirs $(TARGET) sassparser

debug: CXXFLAG=-std=c++17 -O0 -g -Wall
debug: DEBUG=debug
debug: makedirs $(TARGET)

makedirs:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi;
	
	if [ ! -d $(BIULD_DIR) ]; then mkdir -p $(shell dirname $(OBJS)); fi;

sassparser:
	source $(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/setup_environment.sh $(DEBUG) && make -C $(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/ -j

$(TARGET): LIB_DIRS+=-L$(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/lib/$(DEBUG) -Wl,-rpath=$(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/lib/$(DEBUG)
$(TARGET): $(OBJS) $(HEADS) sassparser
	$(CXX) $(LIB_DIRS) $(CXXFLAG) -o $(BIN_DIR)/$(TARGET) $(OBJS) $(INCLUDES) $(LIBS)

$(BIULD_DIR)/%.o: %.cc
	$(CXX) $(CXXFLAG) -c $< -o $@ $(INCLUDES)

clean:
	rm -f $(BIN_DIR)/$(TARGET)
	rm -f $(OBJS)
	rm -r $(BIULD_DIR)
	make -C $(GCOM_DIR)/third-party/accel-sim/SASS_parser/gpu-simulator/ clean
# clean:
# 	rm -rf $(BIN_DIR) $(BIULD_DIR)
# 	make -C third-party/accel-sim/SASS_parser/gpu-simulator/ clean