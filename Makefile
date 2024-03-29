# * URL: http://www.partow.net/programming/makefile/index.html *
# make all
# make clean
# make program
# make build
# make release
# make debug

GIT_COMMIT  := $(shell git rev-parse --verify HEAD)
GIT_DATE    := $(firstword $(shell git --no-pager show --date=iso-strict --format="%ad" --name-only))
BUILD_DATE  := $(shell date --iso=seconds)

CXX      := -g++
CXXFLAGS := -std=c++11 -pedantic-errors -Wall -Wextra -Werror -pthread \
			-DGIT_COMMIT=\"'$(GIT_COMMIT)'\"\
 			-DGIT_DATE=\"'$(GIT_DATE)'\"\
 			-DBUILD_DATE=\"'$(BUILD_DATE)'\"
LDFLAGS  := -L/usr/lib
#LDLIBS   := -lm -lstdc++fs -lboost_system -lboost_filesystem -lboost_chrono
LDLIBS   := -lm -lstdc++fs
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/bin
TARGET   := malaria_model
INCLUDE  := -Iinclude/malaria_model/ -Iinclude/third_party/
SRC      :=                      \
	$(wildcard src/util/*.cpp) \
	$(wildcard src/common/*.cpp) \
	$(wildcard src/human/*.cpp) \
	$(wildcard src/*.cpp)         \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

TARGET_TEST  := test
INCLUDE_TEST := $(INCLUDE) -Itest/third_party/catch2
SRC_TEST 	 := $(wildcard test/*.cpp)
OBJECTS_TEST := $(filter-out $(OBJ_DIR)/src/$(TARGET).o, $(OBJECTS)) $(SRC_TEST:test/%.cpp=$(OBJ_DIR)/test/%.o)


.DELETE_ON_ERROR:

.PHONY: all
all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LDLIBS)


.PHONY: test

test: build $(APP_DIR)/$(TARGET_TEST)

$(OBJ_DIR)/test/%.o: test/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET_TEST): $(OBJECTS_TEST)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE_TEST) $(LDFLAGS) -o $(APP_DIR)/$(TARGET_TEST) $(OBJECTS_TEST) $(LDLIBS)


.PHONY: build clean
build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	@rm -r $(OBJ_DIR)/*
	@rm -r $(APP_DIR)/*	
