# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I./include -I./include/glo_stiff_system -I./include/grid_generation -I./include/gradient_recovery

# Debugging and Sanitizer flags
DEBUG_FLAGS := -g -O0
ASAN_FLAGS := -fsanitize=address -fno-omit-frame-pointer

# Directories
SRC_DIR := ./scr
OBJ_DIR := ./obj
BIN_DIR := ./bin

# Output binary
TARGET := $(BIN_DIR)/main

# Test binary
TEST_TARGET := $(BIN_DIR)/test

# Source and object files
SRCS := $(shell find $(SRC_DIR) -name '*.cpp') test.cpp
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# Enable debug or sanitizer flags if set
ifeq ($(DEBUG), 1)
  CXXFLAGS += $(DEBUG_FLAGS)
endif

ifeq ($(SANITIZE), 1)
  CXXFLAGS += $(ASAN_FLAGS)
endif

# Default rule
all: $(TARGET)

# Link the final binary
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile each source file into an object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Test target
test: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) test.cpp $^
	@echo "Running tests..."
	$(TEST_TARGET)

# Valgrind target
valgrind: $(TARGET)
	@echo "Running valgrind on $(TARGET)..."
	valgrind --leak-check=full --show-leak-kinds=all $(TARGET)

valgrind-test: $(TEST_TARGET)
	@echo "Running valgrind on $(TEST_TARGET)..."
	valgrind --leak-check=full --show-leak-kinds=all $(TEST_TARGET)

# Clean rule
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean test valgrind valgrind-test
