# --- Paths ---
SRC_DIR     := src
PIC_DIR     := $(SRC_DIR)/pic_method
OBJ_DIR     := obj
INCLUDE_DIR := include

# --- Tools ---
CXX      ?= clang++
CXXFLAGS := -std=c++17 -O2 -Wall -I$(INCLUDE_DIR) -I$(SRC_DIR)
LDFLAGS  := -framework OpenGL -framework GLUT          # macOS
# LDFLAGS := -lglut -lGL -lX11                         # linux

# --- Sources & objects ---
SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(PIC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
TARGET := fluid

# --- Default rule ---
.PHONY: all
all: $(TARGET)

# --- Link step ---
$(TARGET): $(OBJS)
	$(CXX) $^ $(LDFLAGS) -o $@

# --- Compile step (duplicate folder structure inside OBJ_DIR) ---
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Clean ---
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(TARGET)