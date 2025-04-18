# === Directory ===
SRC_DIR := src
INCLUDE_DIR := include

# === Tool ===
CXX := clang++
CXXFLAGS := -std=c++17 -I$(INCLUDE_DIR)
LDFLAGS := -framework OpenGL -framework GLUT

# === Source, object, and target names ===
# SRCS := $(SRC_DIR)/particle.cpp
SRCS := $(SRC_DIR)/main.cpp
# SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:%.cpp=%.o)
TARGET := fluid


# === Default ===
.PHONY: all
all: $(TARGET)

# === Link ===
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# === Compile ===
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# === Clean ===
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)