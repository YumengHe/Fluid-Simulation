# === Directory ===
SRC_DIR := src
INCLUDE_DIR := include

# === Tool ===
CXX := g++
CXXFLAGS := -std=c++17 -I$(INCLUDE_DIR)

# === Source, object, and target names ===
SRCS := $(SRC_DIR)/particle.cpp
OBJS := $(SRCS:%.cpp=%.o)
TARGET := particle

# === Default ===
.PHONY: all
all: $(TARGET)

# === Link ===
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@

# === Compile ===
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# === Clean ===
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)