# Fluid-Simulation
We compare grid based method (Stable Fluids) with particle based method (Smoothed Particle Hydrodynamics)

## Sample Usage

### Step 1: Install Dependencies(macOS)
Make sure Xcode Command Line Tools are installed (includes OpenGL and GLUT frameworks):

```bash
xcode-select --install
```
### Step 2: Compile:
```bash
clang++ src/main.cpp -framework OpenGL -framework GLUT -o app
```
### Step 3: Run:
```bash
./app 
```
## Directory Structure
```bash
├── src
│   ├── main.cpp
│   ├── input.h ── input.cpp
│   ├── grid.h ── grid.cpp
│   └── particle.h ── particle.cpp
├── include
│   └── Eigen
├── shaders
│   └── 
├── LICENSE
├── README.md
└── .gitignore
```