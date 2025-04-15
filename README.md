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
clang++ ./src/main.cpp -framework OpenGL -framework GLUT -o fluid
```
### Step 3: Run:
```bash
./fluid
```
### Particle Initialization
To generate initial particle configurations, use the included script:
```bash
python particle_init.py
```
This will generate a particles.json file with a grid of particles placed at a specified region within the simulation domain.

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
├── particle_init.py
└── .gitignore
```