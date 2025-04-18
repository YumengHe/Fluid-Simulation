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
make
```
<!-- ```bash
g++ -std=c++14 src/*.cpp -Iinclude -o fluid -framework OpenGL -framework GLUT
``` -->
Clean Build Files
```bash
make clean
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
### Interaction Features
	•	Mouse Click: Position Output
	•	Esc Key: Exit Program
## Directory Structure
```bash
├── src
│   ├── main.cpp
│   ├── input.h ── input.cpp
│   ├── grid.h ── grid.cpp
│   └── particle.h ── particle.cpp
├── include
│   ├── Eigen
│   └── json
├── shaders
│   └── 
├── LICENSE
├── README.md
├── particle_init.py
└── .gitignore
```

## Grid
Based on stable fluid method [[1]](#1)[[2]](#2).

## References
<a id="1">[1]</a> 
Stam, J. & Alias wavefront. (1999b). Stable fluids. In Alias Wavefront (p. 121) [Journal-article]. Alias wavefront. https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf

<a id="1">[2]</a> 
Stam, J. (2003, March). Real-time fluid dynamics for games. In Proceedings of the game developer conference (Vol. 18, No. 11). https://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf 


### Particles
```bash
@inproceedings{clavet2005viscoelastic,
  title={Particle-based viscoelastic fluid simulation},
  author={Clavet, Simon and Beaudoin, Philippe and Poulin, Pierre},
  booktitle={Proceedings of the 2005 ACM SIGGRAPH/Eurographics Symposium on Computer Animation},
  pages={219--228},
  year={2005}
}
```