# Fluid-Simulation
We compare grid based method (Stable Fluids), particle based method (Smoothed Particle Hydrodynamics), particle-in-cell method and APIC.

![fluid simualtion methods](/images/fluid_simulation_methods.png)

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
Clean Build Files
```bash
make clean
```
### Step 3: Run:
Run Grid/Particle Demo
```bash
./fluid ./fluids/grid_or_particle_file_name
```
Run PIC/FLIP Demo
```bash
make picflip_demo
./picflip_demo
```
Run APIC Demo
```bash
make apic_demo
./apic_demo
```
### Grid File Generation
```bash
cd fluids
g++ -o grid_generator gridGenerator.cpp
./grid_generator template.grid
```
### Particle File Generation
```bash
cd fluids
g++ -o particle_gen particleGenerator.cpp
./particle_gen
```
### Interaction Features
	•	Mouse Click: Position Output
	•	Esc Key: Exit Program
## Directory Structure
```bash
├── src
│   ├── main.cpp
│   ├── input.h ── input.cpp
│   ├── grid.h ── grid.cpp
│   ├── particle.h ── particle.cpp
│   ├── pic_method
│   │   ├──grid_pic.h
│   │   ├──particle_pic.h
│   │   ├──pic.h ── pic.cpp
│   │   ├──apic.h ── apic.cpp
│   │   ├──picflip.h ── picflip.cpp
│   │   └──picflip_main.cpp
│   └── apic_method
│       ├──fluidsim.h
│       ├──fluidsim.cpp
│       └──main_apic.cpp
├── include
│   └── Eigen
├── fluids
│   ├── gridGenerator.cpp
│   └── template.grid
│   └── particleGenerator.cpp
│   └── test.par
├── images
├── LICENSE
├── README.md
└── .gitignore
```

## Grid
Based on stable fluid method [[1]](#1)[[2]](#2).
```bash
./fluid ./fluids/grid/10.grid
./fluid ./fluids/grid/15.grid
./fluid ./fluids/grid/20.grid
./fluid ./fluids/grid/50.grid
```

## Particle
Based on SPH.
```bash
./fluid
./fluid ./fluids/grid/template.grid
./fluid ./fluids/particle/400particles.par
./fluid ./fluids/400particles.par
```
## PIC/FLIP
Based on PIC and FLIP method [[3]](#3)
```bash
make picflip_demo
./picflip_demo
```


## APIC
Based on PIC and APIC method [[4]](#4).

## References
<a id="1">[1]</a> 
Stam, J. & Alias wavefront. (1999b). Stable fluids. In Alias Wavefront (p. 121) [Journal-article]. Alias wavefront. https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf

<a id="1">[2]</a> 
Stam, J. (2003, March). Real-time fluid dynamics for games. In Proceedings of the game developer conference (Vol. 18, No. 11). https://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf 

<a id="3">[3]</a> 
Robert Bridson. (2015). Fluid Simulation for Computer Graphics (2nd ed., Chapter 7.6). CRC Press. [Book]. 

<a id="4">[4]</a>
Jiang, C., Schroeder, C., Selle, A., Teran, J., & Stomakhin, A. (2015). The affine particle-in-cell method. ACM Transactions on Graphics, 34(4), Article 51. https://doi.org/10.1145/2766996


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