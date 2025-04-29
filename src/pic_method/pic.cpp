#include "pic.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;

std::vector<Particle_PIC> pic_particles;
Grid_PIC pic_grid(0, 0, 0.0f);  // 初始化为默认值
float dt_pic = 0.1f; 


// Initialize particles in a grid
// particles: The particle container
// num_w: Number of grid cells in x direction
// num_h: Number of grid cells in y direction
// dx: grid cell size
void init_particle(std::vector<Particle_PIC>& particles,int num_w,int num_h,float dx){
    //Clear existing particles
    particles.clear();

    for (int i=1;i<num_w-1;i++){
        for (int j=1;j<num_h-1;j++){
            float x=(i+0.5f)*dx;
            float y=(j+0.5f)*dx;
            Particle_PIC particle(x,y,0.0f,0.0f);
            particles.push_back(particle);            
        }
    }
}


// Reset all grid velocities and masses to zero
// grid: The MAC grid
void reset_grid(Grid_PIC& grid){
    for (int i=0;i<=grid.grid_width;i++){
        for (int j=0;j<grid.grid_height;j++){
            grid.g_velocity_x[i][j]=0.0f;
            grid.g_mass_x[i][j]=0.0f;
        }
    }
    for (int i=0;i<grid.grid_width;i++){
        for (int j=0;j<=grid.grid_height;j++){
            grid.g_velocity_y[i][j]=0.0f;
            grid.g_mass_y[i][j]=0.0f;
        }
    }
}

// helper function
// Transfer a single particle's velocity component to nearby grid cells using B-spline
// g_velocity: 2D grid storing velocity
// g_mass: 2D grid storing mass weights
// p_x: x-coordinate of particle
// p_y: y-coordinate of particle
// velocity_component: particle velocity in either x or y
// offset_x: grid offset in x
// offset_y: grid offset in y
// grid_width: width of the g_velocity grid
// grid_height: height of the g_velocity grid
// dx: grid cell size
void transfer_velocity(std::vector<std::vector<float>>& g_velocity,std::vector<std::vector<float>>& g_mass,
    float p_x, float p_y,float velocity_component,float offset_x, float offset_y,
    int grid_width, int grid_height,float dx){
        float g_x=p_x/dx-offset_x;
        float g_y=p_y/dx-offset_y;

        int int_x=static_cast<int>(g_x);
        int int_y=static_cast<int>(g_y);

        float x_f=g_x-int_x;
        float y_f=g_y-int_y;

    
        //B-spline kernels
        float w_x[3]={0.5f*(1.5f-x_f)*(1.5f-x_f),0.75f-(x_f-1.0f)*(x_f-1.0f),0.5f*(x_f-0.5f)*(x_f-0.5f)};
        float w_y[3]={0.5f*(1.5f-y_f)*(1.5f-y_f),0.75f-(y_f-1.0f)*(y_f-1.0f),0.5f*(y_f-0.5f)*(y_f-0.5f)};

        //3*3 grid
        for (int i=0;i<3;i++){
            int x_i=int_x+i;
            if (x_i<0||x_i>=(int)g_velocity.size()){
                continue;
            }

            for (int j=0;j<3;j++){
                int y_i=int_y+j;
                if (y_i<0||y_i>=(int)g_velocity[0].size()){
                    continue;
                }

                float weight=w_x[i]*w_y[j];

                g_velocity[x_i][y_i]+=velocity_component*weight;
                g_mass[x_i][y_i]+=weight;

            }
        }
    }

// Transfer all particles' velocities to the grid using B-spline interpolation and normalize
// particles: array containing PIC particles
// grid: a MAC grid
void p2g(const std::vector<Particle_PIC>& particles,Grid_PIC& grid){
    for (auto& p:particles){
        float p_x=p.x;
        float p_y=p.y;

        // velocity of x direction
        transfer_velocity(grid.g_velocity_x,grid.g_mass_x,p_x,p_y,p.velocity_x,0.5f,0.0f,grid.grid_width+1,grid.grid_height,grid.grid_dx);
        // velocity of y direction
        transfer_velocity(grid.g_velocity_y,grid.g_mass_y,p_x,p_y,p.velocity_y,0.0f,0.5f,grid.grid_width,grid.grid_height+1,grid.grid_dx);
    }
    //normalization
    for (int i=0;i<=grid.grid_width;i++){
        for (int j=0;j<grid.grid_height;j++){
            if (grid.g_mass_x[i][j]>0.0f){
                grid.g_velocity_x[i][j]/=grid.g_mass_x[i][j];
            }
        }
    }

    for (int i=0;i<grid.grid_width;i++){
        for (int j=0;j<=grid.grid_height;j++){
            if (grid.g_mass_y[i][j]>0.0f){
                grid.g_velocity_y[i][j]/=grid.g_mass_y[i][j];
            }
        }
    }
}

// Apply gravity
// grid: The grid storing velocities
// dt: Timestep
void apply_gravity(Grid_PIC& grid,float dt){
    for (int i=0;i<grid.grid_width;i++){
        for (int j=0;j<=grid.grid_height;j++){
            grid.g_velocity_y[i][j]+=-9.8f*dt;
        }
    }
}



// helper function 
// Convert 2D grid indices to 1D index for sparse matrix representation
// i: Row index
// j: Column index
// height: Height of the 2D grid
int transfer_index(int i,int j, int height){
    return i*height +j;
}



// Solve the Poisson equation ∇²p = ∇·u using a CG solver
// grid: MAC grid storing velocities and pressure
// dt: Time step
void solve_pressure(Grid_PIC& grid,float dt){
    using SPM=Eigen::SparseMatrix<float>;
    using Triplet =Eigen::Triplet<float>; //A element(row,column,value)
    using VEC= Eigen::VectorXf;

    int w =grid.grid_width;
    int h =grid.grid_height;
    int N=w*h;
    float dx=grid.grid_dx;

    // Construct sparse matrix A
    SPM A(N,N);
    std::vector<Triplet> coef;
    for (int i=1;i<w-1;i++){
        for (int j=1;j<h-1;j++){
            int index=transfer_index(i,j,h);
            coef.emplace_back(index,index,4.0f);
            coef.emplace_back(index,transfer_index(i-1,j,h),-1.0f);//left
            coef.emplace_back(index,transfer_index(i+1,j,h),-1.0f);//right
            coef.emplace_back(index,transfer_index(i,j-1,h),-1.0f);//bottom
            coef.emplace_back(index,transfer_index(i,j+1,h),-1.0f);//top
        }
    }
    A.setFromTriplets(coef.begin(),coef.end());// Load the Triplet list into the sparse matrix

    // Construct the vector on the right,∇·u
    VEC B(N);
    for (int i=1;i<w-1;i++){
        for (int j=1;j<h-1;j++){
            //right-left/dx+top-bottom/dy
            float div=(grid.g_velocity_x[i+1][j]-grid.g_velocity_x[i][j])/dx+(grid.g_velocity_y[i][j+1]-grid.g_velocity_y[i][j])/dx;
            B[transfer_index(i,j,h)]=div;
        }
    }

    // A * x = B
    Eigen::ConjugateGradient<SPM,Eigen::Lower|Eigen::Upper> solver;
    solver.compute(A);
    VEC x=solver.solve(B);
    if (solver.info() != Eigen::Success) {
        // std::cerr << "⚠️ Pressure solve failed! Setting pressure to zero.\n";
        x.setZero();  // 万一solver失败，就不要让x是垃圾数，而是直接清零
    }

    // Fill the solution into the perssure grid
    for (int i=1;i<w-1;i++){
        for (int j=1;j<h-1;j++){
            grid.g_pressure[i][j]=x[transfer_index(i,j,h)];
        }
    }

    // Correct the velocity field
    for (int i=1;i<w;i++){
        for (int j=0;j<h;j++){
            grid.g_velocity_x[i][j]-=dt*(grid.g_pressure[i][j]-grid.g_pressure[i-1][j])/dx;
        }
    }
    for (int i=0;i<w;i++){
        for (int j=1;j<h;j++){
            grid.g_velocity_y[i][j]-=dt*(grid.g_pressure[i][j]-grid.g_pressure[i][j-1])/dx;
        }
    }
}



// helper function 
// Interpolate grid velocity at the given particle position using B-spline
// g_velocity: 2D grid storing velocity
// p_x: x-coordinate of particle
// p_y: y-coordinate of particle
// offset_x: Grid offset(x-direction)
// offset_y: Grid offset(y-direction)
// dx: grid cell size
// grid_width: width of g_velocity grid
// grid_height: height of g_velocity grid
float interp_velocity(const std::vector<std::vector<float>>& g_velocity,float p_x, float p_y,
                     float offset_x, float offset_y,float dx,int grid_width, int grid_height) {
    float g_x= p_x/dx - offset_x;
    float g_y= p_y/dx - offset_y;

    int int_x= static_cast<int>(g_x);
    int int_y= static_cast<int>(g_y);

    float x_f= g_x - int_x;
    float y_f= g_y - int_y;

    //B-spline kernels
    float w_x[3]={0.5f*(1.5f-x_f)*(1.5f-x_f),0.75f-(x_f-1.0f)*(x_f-1.0f),0.5f*(x_f-0.5f)*(x_f-0.5f)};
    float w_y[3]={0.5f*(1.5f-y_f)*(1.5f-y_f),0.75f-(y_f-1.0f)*(y_f-1.0f),0.5f*(y_f-0.5f)*(y_f-0.5f)};

    float result = 0.0f;
    for (int i = 0; i < 3; ++i) {
        int x_i = int_x + i;
        if (x_i < 0|| x_i >= (int)g_velocity.size()){
            continue;
        }

        for (int j = 0; j < 3; ++j) {
            int y_i = int_y + j;
            if (y_i < 0|| y_i >= (int)g_velocity[0].size()){
                continue;
            }

            float weight = w_x[i]*w_y[j];
            result+= g_velocity[x_i][y_i]*weight;
        }
    }
    return result;
}

// Update particles'velocities by interpolating from the grid
// particles: Array of particles to update
// grid: The MAC grid 
void g2p(std::vector<Particle_PIC>& particles, const Grid_PIC& grid){
    for (auto&p:particles){
        float p_x=p.x;
        float p_y=p.y;

        p.velocity_x=interp_velocity(grid.g_velocity_x,p_x,p_y,0.5f,0.0f,grid.grid_dx,grid.grid_width+1,grid.grid_height);
        p.velocity_y=interp_velocity(grid.g_velocity_y,p_x,p_y,0.0f,0.5f,grid.grid_dx,grid.grid_width,grid.grid_height+1);
    }
}

// Move particles according to their velocities and handle boundary
// particles: Array of particles to update
// dt: Time step
// grid: Grid for clamping positions
void advect(std::vector<Particle_PIC>& particles, float dt,const Grid_PIC& grid){
    float dx=grid.grid_dx;

    //boundary restrictions
    float min_x=dx;
    float max_x=dx*grid.grid_width-dx;
    float min_y=dx;
    float max_y=dx*grid.grid_height-dx;
    
    for (auto& p:particles){
        p.x+=p.velocity_x*dt;
        p.y+=p.velocity_y*dt;

        //After hitting the boundary, speed directly becomes 0
        if(p.x<min_x){
            p.x=min_x;
            p.velocity_x=0.0f;
        }
        if(p.y<min_y){
            p.y=min_y;
            p.velocity_y=0.0f;
        }
        if(p.x>max_x){
            p.x=max_x;
            p.velocity_x=0.0f;
        }
        if(p.y>max_y){
            p.y=max_y;
            p.velocity_y=0.0f;
        }

    }
}

// Perform a full simulation step:reset grid, P2G, gravity, pressure solve, G2P, and advection.
// particles: Particle array in simulation
// grid: MAC grid storing velocity and pressure
// dt_pic: Time step size
void simul_step(std::vector<Particle_PIC>& particles,Grid_PIC& grid, float dt_pic){
    reset_grid(grid);
    p2g(particles,grid);
    apply_gravity(grid,dt_pic);
    solve_pressure(grid,dt_pic);
    g2p(particles,grid);
    advect(particles,dt_pic,grid);
}