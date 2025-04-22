#include "pic.h"

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

//helper function
//p_x:
//p_y:
//velocity_component:
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


void p2g(const std::vector<Particle_PIC>& particles,Grid_PIC& grid){
    for (auto& p:particles){
        float p_x=p.x;
        float p_y=p.y;

        // velocity of x direction
        transfer_velocity(grid.g_velocity_x,grid.g_mass_x,p_x,p_y,p.velocity_x,0.5f,0.0f,grid.grid_width+1,grid.grid_height,grid.grid_dx);
        //velocity of y direction
        transfer_velocity(grid.g_velocity_y,grid.g_mass_y,p_x,p_y,p.velocity_y,0.0f,0.5f,grid.grid_width,grid.grid_height+1,grid.grid_dx);

        //归一化
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
}



void apply_gravity(Grid_PIC& grid,float dt){
    for (int i=0;i<grid.grid_width;i++){
        for (int j=0;j<=grid.grid_height;j++){
            grid.g_velocity_y[i][j]+=-9.8f*dt;
        }
    }
}

//未完成
void solve_pressure(Grid_PIC& grid,float dt);




void g2p(std::vector<Particle_PIC>& particles, const Grid_PIC& grid);






void advect(std::vector<Particle_PIC>& particles, float dt);

