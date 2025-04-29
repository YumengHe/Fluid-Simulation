#include <iostream>
#include <fstream>
#include <vector>
#include "../src/constants.h"


/**
 * @brief 在网格模式下生成APIC粒子
 * @param dt 时间步长
 * @param grid_width 网格宽度
 * @param grid_height 网格高度
 * @param grid_dx 网格大小
 * @param rows 粒子网格的行数
 * @param cols 粒子网格的列数
 * @param spacing 相邻粒子之间的距离
 * @param outputFile 粒子数据的输出文件路径
 */
void generateAPICParticles(
    float dt,
    int grid_width,
    int grid_height,
    int rows,
    int cols,
    float spacing,
    const std::string& outputFile = "particles.apic"
) {
    std::ofstream file(outputFile);
    
    // 写入时间步长
    file << dt << "\n";
    
    // 写入网格信息
    file << grid_width << " " << grid_height << "\n";
    
    // 写入粒子部分标记
    
    // 计算起始位置（在视图中居中）
    float grid_dx = VIEW_WIDTH / grid_width;
    float startX = (grid_width * grid_dx - (cols - 1) * spacing) / 2.0f;
    float startY = (grid_height * grid_dx - (rows - 1) * spacing) / 2.0f;

    // 计算起始位置（从左下角开始）
    // float startX = grid_dx;  // 从第一个网格开始
    // float startY = grid_dx;  // 从第一个网格开始
    
    // 在网格模式下生成粒子
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // 计算粒子位置
            float x = startX + j * spacing;
            float y = startY + (rows - 1 - i) * spacing;
            
            // 初始速度设为0
            float vx = 0.0f;
            float vy = 0.0f;
            
            // 初始化仿射矩阵为零矩阵（2x2矩阵，按行存储）
            float B11 = 0.0f, B12 = 0.0f;
            float B21 = 0.0f, B22 = 0.0f;
            
            // 写入粒子数据：位置、速度、仿射矩阵
            file << x << " " << y << " "          // 位置
                 << vx << " " << vy << " "        // 速度
                 << B11 << " " << B12 << " "      // 仿射矩阵第一行
                 << B21 << " " << B22 << "\n";    // 仿射矩阵第二行
        }
    }
    
    file.close();
    
    // 打印生成摘要
    std::cout << "✅ 已生成APIC配置文件到 " << outputFile << std::endl;
    std::cout << "   网格大小: " << grid_width << "x" << grid_height << std::endl;
    std::cout << "   粒子数量: " << rows * cols << std::endl;
}

int main() {
    // 示例用法
    generateAPICParticles(
        0.01f,      // dt: 时间步长
        50,         // number of grids in x direction
        50,         // number of grids in y direction
        20,          // rows: 粒子行数
        20,          // cols: 粒子列数
        10.0f,       // spacing: 粒子间距（设为2倍grid_dx以便观察）
        "100center.apic" // 输出文件名
    );
    return 0;
}