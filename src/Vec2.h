#ifndef VEC2_H
#define VEC2_H

#include <cmath>

class Vec2 {
public:
    float x, y;
    Vec2() : x(0.0f), y(0.0f) {}
    Vec2(float x_, float y_) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& other) const { return Vec2(x + other.x, y + other.y); }
    Vec2& operator+=(const Vec2& other) { x += other.x; y += other.y; return *this; }

    Vec2 operator-(const Vec2& other) const { return Vec2(x - other.x, y - other.y); }
    Vec2& operator-=(const Vec2& other) { x -= other.x; y -= other.y; return *this; }

    Vec2 operator*(float s) const { return Vec2(x * s, y * s); }
    Vec2& operator*=(float s) { x *= s; y *= s; return *this; }

    float dot(const Vec2& other) const { return x * other.x + y * other.y; }
    float length() const { return std::sqrt(x * x + y * y); }
    float lengthSquared() const { return x * x + y * y; }

    Vec2 normalized() const {
        float len = length();
        return len > 1e-6f ? Vec2(x / len, y / len) : Vec2(0.0f, 0.0f);
    }
};

#endif // VEC2_H
