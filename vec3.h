#pragma once
#include <iostream>
#include <cmath>

class vec3
{
private:
    double m_vec[3];

public:
    vec3(); // Create a zero vector
    vec3(double x, double y, double z);
    bool operator==(vec3 rhs);
    vec3 operator+(vec3 rhs);
    vec3 &operator+=(vec3 rhs);
    vec3 operator-(vec3 rhs);
    vec3 &operator-=(vec3 rhs);
    vec3 operator*(vec3 rhs);
    vec3 &operator*=(vec3 rhs);
    vec3 operator/(vec3 rhs);
    vec3 &operator/=(vec3 rhs);
    vec3 operator+(double scalar);
    vec3 &operator+=(double scalar);
    vec3 operator-(double scalar);
    vec3 &operator-=(double scalar);
    vec3 operator*(double scalar);
    vec3 &operator*=(double scalar);
    vec3 operator/(double scalar);
    vec3 &operator/=(double scalar);
    inline vec3 operator-() { return vec3(-m_vec[0], -m_vec[1], -m_vec[2]); }
    void add(vec3 &rhs) {
        m_vec[0] += rhs.x();
        m_vec[1] += rhs.y();
        m_vec[2] += rhs.z();
    }
    void addAndMultiply(vec3 &rhs, double scalar) {
        m_vec[0] += rhs.x()*scalar;
        m_vec[1] += rhs.y()*scalar;
        m_vec[2] += rhs.z()*scalar;
    }
    vec3 cross(vec3 &rhs);
    double dot(vec3 &rhs);
    void normalize();
    void setToZero();
    void randomGaussian(double mean, double standardDeviation);
    void randomUniform(double min, double max);
    void set(double x, double y, double z);
    inline double x() const { return m_vec[0]; }
    inline double y() const { return m_vec[1]; }
    inline double z() const { return m_vec[2]; }
    inline double &operator[](int index) { return m_vec[index]; }
    inline double operator[](int index) const { return m_vec[index]; }
    inline double lengthSquared() { return m_vec[0]*m_vec[0] + m_vec[1]*m_vec[1] + m_vec[2]*m_vec[2]; }
    inline double length() { return sqrt(lengthSquared()); }
    inline void subtract(const vec3 &v1, const vec3 &v2) { m_vec[0] = v1[0] - v2[0]; m_vec[1] = v1[1] - v2[1]; m_vec[2] = v1[2] - v2[2]; }

private:
    friend std::ostream& operator<<(std::ostream&stream, vec3 vec);

};
