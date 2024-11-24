#pragma once
#include<math.h>
#include<cmath>
#include <algorithm>
namespace vecs {
	//Integer 2D Vector

	struct ivec2 {
		int x, y;
		double len() {
			return pow(pow(x, 2) + pow(y, 2), 0.5);
		}
		void norm() {
			y = y / this->len();
			x = x / this->len();
		}
		struct ivec2& operator+=(const ivec2& rhs) { x += rhs.x; y += rhs.y; return *this; }
		struct ivec2& operator-=(const ivec2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
		struct ivec2& operator*=(const ivec2& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
		struct ivec2& operator*=(const int& rhs) { x *= rhs; y *= rhs; return *this; }
	};
	ivec2 operator+(ivec2 lhs, const ivec2& rhs) { return lhs += rhs; }
	ivec2 operator-(ivec2 lhs, const ivec2& rhs) { return lhs -= rhs; }
	ivec2 operator*(ivec2 lhs, const ivec2& rhs) { return lhs *= rhs; }
	ivec2 operator*(ivec2 lhs, const int k) { return lhs *= k; }
	int dot(ivec2 p1, ivec2 p2) {
		return p1.x * p2.x + p1.y * p2.y;
	}
	double dist(ivec2 o1, ivec2 o2) {
		return pow(pow(o1.x - o2.x, 2) + pow(o1.y - o2.y, 2), 0.5);
	}
	//2D Vector
	struct vec2 {
		double x, y;
		double len() {
			return pow(pow(x, 2) + pow(y, 2), 0.5);
		}
		void norm() {
			y = y / this->len();
			x = x / this->len();
		}
		struct vec2& operator+=(const vec2& rhs) { x += rhs.x; y += rhs.y; return *this; }
		struct vec2& operator-=(const vec2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
		struct vec2& operator*=(const vec2& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
		struct vec2& operator*=(const double& rhs) { x *= rhs; y *= rhs; return *this; }
	};
	vec2 operator+(vec2 lhs, const vec2& rhs) { return lhs += rhs; }
	vec2 operator-(vec2 lhs, const vec2& rhs) { return lhs -= rhs; }
	vec2 operator*(vec2 lhs, const vec2& rhs) { return lhs *= rhs; }
	vec2 operator*(vec2 lhs, const double k) { return lhs *= k; }
	vec2 normalize(vec2 p1, vec2 p2) {
		vec2 sub = p2 - p1;
		sub *= (1.0 / sub.len());
		return sub;
	}
	double dot(vec2 p1, vec2 p2) {
		return p1.x * p2.x + p1.y * p2.y;
	}
	double dist(vec2 o1, vec2 o2) {
		return pow(pow(o1.x - o2.x, 2) + pow(o1.y - o2.y, 2), 0.5);
	}

	// 3D Vector
	struct vec3 {
		double x, y, z;
		vec2 xy = { this->x, this->y };
		vec2 xz = { this->x, this->z };
		vec2 yz = { this->y, this->z };
		vec2 yx = { this->y, this->x };
		vec2 zx = { this->z, this->x };
		vec2 zy = { this->z, this->y };
		double len() {
			return sqrt(x*x + y*y + z*z);
		}
		void norm() {
			double length = this->len();
			x /= length;
			y /= length;
			z /= length;

		}
		struct vec3& operator+=(const vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
		struct vec3& operator-=(const vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
		struct vec3& operator*=(const vec3& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
		struct vec3& operator*=(const double& rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
	};
	vec3 operator+(vec3 lhs, const vec3& rhs) { return lhs += rhs; }
	vec3 operator-(vec3 lhs, const vec3& rhs) { return lhs -= rhs; }
	vec3 operator*(vec3 lhs, const vec3& rhs) { return lhs *= rhs; }
	vec3 operator*(vec3 lhs, const double k) { return lhs *= k; }
	vec3 normalize(vec3 p1, vec3 p2) {
		vec3 sub = p2 - p1;
		sub *= (1.0 / sub.len());
		return sub;
	}
	double dot(vec3 p1, vec3 p2) {
		return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
	}
	double dist(vec3 o1, vec3 o2) {
		return sqrt(pow(o1.x - o2.x, 2) + pow(o1.y - o2.y, 2) + pow(o1.z - o2.z, 2));
	}
	vec3 cross(vec3 p1, vec3 p2) {
		return { p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z,p1.x * p2.y - p1.y * p2.x };
	}
	vec3 v_abs(vec3 p1) {
		p1.x = abs(p1.x);
		p1.y = abs(p1.y);
		p1.z = abs(p1.z);
		return p1;
	}
	vec3 v_min(vec3 p1, double k) {
		p1.x = std::min(p1.x, k);
		p1.y = std::min(p1.y, k);
		p1.z = std::min(p1.z, k);
		return p1;
	}
	vec3 v_max(vec3 p1, double k) {
		p1.x = std::max(p1.x, k);
		p1.y = std::max(p1.y, k);
		p1.z = std::max(p1.z, k);
		return p1;
	}
	vec2 toVec2(ivec2 iv) {
		return { (double)iv.x, (double)iv.y };
	}
	ivec2 toiVec2(vec2 iv) {
		return { (int)iv.x, (int)iv.y };
	}
}