#pragma once
#include "F:\source\repos\raymarch\raymarch\vecs.h"
using namespace vecs;
namespace intersect {
	double triangle(vec3 rpos, vec3 rdir, vec3 t1, vec3 t2, vec3 t3) {
		t1 -= rpos;
		t2 -= rpos;
		t3 -= rpos;
		double t1dist = dot(t1, rdir);
		double t2dist = dot(t2, rdir);
		double t3dist = dot(t3, rdir);
		t1 -= rdir * t1dist;
		t2 -= rdir * t2dist;
		t3 -= rdir * t3dist;
		double w1 = ((t2.y - t3.y)*(-t3.x) + (t3.x - t2.x) * (-t3.y)
	}
}