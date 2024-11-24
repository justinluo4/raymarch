#define _CRT_SECURE_NO_DEPRECATE
#include "F:\source\repos\raymarch\raymarch\vecs.h"
#include "C:\Users\justi\source\repos\physsim_backend\cdraw.h"
#include <algorithm>
#include <cstdlib>
#include <math.h> 
#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <stdarg.h>
#include <windows.h>
#include <ctime>
#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
# define M_PI           3.14159265358979323846
using namespace std;
using namespace vecs;
using namespace draw;
struct ray {
	vec3 pos;
	vec3 dir;
	vec3 color;
	bool complete;
	int march_count;
	double travel_dist;
	double amb_dist;
	double closest_dist;
	double cur_IOR;
	void print_pos() {
		cout  << "position" << endl << pos.x << " " << pos.y << " " << pos.z << endl;
	}
};
struct light_point {
	vec3 pos;
	vec3 color;
	double brightness;
};
struct material {
	vec3 color = { 255, 255, 255 };
	bool reflective;
	bool transparent;

	double IOR;
};
struct sphere_surf : material {
	vec3 pos;
	double radius;
	double DE(vec3 p1) {
		//p1.x = fmod(p1.x, 1);
		//p1.y = fmod(p1.y, 1);
		//p1.z = fmod(p1.z, 1);
		double freq = 5;
		//double offset = sin(sin(sin(sin((p1.x - pos.x) * freq)*2) * 2) * 2) * sin(sin((p1.y - pos.y) * freq))* sin(sin((p1.z - pos.z) * freq));
		//offset *= 0.25;
		return dist(pos , p1) - radius;
	}
};
struct plane_surf : material {
	vec3 v1;
	vec3 v2;
	double y;
	double DE(vec3 p1) {
		vec3 offset_pos = p1;
		offset_pos.y -= y;
		vec3 norm_vec = cross(v1, v2);
		norm_vec.norm();
		return -dot(offset_pos, norm_vec);
	}

};
struct box_surf : material {
	vec3 dims;
	vec3 pos;
	double DE(vec3 p1) {
		
		vec3 corner_diff = dims - v_abs(p1);
		return v_min(corner_diff, 0).len() - min(0, max(max(corner_diff.x, corner_diff.y), corner_diff.z));
	}
};
struct mandelbulb : material {
	double dr = 1.0;
	double r = 0.0;
	int Iterations = 100;
	double Bailout = 5;
	int Power = 8;
	double mandel_DE(vec3 p1) {
		vec3 z = p1;

		for (int i = 0; i < Iterations; i++) {
			r = z.len();
			if (r > Bailout) break;

			// convert to polar coordinates
			double theta = acos(z.z / r);
			double phi = atan(z.y / z.x);
			if (z.x < 0) {
				phi += M_PI;
			}
			dr = pow(r, Power - 1.0) * Power * dr + 1.0;

			// scale and rotate the point
			double zr = pow(r, Power);
			theta = theta * Power;
			phi = phi * Power;

			// convert back to cartesian coordinates
			z = { sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta) };
			z *= zr;
			z += p1;
		}
		return 0.5 * log(r) * r / dr;
	}
};
double clamp(double a, double lb, double ub) {
	if (a < lb) {
		return lb;
	}
	if (a > ub) {
		return ub;
	}
	return a;
}
float smin(float a, float b, float k)
{
	float h = max(k - abs(a - b), 0.0) / k;
	return min(a, b) - h * h * k * (1.0 / 4.0);
}
float smax(float a, float b, float k)
{
	float h = max(k + abs(a + b), 0.0) / k;
	return min(a, -b) + h * h * k * (1.0 / 4.0);
}
double DE_main(vec3 p1, material &mat) {
	mat.reflective = false;
	mat.transparent = false;
	mat.IOR = 1;
	double center_radius = 1;
	vec3 dims = { center_radius, center_radius, center_radius };
	vec3 corner_diff = v_abs(p1) - dims ;
	double d = v_max(corner_diff, 0).len() - min(0, max(max(corner_diff.x, corner_diff.y), corner_diff.z)) - 0.3;
	double d2 = dist(p1, { 0, 1.5, 0 }) - 0.9;
	d = smin(d, d2, 0.8);
	/*
	vec2 pillar_center = { 1, 1 };
	vec2 pillar_bounds = { 0, 3 };
	double vert_pillar_dist = max(0, abs(p1.y - (pillar_bounds.x + pillar_bounds.y) / 2.) - (pillar_bounds.y - pillar_bounds.x) / 2);
	double horiz_pillar_dist = max(dist({ p1.x, p1.z }, pillar_center) - 0.7, 0);
	double d3 = sqrt(vert_pillar_dist * vert_pillar_dist + horiz_pillar_dist * horiz_pillar_dist);
	d = smin(d, d3, 0.3);

	if (abs(p1.x) < 2 && abs(p1.y) < 2 && abs(p1.z) < 2) {
		for (int octave = 1; octave < 10; octave *= 2) {
			double box_size = center_radius * 2 / octave;
			vec3 off = { center_radius, center_radius, center_radius };
			vec3 box_pos = p1 + off;
			box_pos.x = floor(box_pos.x*(octave/2.0))/ (octave / 2.0);
			box_pos.y = floor(box_pos.y* (octave / 2.0))/ (octave / 2.0);
			box_pos.z = floor(box_pos.z* (octave / 2.0))/ (octave / 2.0);
			
			for (int i = 0; i < 8; i++) {
				int j = i;
				vec3 box_corner = box_pos;
				if (j % 2) {
					box_corner.x += box_size;
				}
				j /= 2;
				if (j % 2) {
					box_corner.y += box_size;
				}
				j /= 2;
				if (j % 2) {
					box_corner.z += box_size;
				}
				j /= 2;
				d = max(d, -(dist(box_corner, p1 + off) - box_size *0.6*abs(sin(box_corner.x * 1000 + 2) * sin(box_corner.y * 100 + 3) * sin(box_corner.z * 2000 + 4))));
				
			}


		}
	}
	*/
	return d;
}

struct Marcher {
public:
	double contact_thresh = 2e-4;
	double render_dist = 1000;
	int max_steps = 70;
	int max_incs = 5;
	vec3 cam_pos = { 0, 1, 0 };
	vec2 cam_angle = { 0, 0 };
	vec3 cam_angle_vec = { 0, 0, 0 };
	vec2 cam_angle_rads;
	vec2 fov = { 90, 60 };
	vec2 fov_rads;
	vec2 window_res = { 50, 30 };
	vec3 bg_color = { 5, 5, 5 };
	vector<ray> rays;
	vector<sphere_surf> spheres;
	vector<plane_surf> planes;
	vector<light_point> light_points;
	vector<void (*)(vec3 p1)> dist_estims;
	void update_cam_angle_vec() {
		cam_angle_rads = cam_angle * (M_PI / 180);
		cam_angle_vec.x = cos(cam_angle_rads.x) * cos(cam_angle_rads.y);
		cam_angle_vec.z = sin(cam_angle_rads.x) * cos(cam_angle_rads.y);
		cam_angle_vec.y = sin(cam_angle_rads.y);
	}
	void initialize_rays() {
		update_cam_angle_vec();
		rays.clear();
		cam_angle_rads = cam_angle * (M_PI / 180);
		fov_rads = fov * (M_PI / 180);
		vec3 horiz_vec = { sin(cam_angle_rads.x) , 0, -cos(cam_angle_rads.x) };
		vec3 vert_vec = cross(cam_angle_vec, horiz_vec);
		vec3 horiz_max = horiz_vec * tan(fov_rads.x / 2);
		vec3 vert_max = vert_vec * tan(fov_rads.y / 2);
		for (int i = 0; i < window_res.y; i++) {
			for (int j = 0; j < window_res.x; j++) {
				ray r;
				r.pos = cam_pos;
				r.dir = cam_angle_vec;
				r.dir += horiz_vec * (1 - ((double)2 * j / window_res.x));
				r.dir += vert_vec * (1 - ((double)2 * i / window_res.y));
				r.dir.norm();
				r.complete = false;
				r.march_count = 0;
				r.color = { 0, 0, 0 };
				rays.push_back(r);
			}
		}
	}
	double DE(vec3 v1, material& closest_mat) {
		double min_dist = render_dist;
		double d;

		for (auto surf : spheres) {
			d = surf.DE(v1);
			if (d < min_dist) {
				closest_mat.reflective = surf.reflective;
				closest_mat.transparent = surf.transparent;
				closest_mat.IOR = surf.IOR;
				min_dist = d;
			}
		}
		for (auto surf : planes) {
			d = surf.DE(v1);
			if (d < min_dist) {
				closest_mat.reflective = surf.reflective;
				closest_mat.transparent = surf.transparent;
				closest_mat.IOR = surf.IOR;
				min_dist = d;
			}
		}
		material main_mat;
		d = DE_main(v1, main_mat);
		if (d < min_dist) {
			closest_mat = main_mat;
			min_dist = d;
		}
		return min_dist;
	}
	double DE(vec3 v1) {
		double min_dist = render_dist;
		double d;

		for (auto surf : spheres) {
			d = surf.DE(v1);
			if (d < min_dist) {
				min_dist = d;
			}
		}
		for (auto surf : planes) {
			d = surf.DE(v1);
			if (d < min_dist) {
				min_dist = d;
			}
		}
		material main_mat;
		min_dist = min(min_dist, DE_main(v1, main_mat));
		return min_dist;
	}
	double abs_DE(vec3 v1, material& closest_mat) {
		double min_dist = render_dist;
		double d;

		for (auto surf : spheres) {
			d = abs(surf.DE(v1));
			if (d < min_dist) {
				closest_mat.reflective = surf.reflective;
				closest_mat.transparent = surf.transparent;
				closest_mat.IOR = surf.IOR;
				min_dist = d;
			}
		}
		for (auto surf : planes) {
			d = abs(surf.DE(v1));
			if (d < min_dist) {
				closest_mat.reflective = surf.reflective;
				closest_mat.transparent = surf.transparent;
				closest_mat.IOR = surf.IOR;
				min_dist = d;
			}
		}
		material main_mat;
		d = abs(DE_main(v1, main_mat));
		if (d < min_dist) {
			closest_mat = main_mat;
			min_dist = d;
		}
		return min_dist;
	}
	double abs_DE(vec3 v1) {
		double min_dist = render_dist;
		double d;

		for (auto surf : spheres) {
			d = abs(surf.DE(v1));
			if (d < min_dist) {
				min_dist = d;
			}
		}
		for (auto surf : planes) {
			d = abs(surf.DE(v1));
			if (d < min_dist) {
				min_dist = d;
			}
		}
		material main_mat;
		min_dist = min(min_dist, abs(DE_main(v1, main_mat)));
		return min_dist;
	}
	vec3 calc_normal(vec3 v1) {
		double d_p = 1e-6;
		double cur_dist = DE(v1);
		vec3 plus_x = v1;
		vec3 plus_y = v1;
		vec3 plus_z = v1;
		plus_x.x += d_p;
		plus_y.y += d_p;
		plus_z.z += d_p;
		double dist_x = DE(plus_x) - cur_dist;
		double dist_y = DE(plus_y) - cur_dist;
		double dist_z = DE(plus_z) - cur_dist;
		vec3 norm_v = { dist_x, dist_y, dist_z };
		norm_v.norm();
		//norm_v.y += sin(v1.y * 60) * 0.1;
		//norm_v.norm();
		return norm_v;
	}
	vec3 abs_normal(vec3 v1) {
		double d_p = 1e-6;
		double cur_dist = abs_DE(v1);
		vec3 plus_x = v1;
		vec3 plus_y = v1;
		vec3 plus_z = v1;
		plus_x.x += d_p;
		plus_y.y += d_p;
		plus_z.z += d_p;
		double dist_x = abs(DE(plus_x)) - cur_dist;
		double dist_y = abs(DE(plus_y)) - cur_dist;
		double dist_z = abs(DE(plus_z)) - cur_dist;
		vec3 norm_v = { dist_x, dist_y, dist_z };
		norm_v.norm();

		return norm_v;
	}
	void march_light_ray(ray& cur_ray) {
		double step_len = 0;
		cur_ray.travel_dist = 0;
		double last_dist = 0;
		double second_dist = 0;
		cur_ray.closest_dist = render_dist;
		material mat_sample;
		for (int j = 0; j < max_steps; j++) {
			second_dist = last_dist;
			last_dist = step_len;

			step_len = DE(cur_ray.pos, mat_sample) * 0.98;


			if (abs(step_len) < contact_thresh) {
				cur_ray.complete = true;

			}

			else if (dist(cur_ray.pos, cam_pos) > render_dist) {
				cur_ray.complete = true;
			}
			if (cur_ray.complete) {
				break;
			}

			if (last_dist < step_len && second_dist > last_dist) {
				cur_ray.closest_dist = min(cur_ray.closest_dist, (step_len / last_dist) * pow((step_len - last_dist / 2) * (step_len + last_dist / 2), 0.5));
			}
			cur_ray.pos += cur_ray.dir * step_len;
			cur_ray.travel_dist += step_len;
		}
	}
	vec3 march_cam_ray(ray& cur_ray, int inc_count) {
		double step_len = 0;
		cur_ray.travel_dist = 0;
		double last_dist = 0;
		cur_ray.closest_dist = render_dist;
		cur_ray.cur_IOR = 1;
		material mat;
		double surf_dist;
		while (cur_ray.march_count < max_steps) {
			last_dist = step_len;

			step_len = DE(cur_ray.pos, mat);
			if (inc_count > max_incs) {
				cur_ray.color = bg_color;
				cur_ray.complete = true;
			}
			else if (mat.transparent && abs(step_len) < contact_thresh) {
				if (dot(abs_normal(cur_ray.pos), cur_ray.dir) < 0) {
					vec3 unit = { 1,1, 1 };
					ray reflected = cur_ray;
					reflected.march_count = 0;
					vec3 reflect_color;
					reflected.dir = reflect_ray(reflected);


					ray refracted = cur_ray;
					refracted.march_count = 0;
					vec3 refract_color;
					vec3 abs_norm = abs_normal(refracted.pos);
					double next_IOR = get_IOR(refracted.pos);
					refracted.cur_IOR = next_IOR;
					refracted.pos -= abs_norm * 2 * step_len;
					refracted.dir = refract_ray(refracted, cur_ray.cur_IOR, refracted.cur_IOR);
					//cur_ray.color = unit *255* dot(refracted.dir, abs_normal(refracted.pos));
					double cos_i = -dot(abs_norm, cur_ray.dir);
					double sin_i = sqrt(1 - cos_i * cos_i);
					double sin_t = sin_i * cur_ray.cur_IOR / refracted.cur_IOR;
					double R_s = pow(abs((cur_ray.cur_IOR * cos_i - refracted.cur_IOR * sqrt(1 - sin_t * sin_t)) / (cur_ray.cur_IOR * cos_i + refracted.cur_IOR * sqrt(1 - sin_t * sin_t))), 2);
					double R_p = pow(abs((cur_ray.cur_IOR * sqrt(1 - sin_t * sin_t) - refracted.cur_IOR * cos_i) / (cur_ray.cur_IOR * sqrt(1 - sin_t * sin_t) + refracted.cur_IOR * cos_i)), 2);
					double R_eff = (R_s + R_p) / 2;

					reflect_color = march_cam_ray(reflected, inc_count + 1);
					if (refracted.dir.x != -100) {
						refract_color = march_cam_ray(refracted, inc_count + 1);
						cur_ray.color += reflect_color * R_eff;
						cur_ray.color += refract_color * (1 - R_eff);

						//cur_ray.color = unit*255* R_eff;

					}
					else {
						cur_ray.color = reflect_color;
					}

					cur_ray.complete = true;

				}
			}
			else if (step_len < contact_thresh) {
				if (mat.reflective) {
					//cur_ray.color = { 200, 200, 200 };
					//cur_ray.color *= max(0, dot({ 0, 1, 0 }, reflect_ray(cur_ray)));
					//cur_ray.complete = true;
					if (dot(calc_normal(cur_ray.pos), cur_ray.dir) < 0) {
						cur_ray.dir = reflect_ray(cur_ray);

					}
				}
				else {
					cur_ray.color = get_color(cur_ray, mat);
					cur_ray.complete = true;

				}
			}
			else if (cur_ray.travel_dist > render_dist) {
				cur_ray.color = get_color(cur_ray, mat);
				cur_ray.complete = true;
			}
			if (cur_ray.complete) {
				return cur_ray.color;
			}
			cur_ray.pos += cur_ray.dir * step_len;
			cur_ray.travel_dist += step_len;
			cur_ray.march_count++;
		}
	}
	vec3 reflect_ray(ray cur_ray) {
		vec3 surf_norm = abs_normal(cur_ray.pos);
		//cout << dot(cur_ray.dir, surf_norm) << " ";
		cur_ray.dir -= surf_norm * 2 * dot(cur_ray.dir, surf_norm);
		//cout << dot(cur_ray.dir, surf_norm) << endl;
		return cur_ray.dir;
	}
	vec3 refract_ray(ray cur_ray, double IOR1, double IOR2) {
		vec3 surf_norm = abs_normal(cur_ray.pos);
		double cos1 = dot(surf_norm, cur_ray.dir);
		double sin1 = sqrt(1-cos1*cos1);
		double sin2 = sin1 * IOR1 / IOR2;
		if (sin2 >= 1) {
			return { -100, 0, 0 };
		}
		double adjusted_h = sqrt(pow(IOR2/IOR1, 2) - sin1*sin1);
		vec3 refract_dir = cur_ray.dir;
		refract_dir += surf_norm * (adjusted_h - cos1);
		refract_dir.norm();
		return refract_dir;
	}
	double get_IOR(vec3 p1) {
		material mat;
		DE(p1, mat);
		return mat.IOR;
	}
	void render_frame() {

		for (auto& r : rays) {

			march_cam_ray(r, 0);

		}
	}
	void print_frame() {
		for (int i = 0; i < window_res.y; i++) {
			for (int j = 0; j < window_res.x; j++) {
				char character;
				ray r = rays[i * window_res.x + j];

				if (r.color.x == 255) {
					character = '@';
				}
				else if (r.color.x == 0) {
					character = '.';
				}
				else {
					character = ' ';
				}
				printf("\033[%d;%dH%c", i, j, character);
			}
		}
	}

	vec3 phong_shading(ray cur_ray) {
		vec3 final_light = { 0, 0, 0 };
		for (auto light : light_points) {
			double l_brightness = 0;
			double ambient = 0.15;
			double diffuse = 0;
			double specular = 0;
			ray light_ray;
			light_ray.complete = false;
			light_ray.pos = light.pos;
			light_ray.dir = cur_ray.pos - light.pos;
			light_ray.dir.norm();
			march_light_ray(light_ray);
			vec3 surf_norm = calc_normal(cur_ray.pos);

			if (dist(cur_ray.pos, light_ray.pos) < contact_thresh * 4) {
				diffuse += max(-dot(surf_norm, light_ray.dir), 0.0);
				specular += max((-dot(cur_ray.dir, reflect_ray(light_ray)) - 0.8) * 5, 0.0) / 2;
				
				diffuse *= min(5 * light_ray.closest_dist, 1);
				specular *= min(5 * light_ray.closest_dist, 1);
				

			}
			l_brightness += ambient + diffuse + specular;
			l_brightness /= light_ray.travel_dist * light_ray.travel_dist;
			l_brightness *= light.brightness;
			//l_brightness = min(l_brightness, (255 - light_ray.march_count) / 255);
			final_light += light.color * l_brightness;
		}


		//cout << light_ray.dir.len();
		
		return final_light;
	}
	vec3 get_color(ray r, material mat) {
		/*
		vec3 base_color = { 255, 255, 255 };
		return base_color * ((abs((int)floor(p1.x)) % 2 + abs((int)floor(p1.z) % 2)) % 2);
		*/
		if (dist(r.pos, cam_pos) > render_dist) {
			return bg_color;
		}
		vec3 unit = mat.color;
		vec3 surf_lighting = phong_shading(r);
		double white_tile = ((int)(floor(r.pos.x) + floor(r.pos.y) + floor(r.pos.z)) % 2 == 0 ? 1 : 0.5);
		unit.x *= white_tile;
		unit.y *= white_tile;
		unit.z *= white_tile;
		unit.x *= surf_lighting.x;
		unit.y *= surf_lighting.y;
		unit.z *= surf_lighting.z;
		unit = v_min(unit, 255);
		return unit ;
		
	}
	void add_sphere(vec3 s_pos, double rad, bool reflective, double IOR) {
		sphere_surf surf;
		surf.pos = s_pos;
		surf.radius = rad;
		surf.reflective = reflective;
		if (IOR != 0) {
			surf.transparent = true;

		}
		else {
			surf.transparent = false;
		}
		surf.IOR = IOR;
		spheres.push_back(surf);
	}

	void add_plane(vec3 v1, vec3 v2, double y, bool reflective) {
		plane_surf surf;
		surf.v1 = v1;
		surf.v2 = v2;
		surf.y = y;
		surf.reflective = reflective;
		surf.transparent = false;
		planes.push_back(surf);
	}
	void add_light(vec3 v1, vec3 c1, double brightness) {
		light_point lpoint;
		lpoint.pos = v1;
		lpoint.color = c1;
		lpoint.brightness = brightness;
		light_points.push_back(lpoint);
	}
};

void save_image(Marcher march) {
	freopen("pic.ppm", "w", stdout);
	march.window_res = { 1000,1000 };
	march.max_steps = 500;
	march.max_incs = 7;
	march.initialize_rays();
	march.render_frame();
	cout << "P3 " << march.window_res.x << " " << march.window_res.y << " 255" << endl;
	for (int i = 0; i < march.window_res.y; i++) {
		for (int j = 0; j < march.window_res.x; j++) {
			ray r = march.rays[i * (int)march.window_res.x + j];
			cout << r.color.x << " " << r.color.y << " " << r.color.z << " ";
		}
		cout << endl;
	}
	fclose(stdout);
}


int main(void)
{
	double step_size = 0.4;
	bool quit = false;
	POINT p;
	Marcher march;
	SDL_Event e;
	march.window_res = { 30,30 };
	march.max_steps = 70;
	march.render_dist = 40;
	march.cam_angle = { 45, 0 };
	march.cam_pos = { 4, 1.5, 4 };
	
	march.add_plane({ 1,0,0 }, { 0,0,1 }, -1, false);
	march.add_sphere({ 2,4,2 }, 1, false, 1.1);
	//march.add_sphere({ 1,4,0}, 5, false, 0);
	//march.add_light({ -10, 10, 5 }, { 1, 1,  1 }, 30);
	
	mandelbulb m;
	march.add_light({ 10, 10, 0 }, { 1, 1,  1 }, 40);
	SDL_SetMainReady();
	if (!init(600, 1000))
	{
		printf("Failed to initialize!\n");
	}

	while (!quit) {
		GetCursorPos(&p);
		march.cam_angle.x = p.x / 4;
		march.cam_angle.y = 40-p.y / 6;
		march.initialize_rays();
		march.render_frame();
		//march.print_frame();
		//cout << "frame" << endl;
		while (SDL_PollEvent(&e) != 0)
		{
			switch (e.type)
			{
			case SDL_QUIT:

				quit = true;
				break;

			case SDL_KEYDOWN:
				switch (e.key.keysym.sym)
				{
				case SDLK_ESCAPE:  quit = true; break;
				case SDLK_SPACE:  march.cam_pos.y += step_size; break;
				case SDLK_LSHIFT:  march.cam_pos.y -= step_size; break;
				case SDLK_w:  march.cam_pos.x += cos(march.cam_angle_rads.x) * step_size; march.cam_pos.z += sin(march.cam_angle_rads.x) * step_size; break;
				case SDLK_a:  march.cam_pos.x += sin(march.cam_angle_rads.x) * step_size; march.cam_pos.z -= cos(march.cam_angle_rads.x) * step_size; break;
				case SDLK_s:  march.cam_pos.x -= cos(march.cam_angle_rads.x) * step_size; march.cam_pos.z -= sin(march.cam_angle_rads.x) * step_size; break;
				case SDLK_d:  march.cam_pos.x -= sin(march.cam_angle_rads.x) * step_size; march.cam_pos.z += cos(march.cam_angle_rads.x) * step_size; break;
				case SDLK_z:  march.fov *= 0.7; break;
				case SDLK_x:  march.fov *= 1.0 / 0.7; break;
				case SDLK_p:  cout << march.cam_pos.x << " " << march.cam_pos.y << " " << march.cam_pos.z << " " << march.cam_angle.x << " " << march.cam_angle.y << endl; break;
				case SDLK_o:  save_image(march); break;
				}
				break;
			}
		}
		SDL_RenderClear(gRenderer);
		for (int i = 0; i < march.window_res.y; i++) {
			for (int j = 0; j < march.window_res.x; j++) {
				char character;
				ray r = march.rays[i * march.window_res.x + j];

				SDL_Rect rect;
				rect.x = 2*j;
				rect.y = 2*i;
				rect.w = 2;
				rect.h = 2;
				SDL_SetRenderDrawColor(gRenderer, min(255, r.color.x), min(255, r.color.y), min(255, r.color.z), 0xFF);
				SDL_RenderFillRect(gRenderer, &rect);
				SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
			}
		}
		SDL_RenderPresent(gRenderer);
	}
	
	return 0;
}