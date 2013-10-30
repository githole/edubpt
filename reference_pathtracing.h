#ifndef _REFERENCE_PATHTRACING_
#define _REFERENCE_PATHTRACING_

#include "hdr.h"
#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"

namespace edubpt {

// L(xi <- xi+1) を計算する
Color radiance_no_recursion(const Ray &ray, Random *rnd) {
	Color L;
	Ray now_ray = ray;
	Color weight(1, 1, 1);
	double pdf = 1.0;

	for (;;) {
		Intersection intersection;
		if (!intersect_scene(now_ray, &intersection)) {
			L = L + multiply(weight, BackgroundColor) / pdf;
			break;
		}
		// ロシアンルーレット
		const Sphere &now_object = spheres[intersection.object_id];
		const double russian_roulette_probability = russian_roulette(now_object); // using only local information of the vertex
		if (rnd->next01() >= russian_roulette_probability) {
			break;
		}
		pdf *= russian_roulette_probability;

		// 影響加算
		L = L + multiply(weight, now_object.emission) / pdf;
	
		const Hitpoint &hitpoint = intersection.hitpoint;
		const Vec orienting_normal = dot(hitpoint.normal , now_ray.dir) < 0.0 ? hitpoint.normal: (-1.0 * hitpoint.normal); // 交差位置の法線（物体からのレイの入出を考慮）

		switch (now_object.reflection_type) {
		// 完全拡散面
		case REFLECTION_TYPE_DIFFUSE: {
			double pdf_omega;
			const Vec dir = sample_hemisphere_cos_term(orienting_normal, rnd, &pdf_omega);
			now_ray = Ray(hitpoint.position, dir);
			
			const double cos_term = dot(orienting_normal, dir);
			weight = multiply(weight, (now_object.color / kPI) * cos_term); // LambertianのBRDF = ρ/π と cos項
			pdf    *= pdf_omega;
		} break;
		// 完全鏡面
		case REFLECTION_TYPE_SPECULAR: {
			// 完全鏡面なのでレイの反射方向は決定的。
			now_ray = Ray(hitpoint.position, reflection_vector(now_ray.dir, hitpoint.normal));
			weight = multiply(weight, now_object.color);
			pdf    *= 1.0;
		} break;
		// 屈折率kIorのガラス
		case REFLECTION_TYPE_REFRACTION: {
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

			Vec reflection_dir, refraction_dir;
			double fresnel_reflectance, fresnel_transmittance;
			if (check_refraction(into, hitpoint.position, now_ray.dir, hitpoint.normal, orienting_normal, &reflection_dir, &refraction_dir, &fresnel_reflectance, &fresnel_transmittance)) {
				// 屈折 + 反射
				if (rnd->next01() < reflection_probability) { // 反射
					now_ray = Ray(hitpoint.position, reflection_dir);
					weight  = fresnel_reflectance * multiply(weight, now_object.color);
					pdf    *= reflection_probability;
				} else { // 屈折
					// レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
					const double nnt2 = 
						pow(into ?
						refractive_index_of_vaccum / refractive_index_of_object :
						refractive_index_of_object / refractive_index_of_vaccum, 2.0); 
					now_ray = Ray(hitpoint.position, refraction_dir);
					weight  = nnt2 * fresnel_transmittance * multiply(weight, now_object.color);
					pdf    *= 1.0 - reflection_probability;
				}
			} else {
				// 全反射
				now_ray = Ray(hitpoint.position, reflection_dir);
				weight  = multiply(weight, now_object.color);
				pdf    *= 1.0;
				break;
			}
		} break;
		}
	}

	return L;
}

};

#endif