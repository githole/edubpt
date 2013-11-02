﻿#ifndef _PATHTRACING_H_
#define _PATHTRACING_H_

#include <vector>

#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"
#include "vertex.h"

namespace edubpt {
	

struct PathtracingResult {
	Color value;
	int imagebuffer_x, imagebuffer_y;
	bool is_light_hit;

	PathtracingResult(const Color &value, const int image_x, const int image_y, const bool is_light_hit) :
	value(value), imagebuffer_x(imagebuffer_x), imagebuffer_y(imagebuffer_y), is_light_hit(is_light_hit) {}
};

PathtracingResult generate_vertices_by_pathtracing(const Camera &camera, const int imagebuffer_x, const int imagebuffer_y, Random *rnd, std::vector<Vertex> &vertices) {
	Vec position_on_imagesensor, position_on_objectplane, position_on_lens;
	// イメージセンサとレンズ上から最初の点をサンプリング
	double P_Image, P_lens;
	camera.sample_points(imagebuffer_x, imagebuffer_y, rnd, &position_on_imagesensor, &position_on_objectplane, &position_on_lens, &P_Image, &P_lens);

	double total_pdf_A = P_lens;
	Color MC_throughput(1, 1, 1); // PTの場合、これはインポータンス値（brdfとかのGとかの積）
	
	// レンズ上の頂点（x0）を頂点リストに追加
	vertices.push_back(Vertex(position_on_lens, camera.normal_on_lens, camera.normal_on_lens, -1, Vertex::OBJECT_TYPE_LENS, total_pdf_A, MC_throughput));
	
	Ray now_ray(position_on_lens, normalize(position_on_objectplane - position_on_lens));
	double now_sampled_pdf_omega = 1.0;
	Vec previous_normal = camera.normal_on_lens;

	for (int now_vertex_id = 1;; ++now_vertex_id) {
		Intersection intersection;
		if (!intersect_scene(now_ray, &intersection))
			break;
		
		const Sphere &now_object = spheres[intersection.object_id];
		const Hitpoint &hitpoint = intersection.hitpoint;
		const Vec orienting_normal = dot(hitpoint.normal , now_ray.dir) < 0.0 ? hitpoint.normal: (-1.0 * hitpoint.normal); // 交差位置の法線（物体からのレイの入出を考慮）
		const double russian_roulette_probability = russian_roulette(now_object);
		
		// ロシアンルーレット
		if (rnd->next01() >= russian_roulette_probability) {
			break;
		}
		// 新しい頂点がサンプリングされた
		total_pdf_A *= russian_roulette_probability;
		
		const Vec between = now_ray.org - hitpoint.position;
		if (now_vertex_id == 1) 
		{
			// x1だけは確率密度の計算が特殊になる
			const Vec x0_xI = position_on_imagesensor - position_on_lens;
			const Vec x0_xV = position_on_objectplane - position_on_lens;
			const Vec x0_x1 = hitpoint.position - position_on_lens;
			const double PA_x1 = camera.P_Image_to_PA_x1(P_Image, x0_xV, x0_x1, orienting_normal);
			total_pdf_A *= PA_x1;

			// さらに、MC_throughputの計算も若干特殊になる（カメラ周りの処理）
			MC_throughput = camera.W_dash(x0_xV, x0_xI, x0_x1) * MC_throughput;
		} else {
			// 新しい頂点を得たので、今までの頂点をサンプリングする確率の総計を出す
			const double now_sampled_pdf_A = now_sampled_pdf_omega * (dot(normalize(between), orienting_normal) / between.length_squared());
			total_pdf_A *= now_sampled_pdf_A;
		}
		// ジオメトリターム
		const double G =  dot(normalize(between), orienting_normal) * dot(normalize(-1.0 * between), previous_normal) / between.length_squared();
		MC_throughput = G * MC_throughput;

		// 光源にヒットしたらそこで追跡終了（光源の反射率がゼロであることを仮定しており、これ以上の追跡は全てMC_throughputがゼロになるため。光源にも反射率があるならこの限りではない）
		if (now_object.emission.length_squared() > 0) {
			vertices.push_back(Vertex(hitpoint.position, orienting_normal, hitpoint.normal, intersection.object_id, Vertex::OBJECT_TYPE_LIGHT, total_pdf_A, MC_throughput));
			const Color result = multiply(MC_throughput, now_object.emission) / total_pdf_A;
			return PathtracingResult(result, imagebuffer_x, imagebuffer_y, true);
		}

		// 新しい頂点を頂点リストに追加する
		vertices.push_back(Vertex(hitpoint.position, orienting_normal, hitpoint.normal, intersection.object_id, 
				now_object.reflection_type == REFLECTION_TYPE_DIFFUSE ? Vertex::OBJECT_TYPE_DIFFUSE :Vertex::OBJECT_TYPE_SPECULAR, 
				total_pdf_A, MC_throughput));

		switch (now_object.reflection_type) {
		// 完全拡散面
		case REFLECTION_TYPE_DIFFUSE: {
			const Vec dir = sample_hemisphere_cos_term(orienting_normal, rnd, &now_sampled_pdf_omega);
			now_ray = Ray(hitpoint.position, dir);
			MC_throughput = multiply(now_object.color / kPI, MC_throughput);
		} break;
		// 完全鏡面
		case REFLECTION_TYPE_SPECULAR: {
			// 完全鏡面なのでレイの反射方向は決定的。
			now_sampled_pdf_omega = 1.0;
			now_ray = Ray(hitpoint.position, reflection_vector(now_ray.dir, hitpoint.normal));

			// スループット調整
			// Lighttracingと同様、sample_pdg_omegaのdiracとキャンセルするためにこのような形になる
			MC_throughput = multiply(now_object.color / dot(normalize(between), orienting_normal), MC_throughput);
		} break;
		// ガラス
		case REFLECTION_TYPE_REFRACTION: {
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

			Vec reflection_dir, refraction_dir;
			double fresnel_reflectance, fresnel_transmittance;
			if (check_refraction(into, hitpoint.position, now_ray.dir, hitpoint.normal, orienting_normal, &reflection_dir, &refraction_dir, &fresnel_reflectance, &fresnel_transmittance)) {
				// 屈折 + 反射
				if (rnd->next01() < reflection_probability) { // 反射
					now_sampled_pdf_omega = 1.0;
					now_ray = Ray(hitpoint.position, reflection_dir);
					MC_throughput = fresnel_reflectance * multiply(now_object.color / dot(normalize(between), orienting_normal), MC_throughput);
					total_pdf_A *= reflection_probability;
				} else { // 屈折
					const double nnt2 = 
						pow(into ?
						refractive_index_of_vaccum / refractive_index_of_object :
						refractive_index_of_object / refractive_index_of_vaccum, 2.0); 
					now_sampled_pdf_omega = 1.0;
					now_ray = Ray(hitpoint.position, refraction_dir);
					MC_throughput = nnt2 * fresnel_transmittance * multiply(now_object.color / dot(normalize(between), orienting_normal), MC_throughput);
					total_pdf_A *= 1.0 - reflection_probability;
				}
			} else {
				// 全反射
				now_sampled_pdf_omega = 1.0;
				now_ray = Ray(hitpoint.position, reflection_dir);
				MC_throughput = multiply(now_object.color / dot(normalize(between), orienting_normal), MC_throughput);
				break;
			}
		} break;
		}

		previous_normal = orienting_normal;
	}

	return PathtracingResult(Color(), 0, 0, false);
}


};

#endif