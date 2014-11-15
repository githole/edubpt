#ifndef _BIDIRECTIONAL_PATHTRACING_H_
#define _BIDIRECTIONAL_PATHTRACING_H_

#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"
#include "vertex.h"

#include "pathtracing.h"
#include "lighttracing.h"

#include <assert.h>
#include <iostream>

namespace edubpt {

// 頂点fromから頂点nextをサンプリングしたとするとき、面積測度に関するサンプリング確率密度を計算する。
// 今回はLighttracingもPathtracingも同じ戦略でサンプリングを行うため一つの関数で良い。
// サンプリング戦略が異なるなら、それぞれのための関数を用意する必要がある。
double calc_pdf_A(const Camera &camera, const std::vector<const Vertex*> &vs, const int prev_from_idx, const int from_idx, const int next_idx) {
	const Vertex &from_vertex = *vs[from_idx];
	const Vertex &next_vertex = *vs[next_idx];
	Vertex const *prev_from_vertex = NULL;
	if (0 <= prev_from_idx && prev_from_idx < vs.size())
		prev_from_vertex = vs[prev_from_idx];

	const Vec to = next_vertex.position - from_vertex.position;
	const Vec normalized_to = normalize(to);
	double pdf = 0.0;

	// 頂点fromのオブジェクトの種類によって、その次の頂点のサンプリング確率密度の計算方法が変わる。
	switch (from_vertex.type) {
	case Vertex::OBJECT_TYPE_LIGHT:
	case Vertex::OBJECT_TYPE_DIFFUSE:
		pdf = sample_hemisphere_pdf_omega(from_vertex.orienting_normal, normalized_to);
		break;
	case Vertex::OBJECT_TYPE_LENS:
		{
			// レンズ上の点からシーン上の点をサンプリングするときの面積測度に関する確率密度を計算。
			// イメージセンサ上の点のサンプリング確率密度を元に変換する。
			const Ray test_ray(next_vertex.position, -1.0 * normalized_to);
			Vec position_on_lens, position_on_objectplane, position_on_imagesensor, uv_on_imagesensor;
			const double lens_t = camera.intersect_lens(test_ray, &position_on_lens, &position_on_objectplane, &position_on_imagesensor, &uv_on_imagesensor);
			if (kEPS < lens_t) {
				const Vec x0_xI = position_on_imagesensor - position_on_lens;
				const Vec x0_xV = position_on_objectplane - position_on_lens;
				const Vec x0_x1 = test_ray.org - position_on_lens;
					
				const double P_Image = 1.0 / (camera.imagesensor_width * camera.imagesensor_height);
				const double PA_x1 = camera.P_Image_to_PA_x1(P_Image, x0_xV, x0_x1, next_vertex.orienting_normal);
				return PA_x1;
			} 
			return 0.0;
		}
		break;
	case Vertex::OBJECT_TYPE_SPECULAR:
		if (spheres[from_vertex.objectID].reflection_type == REFLECTION_TYPE_GLASS) {			
			// 始点がガラス上に存在した場合、次の点が屈折によってサンプリングされたのか、反射によってサンプリングされたかによって
			// 次の頂点のサンプリング確率が変わるので、その判定を行う。
			if (prev_from_vertex != NULL) {
				// from頂点のひとつ前の頂点に基づいて、物体に入り込んでいるのか、それとも出て行くのかを判定する。
				const Vec into_from_vertex_dir = normalize(from_vertex.position - prev_from_vertex->position);
				const bool into = dot(into_from_vertex_dir, from_vertex.object_normal) < 0.0; // prev_from_vertexからfrom_vertexへのベクトルが、スペキュラ物体に入るのか、それとも出るのか
				const Vec from_new_orienting_normal = into ? from_vertex.object_normal : -1.0 * from_vertex.object_normal;
				
				// 全反射しているのか、していないのかを判定
				Vec reflection_dir, refraction_dir;
				double fresnel_reflectance, fresnel_transmittance;
				if (check_refraction(
					into, from_vertex.position, into_from_vertex_dir, from_vertex.object_normal, from_new_orienting_normal,
					&reflection_dir, &refraction_dir, &fresnel_reflectance, &fresnel_transmittance)) {
					// 反射 or 屈折
					pdf = dot(from_new_orienting_normal, normalized_to) > 0.0 ? reflection_probability : 1.0 - reflection_probability;
				} else {
					// 全反射
					pdf = 1.0;
				}
			}
		} else {
			// 鏡面
			pdf = 1.0;
		}
		break;
	}

	// 次の頂点の法線を、現在の頂点からの方向ベクトルに基づいて改めて求める
	const Vec next_new_orienting_normal = dot(to, next_vertex.object_normal) < 0.0 ? next_vertex.object_normal : -1.0 * next_vertex.object_normal;
	// 立体角測度に関する確率密度を面積測度に関する確率密度に変換
	return pdf * (dot(-1.0 * normalized_to, next_new_orienting_normal) / to.length_squared());
}

double calc_mis_weight(const Camera &camera, const double total_pdf_A, const std::vector<Vertex> &eye_vs, const std::vector<Vertex> &light_vs, const int num_eye_vertex, const int num_light_vertex) {
	std::vector<const Vertex*> vs(num_light_vertex + num_eye_vertex);
	std::vector<double> pi1_pi(num_eye_vertex + num_light_vertex);
	const double PA_y0 = sample_sphere_pdf_A(spheres[LightID].radius); // 光源上のサンプリング確率
	const double PA_x0 = camera.sampling_pdf_on_lens(); // レンズ上のサンプリング確率

	// 頂点を一列に並べる。
	// vs[0] = y0, vs[1] = y1, ... vs[k-1] = x1, vs[k] = x0
	const int k = num_eye_vertex + num_light_vertex - 1;
	for (int i = 0; i < num_light_vertex; ++i)
		vs[i] = &light_vs[i];
	for (int i = num_eye_vertex - 1; i >= 0; --i)
		vs[num_light_vertex + num_eye_vertex - 1 - i] = &eye_vs[i];

	// ロシアンルーレットの確率を忘れずに入れる。
	pi1_pi[0] = PA_y0 / (calc_pdf_A(camera, vs, 2, 1, 0) * russian_roulette(spheres[vs[0]->objectID]));
	for (int i = 1; i < k; ++i) {
		const double a = calc_pdf_A(camera, vs, i - 2, i - 1, i);
		const double b = calc_pdf_A(camera, vs, i + 2, i + 1, i);
		pi1_pi[i] = a / b;
	}
	pi1_pi[k] = (calc_pdf_A(camera, vs, k - 2, k - 1, k) * russian_roulette(spheres[vs[k]->objectID])) / PA_x0;

	// pを求める
	std::vector<double> p(num_eye_vertex + num_light_vertex + 1);
	p[num_light_vertex] = total_pdf_A;
	for (int i = num_light_vertex; i <= k; ++i) {
		p[i + 1] = p[i] * pi1_pi[i];
	}
	for (int i = num_light_vertex - 1; i >= 0; --i) {
		p[i] = p[i + 1] / pi1_pi[i];
	}

	// Specular処理
	for (int i = 0; i < vs.size(); ++i) {
		if (vs[i]->type == Vertex::OBJECT_TYPE_SPECULAR) {
			p[i] = 0.0;
			p[i + 1] = 0.0;
		}
	}

	// Power-heuristic
	double mis_weight = 0.0;
	for (int i = 0; i < p.size(); ++i) {
		const double v = p[i] / p[num_light_vertex];
		mis_weight += v * v; // beta = 2
	}
	return 1.0 / mis_weight;
}

struct BidirectionalPathtracingResult {
	struct Sample {
		int imagebuffer_x, imagebuffer_y;
		Color value;
		bool start_from_pixel;
		
		Sample(const int imagebuffer_x, const int imagebuffer_y, const Color &value, const bool start_from_pixel) :
		imagebuffer_x(imagebuffer_x), imagebuffer_y(imagebuffer_y), value(value), start_from_pixel(start_from_pixel) {}
	};
	std::vector<Sample> samples;
};
BidirectionalPathtracingResult bidirectional_pathtracing(const Camera &camera, const int imagebuffer_x, const int imagebuffer_y, Random *rnd) {
	BidirectionalPathtracingResult bpt_result;

	std::vector<Vertex> eye_vs, light_vs;
	const PathtracingResult  pt_result = generate_vertices_by_pathtracing(camera, imagebuffer_x, imagebuffer_y, rnd, eye_vs);
	const LighttracingResult lt_result = generate_vertices_by_lighttracing(camera, rnd, light_vs);

	// num_light_vertex == 0のとき、カメラ側からのパスが光源にヒットしていれば、それをそのまま使用する。
	if (pt_result.is_light_hit) {
		const double mis_weight = calc_mis_weight(camera, eye_vs[eye_vs.size()-1].total_pdf_A, eye_vs, light_vs, (const int)eye_vs.size(), 0);
		const Color result = mis_weight * pt_result.value;
		bpt_result.samples.push_back(BidirectionalPathtracingResult::Sample(imagebuffer_x, imagebuffer_y, result, true));
	}
	// num_eye_vertex == 0のとき、光源側からのパスがレンズにヒットしていれば、それをそのまま使用する。
	if (lt_result.is_lens_hit) {
		const double mis_weight = calc_mis_weight(camera, light_vs[light_vs.size()-1].total_pdf_A, eye_vs, light_vs, 0, (const int)light_vs.size());
		const int lx = lt_result.imagebuffer_x, ly = lt_result.imagebuffer_y;
		const Color result =  mis_weight * lt_result.value;
		bpt_result.samples.push_back(BidirectionalPathtracingResult::Sample(lx, ly, result, false));
	}

	// 各頂点間を接続する。
	// eyeサブパス:カメラ側のパスの頂点の一部を使って作られるパス
	// lightサブパス:光源側のパスの頂点の一部を使って作られるパス
	for (int num_eye_vertex = 1; num_eye_vertex <= eye_vs.size(); ++num_eye_vertex) {
		for (int num_light_vertex = 1; num_light_vertex <= light_vs.size(); ++num_light_vertex) {
			int target_x = imagebuffer_x, target_y = imagebuffer_y;
			const Vertex &eye_end = eye_vs[num_eye_vertex - 1];
			const Vertex &light_end = light_vs[num_light_vertex - 1];

			// トータルの確率密度計算
			const double total_pdf_A = eye_end.total_pdf_A * light_end.total_pdf_A;
			if (total_pdf_A == 0)
				continue;
					
			// MCスループット
			Color eye_throughput = eye_end.throughput;
			Color light_throughput = light_end.throughput;
			// 頂点を接続することで新しく導入される項（基本的にBRDFとG項）
			Color connection_throughput(1, 1, 1);
			
			// num_light_vertex == 1のとき、非完全拡散光源の場合は相手の頂点の位置次第でMCスループットが変化するため改めて光源からの放射輝度値を計算する。
			// 今回は完全拡散光源なので単純にemissionの値を入れる。
			if (num_light_vertex == 1) {
				light_throughput = spheres[light_vs[0].objectID].emission;
			}

			// 端点間でレイトレーシングしておく。
			Intersection intersection;
			const Vec light_end_to_eye_end = eye_end.position - light_end.position;
			const Ray test_ray = Ray(light_end.position, normalize(light_end_to_eye_end));
			intersect_scene(test_ray, &intersection);
					
			// eyeサブパスのスループット調整のための重み計算（主にBRDFの考慮のため）
			// eyeサブパスの端点の材質の種類によって変わる。
			switch (eye_end.type) {
			case Vertex::OBJECT_TYPE_DIFFUSE:
				connection_throughput = multiply(connection_throughput, spheres[eye_end.objectID].color / kPI);
				// 端点同士が別の物体で遮蔽されるかどうかを判定する。遮蔽されていたら処理終わり。
				if ((intersection.hitpoint.position - eye_end.position).length() >= kEPS) {
					continue;
				}
				break;
			case Vertex::OBJECT_TYPE_LENS:
				// ここに入るのはnum_eye_vertex == 1のときだけ
				{
					// lightサブパスを直接レンズにつなげるパターン		
					Vec position_on_lens, position_on_objectplane, position_on_imagesensor, uv_on_imagesensor;
					const double lens_t = camera.intersect_lens(test_ray, &position_on_lens, &position_on_objectplane, &position_on_imagesensor, &uv_on_imagesensor);		
					if (kEPS < lens_t && lens_t < intersection.hitpoint.distance) {
						// レイがレンズにヒット＆イメージセンサにヒット
						const Vec x0_xI = position_on_imagesensor - position_on_lens;
						const Vec x0_xV = position_on_objectplane - position_on_lens;
						const Vec x0_x1 = test_ray.org - position_on_lens;

						target_x = (int)uv_on_imagesensor.x;
						target_y = (int)uv_on_imagesensor.y;
						clamp(target_x, 0, camera.image_width_px);
						clamp(target_y, 0, camera.image_height_px);

						connection_throughput = camera.W_dash(x0_xV, x0_xI, x0_x1) * connection_throughput;
					} else {
						// lightサブパスを直接レンズにつなげようとしたが、遮蔽されたりイメージセンサにヒットしなかった場合、終わり。
						continue;
					}
				}
				break;
				// eyeサブパスの端点がスペキュラや光源（反射率0）だった場合は重みがゼロになりパス全体の寄与もゼロになるので、処理終わり。
			case Vertex::OBJECT_TYPE_LIGHT:
			case Vertex::OBJECT_TYPE_SPECULAR:
				continue;
			}
			
			// lightサブパスのスループット調整のための重み計算（主にBRDFの考慮のため）
			// lightサブパスの端点の材質の種類によって変わる。
			switch (light_end.type) {
			case Vertex::OBJECT_TYPE_DIFFUSE:
				connection_throughput = multiply(connection_throughput, spheres[light_end.objectID].color / kPI);
				break;
			case Vertex::OBJECT_TYPE_LIGHT:
				// 光源の反射率0を仮定しているため、ライトトレーシングの時、最初の頂点以外は光源上に頂点生成されず、
				// num_light_vertex == 1以外でここに入ってくることは無い。
				break;
				// lightサブパスの端点がスペキュラやレンズだった場合は重みがゼロになりパス全体の寄与もゼロになるので、処理終わり。
			case Vertex::OBJECT_TYPE_LENS:
			case Vertex::OBJECT_TYPE_SPECULAR:
				continue;
			}

			// 端点間のジオメトリターム
			const double G = (
				std::max(dot(normalize(-1.0*light_end_to_eye_end), eye_end.orienting_normal), 0.0) *
				std::max(dot(normalize(light_end_to_eye_end), light_end.orienting_normal), 0.0)) /
				light_end_to_eye_end.length_squared();
			connection_throughput = G * connection_throughput;

			// MISの重みを計算する。
			const double mis_weight = calc_mis_weight(camera, total_pdf_A, eye_vs, light_vs, num_eye_vertex, num_light_vertex);

			if (isnan(mis_weight)) {
				continue;
			}

			// 最終的なモンテカルロコントリビューション =
			// MIS重み * 頂点を接続することで新しく導入される項 * eyeサブパススループット * lightサブパススループット / パスのサンプリング確率密度の総計
			const Color result = mis_weight * multiply(multiply(connection_throughput, eye_throughput), light_throughput) / total_pdf_A;
			bpt_result.samples.push_back(BidirectionalPathtracingResult::Sample(target_x, target_y, result, num_eye_vertex <= 1 ? false : true));
		}
	}

	return bpt_result;
}

};

#endif
