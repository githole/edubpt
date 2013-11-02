#ifndef _RENDER_H_
#define _RENDER_H_

#include "reference_pathtracing.h"
#include "pathtracing.h"
#include "lighttracing.h"
#include "bidirectional_pathtracing.h"

#include "ppm.h"

#define USE_OPENMP

#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <assert.h>

namespace edubpt {

	
int render_by_refernce_pathtracing(const Camera &camera, const int num_threads) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height);
#ifdef USE_OPENMP
	omp_set_num_threads(num_threads);
#endif
	
	for (int sample = 0; sample < spp; ++sample) {
	std::cerr << "Rendering (" << (sample + 1) << " samples/pixel) " << (100.0 * (sample + 1) / spp) << "%" << std::endl;

#pragma omp parallel for schedule(dynamic, 1)
		for (int y = 0; y < height; ++y) {
			Random random(y * spp + sample);
			for (int x = 0; x < width; ++x) {
				Vec position_on_imagesensor, position_on_objectplane, position_on_lens;
				double P_Image, P_lens;
				camera.sample_points(x, y, &random, &position_on_imagesensor, &position_on_objectplane, &position_on_lens, &P_Image, &P_lens);

				// 方向
				const Vec dir = normalize(position_on_objectplane - position_on_lens);
				// L(x0 <- x1) (= L(xI <- x0) = L(x0 <- xV))
				const Color L = radiance_no_recursion(Ray(position_on_lens, dir), &random);
				if (!is_valid_value(L)) // 異常値は除外する
					continue;

				const Vec x0_xI = position_on_imagesensor - position_on_lens;
				const double coefficient = pow(dot(normalize(-1.0 * camera.imagesensor_dir), normalize(x0_xI)), 2.0) / x0_xI.length_squared(); // センサ上の点の半球積分を、レンズ上の点の積分に変数変換した時に導入される係数

				image_buffer[y * width + x] = image_buffer[y * width + x] +
					camera.W * L * coefficient / (P_Image * P_lens);
			}
		}
		/*
		if (sample % 128 == 0) {
			char buf[256];
			sprintf(buf, "ref_pt_%03d.hdr", sample);
			for (int y = 0; y < height; ++y)
				for (int x = 0; x < width; ++x)
					*hdr.image_ptr(width - x - 1, height - y - 1) = image_buffer[y * width + x] / ((double)width * height * (sample + 1));
			hdr.save(buf);
		}
		*/
	}

	// サンプル数で割る + 左右反転
	std::vector<Color> finale_buffer(width * height);
	for (int y = 0; y < height; ++y)
		for (int x = 0; x < width; ++x)
			finale_buffer[y * width + x] = image_buffer[y * width + (width - x - 1)] / ((double)width * height * spp);
	// 出力
	save_ppm_file("image_ref_pt.ppm", &finale_buffer[0], width, height);

	return 0;
}

int render_by_pathtracing(const Camera &camera, const int num_threads) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height);
#ifdef USE_OPENMP
	omp_set_num_threads(num_threads);
#endif
	
	for (int sample = 0; sample < spp; ++sample) {
		std::cerr << "Rendering (" << (sample + 1) << " samples/pixel) " << (100.0 * (sample + 1) / spp) << "%" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
		for (int y = 0; y < height; ++y) {
			Random random(y * spp + sample);
			for (int x = 0; x < width; ++x) {
				std::vector<Vertex> vertices;
				const PathtracingResult result = generate_vertices_by_pathtracing(camera, x, y, &random, vertices);

				if (result.is_light_hit) {
					if (is_valid_value(result.value)) {
						image_buffer[y * width + x] = image_buffer[y * width + x] + result.value;
					}
				}
			}
		}
		/*
		if (sample % 128 == 0) {
			char buf[256];
			sprintf(buf, "pt_%03d.hdr", sample);
			for (int y = 0; y < height; ++y)
				for (int x = 0; x < width; ++x) {
					// std::cout << image_buffer[y * width + x].x << " ";
					*hdr.image_ptr(width - x - 1, height - y - 1) = image_buffer[y * width + x] / ((double)width * height * (sample + 1));
				}
			hdr.save(buf);
		}
		*/
	}
	
	// サンプル数で割る + 左右反転
	std::vector<Color> finale_buffer(width * height);
	for (int y = 0; y < height; ++y)
		for (int x = 0; x < width; ++x)
			finale_buffer[y * width + x] = image_buffer[y * width + (width - x - 1)] / ((double)width * height * spp);
	// 出力
	save_ppm_file("image_pt.ppm", &finale_buffer[0], width, height);

	return 0;
}

void render_by_lighttracing(const Camera &camera, const int num_threads) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height * num_threads);
#ifdef USE_OPENMP
	omp_set_num_threads(num_threads);
#endif
	
	const long long samples_per_thread = (long long)spp * width * height / num_threads;
	std::cerr << "Average samples/pixel: " << (samples_per_thread * num_threads / (width * height)) << std::endl;
	int thread_id = 0;
#pragma omp parallel private(thread_id)
	{
#ifdef USE_OPENMP
		thread_id = omp_get_thread_num();
#endif
		Random random(thread_id * 32);

		for (long long sample = 0; sample < samples_per_thread; ++sample) {
			if (thread_id == 0 && sample % (width * height) == 0)
				std::cerr << "Rendering (" << ((double)(num_threads * (sample + 1)) / (width * height)) << " samples/pixel) "
				<< (100.0 * (double)((sample + 1)) / samples_per_thread) << "%" << std::endl;

			std::vector<Vertex> vertices;
			LighttracingResult result = generate_vertices_by_lighttracing(camera, &random, vertices);
			if (result.is_lens_hit) {
				if (is_valid_value(result.value)) {
					const int x = result.imagebuffer_x, y = result.imagebuffer_y;
					image_buffer[(thread_id * width * height) + y * width + x] = image_buffer[(thread_id * width * height) + y * width + x] + result.value;
				}
			}

			/*
			if (sample % (16 * width * height) == 0) {
#pragma omp barrier
				if (thread_id == 0) {
					std::cout << sample << " ";
					char buf[256];
					sprintf(buf, "lt_%03d.hdr", sample / (width * height));
					for (int y = 0; y < height; ++y)
						for (int x = 0; x < width; ++x)
							*hdr.image_ptr(x, y) = Color();
					for (int t = 0; t < num_threads; ++t)
						for (int y = 0; y < height; ++y)
							for (int x = 0; x < width; ++x)
								*hdr.image_ptr(width - x - 1, height - y - 1) = *hdr.image_ptr(width - x - 1, height - y - 1) + image_buffer[(t * width * height) + y * width + x] / (sample + 1) / num_threads;
					hdr.save(buf);
				}
#pragma omp barrier
			}
			*/
		}
	}
	
	// サンプル数で割る + 左右反転
	std::vector<Color> final_buffer(width * height);
	for (int t = 0; t < num_threads; ++t)
		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x) {
				final_buffer[y * width + x] = final_buffer[y * width + x] + image_buffer[(t * width * height) + y * width + (width - x - 1)] / ((double)samples_per_thread * num_threads);
			}
	// 出力
	save_ppm_file("image_lt.ppm", &final_buffer[0], width, height);
}


void render_by_bidirectional_pathtracing(const Camera &camera, const int num_threads) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height * num_threads);
#ifdef USE_OPENMP
	omp_set_num_threads(num_threads);
#endif
	
	const long long iteration_per_thread = spp / num_threads;
	std::cerr << "Average bidirectional sampling/pixel: " << (iteration_per_thread * num_threads) << std::endl;
	int thread_id = 0;
#pragma omp parallel private(thread_id)
	{
		thread_id = omp_get_thread_num();
		Random random(thread_id * 32);
		
		for (int iteration = 0; iteration < iteration_per_thread; ++iteration) {
			if (thread_id == 0)
				std::cerr << "Rendering (" << ((double)(num_threads * (iteration + 1))) << " bidirectional sampling/pixel) "
				<< (100.0 * (double)((iteration + 1)) / iteration_per_thread) << "%" << std::endl;

			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					BidirectionalPathtracingResult bpt_result = bidirectional_pathtracing(camera, x, y, &random);

					for (int i = 0; i < bpt_result.data.size(); ++i) {
						const int ix = bpt_result.data[i].imagebuffer_x;
						const int iy = bpt_result.data[i].imagebuffer_y;
						const int idx = (thread_id * width * height) + (iy * width + ix);
						if (is_valid_value(bpt_result.data[i].value))
							image_buffer[idx] = image_buffer[idx] + bpt_result.data[i].value;
					}
				}
			}
			/*
			if (1) {
#pragma omp barrier
				if (thread_id == 0) {
					std::cout << sample << " ";
					char buf[256];
					sprintf(buf, "bpt_%03d.hdr", sample);
					for (int y = 0; y < height; ++y)
						for (int x = 0; x < width; ++x)
							*hdr.image_ptr(x, y) = Color();
					for (int t = 0; t < num_threads; ++t)
						for (int y = 0; y < height; ++y)
							for (int x = 0; x < width; ++x) {
								const int idx = (t * width * height) + (y * width + x);
								*hdr.image_ptr(width - x - 1, height - y - 1) = 
									*hdr.image_ptr(width - x - 1, height - y - 1) + image_buffer[idx]	 / ((double) width * height * (sample + 1)) / num_threads;
							}
					hdr.save(buf);
				}
#pragma omp barrier
			}
			*/
		}
	}
	// サンプル数で割る + 左右反転
	std::vector<Color> final_buffer(width * height);
	for (int t = 0; t < num_threads; ++t)
		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x) {
				final_buffer[y * width + x] = final_buffer[y * width + x] + image_buffer[(t * width * height) + y * width + (width - x - 1)] /
					((double)width * height * iteration_per_thread * num_threads);
			}
	// 出力
	save_ppm_file("image_bpt.ppm", &final_buffer[0], width, height);
}


};

#endif
