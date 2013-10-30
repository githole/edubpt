#ifndef _HDR_H_
#define _HDR_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "material.h"

namespace edubpt {

struct HDRImage {
private:
	struct HDRPixel {
		unsigned char r, g, b, e;
		explicit HDRPixel(const unsigned char r = 0, const unsigned char g = 0, const unsigned char b = 0, const unsigned char e = 0) :
		r(r), g(g), b(b), e(e) {};
		unsigned char get(int idx) {
			switch (idx) {
			case 0: return r;
			case 1: return g;
			case 2: return b;
			case 3: return e;
			} return 0;
		}

		explicit HDRPixel(const Color &color) {
			double d = std::max(color.x, std::max(color.y, color.z));
			if (d <= 1e-32) {
				r = g = b = e = 0;
				return;
			}
			int et;
			const double m = frexp(d, &et); // d = m * 2^e
			d = m * 256.0 / d;

			r = (unsigned char)(color.x * d);
			g = (unsigned char)(color.y * d);
			b = (unsigned char)(color.z * d);
			e = (unsigned char)(et + 128);
		}
	};
	std::vector<Color> image_;
	int width_;
	int height_;
public:

	void clone(std::vector<Color> *target) {
		*target = image_;
	}

	HDRImage(const int width, const int height) :
	width_(width), height_(height) {
		image_.resize(width_ * height_);

	}

	HDRImage() {
	}

	inline int width() {
		return width_;
	}
	inline int height() {
		return height_;
	}

	inline void set(const int x, const int y, const Color &v) {
		image_[y * width_ + x] = v;
	}
	inline Color* image_ptr(const int x, const int y) {
		return &image_[y * width_ + x];
	}


	void load_unsafe(const std::string &filename) {
		FILE *fp = fopen(filename.c_str(), "rb");
		const int BufSize = 4096;
		char buf[BufSize];
		if (fp == NULL) {
			std::cerr << "Error: " << filename << std::endl;
			return;
		}

		bool valid = false;
		enum FileType {
			NONE,
			RLE_RGBE_32,
		} type;

		type = NONE;
		double exposure = 1.0;

		// ヘッダ読み込み
		for (;;) {
			fgets(buf, BufSize, fp);
			std::cout << buf;
			if (buf[0] == '#') {
				if (strcmp(buf, "#?RADIANCE\n") == 0)
					valid = true;
			} else {
				if (strstr(buf, "FORMAT=") == buf) {
					char buf2[BufSize];
					sscanf(buf, "FORMAT=%s", buf2);
					if (strcmp(buf2, "32-bit_rle_rgbe") == 0)
						type = RLE_RGBE_32;

				} else if (strstr(buf, "EXPOSURE=") == buf) {
					sscanf(buf, "FORMAT=%f", &exposure);
				}
			} 
			
			if (buf[0] == '\n')
				break;
		}

		if (!valid) {
			std::cerr << "Invalid HDR File: " << filename << std::endl;
			return;
		}

		// サイズ読み込み
		char buf2[BufSize], buf3[BufSize];
		int width, height;
		fgets(buf, BufSize, fp);
		sscanf(buf, "%s %d %s %d", buf2, &height, buf3, &width);

		// -Y @ +X @ しかあつかわないことにする　
		if (strcmp(buf2, "-Y") != 0 && strcmp(buf3, "+X") != 0) {
			std::cerr << "Invalid HDR File: " << filename << " " << buf2 << " " << buf3 << std::endl;
			return;
		}

		image_.resize(width * height);
		width_ = width;
		height_= height;

		std::vector<unsigned char> tmp_data(width_ * height_ * 4);

		int nowy = 0;
		for (;;) {
			const int now = fgetc(fp);
			if (now == EOF)
				break;
			const int now2 = fgetc(fp);

			if (now != 0x02 || now2 != 0x02) {
				break;
			}

			const int width = (fgetc(fp) << 8) + fgetc(fp);

			int nowx = 0;
			int nowvalue = 0;
			for (;;) {
				if (nowx >= width) {
					nowvalue ++;
					nowx = 0;
					if (nowvalue == 4)
						break;
				}

				const int info = fgetc(fp);
				if (info <= 128) {
					for (int i = 0; i < info; ++i) {
						const int data = fgetc(fp);
						tmp_data[(nowy * width_ + nowx) * 4 + nowvalue] = data;
						nowx ++;
					}
				} else {
					const int num = info - 128;
					const int data = fgetc(fp);
					for (int i = 0; i < num; ++i) {
						tmp_data[(nowy * width_ + nowx) * 4 + nowvalue] = data;
						nowx ++;
					}
				}
			}

			nowy ++;
		}

		// 展開
		for (int y = 0; y < height_; ++y) {
			for (int x = 0; x < width_; ++x) {
				const int e = tmp_data[(y * width_ + x) * 4 + 3];
				image_[(height_ - 1 - y) * width_ + x].x = tmp_data[(y * width_ + x) * 4 + 0] * pow(2, e - 128.0) / 256.0;
				image_[(height_ - 1 - y) * width_ + x].y = tmp_data[(y * width_ + x) * 4 + 1] * pow(2, e - 128.0) / 256.0;
				image_[(height_ - 1 - y) * width_ + x].z = tmp_data[(y * width_ + x) * 4 + 2] * pow(2, e - 128.0) / 256.0;

			}
		}
	}

	void save(const std::string &filename) {

		FILE *fp = fopen(filename.c_str(), "wb");
		if (fp == NULL) {
			std::cerr << "Error: " << filename << std::endl;
			return;
		}
		// .hdrフォーマットに従ってデータを書きだす
		// ヘッダ
		unsigned char ret = 0x0a;
		fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
		fprintf(fp, "# Made with hole's renderer%c", ret);
		fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
		fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

		// 輝度値書き出し
		fprintf(fp, "-Y %d +X %d%c", height_, width_, ret);
		for (int i = height_ - 1; i >= 0; i --) {
			std::vector<HDRPixel> line;
			for (int j = 0; j < width_; j ++) {
				HDRPixel p(image_[j + i * width_]);
				line.push_back(p);
			}
			fprintf(fp, "%c%c", 0x02, 0x02);
			fprintf(fp, "%c%c", (width_ >> 8) & 0xFF, width_ & 0xFF);
			for (int i = 0; i < 4; i ++) {
				for (int cursor = 0; cursor < width_;) {
					const int cursor_move = std::min(127, width_ - cursor);
					fprintf(fp, "%c", cursor_move);
					for (int j = cursor; j < cursor + cursor_move; j ++)
						fprintf(fp, "%c", line[j].get(i));
					cursor += cursor_move;
				}
			}
		}

		fclose(fp);
	}
};

}

#endif