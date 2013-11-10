#ifndef _CAMERA_H_
#define _CAMERA_H_

namespace edubpt {

struct Camera {
	// 画像サイズ
	int image_width_px, image_height_px;
	int samples_per_pixel; // 主にパストレ用
	
	// イメージセンサ位置
	Vec imagesensor_center;
	Vec imagesensor_dir;
	Vec imagesensor_up;
	Vec imagesensor_u, imagesensor_v;
	// イメージセンササイズ
	double imagesensor_width, imagesensor_height;
	double pixel_width, pixel_height;
	
	// レンズパラメータ
	double imagesensor_to_lens_distance;
	double focal_length;
	double lens_radius;

	// オブジェクトプレーン位置
	Vec objectplane_center;
	Vec objectplane_u, objectplane_v;

	// レンズ位置
	Vec lens_center;
	Vec lens_u, lens_v;
	Vec normal_on_lens;

	// イメージセンサの感度
	double W;

	Camera(const int image_width_px, const int image_height_px, const Vec &imagesensor_center, const Vec &imagesensor_dir, const Vec &imagesensor_up, const double imagesensor_size,
		const double imagesensor_to_lens_distance, const double focal_length, const double lens_radius, const double W_scale, const int samples_per_pixel) :
	image_width_px(image_width_px), image_height_px(image_height_px), imagesensor_center(imagesensor_center), imagesensor_dir(imagesensor_dir), imagesensor_up(imagesensor_up), 
		imagesensor_to_lens_distance(imagesensor_to_lens_distance), focal_length(focal_length), lens_radius(lens_radius), samples_per_pixel(samples_per_pixel)
	{
		imagesensor_width = imagesensor_size * image_width_px / image_height_px;
		imagesensor_height= imagesensor_size;
		
		pixel_width = imagesensor_width / image_width_px;
		pixel_height = imagesensor_height / image_height_px;
		
		imagesensor_u = normalize(cross(imagesensor_dir, imagesensor_up)) * imagesensor_width;
		imagesensor_v = normalize(cross(imagesensor_u, imagesensor_dir))  * imagesensor_height;
		
		objectplane_center = imagesensor_center + (focal_length + imagesensor_to_lens_distance) * imagesensor_dir;
		objectplane_u = imagesensor_u;
		objectplane_v = imagesensor_v;

		lens_center = imagesensor_center + imagesensor_to_lens_distance * imagesensor_dir;
		lens_u = lens_radius * normalize(imagesensor_u);
		lens_v = lens_radius * normalize(imagesensor_v);
		normal_on_lens = normalize(imagesensor_dir);

		// W(xI <- x0) センサのセンシティビティは簡単のため定数にしておく（ピクセルのサイズが小さくなったら感度あげる）
		W = W_scale / (pixel_width * pixel_height);
	}
	
	// 幾何的な係数計算 + センサーセンシティビティの項を計算する
	// x1 -> x0への放射輝度が最終的にイメージセンサに与える寄与度
	double W_dash(const Vec &x0_xV, const Vec &x0_xI, const Vec &x0_x1) const {
		return
			W *
			(x0_xV.length_squared() / x0_xI.length_squared()) *
			pow(
			(imagesensor_to_lens_distance * dot(normalize(x0_xI), normalize(-1.0 * imagesensor_dir))) /
			(focal_length * dot(normalize(x0_x1), normalize(imagesensor_dir))), 2.0);
	}

	// イメージセンサ上のサンプリング確率密度（イメージセンサの面積測度に関する確率密度）をシーン上のサンプリング確率密度（面積測度に関する確率密度）に変換する
	double P_Image_to_PA_x1(const double P_Image, const Vec &x0_xV, const Vec &x0_x1, const Vec &orienting_normal) const {
		return
			P_Image *
			pow(imagesensor_to_lens_distance / focal_length, 2.0) *
			(x0_xV.length_squared() / x0_x1.length_squared()) *
			(dot(normalize(-1.0 * x0_x1), orienting_normal) / dot(normalize(x0_x1), normalize(imagesensor_dir)));
	}
	
	// レンズとレイの交差判定関数
	double intersect_lens(const Ray &now_ray, Vec *position_on_lens, Vec *position_on_objectplane, Vec *position_on_imagesensor, Vec *uv_on_imagebuffer) const {
		// レンズ平面
		const Vec lens_normal = normalize(imagesensor_dir);
		// オブジェクトプレーン平面
		const Vec objectplane_normal = lens_normal;

		// レンズと判定
		const double lens_t = plane_intersection(lens_normal, lens_center, now_ray);
		if (kEPS < lens_t) {
			*position_on_lens = now_ray.org + lens_t * now_ray.dir; // x0
			if ((*position_on_lens - lens_center).length() < lens_radius && dot(lens_normal, now_ray.dir) <= 0.0)
			{
				// レンズとヒット
				// オブジェクトプレーンとの交差点を計算
				const double objectplane_t = plane_intersection(objectplane_normal, objectplane_center, now_ray);
				*position_on_objectplane = now_ray.org + objectplane_t * now_ray.dir; // xV
				const double u_on_objectplane = dot(*position_on_objectplane - objectplane_center, normalize(objectplane_u)) / objectplane_u.length();
				const double v_on_objectplane = dot(*position_on_objectplane - objectplane_center, normalize(objectplane_v)) / objectplane_v.length();

				const double ratio =  imagesensor_to_lens_distance / focal_length;
				const double u_on_imagesensor = -ratio * u_on_objectplane;
				const double v_on_imagesensor = -ratio * v_on_objectplane;
				*position_on_imagesensor = 
					imagesensor_center +
					u_on_imagesensor * imagesensor_u +
					v_on_imagesensor * imagesensor_v; // xI

				// レイが最終的にイメージセンサ内に入ったかどうかを判定する
				if (u_on_imagesensor >= -0.5 && 0.5 > u_on_imagesensor &&
					v_on_imagesensor >= -0.5 && 0.5 > v_on_imagesensor) {

					// ピクセル座標計算
					uv_on_imagebuffer->x = (u_on_imagesensor + 0.5) * image_width_px;
					uv_on_imagebuffer->y = (v_on_imagesensor + 0.5) * image_height_px;

					return lens_t;
				}
			}
		}

		return -kINF;
	}

	double sampling_pdf_on_lens() const {
		return 1.0 / (kPI * lens_radius * lens_radius);
	}
	// イメージセンサ上とレンズ上にサンプルを生成する
	void sample_points(
		const int imagebuffer_x, const int imagebuffer_y, Random *random, 
		Vec *position_on_imagesensor, Vec *position_on_objectplane, Vec *position_on_lens,
		double *P_Image, double *P_lens) const {
		// イメージセンサー上でサンプル生成
		// 生成するのはimagebuffer_x, yのピクセル内

		// ピクセル内の座標、[0, 1]の範囲
		const double u_on_pixel = random->next01();
		const double v_on_pixel = random->next01();
		// イメージセンサ上の座標、[-0,5, 0.5]の範囲（0,0)がイメージセンサ中央を示す。
		const double u_on_imagesensor = ((imagebuffer_x + u_on_pixel) / image_width_px - 0.5);
		const double v_on_imagesensor = ((imagebuffer_y + v_on_pixel) / image_height_px - 0.5);
		*position_on_imagesensor = 
			imagesensor_center +
			u_on_imagesensor * imagesensor_u +
			v_on_imagesensor * imagesensor_v; // xI
				
		// オブジェクトプレーン上の座標計算
		const double ratio = focal_length / imagesensor_to_lens_distance;
		const double u_on_objectplane = -ratio * u_on_imagesensor;
		const double v_on_objectplane = -ratio * v_on_imagesensor;
		*position_on_objectplane = 
			objectplane_center +
			u_on_objectplane * objectplane_u +
			v_on_objectplane * objectplane_v; // xV

		// レンズ上でサンプル生成 [-1, 1]
		const double r0 = sqrt(random->next01());
		const double r1 = random->next01() * 2.0 * kPI;
		const double u_on_lens = r0 * cos(r1);
		const double v_on_lens = r0 * sin(r1);
		*position_on_lens =
			lens_center +
			u_on_lens * lens_u +
			v_on_lens * lens_v; // x0
		
		*P_Image = 1.0 / (pixel_width * pixel_height); // ピクセル内の一点をサンプリングする確率密度関数（面積測度）
		*P_lens = sampling_pdf_on_lens(); // レンズ上の一点をサンプリングする確率密度関数（面積測度）
	}

private:
	double plane_intersection(const Vec &normal, const Vec &pos, const Ray &ray) const {
		const double pn = dot(pos, normal);
		const double on = dot(ray.org, normal);
		const double dn = dot(ray.dir, normal);

		if (fabs(dn) > kEPS) {
			const double t = (pn - on) / dn;
			return t;
		}
		return -kINF;
	}
};

};

#endif