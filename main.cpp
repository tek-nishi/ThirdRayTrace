//
// レイトレ合宿４向け作品
//
// 1. 何もインストールしていないまっさらなマシン上で動作する
// 2. 実行ファイルを叩いたら自動で始まるように(キーボード、マウスの操作を要求する作りにしない)
// 3. おおよそ30秒毎に、レンダリングの途中経過をbmpかpngで連番(000.png, 001.png, ...) で出力
// 4. ネットワーク越しに何かをやるような動作をさせない
// 5. 5分以内に自動で終了
//

// HDR読み込み
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
// PNG書き出し
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#define GLM_SWIZZLE
#define GLM_META_PROG_HELPERS
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/string_cast.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <picojson.h>
#include <limits>
#include <chrono>


enum {
  COMPONENTS = 3,
};


// 繰り返し
glm::vec3 trans(const glm::vec3& p) {
  return glm::mod(p, 8.0f) - 4.0f;
}


namespace Sphere {

// 距離関数(球)
float distance(const glm::vec3& p) {
  return glm::length(p) - 1.0;
}

// 法線を求める
glm::vec3 normal(const glm::vec3& p) {
  const float d = 0.0001;

  return glm::normalize(glm::vec3(distance(p + glm::vec3(d,0.0,0.0)) -distance(p + glm::vec3(-d,0.0,0.0)),
                                  distance(p + glm::vec3(0.0,d,0.0)) -distance(p + glm::vec3(0.0,-d,0.0)),
                                  distance(p + glm::vec3(0.0,0.0,d)) -distance(p + glm::vec3(0.0,0.0,-d))));
}

}


namespace Cube {

float distance(const glm::vec3& p) {
  glm::vec3 q = glm::abs(trans(p));
  return length(glm::max(q - glm::vec3(0.5f, 0.5f, 0.5f), 0.0f));
}

}


namespace RoundedCube {

float distance(const glm::vec3& p, const glm::vec3& c) {
  glm::vec3 q = glm::abs(p - c);
  return length(glm::max(q - glm::vec3(0.5f, 0.5f, 0.5f), 0.0f)) - 0.1f;
}

}


namespace Torus {

float distance(const glm::vec3& p) {
  glm::vec2 t = glm::vec2(0.75f, 0.25f);
  glm::vec2 r = glm::vec2(length(p.xz()) - t.x, p.y);

  return glm::length(r) - t.y;
}

}


namespace Plane {

float distance(const glm::vec3& p) {
  return glm::dot(p, glm::vec3(0.0f, 1.0f, 0.0f)) + 1.0f;
}

}


namespace Julia {

int Iterations;
float Threshold;
glm::vec4 C;

float distance(const glm::vec3& pos) {
	glm::vec4 p = glm::vec4(pos, 0.0);
	glm::vec4 dp = glm::vec4(1.0, 0.0,0.0,0.0);
	for (int i = 0; i < Iterations; i++) {
		dp = 2.0f * glm::vec4(p.x * dp.x - glm::dot(p.yzw(), dp.yzw()), p.x * dp.yzw() + dp.x * p.yzw() + glm::cross(p.yzw(), dp.yzw()));
		p = glm::vec4(p.x * p.x - glm::dot(p.yzw(), p.yzw()), glm::vec3(2.0f * p.x * p.yzw())) + C;
		float p2 = glm::dot(p, p);
		//orbitTrap = glm::min(orbitTrap, glm::abs(glm::vec4(p.xyz(), p2)));
		if (p2 > Threshold) break;
	}
	float r = glm::length(p);
	return  0.5f * r * glm::log(r) / glm::length(dp);
}

}


namespace Mandelbox {

int Iterations;
int ColorIterations;
float MinRad2;
float Scale;
glm::vec3 RotVector;
float RotAngle;

glm::vec4 scale;
glm::mat3 rot;
float absScalem1;
float AbsScaleRaisedTo1mIters;

void init() {
  scale = glm::vec4(Scale, Scale, Scale, glm::abs(Scale)) / MinRad2;
  rot   = glm::mat3_cast(glm::angleAxis(RotAngle, normalize(RotVector)));

  absScalem1 = glm::abs(Scale - 1.0f);
  AbsScaleRaisedTo1mIters = glm::pow(glm::abs(Scale), float(1 - Iterations));
}

// Compute the distance from `pos` to the Mandelbox.
float distance(const glm::vec3& pos) {
  glm::vec4 p = glm::vec4(pos,1), p0 = p;  // p.w is the distance estimate
	
	for (int i = 0; i < Iterations; i++) {
		//p.xyz*=rot;
		p.xyz() = glm::clamp(p.xyz(), -1.0f, 1.0f) * 2.0f - p.xyz();  // min;max;mad
		float r2 = glm::dot(p.xyz(), p.xyz());
		//if (i<ColorIterations) orbitTrap = min(orbitTrap, abs(vec4(p.xyz,r2)));
		p *= glm::clamp(glm::max(MinRad2 / r2, MinRad2), 0.0f, 1.0f);  // dp3,div,max.sat,mul
		p = p * scale + p0;
    if (r2 > 1000.0f) break;
	}
  
	return (glm::length(p.xyz()) - absScalem1) / p.w - AbsScaleRaisedTo1mIters;
}

}


float getDistance(const glm::vec3& p) {
  float d = Julia::distance(p);
  // d = glm::min(d, Plane::distance(p - glm::vec3(0, -1.5, 0)));
  return d;
}

glm::vec3 getNormal(const glm::vec3& p) {
  const float d = 0.0001f;

  return glm::normalize(glm::vec3(getDistance(p + glm::vec3(d,0.0,0.0)) -getDistance(p + glm::vec3(-d,0.0,0.0)),
                                  getDistance(p + glm::vec3(0.0,d,0.0)) -getDistance(p + glm::vec3(0.0,-d,0.0)),
                                  getDistance(p + glm::vec3(0.0,0.0,d)) -getDistance(p + glm::vec3(0.0,0.0,-d))));
}


int g_shadowIteration = 32;
float g_shadowCoef  = 0.4f;
float g_shadowPower = 16.0f;


float genShadow(const glm::vec3& ro, const glm::vec3& rd) {
  float h = 0.0f;
  float c = 0.0f;
  float r = 1.0f;

  for(int i = 0; i < g_shadowIteration; ++i) {
    h = getDistance(ro + rd * c);
    if(h < 0.0001f) {
      return g_shadowCoef;
    }
    r = glm::min(r, h * g_shadowPower / c);
    c += h;
  }
    
  return glm::mix(g_shadowCoef, 1.0f, r);
}


// 露出計算
// exposure 露出値(マイナス値)
glm::vec3 expose(const glm::vec3& light, const float exposure) {
  return glm::vec3(1.0f) - glm::exp(light * exposure);
}


// IBLからUVへ
// SOURCE:http://www_test.dstorm.co.jp/products/lw8/developer/docs/uv.htm
glm::vec2 vec3ToUV(const glm::vec3& normal) {
  // TIPS:glmはコンストラクタで０初期化してくれる
  glm::vec2 uv;
  
  if (normal.x == 0.0f && normal.z == 0.0f) {
    if (normal.y != 0.0f) uv.y = (normal.y < 0.0f) ? -M_PI_2 : M_PI_2;
    else                  uv.y = 0.0f;
  }
  else {
    if (normal.z == 0.0f) uv.x = (normal.x < 0.0f) ? M_PI_2 : -M_PI_2;
    else                  uv.x = glm::atan(normal.x, normal.z);
    
    float x = glm::length(normal.xz());
    if (x == 0.0f) uv.y = (normal.y < 0.0f) ? -M_PI_2 : M_PI_2;
    else           uv.y = glm::atan(normal.y, x);
  }

  uv.x = 0.5f - uv.x / (M_PI * 2.0f);
  uv.y = 0.5f - uv.y / M_PI;

  return uv;
}


template <typename T>
T getVec(const picojson::value& values) {
  T v;
  const auto& array = values.get<picojson::array>();
  for (size_t i = 0; i < T::components; ++i) {
    v[i] = array[i].get<double>();
  }

  return v;
}


float g_alpha = 0.1;

glm::vec3 g_diffuse;
glm::vec3 g_specular;

float sqr(float x) { return x * x; }

float GGX(float alpha, float cosThetaM) {
  float CosSquared = cosThetaM * cosThetaM;
  float TanSquared = (1.0f - CosSquared) / CosSquared;
  return (1.0f / M_PI) * sqr(alpha / (CosSquared * (alpha*alpha + TanSquared)));
}

// フレネル項
float Fresnel(const float c, const float f0) {
  float sf = glm::sqrt(f0);
  float n = (1.0 + sf) / (1.0 - sf);
  float g = glm::sqrt(n * n + c * c - 1.0);
  float ga = (c * (g + c) - 1.0) * (c * (g + c) - 1.0);
  float gb = (c * (g - c) + 1.0) * (c * (g - c) + 1.0);
  
  return (g - c) * (g - c) / (2.0 * (g + c) + (g + c)) * (1.0 + ga / gb);
}

// Cook-Torrance
// SOURCE:https://spphire9.wordpress.com/2013/03/13/webglでcook-torrance/
// L light vec
// V view vec
// N normal vec
glm::vec3 BRDF(const glm::vec3& L, const glm::vec3& V, const glm::vec3& N,
               const glm::vec3& specular, const glm::vec3& diffuse) {
  // ハーフベクトル
  glm::vec3 H = glm::normalize(L + V);

  float hn = glm::dot(H, N);
  float ln = glm::dot(L, N);
  float lh = glm::dot(L, H);
  float vn = glm::dot(V, N);  

  // Cook-Torrance式鏡面反射
  glm::vec3 f = glm::vec3(Fresnel(lh, specular.x), Fresnel(lh, specular.y), Fresnel(lh, specular.z));
  float d = GGX(g_alpha, hn);
  float t = 2.0f * hn / glm::dot(V, H);
  float g = glm::min(1.0f, glm::min(t * vn, t * ln));
  float m = M_PI * vn * ln;
  glm::vec3 spe = glm::max(f * d * g / m, 0.0f);

  // Lambertシェーディング
  glm::vec3 dif = glm::max(ln, 0.0f) * diffuse;

  return spe + dif;
}


struct Info {
  glm::ivec2 iresolution;
  glm::vec2 resolution;

  glm::vec3 cam_pos;
  glm::vec3 cam_dir;
  glm::vec3 cam_up;
  glm::vec3 cam_side;
  
  float focus;

  glm::vec3 light_dir;

  float* ibl_image;
  int bg_x, bg_y, bg_comp;
  glm::vec2 texture;

  int iterate;
  int nest;
  int path_trace;
};


glm::vec3 IBL(const Info& info, const glm::vec3 v) {
  glm::vec2 uv = glm::clamp(vec3ToUV(v) * info.texture, glm::vec2(0.0f, 0.0f), info.texture);
  int index = (int(uv.x) + int(uv.y) * info.bg_x) * COMPONENTS;
  return glm::vec3(info.ibl_image[index], info.ibl_image[index + 1], info.ibl_image[index + 2]);
}


glm::vec3 radiationVector_uniform(const glm::vec3& w) {
  glm::vec3 u = (glm::abs(w.x) > 0.0001f) ? glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), w))
    : glm::normalize(glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), w));
  glm::vec3 v = glm::cross(w, u);

  float phi = 2.0f * M_PI * glm::linearRand(0.0f, 1.0f);
  float cos_theta = glm::linearRand(0.0f, 1.0f);
  float sin_theta = glm::sqrt(1.0f - cos_theta * cos_theta);
  
  return glm::normalize(glm::vec3(u * glm::cos(phi) * sin_theta
                                  + v * glm::sin(phi) * sin_theta
                                  + w * cos_theta));
}


glm::vec3 ortho(const glm::vec3& v) {
  //  See : http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts
  return glm::abs(v.x) > glm::abs(v.z) ? glm::vec3(-v.y, v.x, 0.0f)  : glm::vec3(0.0f, -v.z, v.y);
}

glm::vec3 getSampleBiased(const glm::vec3& dir, const float power) {
	glm::vec3 o1 = glm::normalize(ortho(dir));
	glm::vec3 o2 = glm::normalize(glm::cross(dir, o1));
	glm::vec2 r = glm::linearRand(glm::vec2(0.0f, 0.0f), glm::vec2(1.0f, 1.0f));
	r.x = r.x * 2.0f * M_PI;
	r.y = glm::pow(r.y, 1.0f / (power + 1.0));
	float oneminus = glm::sqrt(1.0f - r.y * r.y);
	return glm::cos(r.x) * oneminus * o1 + glm::sin(r.x) * oneminus * o2 + r.y * dir;
}

glm::vec3 getSample(const glm::vec3& dir) {
  return getSampleBiased(dir, 1.0f);
}


glm::vec3 pathtrace(const Info& info, const int nest,
                    const glm::vec3& ray_dir, const glm::vec3& ray_origin) {
  // レイマーチで交差点を調べる
  float total_distance = 0.0f;
  glm::vec3 pos_on_ray;
  bool hit = false;
  for(int i = 0; !hit && i < info.iterate; ++i) {
    pos_on_ray = ray_origin + ray_dir * total_distance;
    float d = getDistance(pos_on_ray);
    total_distance += d;

    hit = glm::abs(d) < 0.001f;
  }

  if (hit) {
    glm::vec3 normal = getNormal(pos_on_ray);

    glm::vec3 new_ray_dir = glm::reflect(ray_dir, normal);
    glm::vec3 new_pos_on_ray = pos_on_ray + new_ray_dir * 0.001f;

    float shadow = genShadow(new_pos_on_ray + normal * 0.001f, normal);
    glm::vec3 color = IBL(info, normal);
    
    return color * shadow;
  }
  else {
    // どこにも衝突しなかった
    return IBL(info, ray_dir);
  }
}


void render(const Info& info, std::vector<glm::vec3>& pixel) {
  // 全ピクセルでレイマーチ
  for (int y = 0; y < info.iresolution.y; ++y) {
    glm::vec2 coord;
    coord.y = (info.iresolution.y - 1) - y;
    
    for (int x = 0; x < info.iresolution.x; ++x) {
      coord.x = x;

      glm::vec2 pos = (coord.xy() * 2.0f - info.resolution) / info.resolution.y;
      glm::vec3 ray_dir = glm::normalize(info.cam_side * pos.x + info.cam_up * pos.y + info.cam_dir * info.focus);

      glm::vec3 color = pathtrace(info, 0, ray_dir, info.cam_pos);
      pixel[x + y * info.iresolution.x] = color;
    }
  }
}


double lapTime(const std::chrono::system_clock::time_point begin,
               const std::chrono::system_clock::time_point end) {
  std::chrono::microseconds d = end - begin;

  // マイクロ秒→秒
  return d.count() / 1000000.0f;
}


int main() {
  std::vector<std::chrono::system_clock::time_point> lap_times;
  lap_times.push_back(std::chrono::system_clock::now());

  // セッティング
  picojson::value settings;
  {
    std::ifstream ifs("settings.json");
    ifs >> settings;
  }
  
  Info info;

  // 解像度
  // TIPS:キャストを減らすためint型とfloat型の両方を用意
  info.iresolution = getVec<glm::ivec2>(settings.get("resolution"));
  info.resolution  = glm::vec2(info.iresolution);

  // レンダリング結果
  std::vector<glm::vec3> pixel(info.iresolution.x * info.iresolution.y);

  // カメラ
  glm::vec3 cam_rot = glm::radians(getVec<glm::vec3>(settings.get("cam_rot")));
  glm::mat4 transform = glm::rotate(cam_rot.y, glm::vec3(0.0f, 1.0f, 0.0f)) *
    glm::rotate(cam_rot.x, glm::vec3(1.0f, 0.0f, 0.0f)) *
    glm::rotate(cam_rot.z, glm::vec3(0.0f, 0.0f, 1.0f)) *
    glm::translate(getVec<glm::vec3>(settings.get("cam_pos")));
    
  info.cam_pos  = (transform * glm::vec4{ 0.0f, 0.0f, 0.0f, 1.0f }).xyz();
  info.cam_dir  = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 0.0f }).xyz();
  info.cam_up   = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 0.0f }).xyz();
  info.cam_side = glm::cross(info.cam_dir, info.cam_up);

  // 焦点距離
  info.focus = settings.get("focus").get<double>();

  // 光源
  info.light_dir = getVec<glm::vec3>(settings.get("light_dir"));
  
  // IBL用画像
  info.ibl_image = stbi_loadf(settings.get("bg").get<std::string>().c_str(), &info.bg_x, &info.bg_y, &info.bg_comp, 0);
  info.texture = glm::vec2(info.bg_x - 1, info.bg_y - 1);

  // BRDF
  g_alpha    = settings.get("alpha").get<double>();
  g_diffuse  = getVec<glm::vec3>(settings.get("diffuse"));
  g_specular = getVec<glm::vec3>(settings.get("specular"));

  // 反復数
  info.iterate    = settings.get("iterate").get<double>();
  info.nest       = settings.get("nest").get<double>();
  info.path_trace = settings.get("path_trace").get<double>();

  // 影
  {
    const auto shadow = settings.get("shadow");

    g_shadowIteration = shadow.get("iteration").get<double>();
    g_shadowCoef      = shadow.get("coef").get<double>();
    g_shadowPower     = shadow.get("power").get<double>();
  }
  
  // ジュリア集合
  {
    const auto julia = settings.get("julia");

    Julia::Iterations = julia.get("iterations").get<double>();
    Julia::Threshold = julia.get("threshold").get<double>();
    Julia::C = getVec<glm::vec4>(julia.get("C"));
  }

  // Mandelbox集合
  {
    const auto mandelbox = settings.get("Mandelbox");

    Mandelbox::Iterations      = mandelbox.get("Iterations").get<double>();
    Mandelbox::ColorIterations = mandelbox.get("ColorIterations").get<double>();
    Mandelbox::MinRad2         = mandelbox.get("MinRad2").get<double>();
    Mandelbox::Scale           = mandelbox.get("Scale").get<double>();
    Mandelbox::RotVector       = getVec<glm::vec3>(mandelbox.get("RotVector"));
    Mandelbox::RotAngle        = glm::radians(mandelbox.get("RotAngle").get<double>());

    Mandelbox::init();
  }

  lap_times.push_back(std::chrono::system_clock::now());

  render(info, pixel);
  stbi_image_free(info.ibl_image);

  lap_times.push_back(std::chrono::system_clock::now());

  // レンダリング結果をPNG形式に
  float exposure = settings.get("exposure").get<double>();
  std::vector<glm::u8vec3> bitmap(info.iresolution.x * info.iresolution.y);
  for (size_t i = 0; i < pixel.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(expose(pixel[i], exposure) * 255.0f, 0.0f, 255.0f));
  }

  lap_times.push_back(std::chrono::system_clock::now());

  stbi_write_png("test.png", info.iresolution.x, info.iresolution.y, COMPONENTS, bitmap.data(), info.iresolution.x * COMPONENTS);
  
  lap_times.push_back(std::chrono::system_clock::now());

  std::cout << "total:" << lapTime(lap_times.front(), lap_times.back()) << std::endl;
  for (size_t i = 0; i < lap_times.size() - 1; ++i) {
    std::cout << "lap " << i << ":" << lapTime(lap_times[i], lap_times[i + 1]) << std::endl;
  }
}
