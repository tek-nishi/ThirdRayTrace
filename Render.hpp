#pragma once

//
// レンダラー
//

#include <chrono>
#include <iostream>
#include <limits>

#include "Color.hpp"
#include "Texture.hpp"
#include "Plane.hpp"
#include "Sphere.hpp"
#include "Torus.hpp"
#include "Mandelbulb.hpp"
#include "Mandelbox.hpp"
#include "QuaternionJulia.hpp"


// レンダリングに必要な情報
struct RenderParams {
  picojson::value settings;

  glm::ivec2 iresolution;
  glm::vec2  resolution;

  std::vector<glm::vec3> pixel;

  int render_num;
  bool complete;
};


// FIXME:全てグローバル変数...
namespace Info {

glm::ivec2 iresolution;
glm::vec2  resolution;

glm::vec3 cam_pos;
glm::vec3 cam_dir;
glm::vec3 cam_up;
glm::vec3 cam_side;

glm::vec3 Eye;
glm::vec3 Target;
glm::vec3 Up;

float focus;

float focal_distance;
float lens_radius;

Texture image_diffuse;
Texture image_specular;

float Fog;
glm::vec3 FogColor;

float glow_max;
glm::vec3 glow_color;

float Specular;
float SpecularExp;
float SpecularMax;
glm::vec3 SpecularColor;

glm::vec3 SpotLightColor;
glm::vec3 SpotLightDir;

glm::vec3 CamLightColor;
float CamLightMin;

int iterate;
float min_dist;
float max_dist;
float FudgeFactor;

int AoIteration;
float AoIterValue;
float AoStep;
float AoAttenuation;
float AoPower;
float AoMin;

int shadowIteration;
int shadowMinDist;
float shadowCoef;
float shadowPower;

std::vector<Color> colors;


glm::vec3 SphereOffset;
float SphereScale;

};


// 距離関数の結果に色情報を追加している
struct Result {
  float d;
  glm::vec4 orbitTrap;
  int index;

  Result() = default;
  
  Result(const std::pair<float, glm::vec4>& r, const int i)
    : d(r.first),
      orbitTrap(r.second),
      index(i)
  {}
  
};


Result getDistance(const glm::vec3& p) {

  auto d0 = Sphere::distance(p - Info::SphereOffset);
  Result result{ d0, 0 };

  {
    auto d1 = QuaternionJulia::distance((p - Info::SphereOffset) / Info::SphereScale);
    d1.first *= Info::SphereScale;
    if (d1.first < result.d) result = Result(d1, 1);
  }
  
  {
    auto d1 = Mandelbox::distance(p);
    if (d1.first < result.d) result = Result(d1, 2);
  }
  
  return result;
}

glm::vec3 getNormal(const glm::vec3& p) {
  const float d = 0.0001f;

  return glm::normalize(glm::vec3(getDistance(p + glm::vec3(  d, 0.0, 0.0)).d - getDistance(p + glm::vec3( -d, 0.0, 0.0)).d,
                                  getDistance(p + glm::vec3(0.0,   d, 0.0)).d - getDistance(p + glm::vec3(0.0,  -d, 0.0)).d,
                                  getDistance(p + glm::vec3(0.0, 0.0,   d)).d - getDistance(p + glm::vec3(0.0, 0.0,  -d)).d));
}



float genShadow(const glm::vec3& ro, const glm::vec3& rd) {
  float h = 0.0f;
  float c = 0.0f;
  float r = 1.0f;

  for(int i = 0; i < Info::shadowIteration; ++i) {
    h = getDistance(ro + rd * c).d * Info::FudgeFactor;
    if(h < Info::shadowMinDist) {
      return Info::shadowCoef;
    }
    r = glm::min(r, h * Info::shadowPower / c);
    c += h;
  }
    
  return glm::mix(Info::shadowCoef, 1.0f, r);
}

float calcAO(const glm::vec3& pos, const glm::vec3& nor) {
	float occ = 0.0f;
  float sca = 1.0f;
  for (int i = 0; i < Info::AoIteration; ++i) {
    float hr = Info::AoStep * i / Info::AoIterValue;
    glm::vec3 aopos =  nor * hr + pos;
    float dd = getDistance(aopos).d;
    occ += -(dd - hr) * sca;
    sca *= Info::AoAttenuation;
  }
  return glm::clamp(1.0f - Info::AoPower * occ, Info::AoMin, 1.0f);    
}


glm::vec3 lighting(const glm::vec3& n,
                   const glm::vec3& color_diffuse,
                   const glm::vec3& color_specular,
                   const glm::vec3& pos, const glm::vec3& dir) {
	float nDotL = glm::max( 0.0f, glm::dot(n, Info::SpotLightDir));
  glm::vec3 halfVector = glm::normalize(-dir + Info::SpotLightDir);
	float diffuse = nDotL;
	float ambient = glm::max(Info::CamLightMin, glm::dot(-n, dir));
	float hDotN = glm::max(0.0f, glm::dot(n, halfVector));

	// An attempt at Physcical Based Specular Shading:
	// http://renderwonk.com/publications/s2010-shading-course/
	// (Blinn-Phong with Schickl term and physical normalization)
	float specular = ((Info::SpecularExp + 2.0f) / 8.0f)
    * glm::pow(hDotN, Info::SpecularExp)
    * (Info::SpecularExp + (1.0f - Info::SpecularExp) * glm::pow(1.0 - hDotN, 5.0f))
    * nDotL
    * Info::Specular;
	specular = glm::min(Info::SpecularMax, specular);

	return (Info::SpotLightColor * diffuse + Info::CamLightColor * ambient) * color_diffuse + Info::SpecularColor * specular * color_specular;
}


glm::vec3 trace(const glm::vec3& ray_dir, const glm::vec3& ray_origin, const bool reflect) {
  // レイマーチで交差点を調べる
  float td = 0.0f;
  glm::vec3 pos_on_ray;
  bool hit = false;
  float steps = 0.0f;
  Result result;

  for(int i = 0; !hit && i < Info::iterate; ++i) {
    pos_on_ray = ray_origin + ray_dir * td;
    result = getDistance(pos_on_ray);
    float d = result.d * Info::FudgeFactor;
    td += d;
    steps += 1.0f;

    hit = glm::abs(d) < Info::min_dist;
    if (td > Info::max_dist) {
      steps -= (td - Info::max_dist) / d;
      break;
    }
  }

  float step_factor = glm::clamp(steps / Info::glow_max, 0.0f, 1.0f);
  
  if (hit) {
    glm::vec3 normal = getNormal(pos_on_ray);
    
    auto p = pos_on_ray + normal * Info::min_dist;
    float shadow = calcAO(p, normal)
                 * genShadow(p, Info::SpotLightDir);
    
    glm::vec3 reflect_ray_dir = glm::reflect(ray_dir, normal);

    const auto& mat = Info::colors[result.index];
    glm::vec3 c = mat.get(result.orbitTrap);
    
    glm::vec3 color = lighting(normal,
                               Info::image_diffuse.getPixel(normal) * c,
                               Info::image_diffuse.getPixel(reflect_ray_dir) * c,
                               pos_on_ray, ray_dir);
    color *= shadow;

    // OpenGL GL_EXP2 like fog
    color = glm::mix(color, Info::FogColor, 1.0f - glm::exp(-glm::pow(Info::Fog, 4.0f) * td * td));

    if (reflect && mat.isReflection()) {
      glm::vec3 new_pos_on_ray = pos_on_ray + reflect_ray_dir * Info::min_dist;
      color = glm::mix(color, trace(reflect_ray_dir, new_pos_on_ray, false), mat.getReflectionPower());
    }
    
    return color;
  }
  else {
    // どこにも衝突しなかった
    auto color = reflect ? Info::image_diffuse.getPixel(ray_dir) : Info::image_specular.getPixel(ray_dir);
    return color + Info::glow_color * step_factor;
  }
}


void setupParams(const picojson::value& settings) {
  // カメラ
  {
    Info::cam_pos = getVec<glm::vec3>(settings.get("Eye"));
    auto target   = getVec<glm::vec3>(settings.get("Target"));
    Info::cam_up  = glm::normalize(getVec<glm::vec3>(settings.get("Up")));

    Info::cam_dir  = glm::normalize(target - Info::cam_pos);
    Info::cam_side = glm::cross(Info::cam_dir, Info::cam_up);
  }

  // 焦点距離
  Info::focus = settings.get("focus").get<double>();

  // DOF
  Info::focal_distance = settings.get("focal_distance").get<double>();
  Info::lens_radius    = settings.get("lens_radius").get<double>();

  // ライティング
  {
    Info::Specular      = settings.get("Specular").get<double>();
    Info::SpecularExp   = settings.get("SpecularExp").get<double>();
    Info::SpecularMax   = settings.get("SpecularMax").get<double>();
    Info::SpecularColor = getVec<glm::vec3>(settings.get("SpotLightColor"));
  
    float SpotLight      = settings.get("SpotLight").get<double>();
    Info::SpotLightColor = Info::SpecularColor * SpotLight;
    glm::vec2 dir        = getVec<glm::vec2>(settings.get("SpotLightDir"));
    glm::vec3 spot_dir   = glm::vec3(glm::sin(dir.x * M_PI) * glm::cos(dir.y * M_PI / 2.0),
                                     glm::sin(dir.y * M_PI / 2.0) * glm::sin(dir.x * M_PI),
                                     glm::cos(dir.x * M_PI));
    Info::SpotLightDir   = glm::normalize(spot_dir);

    float CamLight      = settings.get("CamLight").get<double>();
    Info::CamLightColor = getVec<glm::vec3>(settings.get("CamLightColor")) * CamLight;
    Info::CamLightMin   = settings.get("CamLightMin").get<double>();
  }

  Info::Fog      = settings.get("Fog").get<double>();
  Info::FogColor = getVec<glm::vec3>(settings.get("FogColor"));
  
  // 光彩
  {
    Info::glow_max   = settings.get("glow_max").get<double>();
    float glow       = settings.get("glow").get<double>();
    Info::glow_color = getVec<glm::vec3>(settings.get("glow_color")) * glow;
  }

  // 反復数とか
  Info::iterate  = settings.get("iterate").get<double>();
  Info::min_dist = glm::pow(10.0, settings.get("detail").get<double>());
  Info::max_dist = settings.get("max_dist").get<double>();
  Info::FudgeFactor = settings.get("FudgeFactor").get<double>();

  // AO
  {
    const auto& p = settings.get("AO");

    Info::AoIteration   = p.get("iteration").get<double>();
    Info::AoIterValue   = Info::AoIteration - 1;
    Info::AoStep        = p.get("step").get<double>();
    Info::AoAttenuation = p.get("attenuation").get<double>();
    Info::AoPower       = p.get("power").get<double>();
    Info::AoMin         = p.get("min").get<double>();
  }
  
  // 影
  {
    const auto& p = settings.get("shadow");

    Info::shadowIteration = p.get("iteration").get<double>();
    Info::shadowMinDist   = glm::pow(10.0, p.get("detail").get<double>());
    Info::shadowCoef      = p.get("coef").get<double>();
    Info::shadowPower     = p.get("power").get<double>();
  }

  // マテリアルカラー
  {
    const auto& params = settings.get("Colors").get<picojson::array>();
    for (const auto& p : params) {
      Info::colors.emplace_back(p);
    }
  }

  Plane::init(settings.get("Plane"));
  Sphere::init(settings.get("Sphere"));
  Torus::init(settings.get("Torus"));
  Mandelbulb::init(settings.get("Mandelbulb"));
  Mandelbox::init(settings.get("Mandelbox"));
  QuaternionJulia::init(settings.get("QuaternionJulia"));

  // 配置
  {
    // float dir = settings.get("SphereDir").get<double>();
    // Info::SphereOffset = getVec<glm::vec3>(settings.get("SphereOffset")) + Info::cam_pos + Info::cam_dir * dir;
    Info::SphereOffset = getVec<glm::vec3>(settings.get("SphereOffset"));
    Info::SphereScale  = settings.get("SphereScale").get<double>();

    // std::cout << "Center distance: " << glm::length(Info::SphereOffset - Info::cam_pos) << std::endl;
  }
}

void render(std::shared_ptr<RenderParams> params) {
  // TIPS:スレッド間で共有している値へのアクセスは重たいので
  //      必要な値はコピーしておく
  Info::iresolution = params->iresolution;
  Info::resolution  = params->resolution;

  // 各種セットアップ
  const auto& settings = params->settings;
  setupParams(settings);
  
  // IBL用画像
  Info::image_diffuse  = Texture(settings.get("image_diffuse").get<std::string>());
  Info::image_specular = Texture(settings.get("image_specular").get<std::string>());

  auto begin_time = std::chrono::system_clock::now();
  
  // レンダリング回数
  int render_iterate = settings.get("render_iterate").get<double>();
  params->render_num = 0;

#pragma omp parallel for
  for (int i = 0; i < render_iterate; ++i) {
    // レンダリング回数から、新しい色の影響力を決める
    // TIPS:ループの内側で変化しない値なので、ここで定義
    float d = 1.0f / float(i + 1);
    
    // 全ピクセルでレイマーチ
    for (int y = 0; y < Info::iresolution.y; ++y) {
      glm::vec2 coord;
      coord.y = (Info::iresolution.y - 1) - y;
    
      for (int x = 0; x < Info::iresolution.x; ++x) {
        coord.x = x;

        glm::vec2 pos = (coord.xy() * 2.0f - Info::resolution) / Info::resolution.y;
        glm::vec3 ray_dir = glm::normalize(Info::cam_side * pos.x + Info::cam_up * pos.y + Info::cam_dir * Info::focus);

        // レンズの屈折をシミュレーション(被写界深度)
        // SOURCE:https://github.com/githole/simple-pathtracer/tree/simple-pathtracer-DOF

        // フォーカスが合う位置
        float ft = glm::abs(Info::focal_distance / glm::dot(Info::cam_dir, ray_dir));
        glm::vec3 focus_pos = Info::cam_pos + ray_dir * ft;

        // 適当に決めたレンズの通過位置とフォーカスが合う位置からRayを作り直す(屈折効果)
        glm::vec2 lens = glm::gaussRand(glm::vec2(-1), glm::vec2(1)) * Info::lens_radius;
        glm::vec3 ray_origin = glm::vec3(Info::cam_side * lens.x + Info::cam_up * lens.y) + Info::cam_pos;
        ray_dir = glm::normalize(focus_pos - ray_origin);

        int offset = x + y * Info::iresolution.x;
        params->pixel[offset] = glm::mix(params->pixel[offset], trace(ray_dir, ray_origin, true), d);
      }
    }
    
    params->render_num += 1;
  }
  
  // レンダリング完了!!
  params->complete = true;

  auto current_time = std::chrono::system_clock::now();
  std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - begin_time);
  std::cout << "Render time "
            << duration.count()
            << " sec."
            << std::endl;
}
