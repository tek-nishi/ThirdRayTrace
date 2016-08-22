#pragma once

//
// レンダラー
//

#include <limits>

#include "Texture.hpp"
#include "Torus.hpp"
#include "Mandelbulb.hpp"
#include "Mandelbox.hpp"


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
  
float focus;

float focal_distance;
float lens_radius;

Texture image_diffuse;
Texture image_specular;

float glow_max;
glm::vec3 glow_color;

float Specular;
float SpecularExp;
float SpecularMax;

glm::vec3 SpotLightColor;
glm::vec3 SpotLightDir;

glm::vec3 CamLightColor;
float CamLightMin;

int iterate;
float min_dist;
float max_dist;

int shadowIteration = 32;
int shadowMinDist;
float shadowCoef    = 0.4f;
float shadowPower   = 16.0f;

};


float getDistance(const glm::vec3& p) {
  float d = Mandelbox::distance(p);
  return d;
}

glm::vec3 getNormal(const glm::vec3& p) {
  const float d = 0.0001f;

  return glm::normalize(glm::vec3(getDistance(p + glm::vec3(d,0.0,0.0)) -getDistance(p + glm::vec3(-d,0.0,0.0)),
                                  getDistance(p + glm::vec3(0.0,d,0.0)) -getDistance(p + glm::vec3(0.0,-d,0.0)),
                                  getDistance(p + glm::vec3(0.0,0.0,d)) -getDistance(p + glm::vec3(0.0,0.0,-d))));
}



float genShadow(const glm::vec3& ro, const glm::vec3& rd) {
  float h = 0.0f;
  float c = 0.0f;
  float r = 1.0f;

  for(int i = 0; i < Info::shadowIteration; ++i) {
    h = getDistance(ro + rd * c);
    if(h < Info::shadowMinDist) {
      return Info::shadowCoef;
    }
    r = glm::min(r, h * Info::shadowPower / c);
    c += h;
  }
    
  return glm::mix(Info::shadowCoef, 1.0f, r);
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

	return (Info::SpotLightColor * diffuse + Info::CamLightColor * ambient) * color_diffuse + specular * color_specular;
}


glm::vec3 trace(const glm::vec3& ray_dir, const glm::vec3& ray_origin) {
  // レイマーチで交差点を調べる
  float td = 0.0f;
  glm::vec3 pos_on_ray;
  bool hit = false;
  float steps = 0.0f;
  for(int i = 0; !hit && i < Info::iterate; ++i) {
    pos_on_ray = ray_origin + ray_dir * td;
    float d = getDistance(pos_on_ray);
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

    // glm::vec3 new_ray_dir = glm::reflect(ray_dir, normal);
    // glm::vec3 new_pos_on_ray = pos_on_ray + new_ray_dir * 0.001f;

    float shadow = genShadow(pos_on_ray + normal * Info::min_dist, normal);
    // glm::vec3 color = IBL(normal);
    glm::vec3 color = lighting(normal,
                               Info::image_diffuse.getPixel(normal), Info::image_specular.getPixel(normal),
                               pos_on_ray, ray_dir);
    
    return color * shadow;
  }
  else {
    // どこにも衝突しなかった
    // return Info::glow_color * step_factor;
    return Info::image_diffuse.getPixel(ray_dir) + Info::glow_color * step_factor;
  }
}


void setupParams(const picojson::value& settings) {
  // カメラ
  glm::vec3 cam_rot = glm::radians(getVec<glm::vec3>(settings.get("cam_rot")));
  glm::mat4 transform = glm::rotate(cam_rot.y, glm::vec3(0.0f, 1.0f, 0.0f)) *
    glm::rotate(cam_rot.x, glm::vec3(1.0f, 0.0f, 0.0f)) *
    glm::rotate(cam_rot.z, glm::vec3(0.0f, 0.0f, 1.0f)) *
    glm::translate(getVec<glm::vec3>(settings.get("cam_pos")));
    
  Info::cam_pos  = (transform * glm::vec4{ 0.0f, 0.0f, 0.0f, 1.0f }).xyz();
  Info::cam_dir  = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 0.0f }).xyz();
  Info::cam_up   = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 0.0f }).xyz();
  Info::cam_side = glm::cross(Info::cam_dir, Info::cam_up);

  // 焦点距離
  Info::focus = settings.get("focus").get<double>();

  // DOF
  Info::focal_distance = settings.get("focal_distance").get<double>();
  Info::lens_radius    = settings.get("lens_radius").get<double>();

  // 光彩
  {
    Info::glow_max   = settings.get("glow_max").get<double>();
    float glow       = settings.get("glow").get<double>();
    Info::glow_color = getVec<glm::vec3>(settings.get("glow_color")) * glow;
  }

  // ライティング
  {
    Info::Specular    = settings.get("Specular").get<double>();
    Info::SpecularExp = settings.get("SpecularExp").get<double>();
    Info::SpecularMax = settings.get("SpecularMax").get<double>();
  
    float SpotLight      = settings.get("SpotLight").get<double>();
    Info::SpotLightColor = getVec<glm::vec3>(settings.get("SpotLightColor")) * SpotLight;
    glm::vec2 dir        = getVec<glm::vec2>(settings.get("SpotLightDir"));
    glm::vec3 spot_dir   = glm::vec3(glm::sin(dir.x * M_PI) * glm::cos(dir.y * M_PI / 2.0),
                                     glm::sin(dir.y * M_PI / 2.0) * glm::sin(dir.x * M_PI),
                                     glm::cos(dir.x * M_PI));
    Info::SpotLightDir   = glm::normalize(spot_dir);

    float CamLight      = settings.get("CamLight").get<double>();
    Info::CamLightColor = getVec<glm::vec3>(settings.get("CamLightColor")) * CamLight;
    Info::CamLightMin   = settings.get("CamLightMin").get<double>();
  }

  // 反復数
  Info::iterate  = settings.get("iterate").get<double>();
  Info::min_dist = glm::pow(10.0, settings.get("detail").get<double>());
  Info::max_dist = settings.get("max_dist").get<double>();

  // 影
  {
    const auto& p = settings.get("shadow");

    Info::shadowIteration = p.get("iteration").get<double>();
    Info::shadowMinDist   = glm::pow(10.0, p.get("detail").get<double>());
    Info::shadowCoef      = p.get("coef").get<double>();
    Info::shadowPower     = p.get("power").get<double>();
  }

  Torus::init(settings.get("Torus"));
  Mandelbulb::init(settings.get("Mandelbulb"));
  Mandelbox::init(settings.get("Mandelbox"));
}

void render(const std::shared_ptr<RenderParams>& params) {
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
  
  // レンダリング回数
  // TIPS:-1で無限ループw
  int render_iterate = settings.get("render_iterate").get<double>();
  params->render_num = 0;
  for (int i = 0; i != render_iterate; ++i) {
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
        params->pixel[offset] = glm::mix(params->pixel[offset], trace(ray_dir, ray_origin), d);
      }
    }
    
    params->render_num += 1;
  }
  
  // レンダリング完了!!
  params->complete = true;
}
