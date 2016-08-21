#pragma once

//
// レンダラー
//

// HDR読み込み
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <limits>

#include "Torus.hpp"
#include "Mandelbulb.hpp"


// レンダリングに必要な情報
struct RenderParams {
  picojson::value settings;

  glm::ivec2 iresolution;
  glm::vec2  resolution;

  std::vector<glm::vec3> pixel;
  
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

float* ibl_image;
int bg_x, bg_y, bg_comp;
glm::vec2 texture;

int iterate;
float min_dist;

int shadowIteration = 32;
int shadowMinDist;
float shadowCoef    = 0.4f;
float shadowPower   = 16.0f;

};


float getDistance(const glm::vec3& p) {
  float d = Mandelbulb::distance(p);
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


// 正規化ベクトルから画素を取り出す
glm::vec3 IBL(const glm::vec3& v) {
  glm::vec2 uv = glm::clamp(vec3ToUV(v) * Info::texture, glm::vec2(0.0f, 0.0f), Info::texture);
  int index = (int(uv.x) + int(uv.y) * Info::bg_x) * COMPONENTS;
  return glm::vec3(Info::ibl_image[index], Info::ibl_image[index + 1], Info::ibl_image[index + 2]);
}


glm::vec3 trace(const glm::vec3& ray_dir, const glm::vec3& ray_origin) {
  // レイマーチで交差点を調べる
  float td = 0.0f;
  glm::vec3 pos_on_ray;
  bool hit = false;
  for(int i = 0; !hit && i < Info::iterate; ++i) {
    pos_on_ray = ray_origin + ray_dir * td;
    float d = getDistance(pos_on_ray);
    td += d;

    hit = glm::abs(d) < Info::min_dist;
  }

  if (hit) {
    glm::vec3 normal = getNormal(pos_on_ray);

    // glm::vec3 new_ray_dir = glm::reflect(ray_dir, normal);
    // glm::vec3 new_pos_on_ray = pos_on_ray + new_ray_dir * 0.001f;

    float shadow = genShadow(pos_on_ray + normal * Info::min_dist, normal);
    glm::vec3 color = IBL(normal);
    // glm::vec3 color = glm::vec3(1);
    
    return color * shadow;
  }
  else {
    // どこにも衝突しなかった
    return IBL(ray_dir);
  }
}


void render(const std::shared_ptr<RenderParams>& params) {
  // TIPS:スレッド間で共有している値へのアクセスは重たいので
  //      必要な値はコピーしておく
  Info::iresolution = params->iresolution;
  Info::resolution  = params->resolution;
  
  // カメラ
  glm::vec3 cam_rot = glm::radians(getVec<glm::vec3>(params->settings.get("cam_rot")));
  glm::mat4 transform = glm::rotate(cam_rot.y, glm::vec3(0.0f, 1.0f, 0.0f)) *
    glm::rotate(cam_rot.x, glm::vec3(1.0f, 0.0f, 0.0f)) *
    glm::rotate(cam_rot.z, glm::vec3(0.0f, 0.0f, 1.0f)) *
    glm::translate(getVec<glm::vec3>(params->settings.get("cam_pos")));
    
  Info::cam_pos  = (transform * glm::vec4{ 0.0f, 0.0f, 0.0f, 1.0f }).xyz();
  Info::cam_dir  = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 0.0f }).xyz();
  Info::cam_up   = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 0.0f }).xyz();
  Info::cam_side = glm::cross(Info::cam_dir, Info::cam_up);

  // 焦点距離
  Info::focus = params->settings.get("focus").get<double>();

  // DOF
  Info::focal_distance = params->settings.get("focal_distance").get<double>();
  Info::lens_radius    = params->settings.get("lens_radius").get<double>();

  // 反復数
  Info::iterate  = params->settings.get("iterate").get<double>();
  Info::min_dist = glm::pow(10.0, params->settings.get("detail").get<double>());

  // 影
  {
    const auto& p = params->settings.get("shadow");

    Info::shadowIteration = p.get("iteration").get<double>();
    Info::shadowMinDist   = glm::pow(10.0, p.get("detail").get<double>());
    Info::shadowCoef      = p.get("coef").get<double>();
    Info::shadowPower     = p.get("power").get<double>();
  }


  // Torus
  {
    const auto& p = params->settings.get("Torus");
    Torus::init(p);
  }
  
  // Mandelbulb集合
  {
    const auto& p = params->settings.get("Mandelbulb");
    Mandelbulb::init(p);
  }

  
  // IBL用画像
  Info::ibl_image = stbi_loadf(params->settings.get("bg").get<std::string>().c_str(), &Info::bg_x, &Info::bg_y, &Info::bg_comp, 0);
  Info::texture = glm::vec2(Info::bg_x - 1, Info::bg_y - 1);


  // 無限ループw
  for (int i = 0; ; ++i) {
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
        // レンダリング回数から、新しい色の影響力を決める
        float d = 1.0f / float(i + 1);
        params->pixel[offset] = glm::mix(params->pixel[offset], trace(ray_dir, ray_origin), d);
      }
    }
  }

  
  // レンダリング完了!!
  params->complete = true;
  
  stbi_image_free(Info::ibl_image);
}
