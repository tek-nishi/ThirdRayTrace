#pragma once

//
// レンダラー
//

// HDR読み込み
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <limits>

#include "Mandelbulb.hpp"


// レンダリングに必要な情報
struct RenderParams {
  picojson::value settings;

  glm::ivec2 iresolution;
  glm::vec2  resolution;

  std::vector<glm::vec3> pixel;
  std::vector<glm::vec3> normal;
  std::vector<float> depth;
  
  bool complete;
};


struct Info {
  glm::ivec2 iresolution;
  glm::vec2  resolution;

  glm::vec3 cam_pos;
  glm::vec3 cam_dir;
  glm::vec3 cam_up;
  glm::vec3 cam_side;
  
  float focus;

  float* ibl_image;
  int bg_x, bg_y, bg_comp;
  glm::vec2 texture;

  int iterate;
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


int g_shadowIteration = 32;
float g_shadowCoef    = 0.4f;
float g_shadowPower   = 16.0f;

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


// 正規化ベクトルから画素を取り出す
glm::vec3 IBL(const Info& info, const glm::vec3& v) {
  glm::vec2 uv = glm::clamp(vec3ToUV(v) * info.texture, glm::vec2(0.0f, 0.0f), info.texture);
  int index = (int(uv.x) + int(uv.y) * info.bg_x) * COMPONENTS;
  return glm::vec3(info.ibl_image[index], info.ibl_image[index + 1], info.ibl_image[index + 2]);
}


std::tuple<glm::vec3, glm::vec3, float> trace(const Info& info, const glm::vec3& ray_dir, const glm::vec3& ray_origin) {
  // レイマーチで交差点を調べる
  float td = 0.0f;
  glm::vec3 pos_on_ray;
  bool hit = false;
  for(int i = 0; !hit && i < info.iterate; ++i) {
    pos_on_ray = ray_origin + ray_dir * td;
    float d = getDistance(pos_on_ray);
    td += d;

    hit = glm::abs(d) < 0.001f;
  }

  if (hit) {
    glm::vec3 normal = getNormal(pos_on_ray);

    glm::vec3 new_ray_dir = glm::reflect(ray_dir, normal);
    // glm::vec3 new_pos_on_ray = pos_on_ray + new_ray_dir * 0.001f;

    float shadow = genShadow(pos_on_ray + normal * 0.001f, normal);
    glm::vec3 color = IBL(info, normal);
    // glm::vec3 color = glm::vec3(1);
    
    return std::make_tuple(color * shadow, normal, td);
  }
  else {
    // どこにも衝突しなかった
    return std::make_tuple(IBL(info, ray_dir), glm::vec3(0.0f), std::numeric_limits<float>::max());
  }
}


void render(const std::shared_ptr<RenderParams>& params) {
  Info info;

  // TIPS:スレッド間で共有している値へのアクセスは重たいので
  //      必要な値はコピーしておく
  info.iresolution = params->iresolution;
  info.resolution  = params->resolution;
  
  // カメラ
  glm::vec3 cam_rot = glm::radians(getVec<glm::vec3>(params->settings.get("cam_rot")));
  glm::mat4 transform = glm::rotate(cam_rot.y, glm::vec3(0.0f, 1.0f, 0.0f)) *
    glm::rotate(cam_rot.x, glm::vec3(1.0f, 0.0f, 0.0f)) *
    glm::rotate(cam_rot.z, glm::vec3(0.0f, 0.0f, 1.0f)) *
    glm::translate(getVec<glm::vec3>(params->settings.get("cam_pos")));
    
  info.cam_pos  = (transform * glm::vec4{ 0.0f, 0.0f, 0.0f, 1.0f }).xyz();
  info.cam_dir  = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 0.0f }).xyz();
  info.cam_up   = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 0.0f }).xyz();
  info.cam_side = glm::cross(info.cam_dir, info.cam_up);

  // 焦点距離
  info.focus = params->settings.get("focus").get<double>();

  // 反復数
  info.iterate = params->settings.get("iterate").get<double>();

  // 影
  {
    const auto& p = params->settings.get("shadow");

    g_shadowIteration = p.get("iteration").get<double>();
    g_shadowCoef      = p.get("coef").get<double>();
    g_shadowPower     = p.get("power").get<double>();
  }

  
  // Mandelbulb集合
  {
    const auto& p = params->settings.get("Mandelbulb");
    Mandelbulb::init(p);
  }

  
  // IBL用画像
  info.ibl_image = stbi_loadf(params->settings.get("bg").get<std::string>().c_str(), &info.bg_x, &info.bg_y, &info.bg_comp, 0);
  info.texture = glm::vec2(info.bg_x - 1, info.bg_y - 1);

  
  // 全ピクセルでレイマーチ
  for (int y = 0; y < info.iresolution.y; ++y) {
    glm::vec2 coord;
    coord.y = (info.iresolution.y - 1) - y;
    
    for (int x = 0; x < info.iresolution.x; ++x) {
      coord.x = x;

      glm::vec2 pos = (coord.xy() * 2.0f - info.resolution) / info.resolution.y;
      glm::vec3 ray_dir = glm::normalize(info.cam_side * pos.x + info.cam_up * pos.y + info.cam_dir * info.focus);
      
      auto results = trace(info, ray_dir, info.cam_pos);
      
      size_t offset = x + y * info.iresolution.x;
      params->pixel[offset]  = std::get<0>(results);
      params->normal[offset] = std::get<1>(results);
      params->depth[offset]  = std::get<2>(results);
    }
  }

  
  // レンダリング完了!!
  params->complete = true;
  
  stbi_image_free(info.ibl_image);
}
