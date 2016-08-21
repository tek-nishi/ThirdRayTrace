#pragma once

//
// 雑多な処理
//

#define GLM_SWIZZLE
#define GLM_META_PROG_HELPERS
#include <glm/glm.hpp>

// 画像書き出し
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>


enum {
  COMPONENTS = 3
};


// picojson::array→glm::vec系
template <typename T>
T getVec(const picojson::value& values) {
  T v;
  const auto& array = values.get<picojson::array>();
  for (size_t i = 0; i < T::components; ++i) {
    v[i] = array[i].get<double>();
  }

  return v;
}


// 途中結果を書き出す
void writeProgressImage(const std::string& path,
                        const std::vector<glm::vec3>& pixel, const int width, const int height) {
  std::vector<glm::u8vec3> bitmap(width * height);
  for (size_t i = 0; i < pixel.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(pixel[i] * 255.0f, 0.0f, 255.0f));
  }
  
  stbi_write_bmp(path.c_str(), width, height, COMPONENTS, bitmap.data());
}

void writeDepthImage(const std::string& path,
                     const std::vector<float>& depth, const int width, const int height) {
  std::vector<glm::u8vec3> bitmap(width * height);
  for (size_t i = 0; i < depth.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(depth[i] * 50.0f, 0.0f, 255.0f));
  }
  
  stbi_write_bmp(path.c_str(), width, height, COMPONENTS, bitmap.data());
}

// 最終結果を書き出す
void writeFinalImage(const std::string& path,
                     const std::vector<glm::vec3>& pixel, const int width, const int height) {
  std::vector<glm::u8vec3> bitmap(width * height);
  for (size_t i = 0; i < pixel.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(pixel[i] * 255.0f, 0.0f, 255.0f));
  }
  
  stbi_write_bmp(path.c_str(), width, height, COMPONENTS, bitmap.data());
}


// vec3からUV値へ
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
