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


// 露出計算
// exposure 露出値(マイナス値)
glm::vec3 expose(const glm::vec3& light, const float exposure) {
  return glm::vec3(1.0f) - glm::exp(light * exposure);
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

// 最終結果を書き出す
void writeFinalImage(const std::string& path,
                     const std::vector<glm::vec3>& pixel, const int width, const int height,
                     const float exposure) {
  std::vector<glm::u8vec3> bitmap(width * height);
  for (size_t i = 0; i < pixel.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(expose(pixel[i], exposure) * 255.0f, 0.0f, 255.0f));
  }
  
  stbi_write_bmp(path.c_str(), width, height, COMPONENTS, bitmap.data());
}
