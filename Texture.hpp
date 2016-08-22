#pragma once

//
// HDR画像
//

// HDR読み込み
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include "Misc.hpp"


class Texture {
  int width, height;
  glm::vec2 size;
  std::vector<glm::vec3> pixel;
  

public:
  Texture() = default;
  
  Texture(const std::string& path) {
    int comp;
    float* image = stbi_loadf(path.c_str(), &width, &height, &comp, 0);
    size = glm::vec2(width - 1, height - 1);

    pixel.resize(width * height);

    for (size_t i = 0; i < pixel.size(); ++i) {
      pixel[i] = glm::vec3(image[i * comp], image[i * comp + 1], image[i * comp + 2]);
    }
    
    stbi_image_free(image);
  }


  const glm::vec3& getPixel(const glm::vec3& v) const {
    glm::vec2 uv = glm::clamp(vec3ToUV(v) * size, glm::vec2(0.0f), size);
    int index = (int(uv.x) + int(uv.y) * width);
    return pixel[index];
  }
};

