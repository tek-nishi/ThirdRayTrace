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
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include <vector>
#include <iostream>


enum {
  WIDTH      = 800,
  HEIGHT     = 600,
  COMPONENTS = 3,
};


// 繰り返し
glm::vec3 trans(const glm::vec3& p) {
  return glm::mod(p, 8.0f) - 4.0f;
}


namespace Sphere {

// 距離関数(球)
float distance(const glm::vec3& p) {
  return glm::length(trans(p)) - 1.0;
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


float getDistance(const glm::vec3& p) {
  return std::min(Torus::distance(p),
                  std::min(RoundedCube::distance(p, glm::vec3(-0.5f, 1.2f, 0.0f)),
                                                 Plane::distance(p)));
}

glm::vec3 getNormal(const glm::vec3& p) {
  const float d = 0.0001f;

  return glm::normalize(glm::vec3(getDistance(p + glm::vec3(d,0.0,0.0)) -getDistance(p + glm::vec3(-d,0.0,0.0)),
                                  getDistance(p + glm::vec3(0.0,d,0.0)) -getDistance(p + glm::vec3(0.0,-d,0.0)),
                                  getDistance(p + glm::vec3(0.0,0.0,d)) -getDistance(p + glm::vec3(0.0,0.0,-d))));
}


float genShadow(const glm::vec3& ro, const glm::vec3& rd) {
    float h = 0.0f;
    float c = 0.001f;
    float r = 1.0f;
    float shadowCoef = 0.5f;

    for(float t = 0.0f; t < 50.0f; t += 1.0f) {
        h = getDistance(ro + rd * c);
        if(h < 0.001f) {
            return shadowCoef;
        }
        r = glm::min(r, h * 4.0f / c);
        c += h;
    }
    
    return 1.0f - shadowCoef + r * shadowCoef;
}


// 露出計算
// exposure 露出値(マイナス値)
glm::vec3 expose(const glm::vec3& light, const float exposure) {
  return glm::vec3(1.0f) - glm::exp(light * exposure);
}


// IBLからUVへ
// SOURCE:http://www_test.dstorm.co.jp/products/lw8/developer/docs/uv.htm
glm::vec2 vec3ToUV(const glm::vec3& normal) {
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


int main() {
  std::vector<glm::vec3> pixel(WIDTH * HEIGHT);

  // カメラ
  glm::mat4 transform = glm::rotate(-0.4f, glm::vec3(0.0f, 1.0f, 0.0f)) *
                        glm::rotate(-0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
    
  glm::vec3 cam_pos{ -2.0f, 2.5f, 5.0f };
  glm::vec3 cam_dir = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 1.0f }).xyz();
  glm::vec3 cam_up  = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 1.0f }).xyz();
  glm::vec3 cam_side = glm::cross(cam_dir, cam_up);
  // 焦点距離
  float focus = 2.5f;

  // 光源
  glm::vec3 light_dir{ -0.577f, 0.577f, 0.577f };
  
  // 画面解像度
  glm::vec2 resolution{ float(WIDTH), float(HEIGHT) };

  // IBL用画像
  int bg_x, bg_y, bg_comp;
  float* ibl_image = stbi_loadf("bg.hdr", &bg_x, &bg_y, &bg_comp, 0);

  glm::vec2 texture(bg_x - 1, bg_y - 1);
  
  // 全ピクセルでレイマーチ
  for (int y = 0; y < HEIGHT; ++y) {
    glm::vec2 coord;
    coord.y = (HEIGHT - 1) - y;
    
    for (int x = 0; x < WIDTH; ++x) {
      coord.x = x;

      glm::vec2 pos = (coord.xy() * 2.0f - resolution) / resolution.y;
      glm::vec3 ray_dir = glm::normalize(cam_side * pos.x + cam_up * pos.y + cam_dir * focus);

      float t = 0.0;
      float d;
      glm::vec3 pos_on_ray = cam_pos;

      // レイマーチで交差点を調べる
      for(int i = 0; i < 256; ++i) {
        d = getDistance(pos_on_ray);
        t += d;
        pos_on_ray = cam_pos + t * ray_dir;
      }

      float shadow = 1.0f;
      glm::vec3 color(0.0f);

      if(glm::abs(d) < 0.001) {
        glm::vec3 normal = getNormal(pos_on_ray);
        float diff = glm::clamp(glm::dot(light_dir, normal), 0.1f, 1.0f);

        // よくある環境光とスペキュラーの計算
        glm::vec3 halfLE = glm::normalize(light_dir - ray_dir);
        float spec = glm::pow(glm::clamp(glm::dot(halfLE, normal), 0.0f, 1.0f), 50.0f);

        // 影の計算
        shadow = genShadow(pos_on_ray + normal * 0.001f, light_dir);

        // IBL
        glm::vec2 uv = glm::clamp(vec3ToUV(normal) * texture, glm::vec2(0.0f, 0.0f), texture);
        int index = (int(uv.x) + int(uv.y) * bg_x) * 3;
        glm::vec3 ibl(ibl_image[index], ibl_image[index + 1], ibl_image[index + 2]);
        color = ibl * glm::vec3(diff) + glm::vec3(spec);
      }
      else {
        glm::vec2 uv = glm::clamp(vec3ToUV(ray_dir) * texture, glm::vec2(0.0f, 0.0f), texture);
        int index = (int(uv.x) + int(uv.y) * bg_x) * 3;

        glm::vec3 ibl(ibl_image[index], ibl_image[index + 1], ibl_image[index + 2]);
        color = ibl;
      }
      
      pixel[x + y * WIDTH] = color * glm::max(0.5f, shadow);
    }
  }

  
  // レンダリング結果をPNG形式に
  std::vector<glm::u8vec3> bitmap(WIDTH * HEIGHT);
  for (size_t i = 0; i < pixel.size(); ++i) {
    bitmap[i] = glm::u8vec3(glm::clamp(expose(pixel[i], -3.0f) * 255.0f, 0.0f, 255.0f));
  }
  
  stbi_write_png("test.png", WIDTH, HEIGHT, 3, bitmap.data(), WIDTH * COMPONENTS);
}
