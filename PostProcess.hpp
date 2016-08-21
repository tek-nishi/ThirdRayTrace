#pragma once

//
// ポストプロセス
//   被写界深度とか
//


// 被写界深度
namespace DOF {


float focal_plane;

// SOURCE:http://d.hatena.ne.jp/colorcle/20090628/1246208145
std::vector<glm::vec3> gaussian_filter(const std::vector<glm::vec3>& pixel, const int width, const int height, float sig) {
  std::vector<glm::vec3> tmp(pixel.size());
  std::vector<glm::vec3> result(pixel.size());
  
  int Wm = int(glm::ceil(3.0f * sig) * 2 + 1);  // 窓幅
  int Rm = (Wm - 1) / 2;                        // 窓半径
    
  // フィルタ
  std::vector<float> msk(Wm);
    
  sig = 2 * sig * sig;
  float div = glm::sqrt(sig * M_PI);
    
  //フィルタの作成
  for (int x = 0; x < Wm; ++x) {
    int p = (x - Rm) * (x - Rm);
    msk[x] = glm::exp(-p / sig) / div;
  }

  // 垂直方向
  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      glm::vec3 sum;
      for (int i = 0; i < Wm; ++i) {
        int p = glm::clamp(y + i - Rm, 0, height - 1);
        sum += msk[i] * pixel[x + p * width];
      }
      tmp[x + y * width] = sum;
    }
  }
  
  //水平方向
  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      glm::vec3 sum;
      for (int i = 0; i < Wm; ++i) {
        int p = glm::clamp(x + i - Rm, 0, width - 1);
        sum += msk[i] * tmp[p + y * width];
      }
      result[x + y * width] = sum;
    }
  }  
  
  return result;
}


std::vector<glm::vec3> exec(const std::vector<glm::vec3>& pixel, const int width, const int height,
                            const std::vector<float>& depth) {

  return gaussian_filter(pixel, width, height, 5.0f);
  
}


}
