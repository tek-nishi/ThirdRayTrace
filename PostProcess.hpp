#pragma once

//
// ポストプロセス
//


namespace Expose {

float value;

// 露出計算
// exposure 露出値(マイナス値)
std::vector<glm::vec3> process(const std::vector<glm::vec3>& in) {
  std::vector<glm::vec3> out(in.size());
  
  for (size_t i = 0; i < in.size(); ++i) {
    out[i] = glm::vec3(1.0f) - glm::exp(in[i] * value);
  }
    
  return out;
}

}


namespace GaussianFilter {

std::vector<glm::vec3> process(const std::vector<glm::vec3>& in,
                               const int width, const int height, float sig) {
  std::vector<glm::vec3> tmp(in.size());
  std::vector<glm::vec3> result(in.size());
  
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
        sum += msk[i] * in[x + p * width];
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

}


namespace DOF {

// 当該ピクセルをぼかしたのを取得
glm::vec3 bokeh(const int pixel_x, const int pixel_y, float sig,
                const std::vector<glm::vec3>& in, const int width, const int height) {
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

  std::vector<glm::vec3> tmp(Wm);

  // 垂直方向
  for (int x = 0; x < Wm; ++x) {
    glm::vec3 sum;
    int ofs_x = glm::clamp(pixel_x + x - Rm, 0, width - 1);
    for (int i = 0; i < Wm; ++i) {
      int p = glm::clamp(pixel_y + i - Rm, 0, height - 1);
      sum += msk[i] * in[ofs_x + p * width];
    }
    tmp[x] = sum;
  }
  
  //水平方向
  glm::vec3 sum;
  for (int i = 0; i < Wm; ++i) {
    sum += msk[i] * tmp[i];
  }
  
  return sum;
}



std::vector<glm::vec3> process(const std::vector<glm::vec3>& in,
                               const std::vector<float>& depth,
                               const int width, const int height) {
  std::vector<glm::vec3> pixel(in.size());

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int offset = x + y * width;
      
      float sig = glm::clamp(glm::abs(depth[offset] - 1.85f) * 2.0f, 0.5f, 10.0f);
      pixel[offset] = bokeh(x, y, sig, in, width, height);
    }
  }

  return pixel;
}

}


namespace Bloom {

float threshold;
float sig;
float intensity;


std::vector<glm::vec3> process(const std::vector<glm::vec3>& in,
                               const int width, const int height) {
  // 高輝度を抽出
  std::vector<glm::vec3> hi_brightness(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    hi_brightness[i] = in[i] * glm::step(threshold, in[i]) * intensity;
  }

  // ガウシアンフィルタでぼかす
  auto result = GaussianFilter::process(hi_brightness, width, height, sig);

  // 合成
  for (size_t i = 0; i < in.size(); ++i) {
    result[i] += in[i];
  }

  return result;
}

}
