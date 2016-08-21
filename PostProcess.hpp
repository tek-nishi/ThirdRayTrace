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
