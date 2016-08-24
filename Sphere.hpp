#pragma once

// 
// 球体
//

namespace Sphere {

float radius;


void init(const picojson::value& params) {
  radius = params.get("radius").get<double>();
}

std::pair<float, glm::vec4> distance(const glm::vec3& p) {
  return std::make_pair(glm::length(p) - radius, glm::vec4());
}

}
