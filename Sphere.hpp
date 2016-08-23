#pragma once

// 
// 球体
//

namespace Sphere {

float radius;


void init(const picojson::value& params) {
  radius = params.get("radius").get<double>();
}

float distance(const glm::vec3& p) {
  return glm::length(p) - radius;
}

}
