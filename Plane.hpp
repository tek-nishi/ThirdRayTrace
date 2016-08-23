#pragma once

//
// 無限平面
//


namespace Plane {

float height;
glm::vec3 normal;


void init(const picojson::value& params) {
  height = params.get("height").get<double>();
  normal = glm::normalize(getVec<glm::vec3>(params.get("normal")));
}

float distance(const glm::vec3& p) {
  return glm::dot(p, glm::vec3(normal)) + height;
}

}
