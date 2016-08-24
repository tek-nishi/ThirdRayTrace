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

std::pair<float, glm::vec4> distance(const glm::vec3& p) {
  return std::make_pair(glm::dot(p, glm::vec3(normal)) + height, glm::vec4());
}

}
