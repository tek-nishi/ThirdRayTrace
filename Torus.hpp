#pragma once


namespace Torus {

float outer_radius;
float inner_radius;


void init(const picojson::value& params) {
  outer_radius = params.get("outer_radius").get<double>();
  inner_radius = params.get("inner_radius").get<double>();
}

std::pair<float, glm::vec4> distance(const glm::vec3& p) {
  glm::vec2 t = glm::vec2(outer_radius, inner_radius);
  glm::vec2 r = glm::vec2(glm::length(p.xz()) - t.x, p.y);

  return std::make_pair(glm::length(r) - t.y, glm::vec4());
}

}
