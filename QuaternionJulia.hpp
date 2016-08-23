#pragma once

//
// Quaternion Julia 集合
//

namespace QuaternionJulia {

int Iterations;
float Threshold;
glm::vec4 C;


void init(const picojson::value& params) {
  Iterations = params.get("iterations").get<double>();
  Threshold  = params.get("threshold").get<double>();
  C          = getVec<glm::vec4>(params.get("C"));
}

float distance(const glm::vec3& pos) {
	glm::vec4 p = glm::vec4(pos, 0.0);
	glm::vec4 dp = glm::vec4(1.0, 0.0,0.0,0.0);
	for (int i = 0; i < Iterations; i++) {
		dp = 2.0f * glm::vec4(p.x * dp.x - glm::dot(p.yzw(), dp.yzw()), p.x * dp.yzw() + dp.x * p.yzw() + glm::cross(p.yzw(), dp.yzw()));
		p = glm::vec4(p.x * p.x - glm::dot(p.yzw(), p.yzw()), glm::vec3(2.0f * p.x * p.yzw())) + C;
		float p2 = glm::dot(p, p);
		//orbitTrap = glm::min(orbitTrap, glm::abs(glm::vec4(p.xyz(), p2)));
		if (p2 > Threshold) break;
	}
	float r = glm::length(p);
	return  0.5f * r * glm::log(r) / glm::length(dp);
}

}
