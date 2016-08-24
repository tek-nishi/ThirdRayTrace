#pragma once

//
// Mandelbox
//


namespace Mandelbox {

int Iterations;
int ColorIterations;
float MinRad2;
float Scale;
glm::vec3 RotVector;
float RotAngle;
float Threshold;

glm::vec4 scale;
glm::mat3 rot;
float absScalem1;
float AbsScaleRaisedTo1mIters;


void init(const picojson::value& params) {
  Iterations      = params.get("Iterations").get<double>();
  ColorIterations = params.get("ColorIterations").get<double>();
  MinRad2         = params.get("MinRad2").get<double>();
  Scale           = params.get("Scale").get<double>();
  RotVector       = getVec<glm::vec3>(params.get("RotVector"));
  RotAngle        = glm::radians(params.get("RotAngle").get<double>());
  Threshold       = params.get("Threshold").get<double>();
  
  scale = glm::vec4(Scale, Scale, Scale, glm::abs(Scale)) / MinRad2;
  rot   = glm::mat3_cast(glm::angleAxis(RotAngle, normalize(RotVector)));

  absScalem1 = glm::abs(Scale - 1.0f);
  AbsScaleRaisedTo1mIters = glm::pow(glm::abs(Scale), float(1 - Iterations));
}

// Compute the distance from `pos` to the Mandelbox.
float distance(const glm::vec3& pos) {
  glm::vec4 p = glm::vec4(pos,1), p0 = p;  // p.w is the distance estimate
	
	for (int i = 0; i < Iterations; i++) {
    glm::vec3 pp(p.xyz());
    pp = rot * pp;
    pp = glm::clamp(pp, -1.0f, 1.0f) * 2.0f - pp;  // min;max;mad
		float r2 = glm::dot(pp, pp);
		//if (i<ColorIterations) orbitTrap = min(orbitTrap, abs(vec4(p.xyz,r2)));

    // FIXME:glmだと p.xyz = pp と書けない
    p.x = pp.x;
    p.y = pp.y;
    p.z = pp.z;
    
		p *= glm::clamp(glm::max(MinRad2 / r2, MinRad2), 0.0f, 1.0f);  // dp3,div,max.sat,mul
		p = p * scale + p0;
    if (r2 > Threshold) break;
	}

  return (glm::length(p.xyz()) - absScalem1) / p.w - AbsScaleRaisedTo1mIters;
}

}
