#pragma once

//
// 色計算
//

class Color {
  glm::vec3 BaseColor;
  float OrbitStrength;
  
  glm::vec4 X;
  glm::vec4 Y;
  glm::vec4 Z;
  glm::vec4 R;

  bool CycleColors;
  float Cycles;

  bool Reflection;
  float ReflectionPower;


  static glm::vec3 cycle(const glm::vec3& c, const float s, const float cycles) {
    return glm::vec3(0.5f) +
      glm::vec3(glm::cos(s * cycles + c.x), glm::cos(s * cycles + c.y), glm::cos(s * cycles + c.z)) * 0.5f;
  }
  

public:
  Color(const picojson::value& params)
    : BaseColor(getVec<glm::vec3>(params.get("BaseColor"))),
      OrbitStrength(params.get("OrbitStrength").get<double>()),
      X(getVec<glm::vec4>(params.get("X"))),
      Y(getVec<glm::vec4>(params.get("Y"))),
      Z(getVec<glm::vec4>(params.get("Z"))),
      R(getVec<glm::vec4>(params.get("R"))),
      CycleColors(params.get("CycleColors").get<bool>()),
      Cycles(params.get("Cycles").get<double>()),
      Reflection(params.get("Reflection").get<bool>()),
      ReflectionPower(params.get("ReflectionPower").get<double>())
  {}

  glm::vec3 get(glm::vec4 orbitTrap) const {
    orbitTrap.w = glm::sqrt(orbitTrap.w);

    glm::vec3 orbitColor;
    if (CycleColors) {
      orbitColor = cycle(X.xyz(), orbitTrap.x, Cycles) * X.w * orbitTrap.x +
                   cycle(Y.xyz(), orbitTrap.y, Cycles) * Y.w * orbitTrap.y +
                   cycle(Z.xyz(), orbitTrap.z, Cycles) * Z.w * orbitTrap.z +
                   cycle(R.xyz(), orbitTrap.w, Cycles) * R.w * orbitTrap.w;
    } else {
      orbitColor = X.xyz() * X.w * orbitTrap.x +
                   Y.xyz() * Y.w * orbitTrap.y +
                   Z.xyz() * Z.w * orbitTrap.z +
                   R.xyz() * R.w * orbitTrap.w;
    }

    return glm::mix(BaseColor, 3.0f * orbitColor,  OrbitStrength);
  }

  bool isReflection() const { return Reflection; }
  float getReflectionPower() const { return ReflectionPower; }
  
};
