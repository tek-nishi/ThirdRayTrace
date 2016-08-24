#pragma once

//
// Mandelbulb集合
// 


namespace Mandelbulb {

// 各種設定値
int Iterations;
int ColorIterations;
float Power;
float Bailout;
bool AlternateVersion;
glm::vec3 RotVector;
float RotAngle;
bool Julia;
glm::vec3 JuliaC;

glm::mat3 rot;


void init(const picojson::value& params) {
  Iterations       = params.get("Iterations").get<double>();
  ColorIterations  = params.get("ColorIterations").get<double>();
  Power            = params.get("Power").get<double>();
  Bailout          = params.get("Bailout").get<double>();
  AlternateVersion = params.get("AlternateVersion").get<bool>();
  RotVector        = getVec<glm::vec3>(params.get("RotVector"));
  RotAngle         = glm::radians(params.get("RotAngle").get<double>());
  Julia            = params.get("Julia").get<bool>();
  JuliaC           = getVec<glm::vec3>(params.get("JuliaC"));
  
  rot = glm::mat3_cast(glm::angleAxis(RotAngle, normalize(RotVector)));
}

void powN1(glm::vec3& z, const float r, float& dr) {
	// extract polar coordinates
	float theta = glm::acos(z.z/r);
	float phi = glm::atan(z.y,z.x);
	dr = glm::pow( r, Power-1.0)*Power*dr + 1.0;
	
	// scale and rotate the point
	float zr = glm::pow( r,Power);
	theta = theta*Power;
	phi = phi*Power;
	
	// convert back to cartesian coordinates
	z = zr*glm::vec3(glm::sin(theta)*glm::cos(phi), glm::sin(phi)*glm::sin(theta), glm::cos(theta));
}

void powN2(glm::vec3& z, const float zr0, float& dr) {
	float zo0 = glm::asin( z.z/zr0 );
	float zi0 = glm::atan( z.y,z.x );
	float zr = glm::pow( zr0, Power-1.0 );
	float zo = zo0 * Power;
	float zi = zi0 * Power;
	dr = zr*dr*Power + 1.0;
	zr *= zr0;
	z  = zr*glm::vec3( glm::cos(zo)*glm::cos(zi), glm::cos(zo)*glm::sin(zi), glm::sin(zo) );
}

std::pair<float, glm::vec4> distance(const glm::vec3& pos) {
	glm::vec3 z = pos;
	float r = glm::length(z);
	float dr = 1.0f;
	int i = 0;
  glm::vec4 orbitTrap(10000.0f);
  
	while(r < Bailout && (i < Iterations)) {
		if (AlternateVersion) {
			powN2(z, r, dr);
		} else {
			powN1(z, r, dr);
		}
		z += (Julia ? JuliaC : pos);
		r = glm::length(z);
		z = rot * z;
    if (i < ColorIterations) orbitTrap = glm::min(orbitTrap, glm::abs(glm::vec4(z.x, z.y, z.z, r * r)));
		i++;
	}
//	if ((type==1) && r<Bailout) return 0.0;
	return std::make_pair(0.5f * glm::log(r) * r / dr, orbitTrap);
	/*
	Use this code for some nice intersections (Power=2)
	float a =  max(0.5*log(r)*r/dr, abs(pos.y));
	float b = 1000;
	if (pos.y>0)  b = 0.5*log(r)*r/dr;
	return min(min(a, b),
		max(0.5*log(r)*r/dr, abs(pos.z)));
	*/
}

}
