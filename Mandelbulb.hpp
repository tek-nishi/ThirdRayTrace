#pragma once

//
// Mandelbulb集合
// SOURCE:
// 


namespace Mandelbulb {

// 各種設定値
int Iterations;
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
  Power            = params.get("Power").get<double>();
  Bailout          = params.get("Bailout").get<double>();
  AlternateVersion = params.get("AlternateVersion").get<bool>();
  RotVector        = getVec<glm::vec3>(params.get("RotVector"));
  RotAngle         = glm::radians(params.get("RotAngle").get<double>());
  Julia            = params.get("Julia").get<bool>();
  JuliaC           = getVec<glm::vec3>(params.get("JuliaC"));
  
  rot = glm::mat3_cast(glm::angleAxis(RotAngle, normalize(RotVector)));
}

void powN1(glm::vec3& z, float r, float& dr) {
	// extract polar coordinates
	float theta = glm::acos(z.z/r);
	float phi = glm::atan(z.y,z.x);
	dr =  glm::pow( r, Power-1.0)*Power*dr + 1.0;
	
	// scale and rotate the point
	float zr = glm::pow( r,Power);
	theta = theta*Power;
	phi = phi*Power;
	
	// convert back to cartesian coordinates
	z = zr*glm::vec3(glm::sin(theta)*glm::cos(phi), glm::sin(phi)*glm::sin(theta), glm::cos(theta));
}

void powN2(glm::vec3& z, float zr0, float& dr) {
	float zo0 = glm::asin( z.z/zr0 );
	float zi0 = glm::atan( z.y,z.x );
	float zr = glm::pow( zr0, Power-1.0 );
	float zo = zo0 * Power;
	float zi = zi0 * Power;
	dr = zr*dr*Power + 1.0;
	zr *= zr0;
	z  = zr*glm::vec3( glm::cos(zo)*glm::cos(zi), glm::cos(zo)*glm::sin(zi), glm::sin(zo) );
}

float distance(const glm::vec3& pos) {
	glm::vec3 z=pos;
	float r;
	float dr=1.0;
	int i=0;
	r=glm::length(z);
	while(r<Bailout && (i<Iterations)) {
		if (AlternateVersion) {
			powN2(z,r,dr);
		} else {
			powN1(z,r,dr);
		}
		z+=(Julia ? JuliaC : pos);
		r=glm::length(z);
		z = rot * z;
		i++;
	}
//	if ((type==1) && r<Bailout) return 0.0;
	return 0.5*glm::log(r)*r/dr;
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
