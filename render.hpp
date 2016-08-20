#pragma once

//
// レンダラー
//

// HDR読み込み
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


// レンダリングに必要な情報
struct RenderParams {
  picojson::value settings;

  glm::ivec2 iresolution;
  glm::vec2  resolution;

  std::vector<glm::vec3> pixel;
  bool complete;
};


struct Info {
  glm::ivec2 iresolution;
  glm::vec2  resolution;

  glm::vec3 cam_pos;
  glm::vec3 cam_dir;
  glm::vec3 cam_up;
  glm::vec3 cam_side;
  
  float focus;

  float* ibl_image;
  int bg_x, bg_y, bg_comp;
  glm::vec2 texture;

  int iterate;
};


void render(const std::shared_ptr<RenderParams>& params) {
  Info info;

  // TIPS:スレッド間で共有している値へのアクセスは重たいので
  //      必要な値はコピーしておく
  info.iresolution = params->iresolution;
  info.resolution  = params->resolution;
  
  // カメラ
  glm::vec3 cam_rot = glm::radians(getVec<glm::vec3>(params->settings.get("cam_rot")));
  glm::mat4 transform = glm::rotate(cam_rot.y, glm::vec3(0.0f, 1.0f, 0.0f)) *
    glm::rotate(cam_rot.x, glm::vec3(1.0f, 0.0f, 0.0f)) *
    glm::rotate(cam_rot.z, glm::vec3(0.0f, 0.0f, 1.0f)) *
    glm::translate(getVec<glm::vec3>(params->settings.get("cam_pos")));
    
  info.cam_pos  = (transform * glm::vec4{ 0.0f, 0.0f, 0.0f, 1.0f }).xyz();
  info.cam_dir  = (transform * glm::vec4{ 0.0f, 0.0f, -1.0f, 0.0f }).xyz();
  info.cam_up   = (transform * glm::vec4{ 0.0f, 1.0f, 0.0f, 0.0f }).xyz();
  info.cam_side = glm::cross(info.cam_dir, info.cam_up);

  // 焦点距離
  info.focus = params->settings.get("focus").get<double>();
  
  // IBL用画像
  info.ibl_image = stbi_loadf(params->settings.get("bg").get<std::string>().c_str(), &info.bg_x, &info.bg_y, &info.bg_comp, 0);
  info.texture = glm::vec2(info.bg_x - 1, info.bg_y - 1);


  // レンダリング完了!!
  params->complete = true;
  
  stbi_image_free(info.ibl_image);
}
