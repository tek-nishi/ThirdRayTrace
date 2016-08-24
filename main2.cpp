//
// レイトレ合宿４向け作品
//
// 1. 何もインストールしていないまっさらなマシン上で動作する
// 2. 実行ファイルを叩いたら自動で始まるように(キーボード、マウスの操作を要求する作りにしない)
// 3. おおよそ30秒毎に、レンダリングの途中経過をbmpかpngで連番(000.png, 001.png, ...) で出力
// 4. ネットワーク越しに何かをやるような動作をさせない
// 5. 5分以内に自動で終了
//

#include <chrono>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <thread>
#include <future>

#include <picojson.h>

#define GLM_SWIZZLE
#define GLM_META_PROG_HELPERS
#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>

#include "Misc.hpp"
#include "Render.hpp"
#include "PostProcess.hpp"


int main() {
#ifdef _OPENMP
  std::cout << "Use OpenMP." << std::endl;
#endif
  
  // メインスレッドでは一定時間ごとにレンダリング結果をファイルに書き出す
  // 4分30秒が経過したら最終レンダリング結果を書き出す
  // FIXME:PNGよりBMP書き出しが良さげ(圧縮とかしないので高速)
  // 起動時間を保持
  auto begin_time = std::chrono::system_clock::now();

  // TIPS:値をスレッド間で共有するのでスマポを利用
  auto params = std::make_shared<RenderParams>();
  
  // セッティング
  {
    std::ifstream ifs("settings.json");
    ifs >> params->settings;
  }

  // 解像度
  // TIPS:キャストを減らすためint型とfloat型の両方を用意
  params->iresolution = getVec<glm::ivec2>(params->settings.get("resolution"));
  params->resolution  = glm::vec2(params->iresolution);

  // レンダリング結果格納先
  params->pixel.resize(params->iresolution.x * params->iresolution.y);
  params->complete = false;
  
  // 一定間隔で進捗を書き出す時間と、レンダリング総時間
  auto sleep_duration     = std::chrono::seconds(int64_t(params->settings.get("render_interval").get<double>()));
  int64_t render_duration = params->settings.get("render_duration").get<double>();
  auto end_time           = begin_time + std::chrono::seconds(render_duration);

  // レンダリング用スレッド開始
  {
    std::packaged_task<void()> task(std::bind(render, params));
    auto future = task.get_future();
    std::thread render_thread{ std::move(task) };
    render_thread.detach();
  }
  
  size_t index = 0;
  while (1) {
    // 指定時間スリープ
    // 残り時間が少ない時は、時間を短く
    {
      auto remain_time = end_time - std::chrono::system_clock::now();
      std::this_thread::sleep_for(remain_time < sleep_duration ? remain_time
                                                               : sleep_duration);
    }

    // 起動からの経過時間を取得
    auto current_time = std::chrono::system_clock::now();

    // 一定時間経過かレンダリング完了報告でループを終了
    if ((current_time >= end_time) || params->complete) break;
    
    // 進捗を書き出す
    std::ostringstream path;
    path << "Progress" << std::setw(2) << std::setfill('0') << index << ".bmp";
    index += 1;

    writeProgressImage(path.str(), params->pixel, params->iresolution.x, params->iresolution.y);

    std::cout << path.str() << std::endl;
  }

  // ポストプロセス
  // Bloom
  Bloom::threshold = params->settings.get("bloom_threshold").get<double>();
  Bloom::sig       = params->settings.get("bloom_sig").get<double>();
  Bloom::intensity = params->settings.get("bloom_intensity").get<double>();
  auto bloomed_pixel = Bloom::process(params->pixel, params->iresolution.x, params->iresolution.y);
  
  // 露出
  Expose::value = params->settings.get("exposure").get<double>();
  auto pixel = Expose::process(bloomed_pixel);
  
  // 最終結果を書き出す
  writeFinalImage("Result.bmp", pixel, params->iresolution.x, params->iresolution.y);

  // レンダリング時間と回数を表示
  auto current_time = std::chrono::system_clock::now();
  std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - begin_time);
  
  std::cout << "Finish!! "
            << duration.count() << " sec. "
            << params->render_num << " times renderd."
            << std::endl;
}
