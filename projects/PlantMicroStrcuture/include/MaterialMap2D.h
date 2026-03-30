//
// Created by samundra on 12/18/25.
//

#ifndef MATERIALMAP2D_H
#define MATERIALMAP2D_H
// MaterialMap2D.hpp
// Header-only utilities to:
//  1) Load a PNG mask/grayscale image (optionally via stb_image)
//  2) Query a Gauss point (x,y) -> phase value (0..1) using bilinear sampling
//  3) Map phase -> (E, nu) using either:
//       - HARD_THRESHOLD: binary fiber/matrix via threshold
//       - SMOOTH_MIXTURE: continuous mixture using grayscale value
//
// ---------------------------------------------
// PNG loading:
//   Option A (recommended): stb_image
//     In ONE .cpp file:
//       #define STB_IMAGE_IMPLEMENTATION
//       #include "stb_image.h"
//     Then everywhere you want to load PNGs:
//       #define MATERIALMAP2D_USE_STB_IMAGE
//       #include "stb_image.h"   // before this header or in a common include
//       #include "MaterialMap2D.hpp"
//
//   Option B: Provide pixels yourself (no stb_image), and use Image(w,h,c,data).
//
// ---------------------------------------------
// Typical use:
//   MaterialMap2D::Domain2D dom{0.0, Lx, 0.0, Ly};
//   MaterialMap2D::PhaseProps props;
//   props.E_f = Ef; props.nu_f = nuf;
//   props.E_m = Em; props.nu_m = num;
//   props.threshold01 = 0.5;
//   props.mode = MaterialMap2D::MixingMode::HARD_THRESHOLD; // or SMOOTH_MIXTURE
//
//   auto img = MaterialMap2D::Image::LoadPNG("mask.png", /*forceGray=*/false);
//   MaterialMap2D::Map map(dom, props, std::move(img), /*flipY=*/true);
//
//   auto mp = map.query(xg, yg); // -> mp.E, mp.nu, mp.phase01, mp.isFiber

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "stb_image.h"

namespace MaterialMap2D {

struct Domain2D {
  double x0 = 0.0, x1 = 1.0; // physical bounds
  double y0 = 0.0, y1 = 1.0;
};

enum class MixingMode {
  HARD_THRESHOLD, // phase >= threshold -> fiber, else matrix
  SMOOTH_MIXTURE  // E,nu are mixed continuously using phase in [0,1]
};

struct PhaseProps {
  double E_f  = 1.0;   // fiber Young's modulus
  double nu_f = 0.30;  // fiber Poisson ratio
  double E_m  = 1.0;   // matrix Young's modulus
  double nu_m = 0.30;  // matrix Poisson ratio

  double threshold01 = 0.5;                 // used by HARD_THRESHOLD (and for isFiber flag in SMOOTH_MIXTURE)
  MixingMode mode = MixingMode::HARD_THRESHOLD;
};

struct MatProps {
  double E = 0.0;
  double nu = 0.0;
  double phase01 = 0.0; // sampled grayscale/mask intensity in [0,1]
  bool isFiber = false; // classification (hard) or diagnostic (smooth)
};

static inline double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

static inline int clampi(int v, int lo, int hi) {
  return std::max(lo, std::min(hi, v));
}

// Minimal image container: pixels are uint8_t, row-major, interleaved channels.
struct Image {
  int W = 0;
  int H = 0;
  int C = 0; // channels: 1 (gray), 3 (RGB), 4 (RGBA), etc.
  std::vector<uint8_t> data;

  Image() = default;

  Image(int w, int h, int c, std::vector<uint8_t> pixels)
      : W(w), H(h), C(c), data(std::move(pixels)) {
    if (W <= 0 || H <= 0 || C <= 0) throw std::invalid_argument("MaterialMap2D::Image: invalid shape.");
    if ((int)data.size() != W * H * C) throw std::invalid_argument("MaterialMap2D::Image: pixel buffer size mismatch.");
  }

  const uint8_t* ptr() const { return data.data(); }
  uint8_t* ptr() { return data.data(); }

  // Requires stb_image.h included by the user (and STB_IMAGE_IMPLEMENTATION in exactly one .cpp).
  static Image LoadPNG(const std::string& path, bool forceGray = false) {
    int w = 0, h = 0, c = 0;
    int reqC = forceGray ? 1 : 0; // 0 = keep original channels
    unsigned char* p = stbi_load(path.c_str(), &w, &h, &c, reqC);
    if (!p) throw std::runtime_error("MaterialMap2D::Image::LoadPNG: stbi_load failed for: " + path);

    const int outC = (reqC != 0) ? reqC : c;
    std::vector<uint8_t> buf((size_t)w * (size_t)h * (size_t)outC);
    std::copy(p, p + buf.size(), buf.begin());
    stbi_image_free(p);
    return Image(w, h, outC, std::move(buf));
  }

};

// Convert one pixel at (i,j) to intensity in [0,1].
// For grayscale: intensity = gray/255.
// For RGB/RGBA: intensity = (0.299 R + 0.587 G + 0.114 B)/255.
// Alpha is ignored.
static inline double intensity01_at_u8(const Image& img, int i, int j) {
  i = clampi(i, 0, img.W - 1);
  j = clampi(j, 0, img.H - 1);
  const uint8_t* p = img.ptr() + (j * img.W + i) * img.C;

  if (img.C == 1) {
    return p[0] / 255.0;
  } else {
    const double r = p[0] / 255.0;
    const double g = (img.C >= 2 ? p[1] / 255.0 : r);
    const double b = (img.C >= 3 ? p[2] / 255.0 : r);
    return 0.299 * r + 0.587 * g + 0.114 * b;
  }
}

// Bilinear sampling in image coordinates (u,v) where u in [0,W-1], v in [0,H-1].
static inline double sampleBilinear01(const Image& img, double u, double v) {
  u = clampd(u, 0.0, (double)(img.W - 1));
  v = clampd(v, 0.0, (double)(img.H - 1));

  const int i0 = (int)std::floor(u);
  const int j0 = (int)std::floor(v);
  const int i1 = std::min(i0 + 1, img.W - 1);
  const int j1 = std::min(j0 + 1, img.H - 1);

  const double fu = u - i0;
  const double fv = v - j0;

  const double p00 = intensity01_at_u8(img, i0, j0);
  const double p10 = intensity01_at_u8(img, i1, j0);
  const double p01 = intensity01_at_u8(img, i0, j1);
  const double p11 = intensity01_at_u8(img, i1, j1);

  const double p0 = (1.0 - fu) * p00 + fu * p10;
  const double p1 = (1.0 - fu) * p01 + fu * p11;
  return (1.0 - fv) * p0 + fv * p1; // in [0,1]
}

// Map physical (x,y) -> image (u,v).
// flipY=true if your image origin is top-left and physical y increases upward.
static inline void physicalToImageUV(double x, double y,
                                     const Domain2D& dom,
                                     const Image& img,
                                     bool flipY,
                                     double& u, double& v) {
  const double tx = (x - dom.x0) / (dom.x1 - dom.x0);
  const double ty = (y - dom.y0) / (dom.y1 - dom.y0);

  const double cx = clampd(tx, 0.0, 1.0);
  const double cy = clampd(ty, 0.0, 1.0);

  u = cx * (img.W - 1);
  v = (flipY ? (1.0 - cy) : cy) * (img.H - 1);
}

class Map {
 public:
  Map() = default;

  Map(Domain2D dom, PhaseProps props, Image img, bool flipY = true)
      : dom_(dom), props_(props), img_(std::move(img)), flipY_(flipY) {
    if (img_.W <= 0 || img_.H <= 0 || img_.C <= 0) {
      throw std::invalid_argument("MaterialMap2D::Map: invalid image.");
    }
    if (dom_.x1 == dom_.x0 || dom_.y1 == dom_.y0) {
      throw std::invalid_argument("MaterialMap2D::Map: invalid domain extents.");
    }
  }

  // One-liner per Gauss point.
  MatProps query(double x, double y) const {
    double u, v;
    physicalToImageUV(x, y, dom_, img_, flipY_, u, v);

    const double phase = sampleBilinear01(img_, u, v); // in [0,1]

    MatProps out;
    out.phase01 = phase;

    switch (props_.mode) {
      case MixingMode::HARD_THRESHOLD: {
        const bool isFiber = (phase >= props_.threshold01);
        out.isFiber = isFiber;
        out.E  = isFiber ? props_.E_f  : props_.E_m;
        out.nu = isFiber ? props_.nu_f : props_.nu_m;
        break;
      }

      case MixingMode::SMOOTH_MIXTURE: {
        // Continuous mixture based on grayscale value.
        // If you want a sharper transition, replace w with a smoothstep around threshold.
        const double w = clampd(phase, 0.0, 1.0);
        out.isFiber = (w >= props_.threshold01); // diagnostic only
        out.E  = w * props_.E_f  + (1.0 - w) * props_.E_m;
        out.nu = w * props_.nu_f + (1.0 - w) * props_.nu_m;
        break;
      }

      default:
        throw std::runtime_error("MaterialMap2D::Map::query: unknown MixingMode");
    }

    return out;
  }

  // Convenience: sample phase only (0..1)
  double queryPhase01(double x, double y) const {
    double u, v;
    physicalToImageUV(x, y, dom_, img_, flipY_, u, v);
    return sampleBilinear01(img_, u, v);
  }

  const Domain2D& domain() const { return dom_; }
  const PhaseProps& phaseProps() const { return props_; }
  const Image& image() const { return img_; }
  bool flipY() const { return flipY_; }

 private:
  Domain2D dom_{};
  PhaseProps props_{};
  Image img_{};
  bool flipY_ = true;
};

} // namespace MaterialMap2D

#endif //MATERIALMAP2D_H
