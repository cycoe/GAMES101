#ifndef __RAY_TRAING_MICROFACET_HPP__
#define __RAY_TRAING_MICROFACET_HPP__

#include "../Vector.hpp"
#include "../global.hpp"
#include "../Material.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

inline float get_cos_theta(Vector3f const& w) {
    return w.z;
}

inline float get_cos_theta_2(Vector3f const& w) {
    return w.z * w.z;
}

inline float get_sin_theta_2(Vector3f const& w) {
    return std::max(0.f, 1.f - get_cos_theta_2(w));
}

inline float get_sin_theta(Vector3f const& w) {
    return std::sqrt(get_sin_theta_2(w));
}

inline float get_cos_phi(Vector3f const& w) {
    float sin_theta = get_sin_theta(w);
    return (sin_theta == 0) ? 1.f : clamp(-1.f, 1.f, w.x / sin_theta);
}

inline float get_sin_phi(Vector3f const& w) {
    float sin_theta = get_sin_theta(w);
    return (sin_theta == 0) ? 0.f : clamp(-1.f, 1.f, w.y / sin_theta);
}

struct TrowbridgeReitzDistribution {
  TrowbridgeReitzDistribution(float alpha_x, float alpha_y)
    : m_alpha_x(alpha_x), m_alpha_y(alpha_y) { }
  TrowbridgeReitzDistribution(TrowbridgeReitzDistribution const&) = default;
  inline float D(Vector3f const& wh) const;
  inline float G1(Vector3f const& w) const;
  inline float Lambda(Vector3f const& w) const;
  inline float G(Vector3f const& wo, Vector3f const& wi) const;

  float m_alpha_x;
  float m_alpha_y;
};

struct MicrofacetReflection : public Material {
  MicrofacetReflection(Vector3f const &R, float alpha_x, float alpha_y,
                       TrowbridgeReitzDistribution const &distribute, bool sample_visible_area=true)
      : m_R(R), m_alpha_x(alpha_x), m_alpha_y(alpha_y),
        m_distribution(distribute), m_sample_visible_area(sample_visible_area) {}
  // sample a ray by Material properties
  virtual Vector3f sample(const Vector3f &wi) override;
  // given a ray, calculate the PdF of this ray
  virtual float pdf(const Vector3f &wi, const Vector3f &wo) override;
  // given a ray, calculate the contribution of this ray
  virtual Vector3f eval(const Vector3f &coords, const Vector3f &wi, const Vector3f &wo) override;
  virtual MaterialType getType() override { return GXX; };

  Vector3f m_R;
  float m_alpha_x;
  float m_alpha_y;
  TrowbridgeReitzDistribution m_distribution;
  bool m_sample_visible_area;
};

float TrowbridgeReitzDistribution::D(const Vector3f &wh) const {
  float cos_theta = wh.z;
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float e = 1 + tan_theta_2 / (m_alpha_x * m_alpha_y);
    return 1.f / (M_PI * m_alpha_x * m_alpha_y * cos_theta_2 * cos_theta_2 * e * e);
  } else {
    float cos_phi = get_cos_phi(wh);
    float cos_phi_2 = cos_phi * cos_phi;
    float sin_phi_2 = 1 - cos_phi_2;
    float e = 1 + tan_theta_2 * (cos_phi_2 / (m_alpha_x * m_alpha_x) +
                                 sin_phi_2 / (m_alpha_y * m_alpha_y));
    return 1.f / (M_PI * m_alpha_x * m_alpha_y * cos_theta_2 * cos_theta_2 * e * e);
  }
}

float TrowbridgeReitzDistribution::G1(const Vector3f &w) const {
  return 1 / (1 + this->Lambda(w));
}

float TrowbridgeReitzDistribution::Lambda(const Vector3f &w) const {
  float cos_theta = w.z;
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float alpha_tan_theta_2 = m_alpha_x * m_alpha_y * tan_theta_2;
    return (-1 + std::sqrt(1.f + alpha_tan_theta_2)) / 2;
  } else {
    float cos_phi = get_cos_phi(w);
    float cos_phi_2 = cos_phi * cos_phi;
    float sin_phi_2 = 1 - cos_phi_2;
    float alpha_2 = cos_phi_2 * m_alpha_x * m_alpha_x + sin_phi_2 * m_alpha_y * m_alpha_y;
    return (-1 + std::sqrt(1.f + alpha_2 * tan_theta_2)) / 2;
  }
}

float TrowbridgeReitzDistribution::G(const Vector3f &wo, const Vector3f &wi) const {
  return 1 / (1 + Lambda(wo) + Lambda(wi));
  //return G1(wo, N) * G1(wi, N);
}

inline void TrowbridgeReitzSample11(float cos_theta, float r1, float r2, float* slope_x, float* slope_y) {
  /// special case (normal incidence)
  if (cos_theta > .9999) {
    float r = std::sqrt(r1 / (1 - r1));
    float phi = 2 * M_PI * r2;
    *slope_x = r * std::cos(phi);
    *slope_y = r * std::sin(phi);
    return;
  }

  float sin_theta = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
  float tan_theta = sin_theta / cos_theta;
  float a = 1.f / tan_theta;
  float G1 = 2 / (1.f + std::sqrt(1.f + 1.f / (a * a)));

  /// sample slope_x
  float A = 2 * r1 / G1 - 1;
  float tmp = 1.f / (A * A - 1.f);
  if (tmp > 1e10) tmp = 1e10;
  float B = tan_theta;
  float D = std::sqrt(std::max(0.f, B * B * tmp * tmp - (A * A - B * B) * tmp));
  float slope_x_1 = B * tmp - D;
  float slope_x_2 = B * tmp + D;
  *slope_x = (A < 0 || slope_x_2 > 1.f / tan_theta) ? slope_x_1 : slope_x_2;

  /// sample slope_y
  float S;
  if (r2 > 0.5f) {
    S = 1.f;
    r2 = 2.f * (r2 - 0.5f);
  } else {
    S = -1.f;
    r2 = 2.f * (0.5f - r2);
  }
  float z = (r2 * (r2 * (r2 * 0.27385f - 0.73369f) + 0.46341f)) /
            (r2 * (r2 * (r2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
  *slope_y = S * z * std::sqrt(1.f + *slope_x * * slope_y);
}

inline Vector3f TrowbridgeReitzSample(Vector3f const& wi, float alpha_x, float alpha_y, float r1, float r2) {
  /// 1. stretch wi
  Vector3f wi_stretched = Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalized();

  /// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
  float slope_x, slope_y;
  TrowbridgeReitzSample11(wi_stretched.z, r1, r2, &slope_x, &slope_y);

  /// 3. rotate
  float cos_theta = wi_stretched.z;
  float r = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
  float tmp = wi_stretched.x / r * slope_x - wi_stretched.y / r * slope_y;
  slope_y = wi_stretched.y / r * slope_x + wi_stretched.x / r * slope_y;
  slope_x = tmp;

  /// 4. unstretch
  slope_x = alpha_x * slope_x;
  slope_y = alpha_y * slope_y;

  /// 3. compute normal
  return Vector3f(-slope_x, -slope_y, 1.f).normalized();
}

Vector3f MicrofacetReflection::sample(const Vector3f &wi) {
  float r1 = get_random_float();
  float r2 = get_random_float();
  Vector3f wh;
  if (!m_sample_visible_area) {
    float cos_theta = 0.f;
    float phi = 2 * M_PI * r2;
    if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
      float tan_theta_2 = m_alpha_x * m_alpha_y * r1 / (1.f - r1);
      cos_theta = 1 / std::sqrt(1 + tan_theta_2);
    } else {
      /// TODO
    }
    float sin_theta = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
    wh = Vector3f(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    if (wh.z * wi.z < 0.f) wh = -wh;
  } else {
    bool flip = wi.z < 0;
    wh = TrowbridgeReitzSample(flip ? -wi : wi, m_alpha_x, m_alpha_y, r1, r2);
    if (flip) wh = -wh;
  }
  return 2 * dotProduct(wi, wh) * wh - wi;
}

float MicrofacetReflection::pdf(const Vector3f &wi, const Vector3f &wo) {
  Vector3f wh = (wi + wo).normalized();
  if (m_sample_visible_area) {
    return m_distribution.D(wh) * m_distribution.G1(wi) * std::abs(dotProduct(wi, wh)) / (4 * wi.z * wi.z);
  } else {
    return m_distribution.D(wh) * std::abs(wh.z) / (4 * wi.z);
  }
}

Vector3f MicrofacetReflection::eval(const Vector3f &coords, const Vector3f &wi, const Vector3f &wo) {
  float cos_theta_i = std::abs(wi.z);
  float cos_theta_o = std::abs(wo.z);
  Vector3f wh = wi + wo;

  if (cos_theta_i < 0.001f || cos_theta_o <  0.001f) return Vector3f(0.f);
  if (wh.norm() < 0.001f) return Vector3f(0.f);

  wh = wh.normalized();
  float d = m_distribution.D(wh);
  float g = m_distribution.G(wi, wo);
  float f = 0.98f;
  Vector3f R = d * g * f / (4 * cos_theta_i * cos_theta_o) * m_R;

  return R;
}

#endif  // __RAY_TRAING_MICROFACET_HPP__