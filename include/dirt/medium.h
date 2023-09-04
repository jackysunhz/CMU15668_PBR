#pragma once

#include <dirt/fwd.h>
#include <dirt/common.h>
#include <dirt/ray.h>
#include <dirt/perlin.h>

struct MediumInteraction
{
  Vec3f p;
  Vec3f wo;
  const Medium* medium = nullptr;

  MediumInteraction() {};

  MediumInteraction(
    const Vec3f &p,
    const Vec3f &wo,
    const Medium* medium
  )
  : p(p)
  , wo(wo)
  , medium(medium) {};

  bool isValid() { return medium != nullptr; };
};

inline Vec3f SphericalDirection(float sinTheta, float cosTheta, float phi)
{
  return Vec3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}
class PhaseFunction {
  public:
    virtual ~PhaseFunction() = default;

    virtual float p(const Vec3f &wo, const Vec3f &wi) const = 0;

    virtual float sample(const Vec3f &wo, Vec3f &wi, const Vec2f &sample) const = 0;
};

class HenyeyGreenstein : public PhaseFunction
{
  public:
    HenyeyGreenstein(const json & j = json::object());

    float p(const Vec3f &wo, const Vec3f &wi) const override;

    float sample(const Vec3f &wo, Vec3f &wi, const Vec2f &sample) const override;

  private:
    float g = 0;
};

class Medium
{
public:
  virtual ~Medium() = default;

  virtual float Tr(const Ray3f &ray, Sampler &sampler) const = 0;

  virtual float Sample(const Ray3f &ray, Sampler &sampler, MediumInteraction &mi) const = 0;

  virtual float density(const Vec3f &p) const = 0;

  std::shared_ptr<PhaseFunction> phase;
};

class HomogeneousMedium: public Medium
{
public:
  HomogeneousMedium(const json &j = json::object());

  float Tr(const Ray3f &ray, Sampler &sampler) const;

  float Sample(const Ray3f &ray, Sampler &sampler, MediumInteraction &mi) const;

  float density(const Vec3f &p) const;

private:
  float sigma_a, sigma_s, sigma_t;
};

class PerlinMedium: public Medium
{
public:
  PerlinMedium(const json &j = json::object());

  float Tr(const Ray3f &ray, Sampler &sampler) const;

  float Sample(const Ray3f &ray, Sampler &sampler, MediumInteraction &mi) const;

private:
  float density(const Vec3f &p) const;

  // params
  float sigma_a, sigma_s;
  float densityScale = 1.0f;
  float densityOffset = 0.0f;
  Vec3f spatialScale = Vec3f(1.0f, 1.0f, 1.0f);

  // computed
  float invMaxDensity, sigma_t;

  Perlin perlin;
};

// MediumInterface Declarations
struct MediumInterface
{
  MediumInterface() : inside(nullptr), outside(nullptr) {}

  MediumInterface(const std::shared_ptr<const Medium> medium) : inside(medium), outside(medium) {}

  MediumInterface(const std::shared_ptr<const Medium> inside, const std::shared_ptr<const Medium> outside): inside(inside), outside(outside) {}

  bool IsMediumTransition() const { return inside != outside; }

  std::shared_ptr<const Medium> getMedium(const Ray3f ray, const HitInfo &hit) const;

  std::shared_ptr<const Medium> inside, outside;
};

Color3f TrL(const Scene &scene, Sampler &sampler, const Ray3f &ray_);
