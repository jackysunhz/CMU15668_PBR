#pragma once

#include <dirt/fwd.h>
#include <dirt/material.h>
#include <dirt/surface.h>
#include <dirt/scene.h>

float FresnelMoment1(float invEta){
	float eta2 = invEta * invEta;
	float eta3 = eta2 * invEta;
	float eta4 = eta3 * invEta;
	float eta5 = eta4 * invEta;
	if (invEta < 1)
        return 0.45966f - 1.73965f * invEta + 3.37668f * eta2 - 3.904945 * eta3 +
               2.49277f * eta4 - 0.68441f * eta5;
    else
        return -4.61686f + 11.1136f * invEta - 10.4646f * eta2 + 5.11455f * eta3 -
               1.27198f * eta4 + 0.12746f * eta5;
}
float FresnelMoment2(float invEta){
	float eta2 = invEta * invEta;
	float eta3 = eta2 * invEta;
	float eta4 = eta3 * invEta;
	float eta5 = eta4 * invEta;
	if (invEta < 1) {
        return 0.27614f - 0.87350f * invEta + 1.12077f * eta2 - 0.65095f * eta3 +
               0.07883f * eta4 + 0.04860f * eta5;
    } else {
        float r_eta = 1 / invEta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
        return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 +
               458.843f * r_eta + 404.557f * invEta - 189.519f * eta2 +
               54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
    }
}

bool CatmullRomWeights(int size, const float *nodes, float x, int *offset,
                       float *weights) {
    // Return _false_ if _x_ is out of bounds
    if (!(x >= nodes[0] && x <= nodes[size - 1])) return false;

    // Search for the interval _idx_ containing _x_
    int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
    *offset = idx - 1;
    float x0 = nodes[idx], x1 = nodes[idx + 1];

    // Compute the $t$ parameter and powers
    float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;

    // Compute initial node weights $w_1$ and $w_2$
    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;

    // Compute first node weight $w_0$
    if (idx > 0) {
        float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        float w0 = t3 - 2 * t2 + t;
        weights[0] = 0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    // Compute last node weight $w_3$
    if (idx + 2 < size) {
        float w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        float w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0;
    }
    return true;
}

float SampleCatmullRom2D(int size1, int size2, const float *nodes1,
                         const float *nodes2, const float *values,
                         const float *cdf, float alpha, float u, float *fval = nullptr,
                         float *pdf = nullptr) {
    // Determine offset and coefficients for the _alpha_ parameter
    int offset;
    float weights[4];
    if (!CatmullRomWeights(size1, nodes1, alpha, &offset, weights)) return 0;

    // Define a lambda function to interpolate table entries
    auto interpolate = [&](const float *array, int idx) {
        float value = 0;
        for (int i = 0; i < 4; ++i)
            if (weights[i] != 0)
                value += array[(offset + i) * size2 + idx] * weights[i];
        return value;
    };

    // Map _u_ to a spline interval by inverting the interpolated _cdf_
    float maximum = interpolate(cdf, size2 - 1);
    u *= maximum;
    int idx =
        FindInterval(size2, [&](int i) { return interpolate(cdf, i) <= u; });

    // Look up node positions and interpolated function values
    float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
    float x0 = nodes2[idx], x1 = nodes2[idx + 1];
    float width = x1 - x0;
    float d0, d1;

    // Re-scale _u_ using the interpolated _cdf_
    u = (u - interpolate(cdf, idx)) / width;

    // Approximate derivatives using finite differences of the interpolant
    if (idx > 0)
        d0 = width * (f1 - interpolate(values, idx - 1)) /
             (x1 - nodes2[idx - 1]);
    else
        d0 = f1 - f0;
    if (idx + 2 < size2)
        d1 = width * (interpolate(values, idx + 2) - f0) /
             (nodes2[idx + 2] - x0);
    else
        d1 = f1 - f0;

    // Invert definite integral over spline segment and return solution

    // Set initial guess for $t$ by importance sampling a linear interpolant
    float t;
    if (f0 != f1)
        t = (f0 - std::sqrt(std::max((float)0, f0 * f0 + 2 * u * (f1 - f0)))) /
            (f0 - f1);
    else
        t = u / f0;
    float a = 0, b = 1, Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        if (!(t >= a && t <= b)) t = 0.5f * (a + b);

        // Evaluate target function and its derivative in Horner form
        Fhat = t * (f0 +
                    t * (.5f * d0 +
                         t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
                              t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 +
               t * (d0 +
                    t * (-2 * d0 - d1 + 3 * (f1 - f0) +
                         t * (d0 + d1 + 2 * (f0 - f1))));

        // Stop the iteration if converged
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        t -= (Fhat - u) / fhat;
    }

    // Return the sample position and function value
    if (fval) *fval = fhat;
    if (pdf) *pdf = fhat / maximum;
    return x0 + width * t;
}

inline float PhaseHG(float cosTheta, float g) {
    float denom = 1 + g * g + 2 * g * cosTheta;
    return 0.25f/M_PI * (1 - g * g) / (denom * std::sqrt(denom));
}

float IntegrateCatmullRom(int n, const float *x, const float *values,
                          float *cdf) {
    float sum = 0;
    cdf[0] = 0;
    for (int i = 0; i < n - 1; ++i) {
        // Look up $x_i$ and function values of spline segment _i_
        float x0 = x[i], x1 = x[i + 1];
        float f0 = values[i], f1 = values[i + 1];
        float width = x1 - x0;

        // Approximate derivatives using finite differences
        float d0, d1;
        if (i > 0)
            d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
        else
            d0 = f1 - f0;
        if (i + 2 < n)
            d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
        else
            d1 = f1 - f0;

        // Keep a running sum and build a cumulative distribution function
        sum += ((d0 - d1) * (1.f / 12.f) + (f0 + f1) * .5f) * width;
        cdf[i + 1] = sum;
    }
    return sum;
}

float InvertCatmullRom(int n, const float *x, const float *values, float u) {
    // Stop when _u_ is out of bounds
    if (!(u > values[0]))
        return x[0];
    else if (!(u < values[n - 1]))
        return x[n - 1];

    // Map _u_ to a spline interval by inverting _values_
    int i = FindInterval(n, [&](int i) { return values[i] <= u; });

    // Look up $x_i$ and function values of spline segment _i_
    float x0 = x[i], x1 = x[i + 1];
    float f0 = values[i], f1 = values[i + 1];
    float width = x1 - x0;

    // Approximate derivatives using finite differences
    float d0, d1;
    if (i > 0)
        d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
    else
        d0 = f1 - f0;
    if (i + 2 < n)
        d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
    else
        d1 = f1 - f0;

    // Invert the spline interpolant using Newton-Bisection
    float a = 0, b = 1, t = .5f;
    float Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        if (!(t > a && t < b)) t = 0.5f * (a + b);

        // Compute powers of _t_
        float t2 = t * t, t3 = t2 * t;

        // Set _Fhat_ using Equation (8.27)
        Fhat = (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 +
               (t3 - 2 * t2 + t) * d0 + (t3 - t2) * d1;

        // Set _fhat_ using Equation (not present)
        fhat = (6 * t2 - 6 * t) * f0 + (-6 * t2 + 6 * t) * f1 +
               (3 * t2 - 4 * t + 1) * d0 + (3 * t2 - 2 * t) * d1;

        // Stop the iteration if converged
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        t -= (Fhat - u) / fhat;
    }
    return x0 + t * width;
}

struct BSSRDFTable{
    const int nRhoSamples, nRadiusSamples;
    std::unique_ptr<float[]> rhoSamples, radiusSamples;
    std::unique_ptr<float[]> profile;
    std::unique_ptr<float[]> rhoEff;
    std::unique_ptr<float[]> profileCDF;

    BSSRDFTable(int nRhoSamples, int nRadiusSamples)
        : nRhoSamples(nRhoSamples), nRadiusSamples(nRadiusSamples),
        rhoSamples(new float[nRhoSamples]),
        radiusSamples(new float[nRadiusSamples]),
        profile(new float[nRadiusSamples * nRhoSamples]),
        rhoEff(new float[nRhoSamples]),
        profileCDF(new float[nRadiusSamples * nRhoSamples]) {}
    
    inline float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * nRadiusSamples + radiusIndex];
    }
};


class BSSRDF
{
public:
	BSSRDF(float eta, const HitInfo& hit_po, Vec3f wo): eta(eta), hit_po(hit_po), wo(wo) {}

    virtual Color3f S(const HitInfo& hit_pi, const Vec3f& wi) = 0;

    /**
     * @brief Sample the BSSRDF.
     * @param scene The scene.
     * @param u1 A random number.
     * @param u2 A random number pair.
     * @param hit The hit information.
     * @param pdf The probability density function.
     * 
     * @return The BSSRDF value for the 2 points and 2 directions
    */
	virtual Color3f SampleS(const Scene& scene, float u1, const Vec2f& u2, HitInfo& hit, float* pdf) const = 0;

protected:
	float eta;							///< The index of refraction.
    const HitInfo& hit_po;
    Vec3f wo;
};

class SeparableBSSRDF : public BSSRDF
{
public:
	SeparableBSSRDF(float eta, const HitInfo& hit_po, Vec3f wo): BSSRDF(eta, hit_po, wo), ns(hit_po.sn), ss(normalize(hit_po.dpdu)), ts(cross(ns, ss)){}

	Color3f S(const HitInfo& hit_pi, const Vec3f& wi) override{
		float Ft = 1 - FrDielectric(dot(wo, hit_po.sn), 1, eta);
		return Ft * Sp(hit_pi.p) * Sw(wi);
	}

	Color3f Sw(const Vec3f& w) const{
		float c = 1 - 2 * FresnelMoment1(1/eta);
		return Color3f((1 - FrDielectric(dot(w, ns), 1, eta)) / (c * M_PI));
	}

	Color3f Sp(const Vec3f& pi) const{
		return Color3f(Sr(length(pi - hit_po.p)));
	}

	Color3f SampleS(const Scene& scene, float u1, const Vec2f& u2, HitInfo& hit_pi, float* pdf) const override{
        
        Color3f Sp = SampleSp(scene, u1, u2, hit_pi, pdf);
        // if (!Sp.IsBlack()) {
        //     // Initialize material model at sampled surface interaction
        //     si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
        //     si->bsdf->Add(ARENA_ALLOC(arena, SeparableBSSRDFAdapter)(this));
        //     si->wo = Vector3f(si->shading.n);
        // }
        return Sp;
    }
	
    Color3f SampleSp(const Scene& scene, float u1, const Vec2f& u2, HitInfo& hit_pi, float* pdf) const{
        // Choose projection axis for BSSRDF sampling
        Vec3f vx, vy, vz;
        if (u1 < .5f) {
            vx = ss;
            vy = ts;
            vz = ns;
            u1 *= 2;
        } else if (u1 < .75f) {
            // Prepare for sampling rays with respect to _ss_
            vx = ts;
            vy = ns;
            vz = ss;
            u1 = (u1 - .5f) * 4;
        } else {
            // Prepare for sampling rays with respect to _ts_
            vx = ns;
            vy = ss;
            vz = ts;
            u1 = (u1 - .75f) * 4;
        }

        // Choose spectral channel for BSSRDF sampling
        int ch = clamp((int)(u1 * 3), 0, 2);
        u1 = u1 * 3 - ch;

        // Sample BSSRDF profile in polar coordinates
        float r = Sample_Sr(ch, u2[0]);
        if (r < 0) return Color3f(0.f);
        float phi = 2 * M_PI * u2[1];

        // Compute BSSRDF profile bounds and intersection height
        float rMax = Sample_Sr(ch, 0.999f);
        if (r >= rMax) return Color3f(0.f);
        float l = 2 * std::sqrt(rMax * rMax - r * r);

        // Compute BSSRDF sampling ray segment
        HitInfo base;
        base.p = hit_po.p + r * (vx * std::cos(phi) + vy * std::sin(phi)) - l * vz * 0.5f;
        base.t = hit_po.t;
        Vec3f pTarget = base.p + l * vz;

        // Intersect BSSRDF sampling ray against the scene geometry

        // Declare _IntersectionChain_ and linked list
        struct HitChain {
            HitInfo si;
            HitChain *next = nullptr;
        };
        HitChain *chain = new HitChain();

        // Accumulate chain of intersections along ray
        HitChain *ptr = chain;
        int nFound = 0;
        while (true) {
            Ray3f r(base.p, pTarget);
            if (r.d == Vec3f(0, 0, 0) || !scene.intersect(r, ptr->si))
                break;

            base = ptr->si;
            // Append admissible intersection to _IntersectionChain_
            if (ptr->si.mat == hit_po.mat) {
                HitChain *next = new HitChain();
                ptr->next = next;
                ptr = next;
                nFound++;
            }
        }

        // Randomly choose one of several intersections during BSSRDF sampling
        if (nFound == 0) return Color3f(0.0f);
        int selected = clamp((int)(u1 * nFound), 0, nFound - 1);
        while (selected-- > 0) chain = chain->next;
        hit_pi = chain->si;

        // Compute sample PDF and return the spatial BSSRDF term $\Sp$
        *pdf = this->Pdf_Sp(hit_pi) / nFound;
        return this->Sp(hit_pi.p);
    }

	float Pdf_Sp(const HitInfo& hit_pi) const{
        // Express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at
        // $\pto$
        Vec3f d = hit_po.p - hit_pi.p;
        Vec3f dLocal(dot(ss, d), dot(ts, d), dot(ns, d));
        Vec3f nLocal(dot(ss, hit_pi.sn), dot(ts, hit_pi.sn), dot(ns, hit_pi.sn));

        // Compute BSSRDF profile radius under projection along each axis
        float rProj[3] = {std::sqrt(dLocal.y * dLocal.y + dLocal.z * dLocal.z),
                        std::sqrt(dLocal.z * dLocal.z + dLocal.x * dLocal.x),
                        std::sqrt(dLocal.x * dLocal.x + dLocal.y * dLocal.y)};

        // Return combined probability from all BSSRDF sampling strategies
        float pdf = 0, axisProb[3] = {.25f, .25f, .5f};
        float chProb = 1 / 3;
        for (int axis = 0; axis < 3; ++axis)
            for (int ch = 0; ch < 3; ++ch)
                pdf += Pdf_Sr(ch, rProj[axis]) * std::abs(nLocal[axis]) * chProb *
                    axisProb[axis];
        return pdf;
    }
	virtual Color3f Sr(float d) const = 0;
	virtual float Sample_Sr(int ch, float u) const = 0;
    virtual float Pdf_Sr(int ch, float r) const = 0;


private:
	const Vec3f ns;						///< The normal of the surface.
	const Vec3f ss;						///< The tangent of the surface.
	const Vec3f ts;						///< The bitangent of the surface.
	shared_ptr<Material> matrial;				///< The material of the surface.
};

class TabulatedBSSRDF : public SeparableBSSRDF
{
public:
    TabulatedBSSRDF(float eta, 
    const HitInfo& hit_po, 
    Vec3f wo, 
    const Color3f& sigma_a, 
    const Color3f& sigma_s, 
    const BSSRDFTable & table)
    : SeparableBSSRDF(eta, hit_po, wo), table(table) {
        sigma_t = sigma_a + sigma_s;
        rho = (sigma_t != Color3f(0)) ? sigma_s / sigma_t : Color3f(0);
    }

    Color3f Sr(float r) const {
        Color3f Sr(0.f);
        for (int ch = 0; ch < 3; ++ch) {
            // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
            float rOptical = r * sigma_t[ch];

            // Compute spline weights to interpolate BSSRDF on channel _ch_
            int rhoOffset, radiusOffset;
            float rhoWeights[4], radiusWeights[4];
            if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(),
                                rho[ch], &rhoOffset, rhoWeights) ||
                !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                                rOptical, &radiusOffset, radiusWeights))
                continue;

            // Set BSSRDF value _Sr[ch]_ using tensor spline interpolation
            float sr = 0;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    float weight = rhoWeights[i] * radiusWeights[j];
                    if (weight != 0)
                        sr += weight *
                            table.EvalProfile(rhoOffset + i, radiusOffset + j);
                }
            }

            // Cancel marginal PDF factor from tabulated BSSRDF profile
            if (rOptical != 0) sr /= 2 * M_PI * rOptical;
            Sr[ch] = sr;
        }
        // Transform BSSRDF value into world space units
        Sr *= sigma_t * sigma_t;
        return clampColor(Sr);
    }

    float TabulatedBSSRDF::Sample_Sr(int ch, float u) const {
        if (sigma_t[ch] == 0) return -1;
        return SampleCatmullRom2D(table.nRhoSamples, table.nRadiusSamples,
                                table.rhoSamples.get(), table.radiusSamples.get(),
                                table.profile.get(), table.profileCDF.get(),
                                rho[ch], u) /
            sigma_t[ch];
    }

    float TabulatedBSSRDF::Pdf_Sr(int ch, float r) const {
        // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
        float rOptical = r * sigma_t[ch];

        // Compute spline weights to interpolate BSSRDF density on channel _ch_
        int rhoOffset, radiusOffset;
        float rhoWeights[4], radiusWeights[4];
        if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch],
                            &rhoOffset, rhoWeights) ||
            !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                            rOptical, &radiusOffset, radiusWeights))
            return 0.f;

        // Return BSSRDF profile density for channel _ch_
        float sr = 0, rhoEff = 0;
        for (int i = 0; i < 4; ++i) {
            if (rhoWeights[i] == 0) continue;
            rhoEff += table.rhoEff[rhoOffset + i] * rhoWeights[i];
            for (int j = 0; j < 4; ++j) {
                if (radiusWeights[j] == 0) continue;
                sr += table.EvalProfile(rhoOffset + i, radiusOffset + j) *
                    rhoWeights[i] * radiusWeights[j];
            }
        }

        // Cancel marginal PDF factor from tabulated BSSRDF profile
        if (rOptical != 0) sr /= 2 * M_PI * rOptical;
        return std::max((float)0, sr * sigma_t[ch] * sigma_t[ch] / rhoEff);
    }

protected:
    const BSSRDFTable& table;
    Color3f sigma_t;//extinction coefficient
    Color3f rho;//albedo
};

float BeamDiffusionSS(float sigma_s, float sigma_a, float g, float eta,
                      float r);
float BeamDiffusionMS(float sigma_s, float sigma_a, float g, float eta,
                      float r);
void ComputeBeamDiffusionBSSRDF(float g, float eta, BSSRDFTable *t);
void SubsurfaceFromDiffuse(const BSSRDFTable &table, const Color3f &rhoEff,
                           const Color3f &mfp, Color3f *sigma_a,
                           Color3f *sigma_s);






