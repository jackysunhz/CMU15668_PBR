/*
    This file is part of Dirt, the Dartmouth introductory ray tracer.
    Copyright (c) 2017-2019 by Wojciech Jarosz
    Dirt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.
    Dirt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <dirt/integrator.h>
#include <dirt/medium.h>
#include <dirt/scene.h>
#include <dirt/pdf.h>

class VolpathTracerUni: public Integrator
{
public:
    VolpathTracerUni(const json& j = json::object())
    {
        m_maxBounces = j.value("max_bounces", m_maxBounces);
    }

    void attachMedium(const HitInfo &hit, Ray3f &ray) const
    {
        if (hit.mat == nullptr && hit.mi->IsMediumTransition())
            ray.medium = hit.mi->getMedium(ray, hit);
        ray.maxt = length(hit.p - ray.o) / length(ray.d) + Epsilon;
    }

    virtual Color3f Li(const Scene & scene, Sampler &sampler, const Ray3f& ray_) const override
    {
        Ray3f ray(ray_);
        Color3f throughput(1.f);
        Color3f result(0.f);

        int bounces = 0;
        while (bounces <= m_maxBounces)
        {
            HitInfo hit;
            bool foundIntersection = scene.intersect(ray, hit);
            if (foundIntersection) attachMedium(hit, ray);

            MediumInteraction mi;
            if (ray.medium)
                throughput *= ray.medium->Sample(ray, sampler, mi);

            if (mi.isValid())
            {
                // sample direction (from phase)
                Vec3f wi;
                float phasePdf = mi.medium->phase->sample(mi.wo, wi, sampler.next2D());
                throughput *= mi.medium->phase->p(mi.wo, wi) / phasePdf;

                // spawn new ray
                ray = Ray3f(mi.p, wi).withMedium(ray.medium);
                bounces++;
            }
            else
            {
                if (!foundIntersection)
                {
                    result += throughput * scene.background(ray);
                    break;
                }

                if (hit.mat == nullptr)
                {
                    ray = Ray3f(hit.p, normalize(ray.d));
                    continue;
                }

                result += throughput * hit.mat->emitted(ray, hit);

                // sample direction (from BSDF)
                ScatterRecord srec;
                if (!hit.mat->sample(ray.d, hit, sampler.next2D(), srec)) break;

                Ray3f scat(hit.p, srec.scattered);
                if (!srec.isSpecular)
                {
                    float bsdfPdf = hit.mat->pdf(ray.d, scat.d, hit);
                    if (bsdfPdf == 0.0f) break;
                    throughput *= hit.mat->eval(ray.d, scat.d, hit) / bsdfPdf;
                }
                else
                {
                    throughput *= srec.attenuation;
                }

                ray = scat;
                bounces++;
            }

            float lum = luminance(throughput);
            const float rrThreshold = 1.0f;
            if (lum < rrThreshold)
            {
                float q = std::max((float).05, 1.0f - lum);
                if (sampler.next1D() < q) break;
                throughput /= (1 - q);
            }
        }

        return result;
    }

private:
    int m_maxBounces = 64;
};
