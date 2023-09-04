#pragma once

#include <dirt/integrator.h>
#include <dirt/scene.h>
#include <dirt/sampler.h>
#include <dirt/spectrum.h>

typedef SampledSpectrum Spectrum;
class PathTracerSpectral : public Integrator
{
public:
    PathTracerSpectral(const json& j = json::object())
    {
        m_maxBounces = j.value("max_bounces", m_maxBounces);
        SampledSpectrum::Init();
    }

    virtual Color3f Li(const Scene & scene, Sampler &sampler, const Ray3f& ray_) const override
    {
        Spectrum combined;
        for(int i = 0; i < combined.nSamples;i++){
            HitInfo hit;
            Ray3f ray(ray_);
            float curWavelengthSample = (float)i/(float)combined.nSamples;
            ray.lambdaSample = curWavelengthSample;
            //cout << ray.lambdaSample << endl;
            Spectrum throughput(1.f);
            Spectrum result(0.f);
            int depth = 0;
            while (depth++ < m_maxBounces)
            {
                if (!scene.intersect(ray, hit))
                {
                    result += throughput * Spectrum::FromRGB(scene.background(ray));
                    break;
                }

                ScatterRecord srec;
                result += throughput * Spectrum::FromRGB(hit.mat->emitted(ray, hit));      

                Vec2f sample = sampler.next2D();

                if(hit.mat->trackSpectral()){
                    Color3f atten;
                    Ray3f scattered;
                    if(hit.mat->scatter(ray, hit, sample, atten, scattered)){
                        //throughput *= Spectrum::FromRGB(atten);
                        ray = scattered;
                        ray.lambdaSample = curWavelengthSample;
                    }
                    else
                        break;
                }
                else if (hit.mat->sample(ray.d, hit, sample, srec))
                {
                    if (!srec.isSpecular)
                    {
                        Ray3f scat(hit.p, srec.scattered);
                        float pdf = hit.mat->pdf(ray.d, scat.d, hit);
                        if (pdf == 0.0f)
                            break;
                        Color3f value = hit.mat->eval(ray.d, scat.d, hit);
                        throughput *= Spectrum::FromRGB(value) / pdf;
                        ray = scat;
                        ray.lambdaSample = curWavelengthSample;
                    }
                    else
                    {
                        ray = Ray3f(hit.p, srec.scattered);
                        ray.lambdaSample = curWavelengthSample;
                        //throughput *= Spectrum::FromRGB(srec.attenuation);
                    }
                }
                else
                    break;
            }
            combined[i] = result[i];
        }
   
        Color3f ret;//= result.SingleToRGB(ray.lambdaSample);
        combined.ToRGB(ret);
        return ret;
    }

    // virtual Color3f Li(const Scene & scene, Sampler &sampler, const Ray3f& ray_) const override
    // {
    //     HitInfo hit;
    //     Ray3f ray(ray_);
    //     Spectrum throughput(1.f);
    //     Spectrum result(0.f);
    //     int depth = 0;

    //     while (depth++ < m_maxBounces)
    //     {
    //         if (!scene.intersect(ray, hit))
    //         {
    //             result += throughput * Spectrum::FromRGB(scene.background(ray));
    //             break;
    //         }

    //         ScatterRecord srec;
    //         result += throughput * Spectrum::FromRGB(hit.mat->emitted(ray, hit));


    //         // if(hit.mat->isEmissive()){
    //         //     //cout << hit.mat->emitted(ray, hit)<<endl;
    //         //     Color3f tmp;
    //         //     Spectrum::FromRGB(Color3f(0.5f)).ToRGB(tmp);
    //         //     //cout << tmp << endl << endl;
    //         //     //cout << Spectrum::FromRGB(hit.mat->emitted(ray, hit)) << endl;
    //         // }
            

    //         Vec2f sample = sampler.next2D();

    //         if(hit.mat->trackSpectral()){
    //             Color3f atten;
    //             Ray3f scattered;
    //             if(hit.mat->scatter(ray, hit, sample, atten, scattered)){
    //                 //throughput *= Spectrum::FromRGB(atten);
    //                 ray = scattered;
    //             }
    //             else
    //                 break;
    //         }
    //         else if (hit.mat->sample(ray.d, hit, sample, srec))
    //         {
    //             if (!srec.isSpecular)
    //             {
    //                 Ray3f scat(hit.p, srec.scattered);
    //                 float pdf = hit.mat->pdf(ray.d, scat.d, hit);
    //                 if (pdf == 0.0f)
    //                     break;
    //                 Color3f value = hit.mat->eval(ray.d, scat.d, hit);
    //                 throughput *= Spectrum::FromRGB(value) / pdf;
    //                 ray = scat;
    //             }
    //             else
    //             {
    //                 ray = Ray3f(hit.p, srec.scattered);
    //                 //throughput *= Spectrum::FromRGB(srec.attenuation);
    //             }
    //         }
    //         else
    //             break;
    //     }


    //     Color3f ret;//= result.SingleToRGB(ray.lambdaSample);
    //     result.ToRGB(ret);
    //     return ret;
    // }

private:
    int m_maxBounces = 64;
    bool m_recursive = true;
};
