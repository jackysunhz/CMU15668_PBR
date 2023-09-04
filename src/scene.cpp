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

#include <dirt/scene.h>
#include <dirt/progress.h>
#include <dirt/sampler.h>
#include <fstream>

/// Construct a new scene from a json object
Scene::Scene(const json & j)
{
    parseFromJSON(j);
}

/// Read a scene from a json file
Scene::Scene(const string & filename)
{
    // open file
    std::ifstream stream(filename, std::ifstream::in);
    if (!stream.good())
    	throw DirtException("Cannot open file: %s.", filename);

    json j;
    stream >> j;
    parseFromJSON(j);
}

Scene::~Scene()
{
    m_materials.clear();
}


shared_ptr<const Material> Scene::findOrCreateMaterial(const json & jp, const string& key) const
{
    auto it = jp.find(key);
    if (it == jp.end())
        return Material::defaultMaterial();
    
    auto j = it.value();
    if (j.is_string())
    {
        string name = j.get<string>();
        // find a pre-declared material
        auto i = m_materials.find(name);
        if (i != m_materials.end())
            return i->second;
        else
            throw DirtException("Can't find a material with name '%s' here:\n%s", name, jp.dump(4));
    }
    else if (j.is_object())
    {
	    // create a new material
        return parseMaterial(j);
    }
    else
        throw DirtException("Type mismatch: Expecting either a material or material name here:\n%s", jp.dump(4));
}

shared_ptr<const Medium> Scene::findOrCreateMedium(const json & jp, const string& key) const
{
    auto it = jp.find(key);
    if (it == jp.end())
        return nullptr;

    auto j = it.value();
    if (j.is_string())
    {
        string name = j.get<string>();
        // find a pre-declared medium
        auto i = m_media.find(name);
        if (i != m_media.end())
            return i->second;
        else
            throw DirtException("Can't find a medium with name '%s' here:\n%s", name, jp.dump(4));
    }
    else if (j.is_object())
    {
	    // create a new medium
        return parseMedium(j);
    }
    else
        throw DirtException("Type mismatch: Expecting either a medium or medium name here:\n%s", jp.dump(4));
}

shared_ptr<const MediumInterface> Scene::findOrCreateMediumInterface(const json & jp, const string& key) const
{
    auto it = jp.find(key);
    if (it == jp.end())
    {
        return std::make_shared<MediumInterface>();
    }
    std::shared_ptr<const Medium> inside = findOrCreateMedium(jp.at(key), "inside");
    std::shared_ptr<const Medium> outside = findOrCreateMedium(jp.at(key), "outside");
    return std::make_shared<MediumInterface>(inside, outside);
}

// compute the color corresponding to a ray by raytracing
Color3f Scene::recursiveColor(Sampler &sampler, const Ray3f &ray, int depth) const
{
    // Pseudo-code:
    //
	// if scene.intersect:
    //      get emitted color (hint: you can use hit.mat->emitted)
	// 		if depth < MaxDepth and hit_material.scatter(....) is successful:
	//			recursive_color = call this function recursively with the scattered ray and increased depth
	//          return emitted color + attenuation * recursive_color
	//		else
	//			return emitted color;
	// else:
	// 		return background color (hint: look at m_background)
    const int maxDepth = 64;
    HitInfo hit;
    if (intersect(ray, hit))
    {
        Ray3f scattered;
        Color3f attenuation;
        Color3f emitted = hit.mat->emitted(ray, hit);
        Vec2f sample = sampler.next2D();
        if (depth < maxDepth && hit.mat->scatter(ray, hit, sample, attenuation, scattered))
        {
            return emitted + attenuation * recursiveColor(sampler, scattered, depth + 1);
        }
        else
        {
            return emitted;
        }
    }
    else
    {
        return m_background->value(ray);
    }
}

// raytrace an image
Image3f Scene::raytrace() const
{
    // allocate an image of the proper size
    auto image = Image3f(m_camera->resolution().x, m_camera->resolution().y);

    if (m_integrator)
        return integrateImage();

    // Pseudo-code:
    //
        // foreach image row (go over image height)
            // foreach pixel in the row (go over image width)
                // init accumulated color to zero
                // repeat m_imageSamples times:
                    // compute a random point within the pixel (you can just add a random number between 0 and 1
                    //                                          to the pixel coordinate. You can use randf() for this)
                    // compute camera ray
                    // accumulate color raytraced with the ray (by calling recursiveColor)
                // divide color by the number of pixel samples

    // Hint: you can create a Progress object (progress.h) to provide a 
    // progress bar during rendering.

    Progress progress("Rendering", m_camera->resolution().x*m_camera->resolution().y);

    int max_threads = omp_get_max_threads();
    //cout << "max threads: " << max_threads << endl;
    int height_per_thread = ceil(image.height() / max_threads);
    int num_threads = ceil(image.height() / height_per_thread);
    #pragma omp parallel for
    // foreach pixel
    for (int t = 0; t < num_threads; ++t) {
        int start_height = t * height_per_thread;
        int end_height = (t + 1) * height_per_thread;
        int m = start_height * image.width();
        for (int j = start_height; j < end_height && j < image.height(); ++j)
        {
            for (int i = 0; i < image.width(); ++i)
            {
                // init accumulated color
                image(i, j) = Color3f(0.f);

                m_sampler->startPixel();

                // foreach sample
                for (int s = 0; s < m_imageSamples; ++s)
                {
                    // set pixel to the color raytraced with the ray
                    INCREMENT_TRACED_RAYS;
                    Vec2f sample = m_sampler->next2D();
                    image(i, j) += recursiveColor(*m_sampler, m_camera->generateRay(i + sample.x, j + sample.y), 0);
                    m_sampler->startNextPixelSample();
                }
                // scale by the number of samples
                image(i, j) /= float(m_imageSamples);

                #pragma omp critical
                ++progress;
            }
        }
    }
    

	// return the ray-traced image
    return image;
}

Image3f Scene::integrateImage() const
{
    // allocate an image of the proper size
    auto image = Image3f(m_camera->resolution().x, m_camera->resolution().y);

    Progress progress("Rendering", m_camera->resolution().x*m_camera->resolution().y);

    int max_threads = omp_get_max_threads();
    cout << "max threads: " << max_threads << endl;
    omp_set_num_threads(max_threads);

    // int grid_dim = 20;
    // int total_grid = grid_dim * grid_dim;
    // int current_grid = 0;
    // int grid_width = ceil(image.width() / grid_dim);
    // int grid_height = ceil(image.height() / grid_dim);

    // #pragma omp parallel
    // {
    //     while(current_grid < total_grid){
    //         int grid_x = current_grid % grid_dim;
    //         int grid_y = current_grid / grid_dim;
    //         int start_height = grid_y * grid_height;
    //         int end_height = (grid_y + 1) * grid_height;
    //         int start_width = grid_x * grid_width;
    //         int end_width = (grid_x + 1) * grid_width;
    //         int m = start_height * image.width();
    //         for (int j = start_height; j < end_height && j < image.height(); ++j)
    //         {
    //             for (int i = start_width; i < end_width && i < image.width(); ++i)
    //             {
    //                 // init accumulated color
    //                 image(i, j) = Color3f(0.f);

    //                 m_sampler->startPixel();

    //                 // foreach sample
    //                 for (int s = 0; s < m_imageSamples; ++s)
    //                 {
    //                     // set pixel to the color raytraced with the ray
    //                     INCREMENT_TRACED_RAYS;
    //                     Vec2f sample = m_sampler->next2D();
    //                     image(i, j) += m_integrator->Li(*this, *m_sampler, m_camera->generateRay(i + sample.x, j + sample.y));
    //                     m_sampler->startNextPixelSample();
    //                 }
    //                 // scale by the number of samples
    //                 image(i, j) /= m_imageSamples;
    //                 #pragma omp critical
    //                 ++progress;
    //             }
    //         }
    //         #pragma omp atomic
    //         current_grid++;
            
    //     }
    // }

    int height_per_thread = ceil(image.height() / max_threads);
    int num_threads = ceil(image.height() / height_per_thread);
    #pragma omp parallel for
    //for each thread
    for (int t = 0; t < num_threads; ++t) {
        int start_height = t * height_per_thread;
        int end_height = (t + 1) * height_per_thread;
        int m = start_height * image.width();
        //cout << omp_get_thread_num() << endl;
        //for each pixel
        for (int j = start_height; j < end_height && j < image.height(); ++j)
        {
            for (int i = 0; i < image.width(); ++i)
            {
                // init accumulated color
                image(i, j) = Color3f(0.f);

                m_sampler->startPixel();

                // foreach sample
                for (int s = 0; s < m_imageSamples; ++s)
                {
                    // set pixel to the color raytraced with the ray
                    INCREMENT_TRACED_RAYS;
                    Vec2f sample = m_sampler->next2D();
                    Ray3f ray = m_camera->generateRay(i + sample.x, j + sample.y);

                    //sample wavelength here
                    ray.lambdaSample = m_sampler->next1D();

                    image(i, j) += clampColor(m_integrator->Li(*this, *m_sampler, ray));
                    m_sampler->startNextPixelSample();
                }
                // scale by the number of samples
                image(i, j) /= m_imageSamples;
                #pragma omp critical
                ++progress;
            }
        }
    }
	// return the ray-traced image
    return image;
}
