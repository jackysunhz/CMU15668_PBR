#include <dirt/common.h>
#include <dirt/background.h>
#include <dirt/scene.h>
#include <dirt/surface.h>
#include <filesystem/resolver.h>

namespace
{

auto g_defaultBackground = make_shared<ConstantBackground>(Color3f(1.f, 1.f, 1.f));

} // namespace


shared_ptr<Background> Background::defaultBackground()
{
	return g_defaultBackground;
}

ConstantBackground::ConstantBackground(const json & j)
{  
    try
    {
        j.get_to(color);
    }
    catch (...)
    {
        color = j.value("color", color);
    }
}

ImageBackground::ImageBackground(const json & j)
{
    std::string filename;

    try
    {
        filename = j.at("filename").get<std::string>();
    }
    catch (...)
    {
        error("No \"filename\" specified for ImageBackground.\n", j.dump());
        tex.resize(1, 1);
        tex(0, 0) = Color3f(1.f, 1.f, 1.f);
    }

    std::string path = getFileResolver().resolve(filename).str();
    if (!tex.load(path))
    {
        error("Cannot load ImageBackground \"%s\".\n", path);
        tex.resize(1, 1);
        tex(0, 0) = Color3f(1.f, 1.f, 1.f);
    }
}

Color3f ImageBackground::value(const Ray3f & ray) const
{
    Vec3f dir = normalize(ray.d);
    float phi = atan2(dir.y, dir.x);
    float theta = asin(dir.z);
    float u = (phi + M_PI) / (2 * M_PI);
    float v = (theta + M_PI / 2) / M_PI;

    int i = clamp((int)round((1 - u) * (tex.sizeX() - 1)), 0, tex.sizeX() - 1);
    int j = clamp((int)round((1 - v) * (tex.sizeY() - 1)), 0, tex.sizeY() - 1);
    return tex(i, j);
}
