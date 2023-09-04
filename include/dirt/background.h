#pragma once

#include <dirt/fwd.h>
#include <dirt/parser.h>
#include <dirt/image.h>


class Background
{
public:
    /// Return a pointer to a global default background
    static shared_ptr<Background> defaultBackground();

    virtual ~Background() = default;

    virtual Color3f value(const Ray3f & ray) const = 0;
};

class ConstantBackground : public Background
{
public:
    ConstantBackground(const Color3f & c) : color(c) {}
    ConstantBackground(const json & j = json::object());

    Color3f value(const Ray3f & ray) const override { return color; }
    Color3f color;
};

class ImageBackground : public Background
{
public:
    ImageBackground(const json & j = json::object());
    Color3f value(const Ray3f & ray) const override;

    Image3f tex;
};
