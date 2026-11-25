#include "formats/dxf/colors.hpp"

#include <cmath>

namespace cadcore::formats::autocad {

Color aci_color(int i) {
    static const Color fixed[10] = {{0, 0, 0, 255},     {255, 0, 0, 255},   {255, 255, 0, 255}, {0, 255, 0, 255},     {0, 255, 255, 255},
                                    {0, 0, 255, 255},   {255, 0, 255, 255}, {255, 255, 255, 255}, {128, 128, 128, 255}, {192, 192, 192, 255}};
    static const std::uint8_t greys[6] = {51, 80, 105, 130, 190, 255};
    if (i >= 1 && i <= 9) return fixed[i];
    if (i >= 250 && i <= 255) return {greys[i - 250], greys[i - 250], greys[i - 250], 255};
    if (i < 10 || i > 249) return fixed[7];
    static const int levels[5] = {255, 204, 153, 127, 76};
    const int hue = (i - 10) / 10, sub = (i - 10) % 10;
    const int v = levels[sub / 2];
    const int lo = sub % 2 ? v / 2 : 0;  // odd entries: half saturation
    // hue in 15 degree steps; channel = lo + (v - lo) * f, f the HSV ramp of that channel
    const double h = hue * 15.0 / 60.0;  // sector 0..6
    auto ramp = [&](double x) {           // value of a channel whose full-on interval is centred at x
        double d = std::fabs(std::fmod(h - x + 6.0, 6.0));
        if (d > 3) d = 6 - d;
        return d <= 1 ? 1.0 : d >= 2 ? 0.0 : 2 - d;
    };
    auto ch = [&](double f) { return std::uint8_t(lo + int(std::floor((v - lo) * f))); };
    return {ch(ramp(0)), ch(ramp(2)), ch(ramp(4)), 255};
}

}  // namespace cadcore::formats::autocad
