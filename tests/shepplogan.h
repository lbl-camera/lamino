

#include <vector>


#ifndef SHEPPLOGAN__H
#define SHEPPLOGAN__H

std::vector<float> sheppLoganPhantom(int size) {
    std::vector<float> phantom(size * size, 0.0f);

    // Define parameters for the ellipses in the Shepp-Logan phantom
    struct Ellipse {
        float A;      // Amplitude
        float a;      // Semi-major axis
        float b;      // Semi-minor axis
        float x0;     // X center
        float y0;     // Y center
        float phi;    // Rotation angle in degrees
    };

    std::vector<Ellipse> ellipses = {
        { 2.0f, 0.69f, 0.92f, 0.0f, 0.0f, 0.0f },
        {-0.98f, 0.6624f, 0.8740f, 0.0f, -0.0184f, 0.0f },
        {-0.02f, 0.1100f, 0.3100f, 0.22f, 0.0f, -18.0f },
        {-0.02f, 0.1600f, 0.4100f, -0.22f, 0.0f, 18.0f },
        { 0.01f, 0.2100f, 0.2500f, 0.0f, 0.35f, 0.0f },
        { 0.01f, 0.0460f, 0.0460f, 0.0f, 0.1f, 0.0f },
        { 0.01f, 0.0460f, 0.0460f, 0.0f, -0.1f, 0.0f },
        { 0.01f, 0.0460f, 0.0230f, -0.08f, -0.605f, 90.0f },
        { 0.01f, 0.0230f, 0.0230f, 0.06f, -0.605f, 90.0f }
    };

    // Generate the phantom
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            float xp = (2.0f * x / (size - 1)) - 1.0f; // Normalize to [-1, 1]
            float yp = (2.0f * y / (size - 1)) - 1.0f; // Normalize to [-1, 1]      
            for (const auto& e : ellipses) {
                float phi_rad = e.phi * 3.14159265f / 180.0f;
                float cos_phi = cos(phi_rad);
                float sin_phi = sin(phi_rad);

                float x_rot = cos_phi * (xp - e.x0) + sin_phi * (yp - e.y0);
                float y_rot = -sin_phi * (xp - e.x0) + cos_phi * (yp - e.y0);

                if ((x_rot * x_rot) / (e.a * e.a) + (y_rot * y_rot) / (e.b * e.b) <= 1.0f) {
                    phantom[y * size + x] += e.A;
                }
            }
        }
    }
    return phantom;
}


#endif  // SHEPPLOGAN__H
