#include "lightmodel.h"
#include <cmath>
#include <iostream>

// This file contains the phong function you need to fill in.
// This should be the ONLY file you modify.

// Helper function to convert illumination to RGBA, applying some form of tone-mapping (e.g. clamping) in the process
RGBA toRGBA(const glm::vec4 &illumination) {
    // Task 1
    auto to8 = [](float x) -> uint8_t { // make sure it's in [0,1 ], then map to [0, 255]
        x = std::min(1.f, std::max(0.f, x));
        return static_cast<uint8_t>(std::lround(255.f * x));
    };
    return RGBA{to8(illumination.r), to8(illumination.g), to8(illumination.b), 255};
}

// Calculates the RGBA of a pixel from intersection infomation and globally-defined coefficients
RGBA phong(glm::vec3  position,
           glm::vec3  normal,
           glm::vec3  directionToCamera,
           Material  &material,
           std::vector<LightInfo> &lights,
           Sampler   &reflectionSampler) {
    // Normalizing directions
    normal            = glm::normalize(normal);
    directionToCamera = glm::normalize(directionToCamera);

    // Output illumination (we can ignore opacity)
    glm::vec4 illumination(0, 0, 0, 1);

    // Task 3: add the ambient term
    illumination += ka * glm::vec4(glm::vec3(material.ambient), 0.f);

    for (const LightInfo &light : lights) {
        glm::vec3 Ldir = light.pos - position;
        float d = glm::length(Ldir);
        glm::vec3 L = glm::normalize(Ldir);
        // Task 6: compute the attenuation factor
        float atten = 1.f / (c1 + c2 * d + c3 * d * d);

        // Task 4, task 6: add the diffuse term
        // Task 5, task 6: add the specular term
        float nl = glm::max(0.f, glm::dot(normal, L));
        if (nl > 0.f) {
            glm::vec3 diffRGB = glm::vec3(material.diffuse) * glm::vec3(light.color);
            illumination += atten * kd * nl * glm::vec4(diffRGB, 0.f);

            glm::vec3  R  = 2.f * nl * normal - L;                // reflect(-L, n)
            float rv = glm::max(0.f, glm::dot(glm::normalize(R), directionToCamera));

            glm::vec3 specRGB = glm::vec3(material.specular) * glm::vec3(light.color);
            illumination += atten * ks * static_cast<float>(std::pow(rv, static_cast<float>(material.shininess)))
                            * glm::vec4(specRGB, 0.f);
        }
    }


    // Task 7: uncomment the following lines and correct the reflection term.
    //      The following code uses Sampler::getReflection(glm::vec3 startPosition, glm::vec3 lightDirection)
    //      to get the reflection intensity when "recursively raytracing" in some direction from some position

//    glm::vec3 reflectedRay = directionToCamera - 2.f*glm::dot(directionToCamera, normal)*normal; // <-- fix this calculation
//    illumination += kr * reflectionSampler.getReflection(position, reflectedRay);   // <-- no need to edit this after uncommenting
    glm::vec3 I = -directionToCamera;
    glm::vec3 reflectedRay = glm::normalize(I - 2.f * glm::dot(normal, I) * normal);
    illumination += kr * reflectionSampler.getReflection(position, reflectedRay);


    RGBA returnValue = toRGBA(illumination);
    return returnValue;
}
