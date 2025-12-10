// BezierRing.cpp
#include "BezierRing.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/constants.hpp>
#include <random>


std::vector<std::vector<BezierRingFrame>> g_bezierTubeFrames;

// -------------------------------------------------------
// Generic Bernstein utilities
// ------------------------------------------------------

static float binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0f;
    if (k == 0 || k == n) return 1.0f;
    k = std::min(k, n - k);
    double c = 1.0;
    for (int i = 0; i < k; ++i) {
        c *= double(n - i) / double(i + 1);
    }
    return static_cast<float>(c);
}

static void bernsteinBasis(int n, float t, float *B) {
    float s = 1.0f - t;
    for (int i = 0; i <= n; ++i) {
        float coeff = binomial(n, i);
        float ti    = (i == 0) ? 1.0f : std::pow(t,  float(i));
        float sni   = (i == n) ? 1.0f : std::pow(s,  float(n - i));
        B[i] = coeff * ti * sni;
    }
}

static void bernsteinBasisDeriv(int n, float t, float *dB) {
    if (n == 0) {
        dB[0] = 0.0f;
        return;
    }
    float Bnm1[32]; // enough for our small degrees
    bernsteinBasis(n - 1, t, Bnm1);

    for (int i = 0; i <= n; ++i) {
        float left  = (i == 0) ? 0.0f : Bnm1[i - 1];
        float right = (i == n) ? 0.0f : Bnm1[i];
        dB[i] = n * (left - right);
    }
}

BezierRing::BezierRing()
{
    m_param1 = 4;   // stacks (height)
    m_param2 = 16;  // slices (angle)

    initControlNet();  // scalar field s(u,v) ≡ 1 for now
    setVertexData();
}

void BezierRing::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}

// Noiseless control net: s(u,v) ≡ 1 everywhere.
// Later we can perturb m_ctrl[i][j].x around 1 to add organic bumps.
// void BezierRing::initControlNet()
// {
//     const int Nu = CTRL_U;
//     const int Nv = CTRL_V;

//     for (int j = 0; j < Nv; ++j) {
//         for (int i = 0; i < Nu; ++i) {
//             m_ctrl[i][j] = glm::vec3(1.0f, 0.0f, 0.0f); // radius scale = 1
//         }
//     }
// }
void BezierRing::initControlNet()
{
    const int Nu = CTRL_U;
    const int Nv = CTRL_V;

    // Per-instance RNG: each BezierRing gets its own perturbation pattern.
    std::random_device rd;
    std::mt19937 rng(rd());

    // Outward-only radial scale jitter:
    // s = 1.0f + jitter, with jitter in [0, maxJitter].
    const float maxJitter = 0.50f; // up to +20% radius
    std::uniform_real_distribution<float> radialJitter(0.0f, maxJitter);

    for (int j = 0; j < Nv; ++j) {
        for (int i = 0; i < Nu; ++i) {
            float jitter = radialJitter(rng);
            float scale  = 1.0f + jitter;   // strictly outward

            // Store scale in .x; y/z unused for now.
            m_ctrl[i][j] = glm::vec3(scale, 0.0f, 0.0f);
        }
    }
}

// Evaluate scalar scale s(u,v) = Σ B_i(u) B_j(v) w_ij,
// where w_ij = m_ctrl[i][j].x.
float BezierRing::evalScale(float u, float v) const
{
    const int Nu   = CTRL_U;
    const int Nv   = CTRL_V;
    const int degU = Nu - 1;
    const int degV = Nv - 1;

    float Bu[CTRL_U];
    float Bv[CTRL_V];

    bernsteinBasis(degU, u, Bu);
    bernsteinBasis(degV, v, Bv);

    float s = 0.0f;

    for (int i = 0; i < Nu; ++i) {
        for (int j = 0; j < Nv; ++j) {
            s += Bu[i] * Bv[j] * m_ctrl[i][j].x;
        }
    }

    return s;
}

float BezierRing::evalScaleDu(float u, float v) const
{
    const int Nu   = CTRL_U;
    const int Nv   = CTRL_V;
    const int degU = Nu - 1;
    const int degV = Nv - 1;

    float dBu[CTRL_U];
    float Bv[CTRL_V];

    bernsteinBasisDeriv(degU, u, dBu);
    bernsteinBasis(degV, v, Bv);

    float s_u = 0.0f;

    for (int i = 0; i < Nu; ++i) {
        for (int j = 0; j < Nv; ++j) {
            s_u += dBu[i] * Bv[j] * m_ctrl[i][j].x;
        }
    }

    return s_u;
}

float BezierRing::evalScaleDv(float u, float v) const
{
    const int Nu   = CTRL_U;
    const int Nv   = CTRL_V;
    const int degU = Nu - 1;
    const int degV = Nv - 1;

    float Bu[CTRL_U];
    float dBv[CTRL_V];

    bernsteinBasis(degU, u, Bu);
    bernsteinBasisDeriv(degV, v, dBv);

    float s_v = 0.0f;

    for (int i = 0; i < Nu; ++i) {
        for (int j = 0; j < Nv; ++j) {
            s_v += Bu[i] * dBv[j] * m_ctrl[i][j].x;
        }
    }

    return s_v;
}

// Tessellate into triangles with positions + normals.
// u = height (0..1), v = angle (0..1, mapped to 0..2π).
void BezierRing::setVertexData()
{
    // if (m_useExternalFrames && m_frames.size() >= 2) {
    //     const int stacks = static_cast<int>(m_frames.size()) - 1;
    //     const int slices = std::max(3, m_param2);

    //     m_vertexData.clear();

    //     for (int i = 0; i < stacks; ++i) {
    //         const auto& A = m_frames[i];
    //         const auto& B = m_frames[i + 1];

    //         for (int j = 0; j < slices; ++j) {
    //             float v0 = float(j)       / slices;
    //             float v1 = float(j + 1)   / slices;

    //             auto ringPoint = [&](const BezierRingFrame& F, float v) {
    //                 float theta = glm::two_pi<float>() * v;
    //                 return F.center +
    //                        F.radius * (std::cos(theta) * F.normal +
    //                                    std::sin(theta) * F.binormal);
    //             };

    //             auto ringNormal = [&](const BezierRingFrame& F, float v) {
    //                 float theta = glm::two_pi<float>() * v;
    //                 glm::vec3 radial =
    //                     std::cos(theta) * F.normal +
    //                     std::sin(theta) * F.binormal;
    //                 return glm::normalize(radial);
    //             };

    //             // Positions
    //             glm::vec3 p00 = ringPoint(A, v0);
    //             glm::vec3 p01 = ringPoint(A, v1);
    //             glm::vec3 p10 = ringPoint(B, v0);
    //             glm::vec3 p11 = ringPoint(B, v1);

    //             // Smooth per-vertex normals
    //             glm::vec3 n00 = ringNormal(A, v0);
    //             glm::vec3 n01 = ringNormal(A, v1);
    //             glm::vec3 n10 = ringNormal(B, v0);
    //             glm::vec3 n11 = ringNormal(B, v1);

    //             auto emitTriVerts =
    //                 [&](const glm::vec3 &pa, const glm::vec3 &na,
    //                     const glm::vec3 &pb, const glm::vec3 &nb,
    //                     const glm::vec3 &pc, const glm::vec3 &nc) {
    //                     insertVec3(m_vertexData, pa);
    //                     insertVec3(m_vertexData, na);
    //                     insertVec3(m_vertexData, pb);
    //                     insertVec3(m_vertexData, nb);
    //                     insertVec3(m_vertexData, pc);
    //                     insertVec3(m_vertexData, nc);
    //                 };

    //             // Same winding as before, but now with smooth normals
    //             emitTriVerts(p00, n00,  p01, n01,  p11, n11);
    //             emitTriVerts(p00, n00,  p11, n11,  p10, n10);
    //         }
    //     }
    //     return;
    // }

    if (m_useExternalFrames && m_frames.size() >= 2) {
        const int frameCount = static_cast<int>(m_frames.size());

        // --- Tessellation resolution ---
        // param1: total stacks along the tube (u direction)
        // But never below the number of frame spans.
        const int minStacks     = std::max(1, frameCount - 1);
        const int desiredStacks = std::max(1, m_param1);
        const int stacks        = std::max(minStacks, desiredStacks);

        // param2: slices around (v direction), same semantics as canonical.
        const int slices = std::max(3, m_param2);

        m_vertexData.clear();

        // --- Interpolate a frame at global parameter u ∈ [0,1] ---
        auto interpFrameBetween = [](const BezierRingFrame &F0,
                                     const BezierRingFrame &F1,
                                     float t) {
            BezierRingFrame F;
            F.center = (1.0f - t) * F0.center + t * F1.center;
            F.radius = (1.0f - t) * F0.radius + t * F1.radius;

            glm::vec3 n = (1.0f - t) * F0.normal   + t * F1.normal;
            glm::vec3 b = (1.0f - t) * F0.binormal + t * F1.binormal;

            if (glm::length(n) < 1e-6f) n = F0.normal;
            if (glm::length(b) < 1e-6f) b = F0.binormal;

            n = glm::normalize(n);
            b = glm::normalize(b);

            F.normal   = n;
            F.binormal = b;
            return F;
        };

        auto sampleFrame = [&](float u) {
            // u ∈ [0,1] → index in [0, frameCount-1]
            float s = u * float(frameCount - 1);
            int   i0 = static_cast<int>(std::floor(s));
            if (i0 >= frameCount - 1) {
                return m_frames.back();
            }
            float localT = s - float(i0);
            return interpFrameBetween(m_frames[i0], m_frames[i0 + 1], localT);
        };

        // --- Position and normal from frame + Bezier scalar field ---
        auto ringPoint = [&](const BezierRingFrame& F, float u, float v) {
            float theta = glm::two_pi<float>() * v;
            glm::vec3 radial =
                std::cos(theta) * F.normal +
                std::sin(theta) * F.binormal;

            // Bezier scalar field on [0,1]×[0,1]
            float s = evalScale(u, v);

            float r = F.radius * s;
            return F.center + r * radial;
        };

        auto ringNormal = [&](const BezierRingFrame& F, float u, float v) {
            float theta = glm::two_pi<float>() * v;
            glm::vec3 radial =
                std::cos(theta) * F.normal +
                std::sin(theta) * F.binormal;
            return glm::normalize(radial);
        };

        auto emitTriVerts =
            [&](const glm::vec3 &pa, const glm::vec3 &na,
                const glm::vec3 &pb, const glm::vec3 &nb,
                const glm::vec3 &pc, const glm::vec3 &nc) {
                insertVec3(m_vertexData, pa);
                insertVec3(m_vertexData, na);
                insertVec3(m_vertexData, pb);
                insertVec3(m_vertexData, nb);
                insertVec3(m_vertexData, pc);
                insertVec3(m_vertexData, nc);
            };

        // --- Tessellate over (u,v) grid, like canonical code ---
        // u = along the chain, v = around the tube.
        for (int i = 0; i < stacks; ++i) {
            float u0 = float(i)       / float(stacks);
            float u1 = float(i + 1)   / float(stacks);

            BezierRingFrame A = sampleFrame(u0);
            BezierRingFrame B = sampleFrame(u1);

            for (int j = 0; j < slices; ++j) {
                float v0 = float(j)       / float(slices);
                float v1 = float(j + 1)   / float(slices);

                // Positions
                glm::vec3 p00 = ringPoint(A, u0, v0);
                glm::vec3 p01 = ringPoint(A, u0, v1);
                glm::vec3 p10 = ringPoint(B, u1, v0);
                glm::vec3 p11 = ringPoint(B, u1, v1);

                // Normals: radial, like canonical’s outward normals
                glm::vec3 n00 = ringNormal(A, u0, v0);
                glm::vec3 n01 = ringNormal(A, u0, v1);
                glm::vec3 n10 = ringNormal(B, u1, v0);
                glm::vec3 n11 = ringNormal(B, u1, v1);

                // Two tris per quad
                emitTriVerts(p00, n00,  p01, n01,  p11, n11);
                emitTriVerts(p00, n00,  p11, n11,  p10, n10);
            }
        }
        return;
    }

    m_vertexData.clear();

    const float baseRadius = 0.5f;
    const float height     = 1.0f;

    int stacks = std::max(1, m_param1); // u: height
    int slices = std::max(3, m_param2); // v: angle

    auto computePos = [&](float u, float v) {
        float s = evalScale(u, v);   // radial scale
        float y = -0.5f * height + u * height;

        float theta = glm::two_pi<float>() * v;
        float c = std::cos(theta);
        float sn = std::sin(theta);

        float r = baseRadius * s;

        return glm::vec3(r * c, y, r * sn);
    };

    auto computeDuDv = [&](float u, float v,
                           glm::vec3 &du, glm::vec3 &dv) {
        float s   = evalScale(u, v);
        float s_u = evalScaleDu(u, v);
        float s_v = evalScaleDv(u, v);

        float theta      = glm::two_pi<float>() * v;
        float dtheta_dv  = glm::two_pi<float>();
        float c          = std::cos(theta);
        float sn         = std::sin(theta);

        float r   = baseRadius * s;
        float r_u = baseRadius * s_u;
        float r_v = baseRadius * s_v;

        // P(u,v) = (r c, y(u), r sn)
        //   y(u) = -0.5*height + u*height
        //
        // ∂P/∂u:
        //   dr/du = r_u
        //   dy/du = height
        du = glm::vec3(
            r_u * c,
            height,
            r_u * sn
            );

        // ∂P/∂v:
        //   dr/dv = r_v
        //   dθ/dv = 2π
        dv = glm::vec3(
            r_v * c - r * sn * dtheta_dv,
            0.0f,
            r_v * sn + r * c * dtheta_dv
            );
    };

    auto pushVertex = [&](float u, float v) {
        glm::vec3 p, du, dv;
        p = computePos(u, v);
        computeDuDv(u, v, du, dv);

        glm::vec3 n = glm::cross(du, dv);
        float len = glm::length(n);
        if (len < 1e-6f) {
            n = glm::vec3(0.0f, 1.0f, 0.0f);
        } else {
            n /= len;
        }

        // Ensure outward normal (radially)
        glm::vec3 out(p.x, 0.0f, p.z);
        if (glm::dot(n, out) < 0.0f) {
            n = -n;
        }

        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, n);
    };

    auto emitTriUV = [&](float u0, float v0,
                         float u1, float v1,
                         float u2, float v2) {
        glm::vec3 p0 = computePos(u0, v0);
        glm::vec3 p1 = computePos(u1, v1);
        glm::vec3 p2 = computePos(u2, v2);

        glm::vec3 g = glm::cross(p1 - p0, p2 - p0);
        glm::vec3 c = (p0 + p1 + p2) / 3.0f;
        glm::vec3 out(c.x, 0.0f, c.z);

        if (glm::dot(g, out) < 0.0f) {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }

        pushVertex(u0, v0);
        pushVertex(u1, v1);
        pushVertex(u2, v2);
    };

    // Note: u = height (stacks), v = angle (slices).
    for (int i = 0; i < stacks; ++i) {
        float u0 = float(i)       / float(stacks);
        float u1 = float(i + 1)   / float(stacks);

        for (int j = 0; j < slices; ++j) {
            float v0 = float(j)       / float(slices);
            float v1 = float(j + 1)   / float(slices);

            emitTriUV(u0, v0,  u0, v1,  u1, v1);
            emitTriUV(u0, v0,  u1, v1,  u1, v0);
        }
    }
}

void BezierRing::setFrames(const std::vector<BezierRingFrame>& frames)
{
    m_frames = frames;
    m_useExternalFrames = true;
    setVertexData();
}

void BezierRing::updateParams(int param1, int param2)
{
    m_param1 = std::max(1, param1); // stacks (height)
    m_param2 = std::max(3, param2); // slices (angle)

    setVertexData();
}
