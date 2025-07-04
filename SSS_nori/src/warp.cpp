/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    float u = sample.x();
    float v = sample.y();
    return (u <= 0.5f) ? Point2f(u * 2, v * u * 2) : Point2f((1 - u) * 2, (1 - u) * (1 - v) * 2);
}

float Warp::squareToTentPdf(const Point2f &p) {
    if (p.x() < 0 || p.x() > 1 || p.y() < 0 || p.y() > 1) return 0.0f;
    return (p.x() <= 0.5f) ? (p.x() / 2) : ((1 - p.x()) / 2);
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = sqrt(sample[0]); // Radial coordinate
    float theta = 2 * M_PI * sample[1]; // Angular coordinate

    // Convert to Cartesian coordinates
    return Point2f(r * cos(theta), r * sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (p.x() * p.x() + p.y() * p.y() > 1) {
        return 0.0f; // Outside the disk
    }
    // PDF for a uniform disk is constant and equal to 1/? over the disk
    return 1 / M_PI;
}

Point2f Warp::squareToUniformTriangle(const Point2f& sample) {
    float su0 = sqrt(sample[0]);

    return Point2f(1 - su0, sample[1] * su0);
}

float Warp::squareToUniformTrianglePdf(const Point2f& p) {
    if (p.x() < 0 || p.y() < 0 || (p.x() + p.y()) > 1) {
        return 0.0f; // Outside the triangle
    }

    // Area of the triangle defined by (0,0), (1,0), (0,1)
    float area = 0.5f; 

    return 1.0f / area;
}


Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float z = 1 - 2 * sample[0];
    float r = sqrt(std::max(0.0f, 1.0f - z * z));    
    float phi = sample[1] * 2 * M_PI; 
    
    // Return Cartesian coordinates on the sphere
    return Vector3f(r * cos(phi), r * sin(phi), z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    
    return  1 / (4 * M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float z = sample[0]; // z = cos(theta)
    float r = sqrt(std::max(0.0f, 1.0f - z * z));
    float phi = 2 * M_PI * sample[1];
    return Vector3f(r * cos(phi), r * sin(phi), z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (v.z() <= 0) return 0.0f; // Below the hemisphere
    return 1 / (2 * M_PI); // PDF for hemisphere
}

Point2f ConcentricSampleDisk(const Point2f &u) {
    // Map uniform random numbers to the [-1, 1] square
    Point2f uOffset = 2.f * u - Point2f(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x() == 0 && uOffset.y() == 0)
        return Point2f(0, 0);

    // Apply concentric mapping to point
    float theta, r;
    if (abs(uOffset.x()) > abs(uOffset.y())) {
        r = uOffset.x();
        theta = M_PI/4 * (uOffset.y() / uOffset.x());
    } else {
        r = uOffset.y();
        theta = M_PI/2 - M_PI/4 * (uOffset.x() / uOffset.y());
    }

    return r * Point2f(std::cos(theta), std::sin(theta));
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f d = ConcentricSampleDisk(sample);
    
    float z = sqrt(std::max(0.0f , 1.0f - d.x() * d.x() - d.y() * d.y()));
    
    return Vector3f(d.x(),d.y(), z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    
    if (v.z() <= 0) {
        return 0.0f; // PDF is zero for directions not in the upper hemisphere
    }
    
    // Return the PDF based on the cosine of the angle with the normal
    return v.z() / M_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float phi = 2 * M_PI * sample.x();
    
    float tanTheta2 = -alpha * alpha * log(1.0f - sample.y()); // tan^2(theta_h)
    float theta = std::atan(std::sqrt(tanTheta2));

    // Uniformly sample the azimuthal angle
    float sinTheta = std::sin(theta);

    // Convert spherical to Cartesian coordinates
    float x = sinTheta * std::cos(phi);
    float y = sinTheta * std::sin(phi);
    float z = std::cos(theta);

    return Vector3f(x, y, z);
}

float Warp::squareToBeckmannPdf(const Vector3f& m, float alpha) {
    if (m.z() <= 0) return 0.0f;

    // Compute theta_h
    float cosTheta = m.z();
    float tanTheta2 = (1.0f - cosTheta * cosTheta) / (cosTheta * cosTheta);

    // Beckmann distribution PDF
    float D = exp(-tanTheta2 / (alpha * alpha)) / (M_PI * alpha * alpha * std::pow(cosTheta, 4));

    // Return the final PDF using p(?h) = D(?h) * cos(?h)
    return D * cosTheta;
}

Point2f Warp::squareToUniformDiskConcentric(const Point2f &sample) {
    return ConcentricSampleDisk(sample);
}

float Warp::squareToUniformDiskConcentricPdf(const Point2f &p) {
    return (p.x() * p.x() + p.y() * p.y() <= 1.0f) ? 1.0f / M_PI : 0.0f;
}




NORI_NAMESPACE_END
