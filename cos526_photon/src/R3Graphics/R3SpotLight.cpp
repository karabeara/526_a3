/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"
#include<random>
#include<cmath>
#include<chrono>
#include <iostream>
using namespace std;


/* Public variables */

R3SpotLight R3null_spot_light;



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3SpotLight);



/* Public functions */

int 
R3InitSpotLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopSpotLight()
{
}



R3SpotLight::
R3SpotLight(void)
{
}



R3SpotLight::
R3SpotLight(const R3SpotLight& light)
    : R3PointLight(light),
      direction(light.direction),
      dropoffrate(light.dropoffrate),
      cutoffangle(light.cutoffangle)
{
    // Make sure direction is normalized
    this->direction.Normalize();
}



R3SpotLight::
R3SpotLight(const R3Point& position, const R3Vector& direction, const RNRgb& color,
	    RNScalar dropoffrate, RNAngle cutoffangle,
	    RNScalar intensity, RNBoolean active,
            RNScalar ca, RNScalar la, RNScalar qa)
    : R3PointLight(position, color, intensity, active, ca, la, qa),
      direction(direction),
      dropoffrate(dropoffrate),
      cutoffangle(cutoffangle)
{
    // Make sure direction is normalized
    this->direction.Normalize();
}



void R3SpotLight::
SetDirection(const R3Vector& direction)
{
    // Set direction
    this->direction = direction;
    this->direction.Normalize();
}



void R3SpotLight::
SetDropOffRate(RNScalar dropoffrate)
{
    // Set drop off rate
    this->dropoffrate = dropoffrate;
}



void R3SpotLight::
SetCutOffAngle(RNAngle cutoffangle)
{
    // Set cut off angle
    this->cutoffangle = cutoffangle;
}



RNScalar R3SpotLight::
IntensityAtPoint(const R3Point& point) const
{
    // Return intensity at point
    RNScalar I = R3PointLight::IntensityAtPoint(point);
    R3Vector ML = point - Position();
    ML.Normalize();
    RNScalar cos_alpha = ML.Dot(Direction());
    if (cos(cutoffangle) > cos_alpha) return 0.0;
    else return (I * pow(cos_alpha, dropoffrate));
}


// Get ray from light source to given point
const R3Ray R3SpotLight::
LightToPointRay(R3Point point) const
{
    R3Ray ray = R3Ray(R3Point(0, 0, 0), R3Point(1, 1, 1));
    return ray;
}


// Give a randomly sampled array from a point light
const R3Ray R3SpotLight::
RandomlySampledRay(void) const
{
    // Properties of the R3SpotLight
    // R3Vector Direction()
    // RNScalar DropOffRate()
    // RNAngle  CutOffAngle()

    R3Point pt_1 = Position();

    float cos_cut_off_angle = cos( CutOffAngle() );

    // Generating Randomly uniformly-distributed point of the surface of a sphere
    // Code referenced from: http://corysimon.github.io/articles/uniformdistn-on-sphere/
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uniform01(0, 1.0);

    float spotCosine = cos_cut_off_angle;
    R3Ray ray = R3Ray(R3Point(0, 0, 0), R3Point(1, 1, 1));

    while (spotCosine <= cos_cut_off_angle)
    {
        double theta = 2 * M_PI * uniform01(generator);
        double phi = acos(1 - 2 * uniform01(generator));
        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);

        R3Point pt_2 = R3Point(x + pt_1.X(), 
                           y + pt_1.Y(), 
                           z + pt_1.Z());

        ray = R3Ray(pt_1, pt_2);

        spotCosine = Direction().Dot( ray.Vector() );
    }


    // L = normalize( light.position.xyz/light.position.w - v_eyeCoords );
    // if (light.spotCosineCutoff > 0.0) { // the light is a spotlight
    // R3Vector D = normalize( Direction() );  // unit vector!
    // float spotCosine = D.Dot(-L);
    // if (spotCosine >= light.spotCosineCutoff) { 
    //     spotFactor = pow(spotCosine,light.spotExponent);
    // }
    // else { // The point is outside the cone of light from the spotlight.
    //     return vec3(0.0); // The light will add no color to the point.
    // }
    // // Light intensity will be multiplied by spotFactor

    return ray;
}


const RNRgb R3SpotLight::
PowerGivenDistance(R3Point reference_point, R3Vector normal) const
{
    // RNRgb resulting_power; 

    // double x = reference_point.X() - Position().X();
    // double y = reference_point.Y() - Position().Y();
    // double z = reference_point.Z() - Position().Z();

    // // Distance between the two points
    // double d = sqrt( x*x + y*y + z*z );

    // double ca = ConstantAttenuation();
    // double la = LinearAttenuation();
    // double qa = QuadraticAttenuation();

    // // Calculate: (r,g,b) * 1.0/ (ca + la*d + qa*d*d)
    // resulting_power = Color() * ( 1.0 / ( ca + (la*d) + (qa*d*d) ) );

    // return resulting_power;

    return Color() * IntensityAtPoint(reference_point);

    //return Color();
}




void R3SpotLight::
Draw(int i) const
{
    // Draw light
    GLenum index = (GLenum) (GL_LIGHT2 + i);
    if (index > GL_LIGHT7) return;
    GLfloat buffer[4];
    buffer[0] = Intensity() * Color().R();
    buffer[1] = Intensity() * Color().G();
    buffer[2] = Intensity() * Color().B();
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);
    buffer[0] = Position().X();
    buffer[1] = Position().Y();
    buffer[2] = Position().Z();
    buffer[3] = 1.0;
    glLightfv(index, GL_POSITION, buffer);
    buffer[0] = Direction().X();
    buffer[1] = Direction().Y();
    buffer[2] = Direction().Z();
    glLightfv(index, GL_SPOT_DIRECTION, buffer);
    buffer[0] = DropOffRate();
    glLightf(index, GL_SPOT_EXPONENT, buffer[0]);
    buffer[0] = RN_RAD2DEG(CutOffAngle());
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    buffer[0] = ConstantAttenuation();
    buffer[1] = LinearAttenuation();
    buffer[2] = QuadraticAttenuation();
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);
    glEnable(index);
}





