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


// // Give a randomly sampled array from a point light
// R3Ray R3SpotLight::
// RandomlySampledRay(void) const
// {
//     R3Point pt_1 = Position();

//     // Generating Randomly uniformly-distributed point of the surface of a sphere
//     // Code referenced from: http://corysimon.github.io/articles/uniformdistn-on-sphere/
//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     std::mt19937 generator (seed);
//     std::uniform_real_distribution<double> uniform01(0.0, 1.0);

//     double theta = 2 * M_PI * uniform01(generator);
//     double phi = acos(1 - 2 * uniform01(generator));
//     double x = sin(phi) * cos(theta);
//     double y = sin(phi) * sin(theta);
//     double z = cos(phi);

//     cout << "Phi: " << phi << endl;
//     cout << "Theta: " << theta << endl;
//     cout << "x: " << x << endl;
//     cout << "y: " << y << endl;
//     cout << "z: " << z << endl;

//     R3Point pt_2 = R3Point(x, y, z);

//     R3Ray ray = R3Ray(pt_1, pt_2);

//     return ray;
// }


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





