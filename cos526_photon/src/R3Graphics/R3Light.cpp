/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"
#include<random>
#include<cmath>
#include<chrono>
#include <iostream>
using namespace std;


/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Light);



/* Public variables */

RNScalar R3ambient_light_intensity = 0.2;
RNRgb R3ambient_light_color(1.0, 1.0, 1.0);



/* Public functions */

int 
R3InitLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopLight()
{
}



R3Light::
R3Light(const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(1),
      intensity(1),
      color(0.5, 0.5, 0.5),
      id(0)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
R3Light(const R3Light& light, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(light.active),
      intensity(light.intensity),
      color(light.color),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
R3Light(const RNRgb& color,
	RNScalar intensity,
	RNBoolean active,
        const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      active(active),
      intensity(intensity),
      color(color),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Light::
~R3Light(void)
{
    // Unload
    if (id > 0) { fprintf(stderr, "Will not get here because id is not used yet.\n"); }

    // Remove light from scene
    if (scene) scene->RemoveLight(this);

    // Free name
    if (name) free(name);
}



void R3Light::
SetName(const char *name)
{
    // Set name
    if (this->name) free(this->name);
    if (name) this->name = strdup(name);
    else this->name = NULL;
}



void R3Light::
SetActive(RNBoolean active)
{
    // Set active
    this->active = active;
}



void R3Light::
SetIntensity(RNScalar intensity)
{
    // Set intensity
    this->intensity = intensity;
}



void R3Light::
SetColor(const RNRgb& color)
{
    // Set color
    this->color = color;
}

// Get ray from light source to given point
const R3Ray R3Light::
LightToPointRay(R3Point point) const 
{
  return R3Ray( R3Point(0, 0, 0), R3Point(1, 1, 1) );
}


const R3Ray R3Light::
RandomlySampledRay(void) const 
{
  return R3Ray( R3Point(0, 0, 0), R3Point(1, 1, 1) );
}


// Give the power of a surface photon given the distance d from the light source
const RNRgb R3Light::
PowerGivenDistance(R3Point reference_point, R3Vector normal) const
{
  cout << "LIGHT default" << endl;
  return RNRgb(0, 0, 0);
}


