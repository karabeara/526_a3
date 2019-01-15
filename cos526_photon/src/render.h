#include <iostream>
using namespace std;

// Include file for the photon map render code

class Photon {

	public:
  		R3Point position;   // Position of the photon
  		R3Ray   direction;  // The direction that the photon came from with the origin of the ray as the "Start" of the ray
  		RNRgb   power;      // Power of the photon, in separate RGB channels

  		// Constructors
  		Photon(void);
  		Photon(R3Point _position, R3Ray _direction, RNRgb _power);

  		// Get Methods
  		R3Point GetPosition();
  		R3Ray   GetDirection();
  		RNRgb   GetPower();

  		// Set Methods
  		void SetPosition(R3Point newPosition);
  		void SetDirection(R3Ray newDirection);
  		void SetPower(RNRgb newPower);

};

R3Kdtree<Photon *> CreatePhotonMap(R3Scene *scene, int photon_count, bool isCausticMap, RNArray<Photon *>* All_Photons);

R3Kdtree<Photon *> CreateVolumePhotonMap(R3Scene *scene, int photon_count, RNArray<Photon *>* All_Photons);

R2Image *RenderImagePhotonMapping(R3Scene *scene, int width, int height, int photon_count, R3Kdtree<Photon *> Global_Photon_Map, R3Kdtree<Photon *> Caustic_Photon_Map, R3Kdtree<Photon *> Volume_Photon_Map);

R2Image *RenderImageRaytracing(R3Scene *scene, int width, int height, int print_verbose);



inline R3Vector 
Reflect(R3Vector rayVector, R3Vector normal)
{
	R3Vector d = R3Vector(rayVector);
    R3Vector n = R3Vector(normal);
    n.Normalize();
    float d_n = d.Dot(n);
    R3Vector r = d -  (2 * d_n * n);

    return r;
}



inline R3Vector 
Refract(const R3Vector &I, R3Vector &N, float ior) 
{ 
    // Refractions calculation reference: https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel

    float cosi = I.Dot(N); // if I.Dot(N) is negative, we are outside the surface; we want cos(theta) to be positive
    float etai = 1, etat = ior; 
    R3Vector n = N; 
    
    if (cosi < 0) { cosi = -cosi; } 
    else { 
    	float temp = etai;
    	etai = etat;
    	etat = temp;

    	n = -N; 
    } 

    float eta = etai / etat; 
    float k = 1 - eta * eta * (1 - cosi * cosi);  // if k < 0, total internal reflection
    
    return R3Vector(k < 0 ? Reflect(I, n) : eta * I + ( eta * cosi - sqrtf(k) ) * n); 
} 
