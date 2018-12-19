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

// Private variables used by internal functions
//RNArray<Photon *>* All_Photons;
//int photon_count;

R3Kdtree<Photon *> EmitPhotons(R3Scene *scene, int photon_count, RNArray<Photon *>* All_Photons);

R3Kdtree<Photon *> ScatterPhotons(R3Scene *scene, int photon_count, RNArray<Photon *>* All_Photons, R3Kdtree<Photon *> Photon_Map);

RNRgb EstimateRadiance(R3Point reference_point, int closest_points_count, RNArray<Photon *> Closest_Photons, R3Kdtree<Photon*> Photon_Map);

R2Image *RenderImagePhotonMapping(R3Scene *scene, int width, int height, int photon_count, R3Kdtree<Photon *> Photon_Map);

R2Image *RenderImageRaytracing(R3Scene *scene, int width, int height, int print_verbose);

