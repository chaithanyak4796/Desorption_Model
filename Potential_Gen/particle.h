/***************************************************************************************************
****************************************************************************************************

	particle.h - contains the definition for the class "particle"

****************************************************************************************************
***************************************************************************************************/
#include <string>
#include <cmath>
#include "Vector3D.h"

class particle
{
	public:
        particle(std::string spec, double m, double q, unsigned long id); // constructor
        unsigned long getID(); // gets the particle's id number
		void setSpecies(std::string spec); // sets the particle species
		std::string getSpecies(); // gets the particle species
		void setMass(double m); // sets the mass of the particle
		double getMass(); // gets the mass of the particle
		void setCharge(double q); // sets the charge of the particle
		double getCharge(); // gets the charge of the particle
		void computeGamma(); // sets the Lorentz factor of the particle
		double getGamma(); // gets the value of gamma
		void setPosition(double, double, double); // sets the position of the particle
		Vector3D getPosition(); // gets the position of the particle
		void computeMomentum(); // computes the momentum from the velocity, mass, and Lorentz factor
		Vector3D getMomentum(); // gets the momentum of the particle
		void setVelocity(double, double, double); // sets the velocity of the particle
		Vector3D getVelocity(); // gets the velocity of the particle
		void computeEnergy(); // computes the kinetic energy of the particle
		double getEnergy(); // gets the particle's kinetic energy
		void setForce(double, double, double); // sets the force on the particle
		Vector3D getForce(); // gets the force on the particle
		void setBody( int ); // sets the body index of the particle
		int getBody(); // gets the body index of the particle

    private:
        std::string species; // particle species
        unsigned long id_number; // particle id number
	int body; // particle body number
        double mass; // particle mass
        double charge; // particle charge
        double gamma; // particle Lorentz factor
        double kinetic_energy; // particle kinetic energy
        Vector3D position; // particle position
        Vector3D momentum; // particle momentum
        Vector3D velocity; // particle velocity
        Vector3D force; // net force on particle
};

// particle constructor - creates a particle with mass m and charge q
particle::particle(std::string spec, double m, double q, unsigned long id)
{
    species = spec;
	mass = m;
	charge = q;
	id_number = id;
    gamma = 1.0;
}

// getID - gets the id number of the particle
unsigned long particle::getID()
{
    return id_number;
}

// setSpecies - sets the species of the particle
void particle::setSpecies(std::string spec)
{
    species = spec;
}

// getSpecies - gets the species of the particle
std::string particle::getSpecies()
{
    return species;
}

//setMass - sets the mass of the particle
void particle::setMass(double m)
{
    mass = m;
}

// getMass - gets the mass of the particle
double particle::getMass()
{
    return mass;
}

// setCharge - sets the charge of the particle
void particle::setCharge(double q)
{
    charge = q;
}

// getCharge - gets the charge of the particle
double particle::getCharge()
{
    return charge;
}

// computeGamma - computes the Lorentz factor of the particle
void particle::computeGamma()
{
    double velx = velocity.x;
    double vely = velocity.y;
    double velz = velocity.z;
    double denom = sqrt(1 - (velz * velz + vely * vely + velz * velz));

    gamma = 1.0 / denom;
}

// getGamma - gets the particle's Lorentz factor
double particle::getGamma()
{
    return gamma;
}

// setPosition - sets the position of the particle
void particle::setPosition(double posx, double posy, double posz)
{
    position.x = posx;
    position.y = posy;
    position.z = posz;
}

// getPosition - gets the position of the particle
Vector3D particle::getPosition()
{
    return position;
}

// computeMomentum - computes the particle momentum from the velocity, mass, and Lorentz factor
void particle::computeMomentum()
{
    // compute gamma in case the current value is incorrect
    computeGamma();

    // compute the momentum components
    momentum.x = mass * gamma * velocity.x;
    momentum.y = mass * gamma * velocity.y;
    momentum.z = mass * gamma * velocity.z;
}

// getMomentum - gets the momentum of the particle
Vector3D particle::getMomentum()
{
    return momentum;
}

// setVelocity - sets the velocity of the particle
void particle::setVelocity(double velx, double vely, double velz)
{
    velocity.x = velx;
    velocity.y = vely;
    velocity.z = velz;
}

// getVelocity - gets the velocity of the particle and returns it in the form of a 3D vector
Vector3D particle::getVelocity()
{
    return velocity;
}

// computeEnergy - computes the particle kinetic energy from the momentum and velocity
void particle::computeEnergy()
{
    // compute gamma in case the current value is not correct
    computeGamma();

    // compute the kinetic energy
    kinetic_energy = mass * (gamma - 1);
}

// getEnergy - gets the particle's kinetic energy
double particle::getEnergy()
{
    return kinetic_energy;
}

// setForce - sets the force on the particle
void particle::setForce(double fx, double fy, double fz)
{
    force.x = fx;
    force.y = fy;
    force.z = fz;
}

// getForce - gets the force on the particle
Vector3D particle::getForce()
{
    return force;
}

// setBody - sets the body index of the particle
void particle::setBody( int ii )
{
	body = ii;
}

// getBody - gets the body index of the particle
int particle::getBody()
{
	return body;
}
