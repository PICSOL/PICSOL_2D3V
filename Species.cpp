#include "Species.h" 

double generateMaxwell(double mu, double sigma)
{
	/* static constant to ensure double number boundary */
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	/* Box-Muller Method */
	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);
	double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	return z0 * sigma + mu;
}

/* generate particle distribution that is propotional to abs(v)*exp(-v^2/2*v_t^2) */
double generateMaxwellFlux( double sigma)
{
	/* static constant to ensure double number boundary */
	static const double epsilon = std::numeric_limits<double>::min();
	
    /* Inverse transform sampling Method */
	double u1;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);
	double z0 = sqrt(-2.0 * log(u1));
	return z0 * sigma ;
}

/* generate particle distribution that is propotional to alpha*abs(v)*(exp(-v^2/2*v_t^2)+ (1-alpha)*v_b*delta(v - v_b)) */
double generateMaxwellBeamFlux(double alpha, double vt, double vb)
{
	/* static constant to ensure double number boundary */
	static const double epsilon = std::numeric_limits<double>::min();

	/* Inverse transform sampling Method for Maxwell Beam system */
	double u;
	do
	{
		u = rand() * (1.0 / RAND_MAX);
	} while (u <= epsilon);

	double c1 = exp(- vb * vb / 2.0);
	double c2 = 1 + alpha / (1 - alpha) * vb;
	double v1 = (1 - c1) / c2;
	double v2 = 1 - c1 / c2;

	
	if (u > v2)
	{
		return vt * sqrt(-2 * log(c2 * (1 - u)));
	}
	else if(u <= v1)
	{
		return vt * sqrt(-2 * log(1 - c2 * u));
	}
	else
	{
		return vt * vb;
	}
}

 Species::Species(double mass, double charge, double tx, double ty, double tz, double dt)
{
	this->mass = mass;
	this->charge = charge;
	this->specific_charge = charge / mass;
	this->tx = tx;
	this->ty = ty;
	this->tz = tz;
	this->dt = dt;
}

void Species::initialMaxwell(Param &param) 
{   
	/* uniformlly distributed in x-y space, maxwell distributed in v space */
	part.reserve(param.INITIAL_PARTICLE * 2);

	for (int i = 1; i <= param.INITIAL_PARTICLE_X; i++)
	{
		double x = param.LBOUND + (param.RBOUND - param.LBOUND) * i / (param.INITIAL_PARTICLE_X + 1);
		for (int j = 1; j <= param.INITIAL_PARTICLE_Y; j++)
		{
			double y = param.LBOUND + (param.UBOUND - param.DBOUND) * i / (param.INITIAL_PARTICLE_Y + 1);
			double vx = generateMaxwell(0, sqrt(tx / mass));
			double vy = generateMaxwell(0, sqrt(ty / mass));
			double vz = generateMaxwell(0, sqrt(tz / mass));
		    
			part.push_back(Particle<double>(x, y, vx, vy, vz));
		}	
	}
}

double Species::kinetic()
{
	double kinetic = 0;
	for (size_t i = 0; i < part.size(); i++)
	{
		kinetic += 0.5 * mass * (part[i].vx * part[i].vx
			                   + part[i].vy * part[i].vy
			                   + part[i].vz * part[i].vz);
	}
	return kinetic;
}

void Species::addParticle(double x, double y,double vx,double vy, double vz)
{
	part.push_back(Particle<double>(x, y, vx, vy, vz));
}

void Species::remParticle(size_t index)
{
	part.erase(part.begin() + index);
}

void Species::saveParticle(string directory)
{
	ofstream output_file(directory.c_str());
	if (!output_file.is_open())
	{
		cerr << " PARTICLE OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	for (size_t i = 0; i < part.size(); i++)
	{
		output_file.setf(ios::fixed);
		output_file.precision(9);
		output_file << part[i].x << '\t' << part[i].y
			<< part[i].vx << '\t' << part[i].vy
			<< '\t' << part[i].vz << endl;
	}
}

void Species::saveGlobal(string directory)
{
	ofstream output_file(directory.c_str());
	if (!output_file.is_open())
	{
		cerr << " PARTICLE OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	output_file << part.size() << '\t' << kinetic() << endl;
}


