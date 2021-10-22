/***************************************************************************************************
****************************************************************************************************

    genC_lammps.cc - Contains the code for generating a hexagonal graphite Carbon lattice in LAMMPS
                     file format.  These files can be imported easily into LAMMPS.

****************************************************************************************************
***************************************************************************************************/
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "particle.h"
#include <iomanip>
using namespace std;

void readInput(double &lattice_a,
               double &lattice_c,
               int &nx,
               int &ny,
               int &nz,
               int &nz_s);

int main()
{
    ofstream output_file; // file to which the particle data will be written
   ofstream vertex_file; // file to which lattice target area data will be written

    // define some constants
    const double m_norm = 6.023e26 * 1.66057788e-27; // avogadro * amu
    const double C_mass = 12.0 / m_norm; // mass in amu
   

    double pos_temp[3]; // intermediate variable in which positions are stored.
    double pos[3];
    pos[0] = 0; pos[1] = 0; pos[2] = 0;


    vector<particle> particle_data;
    int nx, ny, nz, nz_s; // number of particles in each direction
    double lattice_a, lattice_c; // C lattice constants
    
    int surf_layers = 1;
    
    readInput(lattice_a, lattice_c, nx, ny, nz, nz_s);
           
	cout <<"nx "<< nx << endl;
	cout <<"ny "<< ny << endl;
	cout <<"nz "<< nz << endl;
	cout <<"lattice_a "<< lattice_a << endl;
	cout <<"lattice_c "<< lattice_c << endl;

	/*string name;
	std::stringstream convert;
	convert << lattice_a;
	name = convert.str();
	name = "graphite.pos";
	cout << "name of the output file is "<<name<<endl;*/
 
    // determine the unit lengths in each directions
    double u_len[3];
    u_len[0] = 3.0 * lattice_a;
    u_len[1] = sqrt(3.0) * lattice_a;
    u_len[2] = lattice_c;
    
    // Extra height in the z direction to make it a surface
    double extra_height = 20.0;
    cout << "The extra height is " << extra_height << endl;


    // determine the dimensions of the simulation box
    double box[6];
    box[0] = 0;
    box[1] = u_len[0] * nx;

    box[2] = 0;
    box[3] = u_len[1] * ny;

    box[4] = -u_len[2] * (nz-0.5) - extra_height;
    box[5] = u_len[2] * 0.5 + 1.0*extra_height;
    
    unsigned long n_atoms = 4 * nx * ny * nz;

    // output some statistics to the command line
    cout << "N:\t" << n_atoms << endl;
    cout << "unit length x: " << u_len[0] << "\ty: " << u_len[1] << "\tz: " << u_len[2] << endl;
    cout << "box x: " << box[0] << "\ty: " << box[1] << "\tz: " << box[2] << endl;

    

    // add an inert particle to the system
    unsigned long id_num = 0;

    cout<<"Lattice_c = "<<lattice_c<<endl;
    
   for (int i = 0; i< nx; i++) 
     {
       for (int j = 0; j < ny; ++j)
       {
         for (int k = 0; k < nz; k++)
	  	 {  	 	
	  	 	
	  	 	int height_mod = k%2;
	  	 	// basic particles
	  	 	pos[0] = (0.75 + 3.0*i + height_mod)*lattice_a;
	  	 	pos[1] = (0.25*sqrt(3.0) + sqrt(3.0)*j)*lattice_a;
			//	  	 	pos[2] = -(0.5 + k)*lattice_c;
			pos[2] = 1*(k - (nz-1) )*lattice_c;
			//			cout<<pos[2]<<"  "<<1*(k - (4-1) )<<endl;
	  	        
			particle_data.push_back(particle(to_string(k+1), C_mass, 0, id_num++));
			particle_data.back().setPosition(pos[0], pos[1], pos[2]);
    		
    		// extending in 110
	  	 	pos[0] = (1.25 + 3.0*i + height_mod)*lattice_a;
	  	 	pos[1] = (0.25*3.0*sqrt(3.0) + sqrt(3.0)*j)*lattice_a;
	  	 	//pos[3] = (0.5 + k)*lattice_c;
			pos[2] = (k - nz + 1) *lattice_c;

			particle_data.push_back(particle(to_string(k+1), C_mass, 0, id_num++));

			particle_data.back().setPosition(pos[0], pos[1], pos[2]);
    		
    		// extending in 101
	  	 	pos[0] = (2.25 + 3.0*i - 2.0*height_mod)*lattice_a;
	  	 	pos[1] = (0.25*3.0*sqrt(3.0) + sqrt(3.0)*j)*lattice_a;
	  	 	//pos[3] = (0.5 + k)*lattice_c; 
			pos[2] = (k - nz + 1) *lattice_c;

			particle_data.push_back(particle(to_string(k+1), C_mass, 0, id_num++));				
			particle_data.back().setPosition(pos[0], pos[1], pos[2]);
    		
    		// extending in 011
	  	 	pos[0] = (2.75 + 3.0*i - 2.0*height_mod)*lattice_a;
	  	 	pos[1] = (0.25*sqrt(3.0) + sqrt(3.0)*j)*lattice_a;
			//	  	 	pos[3] = (0.5 + k)*lattice_c;	 
			pos[2] = (k - nz + 1) *lattice_c;
			particle_data.push_back(particle(to_string(k+1), C_mass, 0, id_num++));				
			particle_data.back().setPosition(pos[0], pos[1], pos[2]);

    		
    		
	 	 }
	 
       }
      
     }
      
    // Flagging the surface
    cout <<"Number of surface layers = "<< nz_s << endl;
    //    double z_max = (0.5 + nz-nz_s) * lattice_c;
    double z_max = (1-nz_s) * lattice_c;
    double z_min = (1-nz) * lattice_c;
    double z_mid = (0.5 + 0.5*(nz-1)) * lattice_c;
    double skin[3] = {1.0,1.0,2.0};
    int count_surf = 0;
    
    cout<<"z_max = "<<z_max<<endl;
    
    int n_atoms_final = particle_data.size();
    
    for (int i = 0; i<n_atoms_final; i++ )
    {
		Vector3D position = particle_data[i].getPosition();
		pos[0] = position.x;
		pos[1] = position.y;
		pos[2] = position.z;
		
		if(position.z >= z_max )
		{
			particle_data[i].setSpecies("1");
			count_surf++;
		}
		
		if(position.z == z_min)
		{
			particle_data[i].setSpecies("1");
		}
		
		/*if(position.z > z_max)
		{
			particle_data[i].setSpecies("2");
			count_surf++;
		}
		
		if(position.z < z_min)
		{
			particle_data[i].setSpecies("2");
			count_surf++;
		}*/
		
		/*val = 0;
		for (int j=0; j<3; j++) val += pow((pos[j]-center[j])/(axis_len[j]+skin[j]),2);	
		if(val <= 1)
		{
			particle_data[i].setSpecies("2");
			count_surf++;
		}*/
	}
	
	cout<<"Number of surface atoms = "<<count_surf<<endl;
    
    // open the output file for writing
    output_file.open("graphite_slab.lmp");

    // output the file headers
    output_file << "# Hexagonal lattice structure of dimension " << nx << " x " << ny << " x " << nz << " with a = " << lattice_a << " and c = " << lattice_c << endl;
    output_file << "\t\t\t\t " << n_atoms << " atoms" << endl;
    output_file << "\t\t\t\t "<<nz<<" atom types" << endl << endl;

    // output the box dimensions
    //    output_file << "<box lx=\"" << box[0] << "\" ly=\"" << box[1] << "\" lz=\"" << 2*box[2] << "\"/>" << endl;
    output_file <<box[0]<<" \t " << box[1] << "\t xlo xhi" << endl;    
    output_file <<box[2]<<" \t " << box[3] << "\t ylo yhi" << endl; 
    output_file <<box[4]<<" \t " << box[5] << "\t zlo zhi" << endl; 
    // output the particle species
    output_file << endl << "Atoms" << endl << endl;
    output_file << std::setprecision(7) <<fixed;
    for (unsigned long i = 0; i < particle_data.size(); ++i)
    {
        Vector3D position = particle_data[i].getPosition();
        string species = particle_data[i].getSpecies();
        double charge = particle_data[i].getCharge();
        output_file << i+1 << "\t" << species << "\t" << charge << "\t" << position.x << "\t" << position.y << "\t" << position.z << endl;
    }


    // close the file
    output_file.close();
      
}

void readInput(double &lattice_a,
               double &lattice_c,
               int &nx,
               int &ny,
               int &nz,
               int &nz_s)
{
    ifstream input_file;
    const int num_char = 256;

    // read the input file
    input_file.open("input_lammps.dat");

 	input_file.ignore(num_char, '\n');
 	input_file.ignore(num_char, '\n');
    input_file >> lattice_a;
    input_file.ignore(num_char, '\n');
    input_file >> lattice_c;
    input_file.ignore(num_char, '\n');
    input_file >> nx;
    input_file.ignore(num_char, '\n');
    input_file >> ny;
    input_file.ignore(num_char, '\n');
    input_file >> nz;
    input_file.ignore(num_char, '\n');
    input_file >> nz_s;

    if ( nx <=0 || ny <= 0 || nz<=0 || nz_s<0 )
      {
	cout<<"!!!!!!! Error : Negative # of unit cells not allowed !!!!!!!!!!!!"<<endl;
	exit(0);
      }

    input_file.close();
}

