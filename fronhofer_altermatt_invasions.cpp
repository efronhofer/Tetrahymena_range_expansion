//============================================================================
// Name        : Eco-evolutionary feedbacks during experimental range expansions
// Author      : Emanuel A. Fronhofer and Florian Altermatt
// Date	       : February 2015
//============================================================================

/*
	Copyright (C) 2015  Emanuel A. Fronhofer

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

#include <iostream>
#include <cstdlib>								//standard C library
#include <ctime>								//access system time library
#include <fstream>								//file streaming library
#include <string>								//string library included
#include <sstream>								//string streaming for reading numbers from

#include <vector>
#include <cmath>								//standard math library

#include <gsl/gsl_rng.h>						//gsl random number generator
#include <gsl/gsl_randist.h>					//gsl random distributions
#include <gsl/gsl_statistics.h>					//gsl some statistical methods
#include <gsl/gsl_statistics_double.h> 			//gsl some double statistical methods
#include <gsl/gsl_sort_double.h> 				//gsl sort double arrays

#include <algorithm>

using namespace std;

#include "include/procedures.h"					//procedure simplifications
#include "include/classes.h"					//class definitions

//_____________________________________________________________________________
//------------------------------------------------------------ global variables
unsigned int sim_time;															// actual time in simulation
unsigned int burn_in;															// actual time in simulation
int max_runs;																	// no. of repeats
float mut_sd;																	// variance used for mutations
float mut_rate;																	// mutation rate

float max_prey;																	// habitat capacity
float lambda_null;																// fertility

float mu0;																		// dispersal mortality
float epsilon;																	// random patch extinctions

TPatch world[WORLDDIM_X];														// simulated world

float abs_metapopsize;															// absolute metapopulation size
float occupancy;																// metapopulation occupancy
float rel_emigrants_core;														// relative number of emigrants in the core
float rel_emigrants_margin;														// relative number of emigrants in the core
unsigned int metapopsize_core;													// relative number of emigrants in the core
unsigned int metapopsize_margin;												// relative number of emigrants at the margin
unsigned int pocc_core;															// number of occupied patches in the core
unsigned int pocc_margin;														// number of of occupied patches at the margin
int margin_position;															// position of the range margin
double constant_resource_addition;												// constant resource addition
float lambda0_resources;														// resource growth rate
float capacity_resources;														// resource carrying capacity
float foraging_efficiency_max;													// foraging efficiency

float tot_resource_size;														// resource pop size
int resource_occupancy;
float maxDisp;

//_____________________________________________________________________________
//------------------------------------------------------------------ procedures

//------------------------------------------------------------- read parameters
void readParameters(){
	ifstream parinfile("input/parameters.in");							//parameter input file
	string buffer;
	istringstream is;

	getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sim_time;																		//simulation time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> burn_in;																		//burn in time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_runs;																		//no. of repeats
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_sd;																		//variance used for mutations
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate;																		//mutation rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_prey;																		//maximal amount of prey
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> lambda_null;																	//fertility
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> foraging_efficiency_max;														//max foraging efficiency for trade-off function
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> epsilon;																		//random patch extinction probability
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mu0;																			//dispersal mortality
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> constant_resource_addition;													//constant resource addition
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> lambda0_resources;															//resource growth rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> capacity_resources;															//resource carrying capacity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> maxDisp;																		//maximla dispersal value

	parinfile.close();
}

//------------------------------------------------------- calculate foraging efficiency
float foragingMaxTradeOffFunction(float pForMax){
	return(pForMax * max_prey);
}

//------------------------------------------------------- calculate foraging efficiency
float foragingTradeOffFunction(float pFor){
	return(pFor * foraging_efficiency_max);
}

//--------------------------------------------------------------- calculate dispersal rate
float dispersalTradeOffFunction(float pDisp){
	return(pDisp * maxDisp);
}

//------------------------------------------------------- initialize simulation
void Initialize(){
	// initialize patches and individuals in patches
	for (int x = 0; x < WORLDDIM_X; ++x) {

		// clear the world
		world[x].females.clear();
		world[x].newFemales.clear();

		// initialize only in the initialization area
		if (x < INIT_ROWS) {
			// initialize individuals in this patch

			int capacity = round(capacity_resources/max_prey);

			for (int f = 0; f < capacity; ++f) {
				TIndiv newfemale;
				newfemale.pDisp = ran();
				newfemale.pFor = ran();
				newfemale.pForMax = ran();
				newfemale.pElse = ran();

				// renormalize traits
				float traitsum = newfemale.pDisp + newfemale.pElse + newfemale.pFor + newfemale.pForMax ;
				newfemale.pDisp = newfemale.pDisp / traitsum;
				newfemale.pElse = newfemale.pElse / traitsum;
				newfemale.pFor = newfemale.pFor / traitsum;
				newfemale.pForMax = newfemale.pForMax / traitsum;

				world[x].females.push_back(newfemale);
			}
		}

		// initialize resources
		world[x].resources = capacity_resources;
	}
}

// ------------------------------------------------ analyze population dynamics
void Analyze(unsigned int acttime, int actrun, int actworlddim_x){
	//reset metapopulation size and occupancy
	occupancy = 0;
	// reset margin position
	margin_position = 0;

	// reset respurce pop size counter
	tot_resource_size = 0;
	resource_occupancy = 0;

	unsigned int numberoccupied = 0;
	unsigned int metapopsize = 0;

	for (int x = 0; x < actworlddim_x; ++x) {

		tot_resource_size = tot_resource_size + world[x].resources;
		if(world[x].resources > 0){++resource_occupancy;}

		unsigned int localpopsize = world[x].females.size();
		metapopsize += localpopsize;
		if (localpopsize > 0) {
			++numberoccupied;
			if (margin_position < x) {
				margin_position = x;
			}
		}
	}
	// calculate occupancy
	occupancy = float(numberoccupied) / float(WORLDDIM_X);
	abs_metapopsize = metapopsize;
}

// ---------------------------------------------------- save individual results
void saveResults(int actrun, int acttime){
	// output file: individuals
	stringstream outputindiv_path_stream;
	outputindiv_path_stream << "output/output_spatial_profile_run" << actrun << ".out";
	string outputindiv_path = outputindiv_path_stream.str();
	ofstream outputindiv(outputindiv_path.c_str());

	// headers
	outputindiv << "x" << "    " << "pop_size" << "     "  << "emi_rate" << "    " << "for_eff" << "    " << "max_prey" << "    " << "fert" << "    " << "pElse" << "    " << endl;

	for (int x = 0; x < WORLDDIM_X; ++x) {

		// calculate mean values for traits for this x location
		float mean_disp = 0;
		float mean_for = 0;
		float mean_forMax = 0;
		float mean_pElse = 0;

		for (int f = 0; f < world[x].females.size(); ++f) {
			// calculate sums
			mean_disp = mean_disp + dispersalTradeOffFunction(world[x].females.at(f).pDisp);
			mean_for = mean_for + foragingTradeOffFunction(world[x].females.at(f).pFor);
			mean_forMax = mean_forMax + foragingMaxTradeOffFunction(world[x].females.at(f).pForMax);
			mean_pElse = mean_pElse + world[x].females.at(f).pElse;
		}

		// calculate means
		mean_disp = mean_disp/float(world[x].females.size());
		mean_for = mean_for/float(world[x].females.size());
		mean_forMax = mean_forMax/float(world[x].females.size());
		mean_pElse = mean_pElse/float(world[x].females.size());

		outputindiv << x << "     " << world[x].females.size() << "     " << mean_disp << "     " << mean_for << "     " << mean_forMax << "     " << lambda_null << "     " << mean_pElse << endl;;
	}

	// close indiv output file
	outputindiv.close();
}

// --------------------------------------------- find new patch during dispersal
int findNewPatch(int x, int act_worlddim_x){

	int res;
	int dir;

	// nearest neighbour dispersal (nnd8)
	dir = floor(ran()*2);

	switch (dir) {
	case 0:
		res = x-1;
		break;
	case 1:
		res = x+1;
		break;
	default:
		cout << "Error in NND" << endl;
		break;
	}

	// behaviour at world limits
	if (act_worlddim_x == INIT_ROWS) {
		// torus in x direction only during burnin
		if (res == act_worlddim_x){
			res = 0;
		}
		if (res == -1){
			res = act_worlddim_x-1;
		}
	} else {
		// otherwise reset to worldlimit
		if (res == act_worlddim_x){
			res = act_worlddim_x-1;
		}
		if (res == -1){
			res = 0;
		}
	}

	return(res);
}

// -------------------------------------------------------- dispersal procedure
void Dispersal(int act_worlddim_x){

	// core
	unsigned int no_emigrants_core = 0;
	metapopsize_core = 0;
	pocc_core = 0;
	rel_emigrants_core = 0;
	//margin
	unsigned int no_emigrants_margin = 0;
	metapopsize_margin = 0;
	pocc_margin = 0;
	rel_emigrants_margin = 0;

	for (int x = 0; x < act_worlddim_x; ++x) {
		// counter for metapopsize in core
		if (x < CM_AREA){
			metapopsize_core += world[x].females.size();
			pocc_core += 1;
		}
		if (x >= (margin_position - (CM_AREA-1))){
			metapopsize_margin += world[x].females.size();
			if (world[x].females.size() > 0){ pocc_margin += 1;}
		}

		// start with females
		for (unsigned int f = 0; f < world[x].females.size(); ++f) {
			// should the individual disperse?
			if (ran() < dispersalTradeOffFunction(world[x].females.at(f).pDisp)){
				// this individual will disperse
				// increase counter
				if (x < CM_AREA){++no_emigrants_core;}
				if (x >= (margin_position - (CM_AREA-1))){++no_emigrants_margin;}
				// check whether this emigrant survives the dispersal process
				if (ran() > mu0){
					// find new patch (global dispersal)
					int newPatch = findNewPatch(x, act_worlddim_x);
					// copy disperser into new patch
					TIndiv Disperser = world[x].females.at(f);
					world[newPatch].newFemales.push_back(Disperser);
				}
				// delete emigrant from natal patch
				world[x].females.at(f) = world[x].females.back();
				world[x].females.pop_back();
				// decrease female loop counter
				--f;
			}
		}
	}

	// now that dispersal is over, merge philopatrics and residents
	for (int x = 0; x < act_worlddim_x; ++x) {
		// first copy the females
		for (unsigned int f = 0; f < world[x].newFemales.size(); ++f) {
			world[x].females.push_back(world[x].newFemales.at(f));
		}
		// erase the "old" immigrants from newFemales
		world[x].newFemales.clear();

	}

	rel_emigrants_core = float(no_emigrants_core) / float(metapopsize_core);
	rel_emigrants_margin = float(no_emigrants_margin) / float(metapopsize_margin);

}

// ------------------------------------------------------------------ mutations
float mutate(float allele){
	if(ran()< mut_rate){
		float newallele = allele + gauss(mut_sd);

		// boundary conditions
		if (newallele < 0) {
			newallele = 0;
		}

		return(newallele);

	} else {
		return(allele);
	}
}

// ------------------------------------------------------------ foraging
float foraging(int x, unsigned int f){

	float mean_offspring;

	if (world[x].resources > 0){
		// calculated prey items foraged according to Holling Type II functional response
		float prey_foraged = foragingMaxTradeOffFunction(world[x].females.at(f).pForMax) * world[x].resources / (1/foragingTradeOffFunction(world[x].females.at(f).pFor) + world[x].resources);

		// reduce prey population size
		world[x].resources = world[x].resources - prey_foraged;
		if(world[x].resources < 0){
			cout << "ERROR: negative resources! (in foraging function)" << endl;
			world[x].resources = 0;
		}

		mean_offspring = lambda_null * prey_foraged;

	}else{
		mean_offspring = 0;
		world[x].resources = 0;
	}

	return(mean_offspring);
}


// --------------------------------------------------------------- reproduction
void Reproduction(int act_worlddim_x){
	for (int x = 0; x < act_worlddim_x; ++x) {
		// just to be sure: resize new females and males vectors
		world[x].newFemales.clear();

		// for each patch check whether there are females and males
		if (world[x].females.size() > 0 ) {

			// shuffle the females address vector
			int act_females[world[x].females.size()];
			for (unsigned int f = 0; f < world[x].females.size(); ++f) {
				act_females[f] = f;
			}
			shuffle_int(act_females, world[x].females.size());

			// females choose their mates
			for (unsigned int fcnt = 0; fcnt < world[x].females.size(); ++fcnt) {
				// go through the females vector randomly
				int f = act_females[fcnt];
				// calculate number of offspring
				float lambda_indiv = foraging(x,f);
				int no_offspring = poisson(lambda_indiv);

				// loop over offspring
				for (int o = 0; o < no_offspring; ++o) {
					// females
					// initialize new individual
					TIndiv newOffspring;
					newOffspring.pDisp = mutate(world[x].females.at(f).pDisp);
					newOffspring.pElse = mutate(world[x].females.at(f).pElse);
					newOffspring.pFor = mutate(world[x].females.at(f).pFor);
					newOffspring.pForMax = mutate(world[x].females.at(f).pForMax);

					// renormalize traits
					float traitsum = newOffspring.pDisp + newOffspring.pElse + newOffspring.pFor + newOffspring.pForMax;
					newOffspring.pDisp = newOffspring.pDisp / traitsum;
					newOffspring.pElse = newOffspring.pElse / traitsum;
					newOffspring.pFor = newOffspring.pFor / traitsum;
					newOffspring.pForMax = newOffspring.pForMax / traitsum;

					// add new individual to new females vector
					world[x].newFemales.push_back(newOffspring);

				}
			}
		}
	}
}

// --------------------------------------------------------------- resource reproduction
void Resource_Reproduction(int act_worlddim_x){
	for (int x = 0; x < act_worlddim_x; ++x) {
		// make sure that resource value is positive
		if (world[x].resources < 0) {
			cout << "Error: resources negative! " << endl;
			world[x].resources = 0;
		} else {
			// add constant amount of resources
			world[x].resources = world[x].resources + constant_resource_addition*capacity_resources;
			// make sure that this does not get over one
			if (world[x].resources > capacity_resources){world[x].resources = capacity_resources;}
		}

		// for each patch check whether there are resources
		if (world[x].resources > 0 ) {
			// logistic growth as defined by Beverton & Holt
			world[x].resources = world[x].resources * lambda0_resources * 1/(1+(lambda0_resources-1)/float(capacity_resources)*world[x].resources);
		}
	}
}

// -------------------------------------------------- death of annual organisms
void Death(int act_worlddim_x){
	for (int x = 0; x < act_worlddim_x; ++x) {
		int local_offspring_no = world[x].newFemales.size();
		//cout << local_offspring_no << endl;
		if (local_offspring_no > 0) {
			// now clear adult vectors
			world[x].females.clear();
			// include local patch extinctions
			if (ran() > epsilon){
				// now copy new females into females
				for (unsigned int nf = 0; nf < world[x].newFemales.size(); ++nf) {
					world[x].females.push_back(world[x].newFemales.at(nf));
				}
			}
			// clear new females vector
			world[x].newFemales.clear();
		} else {
			world[x].females.clear();
		}
	}
}

//_____________________________________________________________________________
//------------------------------------------------------------------------ main

int main() {
	// random number generator
	//specify_rng(time(NULL));
	specify_rng(RS);

	//read parameters for all simulation runs
	readParameters();

	// repeat loop
	for (int actrun = 0; actrun < max_runs; ++actrun) {

		// output file: metapop
		stringstream outputmetapop_path_stream;
		outputmetapop_path_stream << "output/output_metapop_run" << actrun << ".out";
		string outputmetapop_path = outputmetapop_path_stream.str();
		ofstream outputmetapop(outputmetapop_path.c_str());

		// outputfile header
		outputmetapop << "time" << "    " << "abs_metapopsize" << "    " << "occupancy"<< "    " << "emirate_core" << "    " << "emirate_margin"<< "    " << "popsize_core" << "    " << "popsize_margin" << "    " << "margin_position" << endl;

		Initialize();

		// time loop
		for (unsigned int acttime = 0; acttime < sim_time; ++acttime) {

			// for burnin phase: restrict world dimension to this area (this is only applicatbe to x dim, as y dim is not changed)
			int act_worlddim_x = 0;

			if (acttime < burn_in) {
				act_worlddim_x = INIT_ROWS;
			} else {
				act_worlddim_x = WORLDDIM_X;
			}

			// natal dispersal
			Dispersal(act_worlddim_x);

			// first reproduction of resources
			Resource_Reproduction(act_worlddim_x);

			// reproduction
			Reproduction(act_worlddim_x);

			// density regulation and death of adults
			Death(act_worlddim_x);

			// analyze metapopulation
			Analyze(acttime, actrun, act_worlddim_x);

			// write metapop results to file
			outputmetapop << acttime << "    " << abs_metapopsize << "    " << occupancy << "    " << rel_emigrants_core << "   " << rel_emigrants_margin << "   " << float(metapopsize_core)/float(pocc_core) << "   " << float(metapopsize_margin)/float(pocc_margin) << "   " << margin_position << endl;
			// write run and time to console
			cout << actrun << "   " << acttime << "   " << rel_emigrants_core << "   " << rel_emigrants_margin << "   " << float(metapopsize_core)/float(pocc_core) << "   " << float(metapopsize_margin)/float(pocc_margin) << "   " << margin_position << endl;

			//end of simulation once the world is full
			if (margin_position == (WORLDDIM_X-1)) {
				saveResults(actrun, acttime);
				break;
			}

		}
		// close metapop output file
		outputmetapop.close();
	}
	cout << "simulation over" << endl;
	return 0;
}
