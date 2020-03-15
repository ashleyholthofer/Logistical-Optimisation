// Optimization
// reads in data from a file and finds the averaged optimized position, and also the closest town or city.

// includes

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <time.h>

using namespace std;

#define PI 3.14159

// keep a count function evaluations
int f_evals = 0;

//define a structure that holds the data
struct info {
	string Place;
	string Type;
	double Population;
	double latitude;
	double longitude;
};

//define a function that returns the Haversine distance
double Haversine(double lat1, double lat2, double lon1, double lon2){
	//do the calculation and return the distance d
	double dLat = lat2 - lat1;
	double dLong = lon2 - lon1;
	double a = pow(sin(dLat / 2), 2) + cos(lat1)*cos(lat2) * pow(sin(dLong / 2), 2);
	double c = 2 * atan2(sqrt(a), sqrt(1 - a));
	double d = 6371 * c;
	f_evals++; //adds a value to the counter for how many 
	return d;
}

//define a function that returns the weighted sum of all the Haversine distances from a point
double cost(double longitude, double latitude, vector <info> data, int n, double p){
	//N is the length of the vector 'data' and p is the total population
	//data is the vector holding all the information
	double value = 0;
	double weight;
	for (int i = 0; i < n; i++){
		weight = ((data[i].Population) / p);
		double Distance = Haversine(latitude, data[i].latitude, longitude, data[i].longitude);
		double WDistance = weight * Distance;
		value += WDistance;
	}
	return value;
}

//define a function to return a random number between two limits, lower and upper
double random_number(double lower, double upper, int n) {
	// n is the amount of bits to split the range into
	double r;
	r = lower + (rand() % (n + 1) * (1. / n) * (upper - lower));
	return r;
}

//main program
int main() {
	//declare vector to hold data
	vector < info > Data;

	//open the file
	ifstream file("GBplaces.csv");

	//define temporary values that are used to convert between data types
	string placet;
	string typet;
	string populationt;
	string latitudet;
	string longitudet;
	//define the holder for the data in a line
	string line;

	//count the total population
	double P = 0;

	int N = 0;

	if (file.is_open()) {
		//file did open
		//fet rid of the first line from the file
		getline(file, line);
		//read in the data
		while (!file.eof()) {
			//not at the end of the file
			getline(file, line);

			//now create a temporary info vector to hold the data
			info temp;
			//find the length of the first bit of string
			int position = line.find(',', 0);
			int position1 = line.find(',', position + 1);
			int position2 = line.find(',', position1 + 1);
			int position3 = line.find(',', position2 + 1);

			//add each string to the assigned temporary variable
			placet = line.substr(0, position);
			typet = line.substr(position + 1, position1 - (position + 1));
			populationt = line.substr(position1 + 1, position2 - position1);
			latitudet = line.substr(position2 + 1, position3 - position2);
			longitudet = line.substr(position3 + 1, line.length() - position3);

			//Now set these values into the info temp array, converting into the correct formats
			temp.Place = placet;
			temp.Type = typet;
			temp.Population = atof(populationt.c_str());
			temp.latitude = (PI / 180) * atof(latitudet.c_str());
			temp.longitude = (PI / 180) * atof(longitudet.c_str());

			//Add the population to the total population
			P += temp.Population;

			//push this temporary floating point two point vector into the Data vector
			Data.push_back(temp);
			N++;
		}

		//close the file
		file.close();

		//define variables
		double Cost, oldCost, newCost, bestCost;
		int dlon = 0, dlat = 0;
		vector < vector < double > > Results;

		srand(time(NULL)); // seeds random number generator

		//set the step taken in lon and lat (resolution if you like)
		double step = 0.0001;

		cout << "Input number of runs \n"; //Ask for how many times the program should run, for averaging later.
		int Runs;
		cin >> Runs;

		for (int j = 0; j < Runs; j++){

			// pick a starting point at random from withing the ranges shown (UK)
			double lon = random_number(-0.074316, 0.022654, 1000);
			double lat = random_number(0.879140, 0.997351, 1000);

			//find the cost of going from this place
			Cost = cost(lon, lat, Data, N, P);

			// write it out
			double deglon = (180 / PI)*lon;
			double deglat = (180 / PI)*lat;

			cout << "Starting Position: " << deglon << ", " << deglat << ", with cost: " << Cost << "\n";

			int iteration = 0;

			// main loop to optimize the latitude and longitude
			do {
				// record the fitness of the value before
				oldCost = Cost;
				//for now this is the current best cost
				bestCost = Cost;

				// moving through the neighbouring points, find a point with a better fitness
				for (int k = -1; k <= 1; k++) {
					for (int l = -1; l <= 1; l++) {
						// remove the point where k = l = 0 as this is the previous position
						if (k == 0 && l == 0) {
						}
						else {
							newCost = cost(lon + (step * k), lat + (step * l), Data, N, P); // value at the neighbouring point
							if (newCost <= bestCost) {
								//The value at i, j  is a better cost than the previous
								//save the step taken here
								dlon = k;
								dlat = l;
								//Save the newCost as the current best cost
								bestCost = newCost;
							}
							else {
								//the value for k, l  doesn't improve the cost, so discard and do the next one
							}
						}
					}
				}
				//from the position (lon, lat) the best cost neighbourgh has been identified.
				//now update evrything to this new position
				lon += step * dlon;
				lat += step * dlat;
				Cost = bestCost;

				deglon = (180 / PI)*lon;
				deglat = (180 / PI)*lat;

				iteration++;

			} while (Cost < oldCost);

			cout << "Longitude: " << deglon << ", Latitude: " << deglat << "\n"; //Write out the optimized result
			vector <double> OneResult; //This is hold this data
			OneResult.push_back(deglon);
			OneResult.push_back(deglat);
			//OneResult is a vector of length 2, (Longitude, Latitude) from this result
			Results.push_back(OneResult); //Now this data is stored in the matrix Results
		}

		cout << "Function evaluations: " << f_evals << "\n"; // write out the total number of function evaluations

		//Now do an average of all the results that were collected
		double averagelon, averagelat, sumlon = 0, sumlat = 0;
		for (int m = 0; m < Runs; m++){
			sumlon += Results[m][0];
			sumlat += Results[m][1];
		}
		averagelon = sumlon / Runs;//average from all the runs of optimization
		averagelat = sumlat / Runs;//same for latitude

		

		//what is this close to?
		double alon, alat, distance;//alon and alat are values of the latitude and longitude from a point in Data
		vector <double> distances;//vector to hold the distance from the result position to each town/city
		//converting into radians
		double averlon = (PI / 180) * averagelon;
		double averlat = (PI / 180) * averagelat;

		for (int n = 0; n < N; n++){
			alon = Data[n].longitude;
			alat = Data[n].latitude;
			distance = Haversine(alat, averlat, alon, averlon);
			distances.push_back(distance);
		} //now have a vector of the distances to each city/town

		//need to find the city/town with the lowest distance from the result distance
		//define some indicies
		int m = 0, o = 0;
		//a variable that records the current best index
		int index;
		//define variables, the 100 is arbitary
		double idistance = 100, newdistance, olddistance;
		do {
			olddistance = idistance;//record the distance from the last index
			for (o = 0; o < N-1; o++){
				if (distances[o] < idistance){
					//the distance from the place at index o is closer, so record this better index
					index = o;
				}
				else{
					//do nothing
				}
			}
			//the distance at this new index is recorded as newdistance
			newdistance = distances[index];
			idistance = newdistance;
			//show the consol the distance from this index
			cout << index << ": " << newdistance << "\n";
		} while (olddistance > newdistance);
		//have identified the index that has the closest distance to the location

		//output the results
		cout << "\n" << "RESULTS:\n";
		cout << "Averaged Longitude: " << averagelon << ", Averaged Latitude: " << averagelat << "\n";
		cout << "This location is close to: " << Data[index].Place << "\n";

	}

	else {
		//something went wrong
		cout << "Could not open the file\n";
		//terminate the program
		exit(1);
	}

	return 0;
}