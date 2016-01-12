#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cstring>
#include <vector>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <chrono>
#ifdef __unix__			//delete this and the following includes to get rid of glut
	#include <GL/freeglut.h>
#elif defined(_WIN32) || defined(WIN32)
	#include <freeglut.h>
#endif


using namespace std;

typedef struct {
	double x, y;	//stores the x,y location of the city
	int city;	//stores the id of the city
} city;

vector<city> cities;	//holds a list, in order, of all the cities
vector<vector<city>> chromosomes;	//holds a 2 dimensional vector of all individual chromosomes
vector<vector<city>> gen_buffer;	//stores the next generation while it is being created
vector<city> elite;
double elite_fitness = 100000.0;	//initially very high for comparison
int numcities = 52;
int population_size;
int max_gen_without_improvement = 10;
int mutation_rate = 10;
int original_mut_rate; //this is just to store the original in a dynamic mutation rate
char dyn_mut = 'y';
int seed = 1;
int crossover_rate = 100;
vector<double> average_fitness_vector;
//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
int fitness[150]; 	//this needs to be a bigger array than the population size
int cross_type = 1; //0 = Uniform Order Crossover   1 = 2-point crossover
int k;	//this is the value of k used in the tournament selection
int bit_probability = 5; //probability between 1-10 that the bitmask for oux will create 1's
ofstream outdata;	//this is for writing to the .csv file

const char* const DELIMITER = " ";

/*
	calculates and returns the cartesian distance between two points on a plane.
	It does not square root the result to save computation time and it is irrelivant.
*/
float getCartDistance(float x1, float y1, float x2, float y2){
	float xdiff = pow((x2 - x1), 2);
	float ydiff = pow((y2 - y1), 2);
	return sqrt(xdiff + ydiff);
}
/*
	Adds and returns the route distance of the indexed chromosome
*/
int evaluateFitness(int index){
	double sum = 0;
	for (int i = 1; i < numcities; ++i) {	//add up distance between cities in chromosome
		sum = sum + (getCartDistance(chromosomes[index][i-1].x, chromosomes[index][i-1].y, chromosomes[index][i].x, chromosomes[index][i].y));
	}
	if (sum < elite_fitness) {
		elite_fitness = sum;		//update the elite fitness int if a better one is found
		elite = chromosomes[index];	//undate the elite chromosome
	}
	return sum;
}

/*
	selects a parent from the current generation using the tournament selection technique,
	based on the chromosomes fitness, and returns the index of the winning chromosome
*/
int tournamentSelection(int k){
	int candidates[k];		//array of k indexes of the k potential parents that will be compared
	int temp; //stores the winner of the tournament selection
	for (int i = 0; i < k; ++i) {
		candidates[i] = rand() % population_size;
	}
	temp = candidates[0];
	for (int i = 1; i < k; ++i) {	//scan through to find winner
		if (fitness[candidates[i]] < fitness[temp]) temp = candidates[i];
	}
	return temp;	//return index of winning chromosome
}
/*
	find and returns a suitable city, searching left to right, to place in an open spance in the child
	from the opposite parent.
*/
city findSuitableCity(vector<city> child, vector<city> parent){
	city temp = parent[0];
	bool city_found = false;
	for (int i = 0; i < numcities; ++i){	//start on the left of the parent
		temp = parent[i];
		city_found = true;
		for (int j = 0; j < numcities; ++j){//see if the city already exists in child
			if (temp.city == child[j].city) {
				city_found = false;		//if it does, city_found = false, and break out of nested loop
				break;
			}
		}
		if (city_found) return temp;	//if city_found is still true, return the city
	}
	cout << "error: suitable city not found" << endl;
	return temp;	//this should never be reached
}
/*
	this method takes two parent vectors, performs a Uniform Order crossover with a bit mask, and then 
	returns the child as a vector.
*/
vector<vector<city>> uox(vector<city> parent1, vector<city> parent2){
	vector<vector<city>> child;
	child.resize(2);
	child[0].resize(numcities); 	//child 1
	child[1].resize(numcities); 	//child 2
	int bitmask[numcities];
	int bit;
	for (int i = 0; i < numcities; ++i)	{//create random bitmask
		bit = rand() % 10 + 1; //generate number between 1 - 10
		if (bit <= bit_probability) bitmask[i] = 1;	//if random num is less than bit_prob, set to 1
		else bitmask[i] = 0;	//else set to 0
	}
	
	for (int i = 0; i < numcities; ++i){	//copy bit mask 1's into child
		if (bitmask[i] == 1) {
			child[0][i] = parent1[i];
			child[1][i] = parent2[i];
		}
	}
	for (int i = 0; i < numcities; ++i){
		if (bitmask[i] == 0){
			child[0][i] = findSuitableCity(child[0], parent2);	//fill in missing citys with opposite parent
			child[1][i] = findSuitableCity(child[1], parent1);	//fill in missing citys with opposite parent
		}
	}
	return child;
}
/*
	accepts a vector and a city, checks if the city exists in the list, returns based on this.
*/
bool isValid(vector<city> list, city in){
	bool already_exists = false;
	for (int i = 0; i < numcities; ++i){
		if (in.city == list[i].city) return false;
	}
	return true;
}
/*
	this method takes two parent vectors, performs a 2 Point crossover, and then 
	returns the child as a vector.
*/
vector<vector<city>> twopoint(vector<city> parent1, vector<city> parent2){
	vector<vector<city>> child;
	child.resize(2);
	child[0].resize(numcities); 	//child 1
	child[1].resize(numcities); 	//child 2
	city fake_city;
	fake_city.city = 1000;
	fake_city.x = 10000.0;
	fake_city.y = 10000.0;
	bool already_exists = false;
	int subset_size = rand() % numcities;
	int subset_start = rand() % (numcities-subset_size);
	
	for (int i = 0; i < numcities; ++i){	//copy parents to children
			child[0][i] = parent1[i];
			child[1][i] = parent2[i];
	}

	for (int i = subset_start; i <= subset_start + subset_size; ++i){
		child[0][i] = fake_city;
		child[1][i] = fake_city;
		if (isValid(child[0], parent2[i])){
			child[0][i] = parent2[i];
		} else {
			child[0][i] = parent1[i];
		}
		if (isValid(child[1], parent1[i])){
			child[1][i] = parent1[i];
		} else {
			child[1][i] = parent2[i];
		}
	}
	return child;
}
/*
	crosses over 2 parent vectors and returns a child vector
*/
vector<vector<city>> crossover(vector<city> parent1, vector<city> parent2, int method){	//0 = UOX  1 = 2-point
	if (method == 0) return uox(parent1, parent2);
	else if (method == 1) return twopoint(parent1, parent2);
}
/*
	This mutation works by taking a section of the chromosome, reversing it, and placing it back in the same
	spot
*/
vector<city> mutate(vector<city> mutatee){
	int subset_size = rand() % (numcities/2);
	int subset_start = rand() % (numcities-subset_size);

	reverse(mutatee.begin()+subset_start, mutatee.begin()+subset_start+subset_size);

	return mutatee;
}
/*
	runs the genetic algorithm using the parameters defined by the user
*/
void GeneticAlgorithm(){
	int tournamentWinner;
	double average_fitness = 0;
	outdata.open("results.csv", ios::app);
	int count = 0;
	double elite_compare = elite_fitness;	//this is to determine if the fitness is improving
	vector<city> parent1;
	vector<city> parent2;
	vector<vector<city>> children;		//2 dimensional because crossover returns 2 children
	children.resize(2);		//big enough to store 2 children
	children[0].resize(numcities);	//big enough to store child1 chromosome
	children[1].resize(numcities);	//big enough to store child2 chromosome
	parent1.resize(numcities);
	parent2.resize(numcities);
	outdata << "Population Size: " << population_size << " , K: " << k << " , Seed: " << seed;
	outdata << " , Max Generations w/out improvement: " << max_gen_without_improvement;
	outdata << " , Mutation Rate: " << mutation_rate << " , Crossover Rate: " << crossover_rate << endl;
	while (count < max_gen_without_improvement){	//run GA while the elite is still improving
		for (int i = 0; i < population_size; ++i){	//evaluate the fitnesses and store in corresponding array
			fitness[i] = evaluateFitness(i);
			average_fitness += fitness[i];
		}
		average_fitness = average_fitness / population_size;
		average_fitness_vector.push_back(average_fitness);

		for (int i = 0; i < population_size-1; i += 2){	//tournament selection until generation buffer is full
			tournamentWinner = tournamentSelection(k); //place index of winner in "tournament winner"
			for (int j = 0; j < numcities; ++j){	//copy chromosome of the winning generation into parent1
				parent1[j] = chromosomes[tournamentWinner][j];
			}
			tournamentWinner = tournamentSelection(k); 
			for (int j = 0; j < numcities; ++j){	//copy chromosome of the winning generation into parent2
				parent2[j] = chromosomes[tournamentWinner][j];
			}
			if (rand() % 100 <= crossover_rate) {	//if the random is less than crossover rate
				children = crossover(parent1, parent2, cross_type);
				gen_buffer[i] = children[0];
				gen_buffer[i+1] = children[1];
			} else {								//else, just place these parents into next gen
				gen_buffer[i] = parent1;
				gen_buffer[i+1] = parent2;
			}
			if (rand() % 100 <= mutation_rate){		//simple mutation
				gen_buffer[i] = mutate(gen_buffer[i]);
			}
			if (rand() % 100 <= mutation_rate){
				gen_buffer[i+1] = mutate(gen_buffer[i+1]);
			}
		}
		for (int i = 0; i < population_size; ++i){	//copy buffer into new generation
			for (int j = 0; j < numcities; ++j){
				chromosomes[i][j] = gen_buffer[i][j];
			}
		}
		chromosomes[population_size-1] = elite;	//put elite at the end
		if (dyn_mut == 'y') mutation_rate = original_mut_rate + (count*1.5);	//dynamically change mutation rate as the number of gens w/out improvement increases
		cout << "elite_fitness: " << elite_fitness << "   mutation rate: " << mutation_rate << "   Average Fitness: " << average_fitness << endl;
		outdata << elite_fitness << ",";
		if (elite_compare == elite_fitness) count++;	//if the elite fitness hasnt changed increase count
		else {
			elite_compare = elite_fitness;	//if elite fitness has changed, set count back to 0
			count = 0;
		}
	}
	outdata << endl;
	for (int i = 0 ; i < average_fitness_vector.size(); ++i){
		outdata << average_fitness_vector[i] << ",";
	}
	outdata << endl;
	outdata << endl;
}
/*
	Makes the initial 50 generations as a 2-d vector
*/
void makeInitPop(int population_size){
	default_random_engine engine(seed);
	chromosomes.resize(population_size);	//set vector to be large enough to hold chromosomes
	gen_buffer.resize(population_size);		//set generation buffer vector to match population size
	for (int i = 0; i < population_size; ++i) {
		chromosomes[i].resize(numcities);	//set nested vector to right length to store cities
		gen_buffer[i].resize(numcities);	//set nested generation buffer to fit chromosome lengths
		for (int j = 0; j < numcities; ++j) {//copy cities into the chromosome
			chromosomes[i][j] = cities[j];
		}
		shuffle(chromosomes[i].begin(),chromosomes[i].end(),engine);//shuffle the chromosome
	}
}
/*
	reads the berlin52.tsp file in
*/
int readFileIn(){
	bool start_copying = false;
	 
	cities.resize(52);
	int index = 0;
	ifstream fin;
	fin.open("berlin52.tsp"); // open a file
	if (!fin.good()) 
		return(1); // exit if file not found
	// read each line of the file
	while (!fin.eof()) {
		if (index > 57) {
			cout << "return" << endl;
			return(0);
		}
	// read an entire line into memory
		char buf[46];
		fin.getline(buf, 46);

		// parse the line into blank-delimited tokens
		int n = 0; // a for-loop index

		// array to store memory addresses of the tokens in buf
		const char* token[3] = {}; // initialize to 0

		// parse the line
		token[0] = strtok(buf, DELIMITER); // first token
		if (token[0]) {// zero if line is blank
		
			for (n = 1; n < 3; n++) {
				token[n] = strtok(0, DELIMITER); // subsequent tokens
				if (!token[n]) break; // no more tokens
			}
		}
		if (index >= 6) start_copying = true;

		if (start_copying){
			cities[index-6].city = atoi(token[0]);
			cities[index-6].x = (float) atoi(token[1]);
			cities[index-6].y = (float) atoi(token[2]);
		}
		index++;
	}
}
/*
	draws the map with glut library, delete this method to get rid of the glut
*/
void makemap(void) {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(255, 255, 255);
    glBegin(GL_LINE_LOOP);
    	for (int i = 0; i < numcities; ++i){
    		glVertex2i(elite[i].x, elite[i].y);
    	}
	glEnd();
	glFlush();
}
int main(int argc, char **argv)
{
	cout << "GA Population Size: ";		//user inputs parameters
	cin >> population_size;
	cout << "Value of K: ";
	cin >> k;
	cout << "Maximum allowed generations without improvement: ";
	cin >> max_gen_without_improvement;
	cout << "Mutation rate (0-100): ";
	cin >> mutation_rate;
	cout << "Dynamic mutation (y / n): ";
	cin >> dyn_mut;
	cout << "Crossover rate (0-100): ";
	cin >> crossover_rate;
	cout << "seed: ";
	cin >> seed;
	srand(seed);
	original_mut_rate = mutation_rate;
	readFileIn();
	makeInitPop(population_size);
	GeneticAlgorithm();		
	cout << endl;				//printElite to standard out and to .csv file
	cout << "Final Fitness: " << elite_fitness << endl;	//print out the final results
	cout << "Final Chromosome: " << endl;
	for (int j = 0; j < numcities; ++j){
		outdata << elite[j].city << ",";
		cout << elite[j].city << " ";
	}	cout << endl;
		outdata << endl;
	glutInit(&argc, argv);				//for no glut, erase from here: ...
	glutCreateWindow("Traveling Salesman");
	gluOrtho2D(0, 1800, 0, 1200);
	glutPostRedisplay();
	glutDisplayFunc(makemap);
	glutMainLoop();					//... to here
	
	return 0;
}
