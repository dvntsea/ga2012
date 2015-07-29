//GA_atom.c   Doug Van Nort    GA object for max


//notes:
//GA object accepts lists of parameters and computes new generations using genetic algorithm operators
// many post messages are left here commented out - for convenience of debugging at a later point


#include "ext.h"
#include "ext_obex.h"	

#include <math.h>
#include <stdlib.h>



//#define MAXRAND 32767.0 
#define MAXRAND 65536.0 
//#define MAXRAND RAND_MAX 



#define ASSIST_INLET  (1)	
#define ASSIST_OUTLET (2)



typedef struct genalg

{
	struct object t_ob;
	void *t_out;
	void *t_out1;
	float t_sum;
	double *allfitness;		//from 0....population-1 contains parent fitness, from population....2*population-1 contains chidren's.  we use one array in case we need to sort list.
	double *fitty;
	double *fittest;	
	double t_mutrate;
	double t_mutdepth;
	int t_crossrate;
	Atom *matrix;			//struct containing members of current population and their children
	double *temparray;
	double *outarray;
	double *children;
	double *max;
	double *min;
	double *ascsort;
	long *sortedindex;
	double diameter;
	long population;
	long dimension;
	long m;
	long n;
	long nfit;
	long childex;
	float repnum;
} Genalg;








void *class;
void genalg_bang(Genalg *x);
//void genalg_float(Genalg *x, float f);	//does not currently accept float or int - but they are ready and waiting
//void genalg_int(Genalg *x, long n);
void *genalg_new(Symbol *s, short ac, Atom *av);
void genalg_free(Genalg *x);
void genalg_assist(Genalg *x, Object *b, long msg, long arg, char *s);

 
void genalg_fitlist(Genalg *x, Symbol *s, short ac, Atom *av); //  input for list of fitness values
void genalg_member(Genalg *x, Symbol *s, short ac, Atom *av);  //stores members of population


void genalg_clear(Genalg *x, Symbol *s, short ac, Atom *av); 

void genalg_outmeth(Genalg *x);		//method for outputting list

void genalg_randomize(Genalg *x);	//randomize population members

void genalg_mutate(Genalg *x,int n);

//void genalg_mutrate(Genalg *x, float f);
void genalg_mutrate(Genalg *x, Symbol *s, short msg, Atom *av);
void genalg_mutdepth(Genalg *x, Symbol *s, short msg, Atom *av);

void genalg_replacemeth(Genalg *x, float f);

void genalg_max(Genalg *x, Symbol *s, short ac, Atom *av);
void genalg_min(Genalg *x, Symbol *s, short ac, Atom *av); 

void genalg_newfit(Genalg *x, Symbol *s, short ac, Atom *av); 

void genalg_cross(Genalg *x);		//crossover operator
void genalg_select(Genalg *x);		//selection operator
void genalg_replace(Genalg *x);		//replacement operator
void genalg_testfit(Genalg *x);		//fitness metric - tested on sucessive generations. overrided when user fitness assigned
void genalg_diam(Genalg *x);		//calculates diameter/diagonal of parameter space
void genalg_fittest(Genalg *x);		//finds the fittest member and stores it in an array


//void genalg_crossrate(Genalg *x, double f);		//crossover rate
double amax(double arg1, double arg2);			//my max function
double amin(double arg1, double arg2);			//my min function
void swap(double *i, double *j);				//my swap function
double fit_dist(Genalg *x, short i, double *u, double *v);	//euclidean distance metric
//void sort_fit(Genalg *x, double *u, long begin, long end);	

//**quicksort methods - used in replacement scheme 2*//
void quickSort(double *a, int lb, int ub);
void insertSort(double *a, int lb, int ub);
int partition(double *a, int lb, int ub);
//****////

void main(void)
{
	setup((t_messlist **)&class, (method)genalg_new,(method)genalg_free, (short)sizeof(Genalg), 0L, A_GIMME, 0);
	
	addbang((method)genalg_bang);		//bang messages calls operators to complete one ga cycle
	
//	addfloat((method)genalg_float);
//	addint((method)genalg_int);

	///****add all the messages ga likes to see***/////
	addmess((method)genalg_fitlist, "fitlist", A_GIMME, 0);
	addmess((method)genalg_member, "member", A_GIMME, 0);
	addmess((method)genalg_clear, "clear", A_GIMME, 0);
	addmess((method)genalg_mutate, "mutate", A_GIMME, 0);
	addmess((method)genalg_mutrate, "mutrate", A_GIMME, 0);
	addmess((method)genalg_mutdepth, "mutdepth", A_GIMME, 0);
	addmess((method)genalg_randomize, "randomize", A_GIMME, 0);
	addmess((method)genalg_replacemeth, "replacemeth", A_GIMME, 0);
	//addmess((method)genalg_crossrate, "crossrate", A_GIMME, 0);
	addmess((method)genalg_max, "max", A_GIMME, 0);
	addmess((method)genalg_min, "min", A_GIMME, 0);
	
	addmess((method)genalg_newfit, "newfit", A_GIMME, 0);
	
	addmess((method)genalg_assist, "assist", A_CANT,0); 
	
	/////***************************////////////
	
	finder_addclass("All Objects", "ga");		//add ga to the list of objects in Max
	
}

void *genalg_new(Symbol *s, short ac, Atom *av)
{
	Genalg *x;
	short i;
	
	x = (Genalg *)newobject(class);			//create our object
	
	
	//******memory allocation*****//////
	
	x->matrix = (Atom *)getbytes((unsigned short)(2*sizeof(x->t_mutrate)*av[0].a_w.w_long*av[1].a_w.w_long));
	
	x->temparray = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[0].a_w.w_long*av[1].a_w.w_long));
	
	x->outarray = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[1].a_w.w_long));
	
	x->children = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[0].a_w.w_long*av[1].a_w.w_long));
	
	x->allfitness = (double *)getbytes((unsigned short)(2*sizeof(x->t_mutrate)*av[0].a_w.w_long));
		
	x->fitty = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*(av[0].a_w.w_long) + 1));
	
	x->fittest = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[1].a_w.w_long));
	
	x->max = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[1].a_w.w_long));
	
	x->min = (double *)getbytes((unsigned short)(sizeof(x->t_mutrate)*av[1].a_w.w_long));
	
	x->ascsort = (double *)getbytes((unsigned short)(2*sizeof(x->t_mutrate)*av[0].a_w.w_long));
	
	x->sortedindex = (long *)getbytes((unsigned short)(2*sizeof(x->nfit)*av[0].a_w.w_long));
	
	////**********************///////////////
	
	
	x->population = av[0].a_w.w_long;	//object accepts values for population size
	x->dimension = av[1].a_w.w_long;	//and dimension
	//post("population size = %d",x->population);		//tell us what they are
	//post("dimension size = %d",x->dimension);
	
	//floatin(x, 0);
	x->t_out = listout(x);	
	x->t_mutrate = 0.0;		//default mutation rate
	x->t_mutdepth = 0.0;		//default mutation depth
	x->t_crossrate = (int)(x->population/2);		//default crossover rate. corresponds to 1.0.
	x->childex = 0; //default child index.  Used to keep track of which children to replace.
	x->repnum = 1.0;  //default replacement method: 1 = children replace parents.

//post("RAND_MAX is %d",RAND_MAX);


//set default max to 1. and min to 0.
short j;

	{
	for (j=0;j < av[1].a_w.w_long;j++)
		{
		x->max[j] = 1.;
		}
	}


	{
	for (j=0;j < av[1].a_w.w_long;j++)
		{
		x->min[j] = 0.;
		}
	}
	
////////////

	
	return (x);		//give us our object
}




void genalg_bang(Genalg *x)			//calls operators for selection, crossover, mutation, and replacement. also tests fitness.
{
short i;
	x->childex = 0;		//our pointer for replacement of previous generation; initialized to location 0.
	for (i=0;i < x->t_crossrate;i++)
		{
		genalg_select(x);
		genalg_cross(x);
		}
	genalg_mutate(x,1);
	genalg_testfit(x);
	genalg_replace(x);
	genalg_outmeth(x);		
}


/*void genalg_float(Genalg *x, float f)
{
//blank for now	 
}

void genalg_int(Genalg *x, long n)
{
//blank for now	 
}	*/


void genalg_max(Genalg *x, Symbol *s, short ac, Atom *av)		//define max in each dimension
{
short j;
if(ac != x->dimension)						//error checking///
		post("error! max list must be of length dimension");
else
	{
	for (j=0;j < ac;j++)
		{
		x->max[j] = av[j].a_w.w_float;
		//post("max[%d] = %f",j,x->max[j]);
		}
	}
	
genalg_diam(x);		//calculate new diameter of parameter space

}
void genalg_min(Genalg *x, Symbol *s, short ac, Atom *av)		//define max in each dimension
{
short j;
if(ac != x->dimension)						//error checking///
		post("error! min list must be of length dimension");
else
	{
	for (j=0;j < ac;j++)
		{
		x->min[j] = av[j].a_w.w_float;
		//post("min[%d] = %f",j,x->min[j]);
		}
	}
	
genalg_diam(x);			//calculate new diameter of parameter space

} 


void genalg_diam(Genalg *x)		//finds maximum distance of parameter space. used to normalize fitness function
{
short j;
double dist = 0;
for (j = 0;j < x->dimension;j++)
	dist += (x->max[j] - x->min[j]) * (x->max[j] - x->min[j]);
	
x->diameter = sqrt(dist);

}



void genalg_fitlist(Genalg *x, Symbol *s, short ac, Atom *av)   // new fitness list input
{
short i;
short j;
double damax;
double temp;
short n;

x->t_sum = 0;					//clear sum
	
	for (i = 0; i < ac; i++)     //running sum 
	{
		if (av->a_type == A_LONG)
			x->t_sum += av[i].a_w.w_long;
		else if (av->a_type == A_FLOAT)
			x->t_sum += av[i].a_w.w_float;
			
	}
	
	for (i = 0; i < ac; i++)     //normalize fitness values 
	{
		x->allfitness[i] = av[i].a_w.w_float / x->t_sum;
			
	}
	
	damax = 0;  
	temp = damax;

	for (i=0;i < x->population;i++)				//find the fittest 
		{
		temp = amax(damax,x->allfitness[i]);
		if(temp != damax)
			{
			x->nfit = i;
			damax = temp;
			}
		}
		
	genalg_fittest(x);	//only time we store a new fittest member is when user inputs new fitness list, or defines new fittest
}

void genalg_fittest(Genalg *x)
{
short j;

	for (j=0;j < x->dimension;j++)			//save copy of the fitttest
		{
		x->fittest[j] = x->matrix[x->nfit*x->dimension + j].a_w.w_float;
		//post("fittest[%d] = %f at %d",j,x->fittest[j],x->nfit);
		}
	 
	 
}

void genalg_newfit(Genalg *x, Symbol *s, short ac, Atom *av)		//input a new fittest member. first index is method type, others are parameters of new fittest
{
short j;

	if(av[0].a_w.w_long == 1)	//newfit scheme is 1, add new fitness member to current population
	{
	int rando;
	//post("newfit method 1");
	rando = (int)((Random()/MAXRAND) * x->population-1);	//randomly select one member of population to replace
	rando = abs(rando);
	//post("rando is %d",rando);
	for (j=1;j <= x->dimension;j++)
		{
		x->fittest[j-1] = av[j].a_w.w_float;	//define new fittest
		//post("new fittest %d is %f inserted in pop at %d",j,x->fittest[j-1],rando);
		x->matrix[rando*x->dimension + j-1].a_w.w_float = x->fittest[j-1];	// place fittest in population at chosen random spot
		}
	genalg_testfit(x);		//update all fitness values
	 }	
	 
	 else if(av[0].a_w.w_long == 2)
	 {
	 //post("newfit method 2");
	 for (j=1;j <= x->dimension;j++)
		{
		x->fittest[j-1] = av[j].a_w.w_float;	//define new fittest
		//post("new fittest %d is %f",j-1,x->fittest[j-1]);
		}
	 genalg_testfit(x);		//update all fitness values
	 }	

}

void genalg_testfit(Genalg *x)		//compute distance from each parent and child to the fittest member
{
short i;
short j;
double dist = 0;
double dist2 = 0;
double sum = 0;
double damax;
double temp;
//double temparray[x->population*x->dimension];


for (i=0;i < x->population;i++)		//need to copy Atom parents into local array of doubles to pass to fit_dist
	{
	for (j=0;j < x->dimension;j++)
		x->temparray[i*x->dimension + j] = x->matrix[i*x->dimension + j].a_w.w_float;
	}


for (i=0;i < x->population;i++)		//test children
	{
	dist = fit_dist(x,i,x->fittest,x->children);
//	post("child %d fitdist = %f",i,dist);
	x->allfitness[x->population + i] = ((x->diameter - dist) / x->diameter) + 0.001;		//our fitness metric: linearly related to distance, in range [0.001,1.0001] (need to avoid 0 fitness)
	sum += x->allfitness[x->population + i];								//running sum of children fitnesses
	//post("childfitness %d = %f",i,x->allfitness[x->population + i]);
	}
	
for (i=0;i < x->population;i++)		//test parents
	{
	dist2 = fit_dist(x,i,x->fittest,x->temparray);
//	post("parent %d fitdist = %f",i,dist2);
	x->allfitness[i] = ((x->diameter - dist2) / x->diameter) + 0.001;
	sum += x->allfitness[i];												//running sum of parent fitnesses added to children's
//	post("parentfitness %d = %f",i,(x->diameter - dist2) / x->diameter);
	}


	

for (i=0;i < 2*x->population;i++)	//normalize new fitness values
	{
	x->allfitness[i] = x->allfitness[i] / sum;
	
	}
	
	for (j=0;j < 2*x->population;j++)	
	{
	for (i=0;i < 2*x->population && i != j;i++)	//compare fitness to each other member's
		{
		while (x->allfitness[j] == x->allfitness[i])	//while two in the list are equal
				x->allfitness[i] = x->allfitness[i] + (Random() / MAXRAND)*.00001;	//slightly perturb one. this makes sorting much easier and doesn't affect selection much
		}
	//post("allfitness %d = %f",j,x->allfitness[j]);
	}
	
	
	
	//****find new fittest*****////
	damax = 0;  
	temp = damax;

	for (i=0;i < 2*x->population;i++)	//find the current fittest, but don't replace it with our stored fittest - which is still our "goal" 
		{
	//	post("allfitness %d = %f",i,x->allfitness[i]);
		temp = amax(damax,x->allfitness[i]);
		if(temp != damax)
			{
			x->nfit = i;
			damax = temp;
			}
		}
	//////////**********///////////
	
}

void genalg_member(Genalg *x, Symbol *s, short ac, Atom *av)	//first index is individual #, others are the parameter values
{
short j;


double damax = x->max[0];  
double temp = damax;
double damin = x->min[0];
for (j=1;j < x->dimension;j++)
		{
		temp = amax(damax,x->max[j]);
		if(temp != damax)
			{
			damax = temp;
			}
		}
for (j=1;j < x->dimension;j++)
		{
		temp = amin(damin,x->min[j]);
		if(temp != damin)
			{
			damin = temp;
			}
		}

if (damax == 0 && damin == 0)				//error checking. if max and min arrays are identically zero...they are not defined.
	post("error! must define min and max for each dimension");
else
	{	
	if(ac != x->dimension + 1)						//error checking///
		post("error! member list must be of length dimension + 1");
	else
		{	
		for(j=0;j < x->dimension;j++)				
			{
			if (x->max[j] < av[j+1].a_w.w_float)							//constrain values within [min,max] for each dimension
				SETFLOAT(x->matrix+(av[0].a_w.w_long*x->dimension + j), x->max[j]);
			else if (av[j+1].a_w.w_float < x->min[j])
				SETFLOAT(x->matrix+(av[0].a_w.w_long*x->dimension + j), x->min[j]);
			else
				SETFLOAT(x->matrix+(av[0].a_w.w_long*x->dimension + j), av[j+1].a_w.w_float);
			}			
		 }
	}
}



void genalg_outmeth(Genalg *x)
{
//short i;
//short j;
int sz = (int)(x->dimension*x->population);	//size of our matrix of members

outlet_list(x->t_out, 0L, sz, x->matrix);	//send out list of new pop members
/*for(i=0; i < x->population ;i++)
	{
	
	post("fitness %d = %f \n value = %f,%f,%f,%f",i,x->allfitness[i],x->matrix[i*x->dimension].a_w.w_float,x->matrix[i*x->dimension + 1].a_w.w_float,x->matrix[i*x->dimension + 2].a_w.w_float,x->matrix[i*x->dimension + 3].a_w.w_float);

	}*/
}



void genalg_outmeth_newfangled(Genalg *x)
{
short i;
short j;
int sz = (int)(x->dimension*x->population);	//size of our matrix of members

if (sz>256)
	{
	
	for (i=0;i < x->population;i++)		//need to copy Atom parents into local array of doubles to pass to fit_dist
	{
	for (j=0;j < x->dimension;j++)
		{
		x->outarray[j] = x->matrix[i*x->dimension + j].a_w.w_float;
		}
	outlet_list(x->t_out, 0L, x->dimension, x->outarray);
	}
	
	}
else
	{
outlet_list(x->t_out, 0L, sz, x->matrix);	//send out list of new pop members
	}
/*for(i=0; i < x->population ;i++)
	{
	
	post("fitness %d = %f \n value = %f,%f,%f,%f",i,x->allfitness[i],x->matrix[i*x->dimension].a_w.w_float,x->matrix[i*x->dimension + 1].a_w.w_float,x->matrix[i*x->dimension + 2].a_w.w_float,x->matrix[i*x->dimension + 3].a_w.w_float);

	}*/
}






void genalg_mutate(Genalg *x, int n)
{
short i;
short j;
float randnum;
float randnum2;

if (n == 1)
{

for(i=0;i < 2*x->population;i++)	//for each member
	{
	randnum = Random() / MAXRAND;	//generate random number
	//post("randum is %f",fabs(randnum));
	if(fabs(randnum) < x->t_mutrate)	//if number is within mutation rate 
		{
		//post("mutate %d,%d",i);
		for(j=0;j < x->dimension; j++)	//then mutate each dimension
			{
			randnum2 = Random() / MAXRAND;		// get second random number for each dimension
			//post("mutate %d,%d = %f",i,j,randnum2);
			
			SETFLOAT(x->matrix+(i*x->dimension + j),x->matrix[i*x->dimension + j].a_w.w_float + (x->max[j]/4)*x->t_mutdepth*randnum2);
			if (x->matrix[i*x->dimension + j].a_w.w_float < x->min[j])
				x->matrix[i*x->dimension + j].a_w.w_float = x->min[j];
			if (x->matrix[i*x->dimension + j].a_w.w_float > x->max[j])
				x->matrix[i*x->dimension + j].a_w.w_float = x->max[j];
			
			}
		}
	}
}

if (n == 2)
{
post("mutation scheme 2 not currently implemented");		//mutation scheme 2 has been deleted for now - but it will be back!
}		
	
}

/*
void genalg_mutrate(Genalg *x, float f)	//mutation rate defines likelihood of mutation for each member
{
post("incoming f = %f",f);
x->t_mutrate = f;
post("mutation rate prescale = %f",x->t_mutrate);
if(x->t_mutrate > 1.0)
	 x->t_mutrate = 1.0;
if(x->t_mutrate < 0.0)
	 x->t_mutrate = 0.0;
	 
post("mutation rate = %f",x->t_mutrate);
}
*/

void genalg_mutrate(Genalg *x, Symbol *s, short msg, Atom *av)	//mutation rate defines likelihood of mutation for each member, and the depth of the mutation if chosen
{

if (av->a_type == A_LONG)
	x->t_mutrate = (double)av->a_w.w_long;
else if (av->a_type == A_FLOAT)
	x->t_mutrate = av->a_w.w_float;

//post("mutation rate = %f",x->t_mutrate);
}


void genalg_mutdepth(Genalg *x, Symbol *s, short msg, Atom *av)	//mutation rate defines likelihood of mutation for each member, and the depth of the mutation if chosen
{

if (av->a_type == A_LONG)
	x->t_mutdepth = (double)av->a_w.w_long;
else if (av->a_type == A_FLOAT)
	x->t_mutdepth = av->a_w.w_float;

//post("mutation depth = %f",x->t_mutdepth);
}





void genalg_cross(Genalg *x)
{
short j;
short i;
///cross scheme 1 - simple one point cross.
long r = (long)((fabs(Random() / MAXRAND) * (x->dimension - 1)) + 0.5);
//post("cross point is %d and cross guys are %d,%d",r,x->n,x->m);
	for(j=0;j <= r;j++)
		{
		x->children[x->childex*x->dimension+j] =  x->matrix[x->n*x->dimension+j].a_w.w_float; 
		x->children[(x->childex+1 % x->population)*x->dimension+j] = x->matrix[x->m*x->dimension+j].a_w.w_float;
		}
	for(j=r+1;j < x->dimension;j++)
		{
		x->children[(x->childex+1 % x->population)*x->dimension+j] =  x->matrix[x->n*x->dimension+j].a_w.w_float; 
		x->children[x->childex*x->dimension+j] = x->matrix[x->m*x->dimension+j].a_w.w_float;
		}
/////end cross scheme 1

		
x->childex = (x->childex + 2) % x->population;  //cycle through pairs of children to replace. updated for next pass.

///*****place children in one matrix with parents for replacement ordering***////
for(i=0;i < x->population;i++)
	{
	for(j=0;j < x->dimension;j++)
		{
		SETFLOAT(x->matrix+((i+x->population)*x->dimension+j), x->children[i*x->dimension+j]);
		//post("child %d,%d in matrix value = %f",i,j,x->matrix[(i+x->population)*x->dimension+j].a_w.w_float);
		}
	}
/////////***********************////////////


/*for(i=0;i < x->population;i++)
	{
	for(j=0;j < x->dimension;j++)
		{
		post("parent %d,%d in matrix value = %f",i,j,x->matrix[i*x->dimension+j].a_w.w_float);
		}
	}*/

	

}

void genalg_select(Genalg *x)
{

short i,n,m;

///selection scheme 1 - prob selection based on fitness
double randy = (fabs(Random() / MAXRAND));

x->fitty[0] = 0;
for(i=1;i <= x->population;i++)
	x->fitty[i] = x->allfitness[i-1] + x->fitty[i-1];

	
for(i=0;i < x->population;i++)	//partition the interval [0,1] with mesh size equal to fitness. Makes prob of selection equivalent to fitness value. so-called "roulette wheel" selection
	{
	if(x->fitty[i] < randy && randy < x->fitty[i+1])
		{
		x->n = i;
		}
	}
	
randy = (fabs(Random() / MAXRAND));

while(x->fitty[x->n] < randy && randy < x->fitty[x->n+1])
	 randy = (fabs(Random() / MAXRAND));


for(i=0;i < x->population;i++)
	{
	if(x->fitty[i] < randy && randy < x->fitty[i+1])
		{
		x->m = i;
		}
	}


///end selection scheme 1 


//selection scheme 2 - a deterministic selection operator that picks best members. not currently implemented
//uses quicksort algorithm
//end selection scheme 2

}

void genalg_randomize(Genalg *x)
{
//shall be method for randomizing population
}

void genalg_replacemeth(Genalg *x, float f)		//number of replacement method to use
{
x->repnum = f;
//post("repnum = %f",x->repnum);

}

void genalg_replace(Genalg *x)
{
short j,i,l,k;

if (x->repnum == 1)	
	{
//replacement scheme 1 - replace parents with children, and corresponding fitness values.
	for (i=0;i < x->population;i++)
		{
		for(j=0;j < x->dimension;j++)
			{
			SETFLOAT(x->matrix+(i*x->dimension+j), x->matrix[(i+x->population)*x->dimension+j].a_w.w_float);
			}
		x->allfitness[i] = x->allfitness[i+x->population];
	   	 }
//end replacement scheme 1
	  }

if (x->repnum == 2)
	{
//replacement scheme 2	- compare fitness of children and parents, replace the weakest members
	for (j=0;j < 2*x->population;j++)
		{
		x->ascsort[j] = x->allfitness[j];
		//post("initial fitness list %d = %f",j,x->ascsort[j]);
		}
      
	
	quickSort(x->ascsort,0,2*x->population-1);	//sort the list of fitness values in ascending order, store in ascsort
	for (i=0;i < 2*x->population;i++)
	//post("sorted fitness list %d = %f",i,x->ascsort[i]);
	
	////////***store indices of ordered fitness values from ascsort***////////////
	for (i=0;i < 2*x->population;i++)
		{
		if (x->ascsort[0] == x->allfitness[i])
			{
			x->sortedindex[0] = i;
		//	post("sorted index 0 = %d",x->sortedindex[0]);
			}
		}

	
	for (j=1;j < 2*x->population;j++)		
		{
			for (i=0;i < 2*x->population;i++)
				{
				if (x->ascsort[j] == x->allfitness[i])
					x->sortedindex[j] = i;
				}
	//	post("sorted index %d = %d",j,x->sortedindex[j]);
		}
	///////////*******************************////////////////////
	
/*	for (i=0;i < 2*x->population;i++)
		{
		for (j=0;j < x->dimension;j++)
		{
			post("test matrix value %d,%d = %f",i,j,x->matrix[i*x->dimension + j].a_w.w_float);
		}*/
		
	/////******keep the x->population most fit members******////////
	for (i=0;i < x->population;i++)
		{
		k = x->sortedindex[2*x->population-1-i];
		//post("sort index value %d = %d",i,k);
		for (j=0;j < x->dimension;j++)
			{
			SETFLOAT(x->matrix+(i*x->dimension + j) , x->matrix[k*x->dimension + j].a_w.w_float);	//move selected members to front of the class
			x->allfitness[i] = x->allfitness[k];															//don't forget to take their fitness along!
		//	post("new generation parent % d,%d= %f from %d with fitness = %f",i,j,x->matrix[i*x->dimension + j].a_w.w_float,k,x->allfitness[i]);
			}
		
		}
	///////////////////////////////
	


//end replacement scheme 2


	}


}

/*void genalg_crossrate(Genalg *x, double f)	//number of iterations of selection/crossover to run.  not being used currently
{
int temp = (int)(x->population/2);		//upper bound is half the population size
x->t_crossrate = (int)(f*temp + 0.5);	//map double [0,1] to int [0,pop/2], with rounding.

if(x->t_crossrate > temp)
	x->t_crossrate = temp;	//in case we round up past maximum

//post("adjusted crossover rate = %d",x->t_crossrate);
}*/

double amax(double arg1, double arg2)
{
double res;
if(arg1 >= arg2)
	res = arg1;
else
	res = arg2;
	
return res;

}

double amin(double arg1, double arg2)
{
double res;
if(arg1 <= arg2)
	res = arg1;
else
	res = arg2;
	
return res;

}

void swap(double *i, double *j)
{
 double t;
 

    t=*i;
    *i=*j;
    *j=t;
    
}

double fit_dist(Genalg *x, short i, double *u, double *v)
{
short j;
double distnce = 0.0;
	for (j=0; j < x->dimension; j++)
		distnce += (u[j] - v[i*x->dimension + j]) * (u[j] - v[i*x->dimension + j]);
distnce = sqrt(distnce);
return(distnce);
}


///**********quicksort methods***********////////////

//int T;          /* type of item to be sorted */
//int tblIndex;   /* type of subscript */

//#define compGT(a,b) (a > b)

void insertSort(double *a, int lb, int ub) {
    double t;
    int i, j;

   /**************************
    *  sort array a[lb..ub]  *
    **************************/
    for (i = lb + 1; i <= ub; i++) {
        t = a[i];

        /* Shift elements down until */
        /* insertion point found.    */
        for (j = i-1; j >= lb && (a[j] > t); j--)
            a[j+1] = a[j];

        /* insert */
        a[j+1] = t;
    }
}

int partition(double *a, int lb, int ub) {
    double t, pivot;
    int i, j, p;

   /*******************************
    *  partition array a[lb..ub]  *
    *******************************/

    /* select pivot and exchange with 1st element */
    p = lb + ((ub - lb)>>1);
    pivot = a[p];
    a[p] = a[lb];

    /* sort lb+1..ub based on pivot */
    i = lb+1;
    j = ub;
    while (1) {
        while (i < j && (pivot > a[i])) i++;
        while (j >= i && (a[j] > pivot)) j--;
        if (i >= j) break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
        j--; i++;
    }

    /* pivot belongs in a[j] */
    a[lb] = a[j];
    a[j] = pivot;

    return j;
}

void quickSort(double *a, int lb, int ub) {
    int m;

   /**************************
    *  sort array a[lb..ub]  *
    **************************/

    while (lb < ub) {

        /* quickly sort short lists */
        if (ub - lb <= 12) {
            insertSort(a, lb, ub);
            return;
        }

        /* partition into two segments */
        m = partition (a, lb, ub);

        /* sort the smallest partition    */
        /* to minimize stack requirements */
        if (m - lb <= ub - m) {
            quickSort(a, lb, m - 1);
            lb = m + 1;
        } else {
            quickSort(a, m + 1, ub);
            ub = m - 1;
        }
    }
}

/////**********end quicksort methods*********//////////////






void genalg_clear(Genalg *x, Symbol *s, short msg, Atom *av)	//calls our memory free method
{
genalg_free(x);
}

void genalg_free(Genalg *x)		//free up memory
{
freebytes(x->matrix, (unsigned short)(2*sizeof(x->t_mutrate)*x->dimension*x->population));
freebytes(x->temparray, (unsigned short)(sizeof(x->t_mutrate)*x->dimension*x->population));
freebytes(x->outarray, (unsigned short)(sizeof(x->t_mutrate)*x->dimension));
freebytes(x->children, (unsigned short)(sizeof(x->t_mutrate)*x->dimension*x->population));
freebytes(x->allfitness, (unsigned short)(2*sizeof(x->t_mutrate)*x->population));
freebytes(x->fitty, (unsigned short)(sizeof(x->t_mutrate)*(x->population + 1)));
freebytes(x->fittest, (unsigned short)(sizeof(x->t_mutrate)*x->dimension));
freebytes(x->max, (unsigned short)(sizeof(x->t_mutrate)*x->dimension));
freebytes(x->min, (unsigned short)(sizeof(x->t_mutrate)*x->dimension));
freebytes(x->ascsort, (unsigned short)(2*sizeof(x->t_mutrate)*x->population));
freebytes(x->sortedindex, (unsigned short)(2*sizeof(x->nfit)*x->population));
}




void genalg_assist(Genalg *x, Object *b, long msg, long arg, char *s)		//assist messages for inlets/outlets
	{
		if (msg == ASSIST_INLET)
		{
		
			
		 sprintf(s, "%s", "various GA commands: mutate, mutrate, member list, etc.");
               
			
		}
		else if (msg == ASSIST_OUTLET)
				{
				sprintf(s, "%s", "(list) New Population Members");
				
				}
	}
