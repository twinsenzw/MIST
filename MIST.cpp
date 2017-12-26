#include <math.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <random>
#include <time.h>
#include <limits>
#include <map>
#include <thread>
/*
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
*/
using namespace::std;


//default_random_engine generator;
//uniform_real_distribution<double> distribution(0.0,1.0);

default_random_engine generator;
normal_distribution<double> distribution(0.0,0.01); //for mu

int random_int_in_range( int first, int last )
{
  /* This function implements the method recommended by the knowing
   * folks at comp.lang.c: http://c-faq.com/lib/randrange.html
   * Returns an integer in [first, last].
   */
  unsigned int N = (last - first <= RAND_MAX)  /* Make sure the algorithm    */
                 ? (last - first + 1U)         /* terminates by keeping N    */
                 : (RAND_MAX + 1U);            /* in rand()'s maximum range. */
  unsigned int x = (RAND_MAX + 1U) / N;
  unsigned int y = x * N;
  unsigned int r;
  do {
    r = rand();
  } while (r >= y);
  return r / x + first;
}

double random_double( int half_open )
{
  /* If half_open is TRUE, returns a random value in [0,1).
   * If half_open is FALSE, returns a random value in [0,1].
   * Remember that the result only has a resolution of RAND_MAX+1 values!      
   */
  return (double)rand() / (RAND_MAX + half_open);
}

double random_double_in_range( double low, double high )
{
  return random_double( 0 ) * (high - low) + low;
}



double log_binomial_ce(int n,int k) // n=dij0+dij1, k=dij0, p=pij
{
	int nk=n-k;
	if(k<nk) k=nk; //the functions runs n-k times, minimize n-k
	double sum=0.0;
	for(int i=1;i<=n-k;i++)
	{
		sum=sum+log(1+k/i);
	}

return sum;
}

class Alignment
{
public:
	vector<string> sequences;
	vector<string> taxaNames;
	vector<int> stateCounts;
	int maxStateCount;
	//state code for the sequences
	vector<vector<int>> counts;
	vector<int> patternWeight;
	vector<int> siteWeights;
	vector<vector<vector<double>>> tipLikelihoods; // #taxa x #sites x #states
	int usingTipLikelihoods = 0;
	vector<vector<int>> sitePatterns; // #patterns x #taxa
	vector<int> patternIndex;
	vector<int> excludedPatterns;
	vector<int> m_nIncluded;
	int isAscertained;

	void extractSitePattern(ifstream &input) //extract site pattern directly from a plain file of strain genomes in 0 and 1
	{
		sitePatterns.clear();
		string line;
		for(int i=0;getline(input,line);++i)
		{
			for(int j=0;j<line.size();j++)
			{
				if(i==0) //first line
				{
					vector<int> empty;
					sitePatterns.push_back(empty);
				}
				
				int nt=line[j]-'0';
				sitePatterns[j].push_back(nt);
			}
		}
	}	

	int getTaxonIndex(string id) //not the fastest
	{
		auto tni=find(taxaNames.begin(),taxaNames.end(),id);
		if(tni==taxaNames.end()) return -1;
		else return distance(taxaNames.begin(), tni);
	}

	int getPatternCount() 
	{
        	return sitePatterns.size();
	}

 	vector<int> getPattern(int patternIndex_) 
	{
        	return sitePatterns[patternIndex_];
	}

	int getPattern(int taxonIndex, int patternIndex_) 
	{
	        return sitePatterns[patternIndex_][taxonIndex];
	}

	vector<double> getTipLikelihoods(int taxonIndex, int patternIndex_) 
	{
		vector<double> empty;
		if (tipLikelihoods.size()==0) return empty;
    		if (taxonIndex >= tipLikelihoods.size() || tipLikelihoods[taxonIndex].size() == 0) 
		{ 
    			return empty; 
    		} 
		else 
		{ 
    			return tipLikelihoods[taxonIndex][patternIndex_];
    		}
    	
	}


	int stateCount=2;
	
	/*
	vector<vector<int>> mapCodeToStateSet;

	vector<int> zero={0};
	vector<int> one={1};

	mapCodeToStateSet.push_back(zero);
	mapCodeToStateSet.push_back(one);
	*/
	vector<int> getStatesForCode(int state) 
	{
		vector<vector<int>> mapCodeToStateSet; //should be class member

		vector<int> zero={0};
		vector<int> one={1};

		mapCodeToStateSet.push_back(zero);
		mapCodeToStateSet.push_back(one);
            	return mapCodeToStateSet[state]; 
	}

	vector<int> getStateSet(int state) 
	{
            	vector<int> stateSet(stateCount);
            	vector<int> stateNumbers = getStatesForCode(state);
            	for (int i : stateNumbers) stateSet[i] = 1;
            	return stateSet;
	} 
};


class Input
{
public:
	int sampleCount;
	int patternCount;

	vector<vector<int>> referenceCoverage; //#sample x #site
	vector<vector<int>> alternativeCoverage; //#sample x #site

	vector<vector<double>> logce;

	void calculateLogCE()
	{
		for(int samplei=0; samplei<sampleCount; samplei++)
		{
			vector<double> samplelogce;
			for(int sitei=0; sitei<patternCount; sitei++)
			{
				double sitelogce=0.0;
				sitelogce=log_binomial_ce(referenceCoverage[samplei][sitei]+alternativeCoverage[samplei][sitei],referenceCoverage[samplei][sitei]);
				//cout<<sitelogce<<"\t"<<referenceCoverage[samplei][sitei]+alternativeCoverage[samplei][sitei]<<"\t"<<referenceCoverage[samplei][sitei]<<endl;
				samplelogce.push_back(sitelogce);
			}
			logce.push_back(samplelogce);
		}
	} 
};

double LogPReconstruction(Input input, vector<vector<double>> strainFrequencies, Alignment alignment) 
{

	double logP=0;
	for(int samplei=0; samplei<input.sampleCount; samplei++)
	{
		for(int sitei=0; sitei<input.patternCount; sitei++)
		{
			double pref=0.0;
			for(int straini=0; straini<strainFrequencies[samplei].size(); straini++)
			{
				pref+=strainFrequencies[samplei][straini]*(1-alignment.sitePatterns[sitei][straini]);
			}
			
			logP=logP+input.logce[samplei][sitei]+input.referenceCoverage[samplei][sitei]*log(pref)+input.alternativeCoverage[samplei][sitei]*log(1-pref);
			
		}
	}
	return logP;
}

class State //build one new state object for each k
{
public:

	//Reconstruction
	Input input;
	int k; //total number of strain
	vector<vector<double>> strainFrequencies; //#sample x #strain
	vector<vector<double>> storedstrainFrequencies;
	Alignment alignment;
	Alignment storedalignment;

	//logP
	double logP;
	double storedlogP;

	double reconstructionLogP;
	double storedreconstructionLogP;

	//alpha 2
	double logalpha2=0;

	//memory for number of sites that have been changed
	double presentStdDev=0;
	int presentX;
	int isSequenceSample=0; //if not sequence sample, no need to update stddev
	int sequenceSampleSteps=0;
	

	void prepare()
	{
		;
	}


/*
	void dirichletStrainSampler()
	{
		
		int samplei=random_int_in_range(0,strainFrequencies.size()-1);
		strainFrequencies[samplei].clear();
		
		r2 = gsl_rng_alloc(gsl_rng_rand48);

		double sumk=0.0;
		double r[k]; //dirichelet parameter
		for(int ki=0; ki<k; ++ki) //ki: strain
		{
			r[ki]=0.0;
			for(int ji=0; ji<alignment.sitePatterns.size(); ji++) //ji: site
			{
				r[ki]=r[ki]+(1-alignment.sitePatterns[ji][ki])*input.referenceCoverage[samplei][ji]+alignment.sitePatterns[ji][ki]*input.alternativeCoverage[samplei][ji];
			}
			sumk=sumk+r[ki];
		}

		for(int ki=0; ki<k; ++ki)
		{
			if(sumk!=0) r[ki]=1.0+5.0*r[ki]/sumk;
			else r[ki]=1.0;
		}
		
		double sampleStrainFrequency[k];
		gsl_ran_dirichlet (r2, k, r, sampleStrainFrequency);
		for(int ki=0; ki<k; ki++)
		{
			strainFrequencies[samplei].push_back(sampleStrainFrequency[ki]);
		}
		
		//calc logalpha2
		double sum=0.0;
		for(int ki=0; ki<k; ++ki)
		{
			sum+=(r[ki]-1)*(log(strainFrequencies[samplei][ki])-log(storedstrainFrequencies[samplei][ki]));
		}
		logalpha2=sum;
	}
*/
	void uniformStrainSampler()
	{
		int straini=random_int_in_range(0,k-1);
		//strainFrequencies[][straini];
		
		vector<int> indicesMinus;
		for (int i=0; i<strainFrequencies.size(); ++i) 
		{
			if(strainFrequencies[i][straini]!=0.0) indicesMinus.push_back(i);
		}
		vector<int> indicesAdd;
		for (int i=0; i<strainFrequencies.size(); ++i) indicesAdd.push_back(i);
		random_shuffle(indicesAdd.begin(),indicesAdd.end());
		
		int nonZero=random_int_in_range(1,strainFrequencies.size());
		indicesAdd.resize(nonZero);

		for (int i=0; i<indicesMinus.size(); i++)
		{
			strainFrequencies[indicesMinus[i]][straini]=0.0;
			double sumFrequency=0.0;
			vector<int> needRecalc;
			for(int j=0; j<k; j++)
			{
				if(j!=straini && strainFrequencies[indicesMinus[i]][j]!=0.0) 
				{
					needRecalc.push_back(j);
					strainFrequencies[indicesMinus[i]][j]=random_int_in_range(1,100);
					sumFrequency+=strainFrequencies[indicesMinus[i]][j];
				}
			}
			for(int j=0; j<needRecalc.size();j++)
			{
				strainFrequencies[indicesMinus[i]][needRecalc[j]]=strainFrequencies[indicesMinus[i]][needRecalc[j]]/sumFrequency;
			}
		}


		for (int i=0; i<indicesAdd.size(); i++)
		{
			
			double sumFrequency=0.0;
			vector<int> needRecalc;
			for(int j=0; j<k; j++)
			{
				if(j==straini || strainFrequencies[indicesAdd[i]][j]!=0.0) 
				{
					needRecalc.push_back(j);
					strainFrequencies[indicesAdd[i]][j]=random_int_in_range(1,100);
					sumFrequency+=strainFrequencies[indicesAdd[i]][j];
				}
			}
			for(int j=0; j<needRecalc.size();j++)
			{
				strainFrequencies[indicesAdd[i]][needRecalc[j]]=strainFrequencies[indicesAdd[i]][needRecalc[j]]/sumFrequency;
			}
		}		
	}

	void uniformStrainSampler2() //directly changes the strain frequencies in a sample
	{
		int samplei=random_int_in_range(0,input.sampleCount-1);

		double sumFrequency=0.0;
		for (int straini=0; straini<k; straini++)
		{
			//strainFrequencies[indicesMinus[i]][straini]=0.0;
			strainFrequencies[samplei][straini]=random_double(0);
			sumFrequency=sumFrequency+strainFrequencies[samplei][straini];
		}
		for (int straini=0; sumFrequency!=0 && straini<k; straini++)
		{
			//strainFrequencies[indicesMinus[i]][straini]=0.0;
			strainFrequencies[samplei][straini]=strainFrequencies[samplei][straini]/sumFrequency;
		}
		
	}

	void uniformStrainSampler2_5() //directly changes the strain frequencies in a sample, specifying x # to be zero
	{
		int samplei=random_int_in_range(0,input.sampleCount-1);

		vector<int> nonZeroIndexes;
		for(int i=0;i<k;i++)
		{
			nonZeroIndexes.push_back(i);
		}
		
		random_shuffle(nonZeroIndexes.begin(),nonZeroIndexes.end());

		int nonZeroCount=random_int_in_range(1,k-1);

		double sumFrequency=0.0;
		for (int straini=0; straini<k; straini++)
		{
			strainFrequencies[samplei][straini]=0.0;
		}

		for (int straini=0; straini<nonZeroCount; straini++)
		{
			//strainFrequencies[indicesMinus[i]][straini]=0.0;
			strainFrequencies[samplei][nonZeroIndexes[straini]]=random_double(0);
			sumFrequency=sumFrequency+strainFrequencies[samplei][nonZeroIndexes[straini]];
		}

		
		for (int straini=0; sumFrequency!=0 && straini<nonZeroCount; straini++)
		{
			//strainFrequencies[indicesMinus[i]][straini]=0.0;
			strainFrequencies[samplei][nonZeroIndexes[straini]]=strainFrequencies[samplei][nonZeroIndexes[straini]]/sumFrequency;
			
		}
		
		
	}
	void uniformStrainSampler3() //fine tune by moving a proportion between two strains in one sample
	{
		int samplei=random_int_in_range(0,input.sampleCount-1);
		int fromi=random_int_in_range(0,k-1);
		int toi=random_int_in_range(0,k-1);
		for(;toi==fromi;toi=random_int_in_range(0,k-1)) ;

		double proportion=random_double(0)*0.1*strainFrequencies[samplei][fromi];
		strainFrequencies[samplei][fromi]=strainFrequencies[samplei][fromi]-proportion;
		strainFrequencies[samplei][toi]=strainFrequencies[samplei][toi]+proportion;
	}

	void bitflip()
	{
		int straini=random_int_in_range(0,k-1); //# strains
		int sitei=random_int_in_range(0,alignment.sitePatterns.size()-1); //# of sites
		alignment.sitePatterns[sitei][straini]=1-alignment.sitePatterns[sitei][straini];
	}

	void sequenceResampler()
	{
		int straini=random_int_in_range(0,k-1); //# strains


		int site1=random_int_in_range(0,alignment.sitePatterns.size()-1); //starting site
		int site2=random_int_in_range(site1,alignment.sitePatterns.size()-1); //ending site
		//for(;site2==site1;site2=random_int_in_range(0,alignment.sitePatterns.size()-1)) ;

		int fromi=max(site1,site2);
		int toi=min(site1,site2);

		for(int sitei=fromi;sitei<=toi;sitei++)
		{
			alignment.sitePatterns[sitei][straini]=random_int_in_range(0,1);
		}
	}

	void columnSorting1() //sort x number of columns
	{
		int x=random_int_in_range(1,alignment.sitePatterns.size());
		//what are the x sites?
		vector<int> siteIndex;
		for(int i=0;i<alignment.sitePatterns.size();i++)
		{
			siteIndex.push_back(i);
		}
		random_shuffle(siteIndex.begin(),siteIndex.end());

		for(int i=0; i<x ; i++)
		{
			random_shuffle(alignment.sitePatterns[siteIndex[i]].begin(),alignment.sitePatterns[siteIndex[i]].end());
		}
		int straini=random_int_in_range(0,k-1); //# strains
	}

	void columnSorting2() //sort 1 column x
	{
		int x=random_int_in_range(0,alignment.sitePatterns.size()-1);
		random_shuffle(alignment.sitePatterns[x].begin(),alignment.sitePatterns[x].end());
		int straini=random_int_in_range(0,k-1); //# strains
	}

	void adaptiveSequenceResampler(int step)
	{	
		double stddev=0;
		if(step<5000)
		{
			stddev=input.patternCount*0.01;
		}
		else
		{
			stddev=presentStdDev;
		}
		normal_distribution<double> Xdistribution(0.0,stddev); //for mu
		presentX=Xdistribution(generator);
		if(presentX<0) presentX=-presentX;

		map<int,vector<int>> flipped; //#strains -> flipped sites
		for(int i=0;i<presentX;i++) //bitflip i random sites
		{
			int strainID=random_int_in_range(0,k-1);
			int siteID=random_int_in_range(0,input.patternCount-1);
			/*
			if(flipped.find(strainID)!=flipped.end())
			{
				for(;find(flipped[strainID].begin(),flipped[strainID].end(),siteID)!=flipped[strainID].end();)
				{
					strainID=random_int_in_range(0,k-1);
					siteID=random_int_in_range(0,input.patternCount-1);
					if(flipped.find(strainID)==flipped.end()) break;
				}
			}*/
			alignment.sitePatterns[siteID][strainID]=1-alignment.sitePatterns[siteID][strainID];
			flipped[strainID].push_back(siteID);
		}

	}

	void computePresentStdDev(int step)
	{
		if(step!=0) presentStdDev=pow((presentStdDev*presentStdDev*(step-1) + presentX * presentX) / step, 0.5);
		else presentStdDev=0;
	}
	void operateLite(ofstream &output, int step, int stepInput) //one parameter at a time
	{
		logalpha2=0;
		double updateDrawProbability=random_double(0);
		if(step<stepInput*0.2)
		{
			/*
			if(updateDrawProbability<=0.125) 
			{
				isSequenceSample=0;
				dirichletStrainSampler(); //update pool
			}

			*/
			if(updateDrawProbability<=0.33)
			{
				isSequenceSample=0;
				uniformStrainSampler2(); //update pool
			}
			
			/*
			else if(updateDrawProbability>0.45 && updateDrawProbability<=0.5) 
			{
				uniformStrainSampler3(); //update pool
				isSequenceSample=0;
			}*/

			
			else if(updateDrawProbability>0.33&&updateDrawProbability<=0.43) //bitflip the alignment
			{
				isSequenceSample=0;
				bitflip();
			}

			
			else if(updateDrawProbability>0.43&&updateDrawProbability<=0.70) //resample a whole strain
			{
				isSequenceSample=0;
				sequenceResampler();
			}

			else if(updateDrawProbability>0.75&&updateDrawProbability<=0.9) //sort columns
			{
				isSequenceSample=0;
				columnSorting1();
			}

			else if(updateDrawProbability>0.9) //sort one column
			{
				isSequenceSample=0;
				columnSorting2();
			}
			/*
			else if(updateDrawProbability>0.75) //resample a whole strain
			{
				isSequenceSample=1;
				sequenceSampleSteps++;
				adaptiveSequenceResampler(sequenceSampleSteps);
			}*/
		}
		else
		{
			/*
			if(updateDrawProbability<=0.125) 
			{
				isSequenceSample=0;
				dirichletStrainSampler(); //update pool
			}*/

			if(updateDrawProbability<=0.15)
			{
				isSequenceSample=0;
				uniformStrainSampler2(); //update pool
			}
			/*
			if(updateDrawProbability>0.25 && updateDrawProbability<=0.5) 
			{
				uniformStrainSampler(); //update pool
				isSequenceSample=0;
			}
			*/
			else if(updateDrawProbability>0.15 && updateDrawProbability <= 0.5) //bitflip the alignment
			{
				isSequenceSample=0;
				bitflip();
			}

			else if(updateDrawProbability>0.5 && updateDrawProbability <= 0.70) //resample a whole strain
			{
				isSequenceSample=0;
				sequenceResampler();
			}
			else if(updateDrawProbability>0.70 && updateDrawProbability <= 0.9) //column sorting
			{
				isSequenceSample=0;
				columnSorting2(); 
			}	
			else if(updateDrawProbability>0.9) //column sorting
			{
				isSequenceSample=0;
				columnSorting1(); 
			}			
		}
	}
	
	void score()
	{
		
		logP=0.0;
			
		reconstructionLogP=LogPReconstruction(input, strainFrequencies, alignment);
		logP=reconstructionLogP+logalpha2;
		//cout<<coalescent.calculateLogP()<<" "<<treeLikelihood.logP<<" "<<reconstructionLogP<<" "<<logalpha2<<endl;
	}

	void store(int step)
	{
		storedalignment=alignment;
		storedlogP=logP;
		storedreconstructionLogP=reconstructionLogP;
		storedstrainFrequencies=strainFrequencies;
		if(isSequenceSample) computePresentStdDev(sequenceSampleSteps);
	}

	void restore(int step)
	{
		
		alignment=storedalignment;
		logP=storedlogP;
		reconstructionLogP=storedreconstructionLogP;
		strainFrequencies=storedstrainFrequencies;
		if(isSequenceSample) 
		{
			presentX=0;
			computePresentStdDev(sequenceSampleSteps);
		}
	}
		

};

class MCMC
{
public:

	State state;
	Input input;

	int chainLength;
	int storeEvery; //store state every x steps
	int burnIn;
	int presentStep=0;
	
	//heated chain parameters
	int isCold=0;
	double temperature;

	void Log(ofstream &output)
	{
		output<<presentStep<<"\t";
		output<<state.logP<<"\t";

		for(int i=0;i<state.alignment.sitePatterns[0].size();++i)
		{
			output<<"{>strain"<<i<<" ";
			for(int j=0;j<state.alignment.sitePatterns.size();++j)
			{
				if(state.alignment.sitePatterns[j][i]) output<<'A';
				else if(!state.alignment.sitePatterns[j][i]) output<<'T';
			}
			output<<"}";
		}

		output<<"\t";

		for(int i=0;i<state.strainFrequencies.size();i++)
		{
			output<<"[";
			int isFirst=1;
			for(int j=0; j<state.strainFrequencies[i].size(); j++)
			{
				if(isFirst!=1) output<<",";
				else isFirst=0;
				output<<state.strainFrequencies[i][j];
			}
			output<<"]";
		}
		output<<endl;

		//output the latest state
		ofstream frequencyout("latestStrainFrequency.txt");
		for(int i=0;i<state.strainFrequencies.size();i++)
		{
			int isFirst=1;
			for(int j=0; j<state.strainFrequencies[i].size(); j++)
			{
				if(isFirst!=1) frequencyout<<"\t";
				else isFirst=0;
				frequencyout<<state.strainFrequencies[i][j];
			}
			frequencyout<<endl;
		}

		ofstream sequenceout("latestStrainSequence.txt");
		sequenceout<<state.alignment.sitePatterns[0].size()<<" "<<state.alignment.sitePatterns.size()<<endl;
			
		for(int i=0;i<state.alignment.sitePatterns[0].size();++i)
		{
			sequenceout<<">strain"<<i<<"\n";
			for(int j=0;j<state.alignment.sitePatterns.size();++j)
			{
				if(state.alignment.sitePatterns[j][i]) sequenceout<<'A';
				else if(!state.alignment.sitePatterns[j][i]) sequenceout<<'T';
			}
			sequenceout<<"\n";
		}			
	}
	
	void InitializeChain(int chainlength, int storeevery, int burnin)
	{
		chainLength=chainlength;
		storeEvery=storeevery;
		burnIn=burnin;
	}


	void InitializeState(int k, int Ne)
	{
		int patternCount=input.patternCount; //# of sites

		state.alignment.sitePatterns.resize(patternCount);
		state.k=k;

		//----Initialize state::alignment----
		
		for(int i=0;i<k;++i)
		{
			ostringstream taxonNameStream;

			//taxa names
			taxonNameStream << i;
			state.alignment.taxaNames.push_back(taxonNameStream.str());

			//patterns
			
			for(int pi=0;pi<patternCount; pi++)
			{
				int ran=random_int_in_range(1,100);
				int pattern;
				if(ran<=50) pattern=0;
				else pattern=1;
				state.alignment.sitePatterns[pi].push_back(pattern);
			}
			
			//ifstream sequenceInput("Pattern_input");
			//state.alignment.extractSitePattern(sequenceInput);


		}

		//if one site is not polymorphic, bitflip a random strain
		
		for(int pi=0;pi<patternCount; pi++)
		{
			int isPoly=0;
			int startPattern=state.alignment.sitePatterns[pi][0];
			for(int stri=0; stri<k; stri++)
			{
				if(state.alignment.sitePatterns[pi][stri]!=startPattern)
				{
					isPoly=1;
					break;
				}
			}
			if(!isPoly)
			{
				
				int ran=random_int_in_range(0,k-1);
				state.alignment.sitePatterns[pi][ran]=1-state.alignment.sitePatterns[pi][ran];
			}
		}		
		
		
		//----Initialize strainfrequencies----
			
		
			state.strainFrequencies.clear();
			vector<double> empty(k,0.0);
			for(int i=0;i<input.sampleCount; i++)
			{
				state.strainFrequencies.push_back(empty);
			}
					
			
			//how many samples have each strain?  //***v2: now every strain has nonzero freq in each sample
			vector<int> NrOfsamples(k);
			for(int stri=0;stri<k; stri++) 
			{
				//NrOfsamples[stri]=(rand()%input.sampleCount)+1;
				NrOfsamples[stri]=input.sampleCount;
			}
			//what samples have each strain?
			vector<int> sampleIndices;
			for(int smpli=0;smpli<input.sampleCount; smpli++)
			{
				sampleIndices.push_back(smpli);
			}
			
			//what frequencies do these samples have each strain?
			for(int stri=0;stri<k; stri++)
			{
				//random_shuffle(sampleIndices.begin(),sampleIndices.end()); //now not needed
				for(int smpli=0;smpli<NrOfsamples[stri]; smpli++)
				{
					state.strainFrequencies[sampleIndices[smpli]][stri]=(double)(random_int_in_range(1,100));
					
				}
			}
			//normalize
			
			for(int smpli=0;smpli<input.sampleCount; smpli++)
			{
				double sum=0.0;
				for(int stri=0; stri<k; stri++) 
				{			
					sum=sum+state.strainFrequencies[smpli][stri];
				}
				for(int stri=0; stri<k; stri++) 
				{			
					state.strainFrequencies[smpli][stri]=state.strainFrequencies[smpli][stri]/sum;
					
				}
			}		
			

		//----score and store state----
		
		state.score();
		
		state.store(presentStep);
	}
	
	void run(ofstream &output,int stepInput)
	{
		for(int step=0;step<chainLength;++step)
		{
			//cout<<"Present step: "<<step<<"..."<<endl;

			presentStep++;

			//propose new state
			state.operateLite(output,presentStep,stepInput);
			//score new state
			state.score();
			//propagate probability
			double propP;
			
			if (std::isinf(state.storedlogP)||std::isnan(state.storedlogP)) propP=1; //if the last state is inf, propagate
			else
			{
				if(std::isinf(state.logP)||std::isnan(state.logP)) propP=-INFINITY; //if the last state is not inf, the next state is, do not propagate
				else propP=min(0.0,temperature*(state.logP-state.storedlogP));
			}
			//propagate?
			
			double prop=random_double(0);
			double logprop=log(prop);
			for(;std::isinf(logprop);)
			{
				prop=random_double(0);
				logprop=log(prop);
			}
			//cout<<"hastings: "<<logprop<<" "<<propP<<endl;
			//cout<<state.treeLikelihood.branchRateModel.mu<<endl;
			if(logprop<=propP)
			{
				state.store(presentStep);
				
			}
			else
			{
				state.restore(presentStep);
				
			}
			//log
			
			if(step%storeEvery==0) Log(output);
		}
	}
		
		
};

void getInput(ifstream &referenceCount, ifstream &alternativeCount, Input &input)
{
	cout<<"Reading inputs..."<<endl;

	string line;
	int patternCount=0;
	int sampleCount=0;
	//reading referenceCounts
	
	vector<int> vempty;
	for(int isFirst=1;getline(referenceCount,line);)
	{
		istringstream linestream(line);
		string field;
		
		for(int samplei=0;getline(linestream,field,'\t');samplei++)
		{
			
			if(isFirst)
			{
				input.referenceCoverage.push_back(vempty);
				sampleCount++;
			}
			int count=stoi(field);
			input.referenceCoverage[samplei].push_back(count);
		}
		patternCount++;
		if(isFirst) isFirst=0;
	}

	//reading alternativeCounts
	for(int isFirst=1;getline(alternativeCount,line);)
	{
		istringstream linestream(line);
		string field;
		
		for(int samplei=0;getline(linestream,field,'\t');samplei++)
		{
			if(isFirst)
			{
				input.alternativeCoverage.push_back(vempty);
			}
			int count=stoi(field);
			input.alternativeCoverage[samplei].push_back(count);
		}
			
		if(isFirst) isFirst=0;
	}

	input.sampleCount=sampleCount;
	cout<<"patternCount:" << patternCount<<endl;
	input.patternCount=patternCount;
}

	
void run_mcmc(MCMC& mcmc, ofstream& out, int stepInput) //for multithreading
{
	mcmc.run(out,stepInput);
}


void MC3(Input input, int kInput, int stepInput)
{
	const int chainNumber=4;
	const double deltaT=0.15;
	const int swapPeriod=2000;
	const int chainLength=stepInput;
	const int swapChainRounds=chainLength/swapPeriod;
	const int storeEvery=2000;
	const int burnIn=stepInput*0.33;
	const int k=kInput;
	const int N=100;

	thread t[chainNumber];

	//initialize the mcmc objects
	cout<<"Initializing "<<chainNumber<<" MCMC chains, with deltaT="<<deltaT<<", swap period="<<swapPeriod<<endl;
	vector<MCMC> mcmc(chainNumber);

	for(int cni=0;cni<chainNumber;cni++)
	{
		mcmc[cni].input=input;
		mcmc[cni].state.input=input;
		
		mcmc[cni].InitializeChain(swapPeriod, storeEvery, burnIn); //each time only run for swapPeriod steps
		
		mcmc[cni].InitializeState(k, N);
		mcmc[cni].temperature=1.0/(1.0+deltaT*(double)cni);	
	}

	mcmc[0].isCold=1;

	//initialize the output files

	// the following write one output for each chain
	
	vector<ofstream> streams;
	for (int i=0; i<chainNumber; i++)
	{
  		string fileName = "log_chain" + to_string(i) + ".txt";
  		ofstream out(fileName.c_str());
		streams.push_back(move(out));
	}

	// the following write output for the cold chain
	ofstream outCold("log_cold.txt");

	int step=0;
	for(int round=0;round<swapChainRounds;round++)
	{
		
		//run MCMC on different threads
		for(int cni=0;cni<chainNumber;cni++)
		{
			if(mcmc[cni].isCold==1) t[cni] = thread(run_mcmc,std::ref(mcmc[cni]),std::ref(outCold),stepInput);
			else t[cni]=thread(run_mcmc,std::ref(mcmc[cni]),std::ref(streams[cni]),stepInput);
		}

		//join
		for(int cni=0;cni<chainNumber;cni++)
		{
			t[cni].join();
		}
		
		//swap
		int swap1=random_int_in_range(0,chainNumber-1);
		int swap2=swap1;
		for(;swap2==swap1;) swap2=random_int_in_range(0,chainNumber-1);
		cout<<"Proposing to swap chains "<<swap1<<" and "<<swap2<<"...";

		//swap or not?
		double logSwapProbability=(mcmc[swap1].temperature-mcmc[swap2].temperature)*(mcmc[swap2].state.logP-mcmc[swap1].state.logP);
		//double sprop=(double) rand()/RAND_MAX;
		double sprop=random_double(0);
		double logsprop=log(sprop);
		for(;std::isinf(logsprop);)
		{
			sprop=random_double(0);
			logsprop=log(sprop);
		}
		cout<<"\t"<<logsprop<<"\t"<<logSwapProbability<<endl;
		if(logsprop<=logSwapProbability)
		{
			cout<<"Proposal accepted"<<endl;
			double tmp=mcmc[swap1].temperature;
			mcmc[swap1].temperature=mcmc[swap2].temperature;
			mcmc[swap2].temperature=tmp;

			if(mcmc[swap1].isCold==1) 
			{
				mcmc[swap1].isCold=0;
				mcmc[swap2].isCold=1;
			}
			else if(mcmc[swap2].isCold==1)
			{
				mcmc[swap2].isCold=0;
				mcmc[swap1].isCold=1;
			}
		}
		else
		{
			cout<<"Proposal rejected"<<endl;
		}

		for(int coldi=0;coldi<chainNumber;coldi++)
		{
			if(mcmc[coldi].isCold==1) 
			{
				cout<<"\tPresent cold chain: "<<coldi<<", LogP = "<<mcmc[coldi].state.logP<<endl;
				break;
			}
		}
		
	}
}


void MC2(Input input,int kInput, int stepInput)
{

	const int chainLength=stepInput;
	const int storeEvery=2000;
	const int burnIn=stepInput*0.33;
	
	const int k=kInput;
	const int N=100;

	ofstream out("log.txt");

	cout<<"Initializing MCMC chain..."<<endl;
	
	MCMC mcmc;

	mcmc.input=input;
	mcmc.state.input=input;
	mcmc.InitializeChain(chainLength, storeEvery, burnIn); 
	mcmc.InitializeState(k, N);
	mcmc.temperature=1.0;
	
	mcmc.run(out,stepInput);
}

int main(int argc, char * argv[]) // 1-k, 2-step, 3-refinput, 4-altinput
{
	string arg1=argv[1];
	int k=stoi(arg1);
	string arg2=argv[2];
	int step=stoi(arg2);
	string arg3=argv[3];
	string arg4=argv[4];
	srand(time(0));

	//=====================load the test example into input=========================

	Input input;
	ifstream referenceInput(arg3.c_str());
	ifstream alternativeInput(arg4.c_str());
	getInput(referenceInput, alternativeInput, input);
	input.calculateLogCE();

	//=====================use single Chain=====================

	MC3(input,k,step);
		
				
}
