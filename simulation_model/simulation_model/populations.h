//this was mostly pointers before, but I'm going to change it to vectors
#pragma once
#include "rand_nums.h"
#include "chi_square.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <vector>
#include <array>


using namespace std;

class LociTracker
{
public:
	double* LociOnChrom;
};//end LociTracker

class ld_info
{
public:
	double d;
	double dprime;
};

class AlleleFrequencies
{
public:
	LociTracker* allele_freq;
};

class Fst
{
public:
	LociTracker* FstValue;
	LociTracker* WeightedFst;
	LociTracker* PValue;
	LociTracker* BS_Pvalue;
	AlleleFrequencies* Adult_AlleleFreq;
	AlleleFrequencies* Progeny_AlleleFreq;
	LociTracker* Locus;
};//end Fst

class Chromosome
{
public:
	int* LociArray;
	double* allelicEffects;
};//end of Class

class Individual
{
public:
	Chromosome* maternal;//create array of 4 chromosomes
	Chromosome* paternal;
	double phenotype;
	double Genotype;
	bool Alive;
	bool IsFemale;
	int MateFound;
};//end Individual


class Population
{
public:
	//all the types
	int PopulationSize;
	int CarryingCapacity;
	int MaximumFecundity;
	int MaximumEncounters;
	int ProgenyNum, NumMales, NumFemales;
	int ProgFem, ProgMale;
	int TotalLociNo;
	int NumSelGen;
	int NumPolymorphicLoci;
	int MOnumPolyLoci;
	int peak_window;
	double GaussianPrefMean;
	string basefilename;
	Individual* Adult;
	Individual* Progeny;
	int* ChosenSexes;
	int* Mothers;
	int* Fathers;
	int* Nonmated;
	Fst F_Marker;
	Fst F_QTL;
	LociTracker* Polymorphisms;
	LociTracker* OPolymorphisms;
	LociTracker* Locations;
	LociTracker* QTLtracker;
	LociTracker* MarkerLoci;
	LociTracker* QAdultHet;
	LociTracker* QMalHet;
	LociTracker* QOffHet;
	LociTracker* QTotHet;
	LociTracker* QHs;
	LociTracker* QHt;
	LociTracker* MalHet;
	LociTracker* OffHet;
	LociTracker* TotHet;
	LociTracker* HetS;
	LociTracker* AdultHet;
	LociTracker* MaleH;
	LociTracker* AMajor;
	LociTracker* OMajor;
	LociTracker* MMajor;
	LociTracker* AFreq1;
	LociTracker* OFreq1;
	LociTracker* MFreq1;
	LociTracker* ALD;
	LociTracker* OLD;
	LociTracker* MLD;
	LociTracker* A_NumAlleles;
	LociTracker* M_NumAlleles;
	LociTracker* O_NumAlleles;
	LociTracker* PopPrimFreq;
	LociTracker* PopMajorAllele;
	LociTracker* PopLD_a;
	LociTracker* Psig;
	LociTracker* GenetVars;
	LociTracker* GenetMeans;
	LociTracker* SampledLoci;
	LociTracker* ANumQTLAlleles;
	LociTracker* ONumQTLAlleles;
	LociTracker* MNumQTLAlleles;
	LociTracker* TNumQTLAlleles;
	LociTracker* MTNumQTLAlleles;
	LociTracker* FstTracker;
	int* SampledAdults;
	int* SampledOffspring;
	double StartingP, pPprime, EffectPop;
	double MutationRate, MutationalVariance;
	double mMaleTrait, mFemaleTrait;
	double males, females;
	double adultRatio, offRatio;
	double AvgLD, CriticalValue;
	int bs_spurious, bs_detected, ci_spurious, ci_detected, bs_spurious95, bs_detected95, fdr_spur, fdr_det, ci_spurious95, ci_detected95;
	double SmoothedCritValue99, SmoothedCritValue98, SmoothedCritValue95, SmoothedCritValue90, SmoothedCritValue80;
	double bsCritValue99, bsCritValue98, bsCritValue95, bsCritValue90, bsCritValue80, FDR95;
	int MaxNumProg;
	int Mated;
	int NumNotMated;
	int NumSampledLoci, TotalSampled;
	int NumChrom, NumQTLs, NumMarkers, NumAlleles;
	int Counter1, Counter2;
	int FemSample, MaleSample, ProgSample, pAdultSample;
	bool popExtinct;
	double meanLD;
	double RecombRate, SigValue;
	double EnvironmentalVariance, EnvStdDev;
	//for writing to files
	int malN;
	int femN;
	int Generations;
	double allphenavg;
	double malphenavg;
	double femphenavg;
	double allDiff;
	double maleDiff;
	double femDiff;
	double Heritability;
	double GenVar;
	double PhenVar;
	double allelefreq;
	int bs_replicates;
	int NumLD;
	bool bootstrap, prior_qtl;
	double avg_pairwise_d, avg_longdist_d;
	vector <double> avg_d;
	vector <double> long_dist_dprime;
	double sigma;//need to make sure I've got enough loci for this

public: //Population Functions

	void SetParameters()
	{
		//Initialize Parameters
		CarryingCapacity = 5000;
		MaximumFecundity = 4;
		MaxNumProg = CarryingCapacity*MaximumFecundity;
		Generations = 200;
		Adult = new Individual[CarryingCapacity];
		Progeny = new Individual[MaxNumProg];
		ChosenSexes = new int[MaxNumProg];
		Mothers = new int[MaxNumProg];
		Fathers = new int[MaxNumProg];
		Nonmated = new int[CarryingCapacity];
		NumChrom = 4;
		NumMarkers = 1000;
		NumQTLs = 2;
		GaussianPrefMean = 4.0;
		NumAlleles = 6;
		MaximumEncounters = 50;
		NumSampledLoci = NumMarkers;
		RecombRate = 0.2;
		MutationRate = 0;
		MutationalVariance = 0;
		TotalLociNo = NumMarkers*NumChrom;
		EnvironmentalVariance = 0;
		EnvStdDev = sqrt(EnvironmentalVariance);
		NumSelGen = 1;
		sgenrand(time(0));
		prior_qtl = false;
		popExtinct = false;
		TotalSampled = NumSampledLoci*NumChrom;
		pAdultSample = 440;//set sample size here
		ProgSample = 160;//if 0, it's set to ProgenyNum in samplepop 
		SigValue = 0.05;
		sigma = 15;//weighted sigma
		bootstrap = true;
		if (bootstrap == true)
			bs_replicates = 10000;
		NumLD = 2000;
		peak_window = 50;//the window that you can detect a QTL near a peak (in either direction)
	}

	void Initialize(bool known_qtls)
	{
		int qtl_list[] = { 444, 82, 393, 218, 699, 477, 603, 694 };
		int s, iii, jjj, lll, k, eee;
		int ch, nq, m;
		double dNumQTLS, AllelicStdDev;
		mMaleTrait = 0;
		mFemaleTrait = 0;
		NumMales = 0;
		NumFemales = 0;
		ProgenyNum = 0;
		dNumQTLS = NumQTLs;
		AllelicStdDev = 0.5;
		double rand;
		double prob, dLoci, dSample;
		dLoci = NumMarkers;
		dSample = NumSampledLoci;
		prob = dSample / dLoci;
		int counter;
		//allocate memory
		for (iii = 0; iii < CarryingCapacity; iii++)
		{
			Adult[iii].maternal = new Chromosome[NumChrom];
			Adult[iii].paternal = new Chromosome[NumChrom];

			for (jjj = 0; jjj < NumChrom; jjj++)
			{
				Adult[iii].maternal[jjj].LociArray = new int[NumMarkers];
				Adult[iii].maternal[jjj].allelicEffects = new double[NumQTLs];

				Adult[iii].paternal[jjj].LociArray = new int[NumMarkers];
				Adult[iii].paternal[jjj].allelicEffects = new double[NumQTLs];
			}

		}
		for (iii = 0; iii < MaxNumProg; iii++)
		{
			Progeny[iii].maternal = new Chromosome[NumChrom];
			Progeny[iii].paternal = new Chromosome[NumChrom];

			for (jjj = 0; jjj < NumChrom; jjj++)
			{
				Progeny[iii].maternal[jjj].LociArray = new int[NumMarkers];
				Progeny[iii].maternal[jjj].allelicEffects = new double[NumQTLs];

				Progeny[iii].paternal[jjj].LociArray = new int[NumMarkers];
				Progeny[iii].paternal[jjj].allelicEffects = new double[NumQTLs];
			}
		}

		F_Marker.BS_Pvalue = new LociTracker[NumChrom];
		F_Marker.FstValue = new LociTracker[NumChrom];
		F_Marker.PValue = new LociTracker[NumChrom];
		F_Marker.WeightedFst = new LociTracker[NumChrom];
		F_Marker.Adult_AlleleFreq = new AlleleFrequencies[NumChrom];
		F_Marker.Progeny_AlleleFreq = new AlleleFrequencies[NumChrom];
		F_Marker.Locus = new LociTracker[NumChrom];

		F_QTL.BS_Pvalue = new LociTracker[NumChrom];
		F_QTL.FstValue = new LociTracker[NumChrom];
		F_QTL.PValue = new LociTracker[NumChrom];
		F_QTL.WeightedFst = new LociTracker[NumChrom];
		F_QTL.Adult_AlleleFreq = new AlleleFrequencies[NumChrom];
		F_QTL.Progeny_AlleleFreq = new AlleleFrequencies[NumChrom];
		F_QTL.Locus = new LociTracker[NumChrom];

		SampledLoci = new LociTracker[NumChrom];
		Polymorphisms = new LociTracker[NumChrom];
		OPolymorphisms = new LociTracker[NumChrom];
		Locations = new LociTracker[NumChrom];
		QTLtracker = new LociTracker[NumChrom];
		MarkerLoci = new LociTracker[NumChrom];
		PopPrimFreq = new LociTracker[NumChrom];
		PopLD_a = new LociTracker[NumChrom];
		PopMajorAllele = new LociTracker[NumChrom];
		QAdultHet = new LociTracker[NumChrom];
		QHs = new LociTracker[NumChrom];
		QHt = new LociTracker[NumChrom];
		QMalHet = new LociTracker[NumChrom];
		QOffHet = new LociTracker[NumChrom];
		QTotHet = new LociTracker[NumChrom];
		MalHet = new LociTracker[NumChrom];
		OffHet = new LociTracker[NumChrom];
		TotHet = new LociTracker[NumChrom];
		HetS = new LociTracker[NumChrom];
		MaleH = new LociTracker[NumChrom];
		AdultHet = new LociTracker[NumChrom];
		AFreq1 = new LociTracker[NumChrom];
		MFreq1 = new LociTracker[NumChrom];
		OFreq1 = new LociTracker[NumChrom];
		AMajor = new LociTracker[NumChrom];
		MMajor = new LociTracker[NumChrom];
		OMajor = new LociTracker[NumChrom];
		ALD = new LociTracker[NumChrom];
		OLD = new LociTracker[NumChrom];
		MLD = new LociTracker[NumChrom];
		A_NumAlleles = new LociTracker[NumChrom];
		M_NumAlleles = new LociTracker[NumChrom];
		O_NumAlleles = new LociTracker[NumChrom];
		Psig = new LociTracker[NumChrom];
		GenetVars = new LociTracker[NumChrom];
		GenetMeans = new LociTracker[NumChrom];
		ANumQTLAlleles = new LociTracker[NumChrom];
		ONumQTLAlleles = new LociTracker[NumChrom];
		MNumQTLAlleles = new LociTracker[NumChrom];
		TNumQTLAlleles = new LociTracker[NumChrom];
		MTNumQTLAlleles = new LociTracker[NumChrom];
		FstTracker = new LociTracker[NumChrom];
		for (s = 0; s < NumChrom; s++){
			FstTracker[s].LociOnChrom = new double[NumSampledLoci];
			SampledLoci[s].LociOnChrom = new double[NumMarkers];
			GenetVars[s].LociOnChrom = new double[NumQTLs];
			Polymorphisms[s].LociOnChrom = new double[NumSampledLoci];
			OPolymorphisms[s].LociOnChrom = new double[NumSampledLoci];
			Locations[s].LociOnChrom = new double[NumQTLs];
			QTLtracker[s].LociOnChrom = new double[NumMarkers];
			MalHet[s].LociOnChrom = new double[NumSampledLoci];
			OffHet[s].LociOnChrom = new double[NumSampledLoci];
			AdultHet[s].LociOnChrom = new double[NumSampledLoci];
			HetS[s].LociOnChrom = new double[NumSampledLoci];
			TotHet[s].LociOnChrom = new double[NumSampledLoci];
			MaleH[s].LociOnChrom = new double[NumSampledLoci];
			A_NumAlleles[s].LociOnChrom = new double[NumSampledLoci];
			M_NumAlleles[s].LociOnChrom = new double[NumSampledLoci];
			O_NumAlleles[s].LociOnChrom = new double[NumSampledLoci];
			ALD[s].LociOnChrom = new double[NumMarkers];
			OLD[s].LociOnChrom = new double[NumMarkers];
			MLD[s].LociOnChrom = new double[NumMarkers];
			PopPrimFreq[s].LociOnChrom = new double[NumMarkers];
			PopLD_a[s].LociOnChrom = new double[NumMarkers];
			PopMajorAllele[s].LociOnChrom = new double[NumMarkers];
			AFreq1[s].LociOnChrom = new double[NumSampledLoci];
			OFreq1[s].LociOnChrom = new double[NumSampledLoci];
			MFreq1[s].LociOnChrom = new double[NumSampledLoci];
			AMajor[s].LociOnChrom = new double[NumSampledLoci];
			OMajor[s].LociOnChrom = new double[NumSampledLoci];
			MMajor[s].LociOnChrom = new double[NumSampledLoci];
			Psig[s].LociOnChrom = new double[NumSampledLoci];
			GenetMeans[s].LociOnChrom = new double[NumQTLs];
			MarkerLoci[s].LociOnChrom = new double[NumMarkers];
			ANumQTLAlleles[s].LociOnChrom = new double[NumQTLs];
			ONumQTLAlleles[s].LociOnChrom = new double[NumQTLs];
			MNumQTLAlleles[s].LociOnChrom = new double[NumQTLs];
			TNumQTLAlleles[s].LociOnChrom = new double[NumQTLs];
			MTNumQTLAlleles[s].LociOnChrom = new double[NumQTLs];
			QAdultHet[s].LociOnChrom = new double[NumQTLs];
			QHs[s].LociOnChrom = new double[NumQTLs];
			QHt[s].LociOnChrom = new double[NumQTLs];
			QMalHet[s].LociOnChrom = new double[NumQTLs];
			QOffHet[s].LociOnChrom = new double[NumQTLs];
			QTotHet[s].LociOnChrom = new double[NumQTLs];
			F_QTL.BS_Pvalue[s].LociOnChrom = new double[NumQTLs];
			F_QTL.FstValue[s].LociOnChrom = new double[NumQTLs];
			F_QTL.PValue[s].LociOnChrom = new double[NumQTLs];
			F_QTL.WeightedFst[s].LociOnChrom = new double[NumQTLs];
			F_QTL.Adult_AlleleFreq[s].allele_freq = new LociTracker[NumQTLs];
			F_QTL.Progeny_AlleleFreq[s].allele_freq = new LociTracker[NumQTLs];
			F_QTL.Locus[s].LociOnChrom = new double[NumQTLs];
			F_Marker.BS_Pvalue[s].LociOnChrom = new double[NumSampledLoci];
			F_Marker.FstValue[s].LociOnChrom = new double[NumSampledLoci];
			F_Marker.PValue[s].LociOnChrom = new double[NumSampledLoci];
			F_Marker.WeightedFst[s].LociOnChrom = new double[NumSampledLoci];
			F_Marker.Locus[s].LociOnChrom = new double[NumSampledLoci];
			F_Marker.Adult_AlleleFreq[s].allele_freq = new LociTracker[NumSampledLoci];
			F_Marker.Progeny_AlleleFreq[s].allele_freq = new LociTracker[NumSampledLoci];
			for (int alleles = 0; alleles < NumSampledLoci; alleles++)
			{
				F_Marker.Adult_AlleleFreq[s].allele_freq[alleles].LociOnChrom = new double[NumAlleles];
				F_Marker.Progeny_AlleleFreq[s].allele_freq[alleles].LociOnChrom = new double[NumAlleles];
			}
			for (int alleles = 0; alleles < NumQTLs; alleles++)
			{
				F_QTL.Adult_AlleleFreq[s].allele_freq[alleles].LociOnChrom = new double[NumAlleles];
				F_QTL.Progeny_AlleleFreq[s].allele_freq[alleles].LociOnChrom = new double[NumAlleles];

			}
		}

		//assign loci and alleles
		int count = 0;
		for (ch = 0; ch < NumChrom; ch++)
		{
			if (known_qtls)
			{
				for (nq = 0; nq < NumQTLs; nq++)
				{
					Locations[ch].LociOnChrom[nq] = qtl_list[count];
					count++;
				}
			}
			else
			{
				for (nq = 0; nq < NumQTLs; nq++)
					Locations[ch].LociOnChrom[nq] = randnum(NumMarkers);
			}

			for (m = 0; m < NumMarkers; m++)
				MarkerLoci[ch].LociOnChrom[m] = randnum(NumMarkers);
		}
		bool QTLfound;
		for (ch = 0; ch < NumChrom; ch++)
		{
			for (k = 0; k < NumMarkers; k++)
			{
				QTLfound = false;
				for (nq = 0; nq < NumQTLs; nq++)
				{
					if (!QTLfound)
					{
						if (Locations[ch].LociOnChrom[nq] == k)
						{//it is a QTL
							QTLtracker[ch].LociOnChrom[k] = 1;
							QTLfound = true;
						}
						else //it's not a QTL
							QTLtracker[ch].LociOnChrom[k] = 0;
					}
				}
			}
		}


		bool found;
		for (ch = 0; ch < NumChrom; ch++)
		{
			avg_d.push_back(0);
			long_dist_dprime.push_back(0);
			counter = 0;
			for (nq = 0; nq < NumMarkers; nq++)
			{
				rand = genrand();
				if (rand < prob && counter < NumSampledLoci)
				{
					SampledLoci[ch].LociOnChrom[nq] = 1;
					counter++;
				}
				else
					SampledLoci[ch].LociOnChrom[nq] = 0;
			}
			if (counter < NumSampledLoci)
			{
				for (int l = 0; l < NumMarkers; l++)
				{
					if (counter < NumSampledLoci)
					{
						found = false;
						for (int comp = 0; comp < counter; comp++)
						{
							if (SampledLoci[ch].LociOnChrom[comp] == 1)
								found = true;
						}
						if (found == false)
						{
							SampledLoci[ch].LociOnChrom[counter] = 1;
							counter++;
						}
					}
				}
			}
		}

		////Sort Sampled Loci
		//double swapholder[1];
		//int last = NumSampledLoci - 1;
		//for (int chrom = 0; chrom < NumChrom; chrom++)
		//{
		//	for (int counter = NumSampledLoci - 1; counter >= 0; counter--)
		//	{
		//		for (int index = 0; index < last; index++)
		//		{
		//			if (SampledLoci[chrom].LociOnChrom[index] > SampledLoci[chrom].LociOnChrom[index + 1])
		//			{//basic swap
		//				swapholder[0] = SampledLoci[chrom].LociOnChrom[index + 1];
		//				//move index to index + 1
		//				SampledLoci[chrom].LociOnChrom[index + 1] = SampledLoci[chrom].LociOnChrom[index];
		//				//move swapholder to index
		//				SampledLoci[chrom].LociOnChrom[index] = swapholder[0];
		//			}
		//		}//end index
		//		last--;
		//	}//end counter
		//}//end chrom
		//
		double* tempallele = new double[NumAlleles];
		for (iii = 0; iii < NumAlleles; iii++)
			tempallele[iii] = randnorm(0, AllelicStdDev);
		for (iii = 0; iii<CarryingCapacity; iii++){
			Adult[iii].phenotype = 0;
			Adult[iii].Genotype = 0;
			for (jjj = 0; jjj<NumChrom; jjj++)
			{
				//Assign allelic effects
				for (eee = 0; eee < NumQTLs; eee++)
				{
					Adult[iii].maternal[jjj].allelicEffects[eee] = tempallele[iii%NumAlleles];
					Adult[iii].paternal[jjj].allelicEffects[eee] = tempallele[iii%NumAlleles];
					Adult[iii].Genotype = Adult[iii].Genotype +
						Adult[iii].maternal[jjj].allelicEffects[eee] +
						Adult[iii].paternal[jjj].allelicEffects[eee];
				}
				Adult[iii].phenotype = Adult[iii].Genotype + randnorm(0, EnvStdDev);

				//Assign genotypes at loci
				//maternal = randnum(NumAlleles);
				//paternal = randnum(NumAlleles);
				for (lll = 0; lll<NumMarkers; lll++)
				{
					Adult[iii].maternal[jjj].LociArray[lll] = iii%NumAlleles;
					Adult[iii].paternal[jjj].LociArray[lll] = iii%NumAlleles;
				}//end of lll loop
			}//jjj loop
			Adult[iii].Alive = true;
			if (genrand() < 0.5){
				Adult[iii].IsFemale = true;
				mFemaleTrait = mFemaleTrait + Adult[iii].phenotype;
				NumFemales++;
			}
			else{
				Adult[iii].IsFemale = false;
				NumMales++;
				mMaleTrait = mMaleTrait + Adult[iii].phenotype;
			}

		}//iii loop
		PopulationSize = CarryingCapacity;
		males = NumMales;
		females = NumFemales;
		mMaleTrait = mMaleTrait / males;
		mFemaleTrait = mFemaleTrait / females;
		delete[] tempallele;
	}//end initialize

	void DeInitialize(void)
	{
		int iii, jjj, sss;
		for (iii = 0; iii < CarryingCapacity; iii++)
		{
			for (jjj = 0; jjj < NumChrom; jjj++)
			{
				delete[] Adult[iii].paternal[jjj].allelicEffects;
				delete[] Adult[iii].paternal[jjj].LociArray;

				delete[] Adult[iii].maternal[jjj].allelicEffects;
				delete[] Adult[iii].maternal[jjj].LociArray;
			}
			delete[] Adult[iii].maternal;
			delete[] Adult[iii].paternal;
		}
		for (iii = 0; iii < MaxNumProg; iii++)
		{
			for (jjj = 0; jjj < NumChrom; jjj++)
			{
				delete[] Progeny[iii].paternal[jjj].allelicEffects;
				delete[] Progeny[iii].paternal[jjj].LociArray;

				delete[] Progeny[iii].maternal[jjj].allelicEffects;
				delete[] Progeny[iii].maternal[jjj].LociArray;
			}
			delete[] Progeny[iii].maternal;
			delete[] Progeny[iii].paternal;
		}
		delete[] Adult;
		delete[] Progeny;
		delete[] ChosenSexes;
		delete[] Fathers;
		delete[] Mothers;
		delete[] Nonmated;
		for (sss = 0; sss < NumChrom; sss++)
		{
			for (int all = 0; all < NumSampledLoci; all++)
			{
				delete[] F_Marker.Adult_AlleleFreq[sss].allele_freq[all].LociOnChrom;
				delete[] F_Marker.Progeny_AlleleFreq[sss].allele_freq[all].LociOnChrom;
			}
			for (int all = 0; all < NumQTLs; all++)
			{
				delete[] F_QTL.Adult_AlleleFreq[sss].allele_freq[all].LociOnChrom;
				delete[] F_QTL.Progeny_AlleleFreq[sss].allele_freq[all].LociOnChrom;
			}
			delete[] FstTracker[sss].LociOnChrom;
			delete[] Locations[sss].LociOnChrom;
			delete[] Polymorphisms[sss].LociOnChrom;
			delete[] OPolymorphisms[sss].LociOnChrom;
			delete[] A_NumAlleles[sss].LociOnChrom;
			delete[] M_NumAlleles[sss].LociOnChrom;
			delete[] O_NumAlleles[sss].LociOnChrom;
			delete[] MalHet[sss].LociOnChrom;
			delete[] OffHet[sss].LociOnChrom;
			delete[] AdultHet[sss].LociOnChrom;
			delete[] TotHet[sss].LociOnChrom;
			delete[] HetS[sss].LociOnChrom;
			delete[] ALD[sss].LociOnChrom;
			delete[] OLD[sss].LociOnChrom;
			delete[] MLD[sss].LociOnChrom;
			delete[] QHs[sss].LociOnChrom;
			delete[] QHt[sss].LociOnChrom;
			delete[] QMalHet[sss].LociOnChrom;
			delete[] QOffHet[sss].LociOnChrom;
			delete[] QAdultHet[sss].LociOnChrom;
			delete[] QTotHet[sss].LociOnChrom;
			delete[] AFreq1[sss].LociOnChrom;
			delete[] OFreq1[sss].LociOnChrom;
			delete[] MFreq1[sss].LociOnChrom;
			delete[] MaleH[sss].LociOnChrom;
			delete[] AMajor[sss].LociOnChrom;
			delete[] OMajor[sss].LociOnChrom;
			delete[] MMajor[sss].LociOnChrom;
			delete[] Psig[sss].LociOnChrom;
			delete[] PopLD_a[sss].LociOnChrom;
			delete[] PopMajorAllele[sss].LociOnChrom;
			delete[] PopPrimFreq[sss].LociOnChrom;
			delete[] QTLtracker[sss].LociOnChrom;
			delete[] GenetVars[sss].LociOnChrom;
			delete[] GenetMeans[sss].LociOnChrom;
			delete[] MarkerLoci[sss].LociOnChrom;
			delete[] SampledLoci[sss].LociOnChrom;
			delete[] ANumQTLAlleles[sss].LociOnChrom;
			delete[] ONumQTLAlleles[sss].LociOnChrom;
			delete[] MNumQTLAlleles[sss].LociOnChrom;
			delete[] TNumQTLAlleles[sss].LociOnChrom;
			delete[] MTNumQTLAlleles[sss].LociOnChrom;
			delete[] F_QTL.FstValue[sss].LociOnChrom;
			delete[] F_QTL.BS_Pvalue[sss].LociOnChrom;
			delete[] F_QTL.PValue[sss].LociOnChrom;
			delete[] F_QTL.WeightedFst[sss].LociOnChrom;
			delete[] F_QTL.Locus[sss].LociOnChrom;
			delete[] F_QTL.Adult_AlleleFreq[sss].allele_freq;
			delete[] F_QTL.Progeny_AlleleFreq[sss].allele_freq;
			delete[] F_Marker.WeightedFst[sss].LociOnChrom;
			delete[] F_Marker.FstValue[sss].LociOnChrom;
			delete[] F_Marker.BS_Pvalue[sss].LociOnChrom;
			delete[] F_Marker.PValue[sss].LociOnChrom;
			delete[] F_Marker.Locus[sss].LociOnChrom;
			delete[] F_Marker.Adult_AlleleFreq[sss].allele_freq;
			delete[] F_Marker.Progeny_AlleleFreq[sss].allele_freq;
		}

		delete[] SampledLoci;
		delete[] Polymorphisms;
		delete[] OPolymorphisms;
		delete[] Locations;
		delete[] QTLtracker;
		delete[] MarkerLoci;
		delete[] PopPrimFreq;
		delete[] PopLD_a;
		delete[] PopMajorAllele;
		delete[] QAdultHet;
		delete[] QHs;
		delete[] QHt;
		delete[] QMalHet;
		delete[] QOffHet;
		delete[] QTotHet;
		delete[] MalHet;
		delete[] OffHet;
		delete[] TotHet;
		delete[] HetS;
		delete[] MaleH;
		delete[] AdultHet;
		delete[] AFreq1;
		delete[] MFreq1;
		delete[] OFreq1;
		delete[] AMajor;
		delete[] MMajor;
		delete[] OMajor;
		delete[] ALD;
		delete[] OLD;
		delete[] MLD;
		delete[] A_NumAlleles;
		delete[] M_NumAlleles;
		delete[] O_NumAlleles;
		delete[] Psig;
		delete[] GenetVars;
		delete[] GenetMeans;
		delete[] ANumQTLAlleles;
		delete[] ONumQTLAlleles;
		delete[] MNumQTLAlleles;
		delete[] TNumQTLAlleles;
		delete[] MTNumQTLAlleles;
		delete[] FstTracker;

		delete[] F_Marker.FstValue;
		delete[] F_Marker.BS_Pvalue;
		delete[] F_Marker.PValue;
		delete[] F_Marker.WeightedFst;
		delete[] F_Marker.Locus;
		delete[] F_Marker.Adult_AlleleFreq;
		delete[] F_Marker.Progeny_AlleleFreq;
		delete[] F_QTL.BS_Pvalue;
		delete[] F_QTL.PValue;
		delete[] F_QTL.FstValue;
		delete[] F_QTL.WeightedFst;
		delete[] F_QTL.Adult_AlleleFreq;
		delete[] F_QTL.Progeny_AlleleFreq;
		delete[] F_QTL.Locus;
		delete[] SampledAdults;
		delete[] SampledOffspring;
	}//end of DeInitialize

	void StandardizeGenotypes()
	{
		//So I need to standardize to a given mean and variance
		int s, ss, sss;
		vector <double> alleles;
		double mean, std_dev;
		double constant;
		int num_alleles = 0;

		////figure out how many alleles there are
		//for (s = 0; s < PopulationSize; s++)
		//{
		//	for (ss = 0; ss < NumChrom; ss ++)
		//	{
		//		for (sss = 0; sss < NumQTLs; sss++)
		//		{
		//			for (a = 0; a < alleles.size(); a++)
		//			{
		//				if (Adult[s].paternal[ss].allelicEffects[sss] != alleles[a])
		//				{
		//					alleles.push_back(Adult[s].Genotype);
		//					num_alleles++;
		//				}
		//				if (Adult[s].maternal[ss].allelicEffects[sss] != alleles[a] 
		//					&& Adult[s].maternal[ss].allelicEffects[sss] != Adult[s].paternal[ss].allelicEffects[sss])
		//				{
		//					alleles.push_back(Adult[s].Genotype);
		//					num_alleles++;
		//				}
		//			}
		//		}
		//	}
		//}

		//calculate mean and std dev
		mean = 0;
		std_dev = 0;
		for (s = 0; s < PopulationSize; s++)
			mean = mean + Adult[s].Genotype;
		mean = mean / PopulationSize;
		for (s = 0; s < PopulationSize; s++)
			std_dev = std_dev + (Adult[s].Genotype - mean)*(Adult[s].Genotype - mean);
		std_dev = sqrt(std_dev / PopulationSize);
		//cout << "Starting Mean: " << mean << '\t' << "Starting Std Dev: " << std_dev << '\t';
		constant = mean / (NumQTLs * 2 * NumChrom);
		//cout << "Constant: " << constant<<'\t';
		//Now standardize the allelic effects
		for (s = 0; s < PopulationSize; s++)
		{
			Adult[s].phenotype = Adult[s].phenotype - Adult[s].Genotype;
			Adult[s].Genotype = 0;
			for (ss = 0; ss < NumChrom; ss++)
			{
				for (sss = 0; sss < NumQTLs; sss++)
				{
					Adult[s].paternal[ss].allelicEffects[sss] =
						(Adult[s].paternal[ss].allelicEffects[sss] - constant) / std_dev;
					Adult[s].maternal[ss].allelicEffects[sss] =
						(Adult[s].maternal[ss].allelicEffects[sss] - constant) / std_dev;
					Adult[s].Genotype = Adult[s].Genotype + Adult[s].paternal[ss].allelicEffects[sss] + Adult[s].maternal[ss].allelicEffects[sss];
				}
			}
			Adult[s].phenotype = Adult[s].Genotype + Adult[s].phenotype;
		}
		mean = 0;
		std_dev = 0;
		for (s = 0; s < PopulationSize; s++)
			mean = mean + Adult[s].Genotype;
		mean = mean / PopulationSize;
		for (s = 0; s < PopulationSize; s++)
			std_dev = std_dev + (Adult[s].Genotype - mean)*(Adult[s].Genotype - mean);
		std_dev = sqrt(std_dev / PopulationSize);
		//cout << "Stdized Mean: " << mean << '\t' << "Stdized Std Dev: " << std_dev << '\n';	
	}//end StandardizeTraits

	void RecombineChromosome(Chromosome &RecombinedChr, Individual &Parent, int WhichChromosome, double ExpectedRecombEvents)
	{
		int RCi, RCj;
		int NumberRecombEvents = 0;
		int SegmentStart[22], SegmentEnd[22];
		int BreakPoint[20];

		if (ExpectedRecombEvents < 6)
			NumberRecombEvents = poissonrand(ExpectedRecombEvents);
		if (ExpectedRecombEvents > 6)
			NumberRecombEvents = positiveroundnorm(ExpectedRecombEvents, sqrt(ExpectedRecombEvents));
		for (RCi = 0; RCi < 20; RCi++)
			BreakPoint[RCi] = NumMarkers + 1;
		bool SegmentMaternal[22];
		int NumberSegments;
		bool StartMaternal;
		if (NumberRecombEvents > 20)
			NumberRecombEvents = 20;

		if (NumberRecombEvents > 0)
		{
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
				BreakPoint[RCi] = randnum(NumMarkers);
			//sort breakpoints
			sort(begin(BreakPoint), end(BreakPoint));

			//first segment maternal or paternal?
			if (genrand() < 0.5)
				StartMaternal = true;
			else
				StartMaternal = false;

			NumberSegments = 1;
			SegmentStart[0] = 0;
			SegmentMaternal[0] = StartMaternal;
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
			{
				SegmentEnd[RCi] = BreakPoint[RCi];
				SegmentStart[RCi + 1] = BreakPoint[RCi];
				if (SegmentMaternal[RCi])
					SegmentMaternal[RCi + 1] = false;
				else
					SegmentMaternal[RCi + 1] = true;
				NumberSegments++;
			}//end RCi
			SegmentEnd[RCi] = NumMarkers;

			//now pass allelic info to recombined chromosome
			for (RCi = 0; RCi < NumberSegments; RCi++)
			{
				if (SegmentMaternal[RCi])
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
						RecombinedChr.LociArray[RCj] = Parent.maternal[WhichChromosome].LociArray[RCj];
				}
				else
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
						RecombinedChr.LociArray[RCj] = Parent.paternal[WhichChromosome].LociArray[RCj];
				}
			}//end RCi

			//now do the QTLs
			for (RCi = 0; RCi < NumQTLs; RCi++)
			{
				for (RCj = 0; RCj < NumberSegments; RCj++)
				{
					if (Locations[WhichChromosome].LociOnChrom[RCi] >= SegmentStart[RCj]
						&& Locations[WhichChromosome].LociOnChrom[RCi] < SegmentEnd[RCj])
					{
						if (SegmentMaternal[RCj])
							RecombinedChr.allelicEffects[RCi] = Parent.maternal[WhichChromosome].allelicEffects[RCi];
						else
							RecombinedChr.allelicEffects[RCi] = Parent.paternal[WhichChromosome].allelicEffects[RCi];
					}
				}
			}//end RCi
		}//end numb recomb events > 0
		else
		{
			//No recombination
			if (genrand() < 0.5)
			{
				for (RCi = 0; RCi < NumMarkers; RCi++)
					RecombinedChr.LociArray[RCi] = Parent.maternal[WhichChromosome].LociArray[RCi];
				for (RCi = 0; RCi < NumQTLs; RCi++)
					RecombinedChr.allelicEffects[RCi] = Parent.maternal[WhichChromosome].allelicEffects[RCi];
			}
			else
			{
				for (RCi = 0; RCi < NumMarkers; RCi++)
					RecombinedChr.LociArray[RCi] = Parent.paternal[WhichChromosome].LociArray[RCi];
				for (RCi = 0; RCi < NumQTLs; RCi++)
					RecombinedChr.allelicEffects[RCi] = Parent.paternal[WhichChromosome].allelicEffects[RCi];
			}
		}//else (no recomb)
	}//end RecombineChromosome

	void Mating(double GaussianPrefVariance)
	{
		int cc;
		int NumProg;
		double dubrand;
		NumProg = 0;
		//Mate Choice:
		int mm, nn, males;
		int Encounters, MaleID, irndnum, Females;
		double MeanFemaleTrait;
		bool MateFound, MalesPresent;
		double MeanMaleTrait, MateProb;
		int MaleIndex, FemaleID;
		int Counter3;
		vector <int> MaleList;
		NumProg = 0;
		Mated = 0;
		Females = 0;
		ProgFem = 0;
		ProgMale = 0;
		NumFemales = 0;
		MalesPresent = false;
		NumMales = 0;
		MeanMaleTrait = 0;
		double SDMaleTrait = 0;
		MeanFemaleTrait = 0;
		Counter1 = 0;
		Counter2 = 0;
		Counter3 = 0;
		//determine mean male trait
		for (males = 0; males < PopulationSize; males++)
		{
			Adult[males].MateFound = 0;
			if (!Adult[males].IsFemale)
			{
				MalesPresent = true;
				MaleList.push_back(males);
				NumMales++;
				MeanMaleTrait = MeanMaleTrait + Adult[males].phenotype;
			}
		} // end of males
		if (NumMales > 0) {
			MeanMaleTrait = MeanMaleTrait / NumMales;
		}
		else {
			MeanMaleTrait = 0;
			popExtinct = true;
		}
		for (males = 0; males < PopulationSize; males++){
			if (!Adult[males].IsFemale)
				SDMaleTrait = SDMaleTrait + (Adult[males].phenotype - MeanMaleTrait)*(Adult[males].phenotype - MeanMaleTrait);
		}
		SDMaleTrait = sqrt(SDMaleTrait / NumMales);
		for (mm = 0; mm < PopulationSize; mm++)
		{
			MateFound = false;
			if (Adult[mm].IsFemale && MalesPresent)
			{
				FemaleID = mm;
				MeanFemaleTrait = MeanFemaleTrait + Adult[mm].phenotype;
				Females++;
				Encounters = 0;
				while (!MateFound && Encounters <= MaximumEncounters)
				{
					irndnum = randnum(NumMales);
					MaleIndex = MaleList[irndnum];
					if (GaussianPrefVariance > 0)
						MateProb = exp(-0.5 * (Adult[MaleIndex].phenotype - GaussianPrefMean)*
						(Adult[MaleIndex].phenotype - GaussianPrefMean) / GaussianPrefVariance);
					else
						MateProb = 1;
					dubrand = genrand();
					if (dubrand < MateProb)
					{
						MateFound = true;
						MaleID = MaleIndex;
						Fathers[Counter2] = MaleID;
						Adult[MaleID].MateFound++;
						Counter2++;
						Mated++;
					}
					Encounters++;
				}//while
				if (MateFound)
				{
					Mothers[Counter1] = mm;
					Counter1++;
					Adult[mm].MateFound++;
					if (NumProg >= MaxNumProg)
						NumProg = MaxNumProg - 1;
					//mother is Parent 1, mm
					//father is parent 2, MateID
					for (nn = 0; nn < MaximumFecundity; nn++)
					{
						Progeny[NumProg].Alive = true;
						//calculate phenotype in mutation once the genotype is determined
						for (cc = 0; cc < NumChrom; cc++)//go through chromosome by chromosome
						{
							RecombineChromosome(Progeny[NumProg].maternal[cc], Adult[FemaleID], cc, RecombRate);
							RecombineChromosome(Progeny[NumProg].paternal[cc], Adult[MaleID], cc, RecombRate);
							////Does recombination occur?
							//randRecomb = genrand();
							//num1 = genrand();
							//num2 = genrand();
							//ChooseLocus = randnum(NumMarkers);
							//if (num1 < 0.5)
							//{//then we're talking about the maternal chrom from the Adult
							//	if (randRecomb < RecombRate)//then recombination occurs
							//	{
							//		for (ccc = 0; ccc < NumMarkers; ccc++)
							//		{
							//			if (ccc < ChooseLocus)//first part is from maternal
							//				Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].maternal[cc].LociArray[ccc];
							//			else//second part is from paternal
							//				Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].paternal[cc].LociArray[ccc];
							//			for (qq = 0; qq < NumQTLs; qq++)
							//			{
							//				if (Locations[cc].LociOnChrom[qq] < ChooseLocus)
							//					Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].maternal[cc].allelicEffects[qq];
							//				else
							//					Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].paternal[cc].allelicEffects[qq];
							//			}//end of QTLS
							//		}//end of locus
							//	}//end of recomb rate
							//	else//no recombination occurs [maternal]
							//	{
							//		for (ccc = 0; ccc < NumMarkers; ccc++)
							//			Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].maternal[cc].LociArray[ccc];
							//		for (qq = 0; qq < NumQTLs; qq++)
							//			Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].maternal[cc].allelicEffects[qq];
							//	}

							//}//end of recombiation for Parent 1
							//else //then Parent1's paternal chromosome is passed to the progeny's maternal chromosome
							//{//num1 > 0.5
							//	if (randRecomb < RecombRate)//then recombination occurs
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++)
							//		{
							//			if (ccc < ChooseLocus)//first part is from paternal
							//				Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].paternal[cc].LociArray[ccc];
							//			else//second part is from maternal
							//				Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].maternal[cc].LociArray[ccc];
							//		}//end of locus		
							//		for (qq = 0; qq < NumQTLs; qq++)
							//		{
							//			if (Locations[cc].LociOnChrom[qq] < ChooseLocus)
							//				Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].paternal[cc].allelicEffects[qq];
							//			else
							//				Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].maternal[cc].allelicEffects[qq];
							//		}
							//	}
							//	else
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++)
							//			Progeny[NumProg].maternal[cc].LociArray[ccc] = Adult[mm].paternal[cc].LociArray[ccc];
							//		for (qq = 0; qq < NumQTLs; qq++)
							//			Progeny[NumProg].maternal[cc].allelicEffects[qq] = Adult[mm].paternal[cc].allelicEffects[qq];
							//	}
							//}//end of num1 >0.5
							////now Parent 2
							//randRecomb = genrand();//each parent gets its own chance at recombination
							//ChooseLocus = randnum(LociNo);
							//if (num2 < 0.5)//then progeny gets maternal chrom (unless recomb occurs)
							//{
							//	if (randRecomb < RecombRate)
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++)
							//		{
							//			if (ccc < ChooseLocus)//first part is from maternal
							//				Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[MateID].maternal[cc].LociArray[ccc];
							//			if (ccc >= ChooseLocus)//second part is from paternal
							//				Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[MateID].paternal[cc].LociArray[ccc];
							//		}//end of locus
							//		for (qq = 0; qq < NumQTLs; qq++)
							//		{
							//			if (Locations[cc].LociOnChrom[qq] < ChooseLocus)
							//				Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[MateID].maternal[cc].allelicEffects[qq];
							//			else
							//				Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[MateID].paternal[cc].allelicEffects[qq];
							//		}
							//	}//end rand recomb
							//	else//no recombination [maternal]
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++)
							//			Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[MateID].maternal[cc].LociArray[ccc];
							//		for (qq = 0; qq < NumQTLs; qq++)
							//			Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[MateID].maternal[cc].allelicEffects[qq];
							//	}//end of no recomb [maternal]
							//}
							//if (num2 >= 0.5) //then the Adult's paternal chromosome is passed to the progeny's paternal chromosome
							//{
							//	if (randRecomb < RecombRate)
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++){
							//			if (ccc < ChooseLocus)//first part is from paternal
							//				Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[MateID].paternal[cc].LociArray[ccc];

							//			if (ccc >= ChooseLocus)//second part is from maternal
							//				Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[MateID].maternal[cc].LociArray[ccc];
							//		}//end of locus
							//		for (qq = 0; qq < NumQTLs; qq++)
							//		{
							//			if (Locations[cc].LociOnChrom[qq] < ChooseLocus)
							//				Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[MateID].paternal[cc].allelicEffects[qq];
							//			else
							//				Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[MateID].maternal[cc].allelicEffects[qq];
							//		}
							//	}//end recombination
							//	else//no recombination occurs [maternal]
							//	{
							//		for (ccc = 0; ccc < LociNo; ccc++)
							//			Progeny[NumProg].paternal[cc].LociArray[ccc] = Adult[mm].paternal[cc].LociArray[ccc];
							//		for (qq = 0; qq < NumQTLs; qq++)
							//			Progeny[NumProg].paternal[cc].allelicEffects[qq] = Adult[mm].paternal[cc].allelicEffects[qq];
							//	}
							//}//end num2 > 0.5
						}//end of chromosome
						if (genrand() < 0.5)
						{
							Progeny[NumProg].IsFemale = true;
							ProgFem++;
						}
						else
						{
							Progeny[NumProg].IsFemale = false;
							ProgMale++;
						}
						NumProg++;
					}//for nn
				}
				if (!MateFound)
				{//Keep track of the females that didn't mate
					Nonmated[Counter3] = mm;
					Counter3++;
				}
			}//if female
		}//end of mm
		ProgenyNum = NumProg;
		NumFemales = Females;
		NumNotMated = Counter3;
		MeanFemaleTrait = MeanFemaleTrait / Females;
	}//end mating

	void Mutation()
	{
		int m, mm, irand, irand2, gg, ggg, irand3, locus;
		double rnd1, rnd2;
		double IndMutationRate;
		double MutSD;
		bool mutated;

		MutSD = sqrt(MutationalVariance);
		IndMutationRate = MutationRate * 2 * TotalLociNo;

		for (m = 0; m < ProgenyNum; m++)
		{
			rnd1 = genrand();
			Progeny[m].phenotype = 0;
			Progeny[m].Genotype = 0;
			mutated = false;
			if (rnd1 < IndMutationRate)
			{
				irand = randnum(NumChrom);
				irand2 = randnum(NumMarkers);
				locus = MarkerLoci[irand].LociOnChrom[irand2];
				rnd2 = genrand();//to choose maternal or paternal
				if (rnd2 < 0.5)//affects maternal chromosome	
				{
					while (!mutated){
						irand3 = randnum(NumAlleles);
						if (!Progeny[m].maternal[irand].LociArray[locus] == irand3)
						{
							Progeny[m].maternal[irand].LociArray[locus] = irand3;
							mutated = true;
						}
					}
					for (mm = 0; mm < NumQTLs; mm++)
					{
						if (Locations[irand].LociOnChrom[mm] == irand2)
							Progeny[m].maternal[irand].allelicEffects[mm] =
							Progeny[m].maternal[irand].allelicEffects[mm] + randnorm(0, MutSD);
					}
				}
				else//affects paternal chromosome
				{
					while (!mutated){
						irand3 = randnum(NumAlleles);
						if (!Progeny[m].paternal[irand].LociArray[locus] == irand3)
						{
							Progeny[m].paternal[irand].LociArray[locus] = irand3;
							mutated = true;
						}
					}
					for (mm = 0; mm < NumQTLs; mm++)
					{
						if (Locations[irand].LociOnChrom[mm] == irand2)
							Progeny[m].paternal[irand].allelicEffects[mm] =
							Progeny[m].paternal[irand].allelicEffects[mm] + randnorm(0, MutSD);
					}
				}
			}//end of if

			for (gg = 0; gg < NumChrom; gg++)
			{
				for (ggg = 0; ggg < NumQTLs; ggg++){
					Progeny[m].Genotype = Progeny[m].Genotype +
						Progeny[m].maternal[gg].allelicEffects[ggg] + Progeny[m].paternal[gg].allelicEffects[ggg];
				}
			}
			Progeny[m].phenotype = Progeny[m].Genotype + randnorm(0, EnvStdDev);
		}//end of m
	}//mutation

	void SelectionOnPhenotypes(double dSelectionStrength)
	{
		int i, ProgAlive;
		double dSurvProb;
		double drnum1;
		double dOptimum;
		double phenSD = 0;
		double phenMean = 0;
		double num;
		int malecount = 0;
		ProgAlive = 0;
		//calc mean male phenotype & std dev
		for (i = 0; i < ProgenyNum; i++){
			if (!Progeny[i].IsFemale){
				phenMean = phenMean + Progeny[i].phenotype;
				malecount++;
			}
		}
		num = malecount;
		phenMean = phenMean / num;
		for (i = 0; i < ProgenyNum; i++){
			if (!Progeny[i].IsFemale)
				phenSD = phenSD + (Progeny[i].phenotype - phenMean)*(Progeny[i].phenotype - phenMean);
		}
		phenSD = sqrt(phenSD / num);

		dOptimum = 0;
		for (i = 0; i < ProgenyNum; i++) {
			if (Progeny[i].IsFemale)
			{
				Progeny[i].Alive = true;
				ProgAlive++;
			}
			else//selection only on males
			{
				if (dSelectionStrength > 0)
					dSurvProb = exp(-1 * (Progeny[i].phenotype - dOptimum)*(Progeny[i].phenotype - dOptimum)
					/ (2 * dSelectionStrength));
				else
					dSurvProb = 1;
				//cout<<dSurvProb<<'\n';
				drnum1 = genrand();
				if (drnum1 < dSurvProb)
				{
					Progeny[i].Alive = true;
					ProgAlive++;
				}
				else
					Progeny[i].Alive = false;
			}
		} // end of i
	}//end selection on phenotypes

	void DensityRegulation()
	{
		int p, pp, ppp;
		int iNumAdultsChosen;
		double CarCapUnfilled, ProgLeft, KeepProb;
		double DRrandnum;
		NumMales = 0;
		NumFemales = 0;
		ProgLeft = 0;
		//count the ones that are still alive
		for (p = 0; p < ProgenyNum; p++)
		{
			if (Progeny[p].Alive)
				ProgLeft++;
		}
		CarCapUnfilled = CarryingCapacity;
		iNumAdultsChosen = 0;
		for (p = 0; p<ProgenyNum; p++)
		{
			if (Progeny[p].Alive)
			{
				if (ProgLeft == 0)
					KeepProb = 0;
				else
					KeepProb = CarCapUnfilled / ProgLeft;
				DRrandnum = genrand();
				if (DRrandnum<KeepProb)
				{//then turn it into an adult
					Adult[iNumAdultsChosen].Alive = true;
					Adult[iNumAdultsChosen].MateFound = 0;
					for (pp = 0; pp < NumChrom; pp++)
					{
						for (ppp = 0; ppp < NumMarkers; ppp++)
						{
							Adult[iNumAdultsChosen].maternal[pp].LociArray[ppp] = Progeny[p].maternal[pp].LociArray[ppp];
							Adult[iNumAdultsChosen].paternal[pp].LociArray[ppp] = Progeny[p].paternal[pp].LociArray[ppp];
						}
						for (ppp = 0; ppp < NumQTLs; ppp++)
						{
							Adult[iNumAdultsChosen].maternal[pp].allelicEffects[ppp] = Progeny[p].maternal[pp].allelicEffects[ppp];
							Adult[iNumAdultsChosen].paternal[pp].allelicEffects[ppp] = Progeny[p].paternal[pp].allelicEffects[ppp];
						}
					}
					Adult[iNumAdultsChosen].phenotype = Progeny[p].phenotype;
					Adult[iNumAdultsChosen].Genotype = Progeny[p].Genotype;
					if (Progeny[p].IsFemale){
						Adult[iNumAdultsChosen].IsFemale = true;
						ChosenSexes[iNumAdultsChosen] = 1;
						NumFemales++;
					}
					else{
						Adult[iNumAdultsChosen].IsFemale = false;
						ChosenSexes[iNumAdultsChosen] = 0;
						NumMales++;
					}
					CarCapUnfilled = CarCapUnfilled - 1;
					iNumAdultsChosen++;
				}//end of if KeepProb
				else
					Progeny[p].Alive = false;
			}//end of if Alive
			ProgLeft = ProgLeft - 1;
		}//end of for p
		PopulationSize = iNumAdultsChosen;
		if (PopulationSize == 0)
			popExtinct = true;
	}//end Density Regulation

	void CalcMeanTraitValues()
	{
		int p;
		double dmales, dfemales;
		NumMales = 0;
		NumFemales = 0;
		mMaleTrait = 0;
		mFemaleTrait = 0;
		for (p = 0; p < PopulationSize; p++)
		{
			if (Adult[p].IsFemale)
			{
				mFemaleTrait = mFemaleTrait + Adult[p].phenotype;
				NumFemales++;
			}
			if (!Adult[p].IsFemale)
			{
				mMaleTrait = mMaleTrait + Adult[p].phenotype;
				NumMales++;
			}
		}
		dmales = NumMales;
		dfemales = NumFemales;
		mMaleTrait = mMaleTrait / dmales;
		mFemaleTrait = mFemaleTrait / dfemales;
	}//end calc mean trait values

	void AdultSumm()
	{
		int e, f, g, adultcount = 0;
		double AGenmean = 0, AGenVar = 0;
		double APhenMean = 0, APhenVar = 0;
		double dNa = PopulationSize;
		double AHeritability;
		double aMaleTraitMean = 0, aFemTraitMean = 0;
		double ANf, ANm;
		double ATraitSD = 0, amTraitSD = 0, afTraitSD = 0;
		double MSmean = 0, fMSmean = 0, mMSmean = 0;
		NumMales = 0;
		NumFemales = 0;
		double MS;
		double Amprime, Fmprime, Mmprime;
		double* StdATrait = new double[PopulationSize];
		double* FstdAtrait, *MstdAtrait;
		double* RelMS = new double[PopulationSize];

		for (e = 0; e < PopulationSize; e++)
		{
			if (Adult[e].Alive)
			{
				AGenmean = AGenmean + Adult[e].Genotype;
				APhenMean = APhenMean + Adult[e].phenotype;
				MSmean = MSmean + Adult[e].MateFound;
				if (Adult[e].IsFemale)
				{
					fMSmean = fMSmean + Adult[e].MateFound;
					aFemTraitMean = aFemTraitMean + Adult[e].phenotype;
					NumFemales++;
				}
				if (!Adult[e].IsFemale)
				{
					mMSmean = mMSmean + Adult[e].MateFound;
					aMaleTraitMean = aMaleTraitMean + Adult[e].phenotype;
					NumMales++;
				}
			}
		}
		AGenmean = AGenmean / dNa;
		APhenMean = APhenMean / dNa;
		MSmean = MSmean / dNa;
		ANf = NumFemales;
		ANm = NumMales;
		aMaleTraitMean = aMaleTraitMean / ANm;
		mMSmean = mMSmean / ANm;
		aFemTraitMean = aFemTraitMean / ANf;
		fMSmean = fMSmean / ANf;
		//cout << "\nmale ms: " << mMSmean << " female ms: " << fMSmean;
		vector <double> Fem_RelMS;
		vector <double> Mal_RelMS;
		int malecount = 0, femalecount = 0;
		for (e = 0; e < PopulationSize; e++)
		{
			if (Adult[e].Alive)
			{
				adultcount++;
				MS = Adult[e].MateFound;
				if (MSmean > 0)
					RelMS[e] = MS / MSmean;
				else
					RelMS[e] = 0;
				ATraitSD = ATraitSD + (Adult[e].phenotype - APhenMean)*(Adult[e].phenotype - APhenMean);
				if (Adult[e].IsFemale){
					afTraitSD = afTraitSD + (Adult[e].phenotype - aFemTraitMean)*(Adult[e].phenotype - aFemTraitMean);
					if (fMSmean > 0)
						Fem_RelMS.push_back(MS / fMSmean);
					else
						Fem_RelMS.push_back(0);
					femalecount++;
				}
				else{
					amTraitSD = amTraitSD + (Adult[e].phenotype - aMaleTraitMean)*(Adult[e].phenotype - aMaleTraitMean);
					if (mMSmean > 0)
						Mal_RelMS.push_back(MS / mMSmean);
					else
						Mal_RelMS.push_back(0);
					malecount++;
				}
			}
		}
		PopulationSize = adultcount;
		ATraitSD = sqrt(ATraitSD / dNa);
		afTraitSD = sqrt(afTraitSD / ANf);
		amTraitSD = sqrt(amTraitSD / ANm);
		FstdAtrait = new double[NumFemales];
		MstdAtrait = new double[NumMales];
		malecount = 0;
		femalecount = 0;
		for (e = 0; e < PopulationSize; e++)
		{
			if (Adult[e].Alive)
			{
				AGenVar = AGenVar + (Adult[e].Genotype - AGenmean)*(Adult[e].Genotype - AGenmean);
				APhenVar = APhenVar + (Adult[e].phenotype - APhenMean)*(Adult[e].phenotype - APhenMean);
				if (ATraitSD > 0)
					StdATrait[e] = (Adult[e].phenotype - APhenMean) / ATraitSD;
				else
					StdATrait[e] = 0;
				if (Adult[e].IsFemale){
					if (afTraitSD > 0)
						FstdAtrait[femalecount] = (Adult[e].phenotype - aFemTraitMean) / afTraitSD;
					else
						FstdAtrait[femalecount] = 0;
					femalecount++;
				}
				else{
					if (amTraitSD > 0)
						MstdAtrait[malecount] = (Adult[e].phenotype - aMaleTraitMean) / amTraitSD;
					else
						MstdAtrait[malecount] = 0;
					malecount++;
				}
			}
		}
		AGenVar = AGenVar / dNa;
		APhenVar = APhenVar / dNa;
		//Calculate heritability
		AHeritability = AGenVar / APhenVar;
		//Calculate mating differential
		//From Jones (2009), eqn 9, m' = cov(std trait value, rel. mating success)
		Amprime = 0;
		Fmprime = 0;
		Mmprime = 0;
		for (e = 0; e < PopulationSize; e++)
		{
			if (Adult[e].Alive)
				Amprime = Amprime + (StdATrait[e] - APhenMean)*(RelMS[e] - MSmean);
		}
		for (f = 0; f < NumFemales; f++)
			Fmprime = Fmprime + (FstdAtrait[f] - aFemTraitMean)*(Fem_RelMS[f] - fMSmean);
		for (g = 0; g < NumMales; g++)
			Mmprime = Mmprime + (MstdAtrait[g] - aMaleTraitMean)*(Mal_RelMS[g] - mMSmean);
		Amprime = Amprime / dNa;
		Fmprime = Fmprime / ANf;
		Mmprime = Mmprime / ANm;

		malN = NumMales;
		femN = NumFemales;
		adultRatio = ANm / ANf;
		allphenavg = APhenMean;
		malphenavg = aMaleTraitMean;
		femphenavg = aFemTraitMean;
		allDiff = Amprime;
		maleDiff = Mmprime;
		femDiff = Fmprime;
		Heritability = AHeritability;
		GenVar = AGenVar;
		PhenVar = APhenVar;
		PopulationSize = dNa;

		delete[] StdATrait;
		delete[] FstdAtrait;
		delete[] MstdAtrait;
		delete[] RelMS;
	}//end AdultSumm

	void OffSumm()
	{
		int e, progcount = 0;
		double* OFitness = new double[ProgenyNum];
		double* OmFitness, *OfFitness;
		double OGenmean = 0, OGenVar = 0;
		double OPhenMean = 0, OPhenVar = 0;
		double dNo = ProgenyNum;
		double OHeritability;
		double oMaleTraitMean = 0, oFemTraitMean = 0;
		double ONf, ONm;
		double OTraitSD = 0, omTraitSD = 0, ofTraitSD = 0;
		int oFemN = 0, oMalN = 0;
		double FitMean = 0, FitMmean = 0, FitFmean = 0, W;
		double RelFit = 0, RelFitM = 0, RelFitF = 0;
		double* FstdOtrait, *MstdOtrait;
		double* StdOTrait = new double[ProgenyNum];
		double Osprime = 0, Omprime = 0, Ofprime = 0;
		double dnum_fem_tot = 0;
		double dnum_mal_tot = 0;
		//Calculate genetic variation
		//Phenotype and Genotype determined in Mutation
		for (e = 0; e < ProgenyNum; e++)
		{
			W = Progeny[e].Alive;
			FitMean = FitMean + W;
			if (Progeny[e].IsFemale)
			{
				FitFmean = FitFmean + W;
				dnum_fem_tot++;
			}
			else
			{
				FitMmean = FitMmean + W;
				dnum_mal_tot++;
			}
			if (Progeny[e].Alive)
			{
				progcount++;
				OGenmean = OGenmean + Progeny[e].Genotype;
				OPhenMean = OPhenMean + Progeny[e].phenotype;

				if (Progeny[e].IsFemale)
				{
					oFemTraitMean = oFemTraitMean + Progeny[e].phenotype;
					oFemN++;
				}
				else
				{
					oMaleTraitMean = oMaleTraitMean + Progeny[e].phenotype;
					oMalN++;
				}
			}
		}
		ONf = oFemN;
		ONm = oMalN;
		dNo = progcount;
		OGenmean = OGenmean / dNo;
		OPhenMean = OPhenMean / dNo;
		FitMean = FitMean / ProgenyNum;
		oFemTraitMean = oFemTraitMean / ONf;
		oMaleTraitMean = oMaleTraitMean / ONm;
		FitMmean = FitMmean / dnum_mal_tot;
		FitFmean = FitFmean / dnum_fem_tot;
		OmFitness = new double[oMalN];
		OfFitness = new double[oFemN];
		int femalecount = 0, malecount = 0;
		for (e = 0; e < ProgenyNum; e++)
		{
			if (Progeny[e].Alive)
			{
				W = Progeny[e].Alive;
				if (FitMean > 0)
					OFitness[e] = W / FitMean;
				else
					OFitness[e] = 0;
				OTraitSD = OTraitSD + (Progeny[e].phenotype - OPhenMean)*(Progeny[e].phenotype - OPhenMean);
				if (Progeny[e].IsFemale){
					ofTraitSD = ofTraitSD + (Progeny[e].phenotype - oFemTraitMean)*(Progeny[e].phenotype - oFemTraitMean);
					if (FitFmean > 0)
						OfFitness[femalecount] = W / FitFmean;
					else
						OfFitness[femalecount] = 0;
					femalecount++;
				}
				else{
					omTraitSD = omTraitSD + (Progeny[e].phenotype - oMaleTraitMean)*(Progeny[e].phenotype - oMaleTraitMean);
					if (FitMmean > 0)
						OmFitness[malecount] = W / FitMmean;
					else
						OmFitness[malecount] = 0;
					malecount++;
				}
			}
		}
		OTraitSD = sqrt(OTraitSD / dNo);
		ofTraitSD = sqrt(ofTraitSD / ONf);
		omTraitSD = sqrt(omTraitSD / ONm);
		FstdOtrait = new double[oFemN];
		MstdOtrait = new double[oMalN];
		malecount = 0;
		femalecount = 0;
		for (e = 0; e < ProgenyNum; e++)
		{
			if (Progeny[e].Alive)
			{
				OGenVar = OGenVar + (Progeny[e].Genotype - OGenmean)*(Progeny[e].Genotype - OGenmean);
				OPhenVar = OPhenVar + (Progeny[e].phenotype - OPhenMean)*(Progeny[e].phenotype - OPhenMean);
				if (OTraitSD > 0)
					StdOTrait[e] = (Progeny[e].phenotype - OPhenMean) / OTraitSD;
				else
					StdOTrait[e] = 0;
				if (Progeny[e].IsFemale)
				{
					if (ofTraitSD > 0)
						FstdOtrait[femalecount] = (Progeny[e].phenotype - oFemTraitMean) / ofTraitSD;
					else
						FstdOtrait[femalecount] = 0;
					femalecount++;
				}
				else
				{
					if (omTraitSD > 0)
						MstdOtrait[malecount] = (Progeny[e].phenotype - oMaleTraitMean) / omTraitSD;
					else
						MstdOtrait[malecount] = 0;
					malecount++;
				}
			}
		}
		OGenVar = OGenVar / dNo;
		OPhenVar = OPhenVar / dNo;
		OHeritability = OGenVar / OPhenVar;
		//Calculate selection differential
		for (e = 0; e < ProgenyNum; e++)
		{
			if (Progeny[e].Alive){
				Osprime = Osprime + (StdOTrait[e] - OPhenMean)*(OFitness[e] - FitMean);
			}
		}
		for (e = 0; e < oFemN; e++)
			Ofprime = Ofprime + (FstdOtrait[e] - oFemTraitMean)*(OfFitness[e] - FitFmean);
		for (e = 0; e < oMalN; e++)
			Omprime = Omprime + (MstdOtrait[e] - oMaleTraitMean)*(OmFitness[e] - FitMmean);
		Osprime = Osprime / dNo;
		Ofprime = Ofprime / ONf;
		Omprime = Omprime / ONm;
		cout << "\nSelection differential: " << Osprime << " Total Num Offspring: " << ProgenyNum << " Num Surviving Offspring:" << dNo << '\n';
		malN = oMalN;
		femN = oFemN;
		offRatio = ONm / ONf;
		allphenavg = OPhenMean;
		malphenavg = oMaleTraitMean;
		femphenavg = oFemTraitMean;
		allDiff = Osprime;
		maleDiff = Omprime;
		femDiff = Ofprime;
		Heritability = OHeritability;
		GenVar = OGenVar;
		PhenVar = OPhenVar;

		delete[] StdOTrait;
		delete[] FstdOtrait;
		delete[] MstdOtrait;
		delete[] OFitness;
		delete[] OfFitness;
		delete[] OmFitness;
	}//end OffSum

	ld_info AdultPopLD(int chromosome_1, int chromosome_2, int allele_1, int allele_2)
	{
		//To do this, need to compare allele frequencies. 
		//D=x11-p1q1
		//x11 = observed freq of major allele in locus A and major allele in locus B
		//p1 = observed freq of major allele in locus A
		//q1 = observed freq of major allele in locus B
		//if D > 0, Dmax = min(p1q1, p2q2)
		//if D < 0, Dmax = min(p1q2, p2q1)
		//D' = D/Dmax
		ld_info result;
		double dNadult = PopulationSize;
		vector <int> num_allele_1;
		vector <int> num_allele_2;
		vector <double> freq_allele_1;
		vector <double> freq_allele_2;
		double **joint_alleles = new double *[NumAlleles];
		double **D = new double *[NumAlleles];
		double **Dmax = new double *[NumAlleles];
		int count_a = 0;
		int f, ff, fff, al, count;
		double Dprime, d_allele_avgs;
		int numA1B1 = 0;

		for (al = 0; al < NumAlleles; al++)
		{
			joint_alleles[al] = new double[NumAlleles];
			D[al] = new double[NumAlleles];
			Dmax[al] = new double[NumAlleles];
		}
		//figure out which loci are polymorphic
		for (al = 0; al < NumAlleles; al++)
		{
			num_allele_1.push_back(0);
			num_allele_2.push_back(0);
			freq_allele_1.push_back(0);
			freq_allele_2.push_back(0);
		}

		int maternal_1, maternal_2, paternal_1, paternal_2;

		for (fff = 0; fff < PopulationSize; fff++)
		{
			maternal_1 = Adult[fff].maternal[chromosome_1].LociArray[allele_1];
			maternal_2 = Adult[fff].maternal[chromosome_2].LociArray[allele_2];
			paternal_1 = Adult[fff].paternal[chromosome_1].LociArray[allele_1];
			paternal_2 = Adult[fff].paternal[chromosome_2].LociArray[allele_2];

			num_allele_1[maternal_1]++;
			num_allele_1[paternal_1]++;
			num_allele_2[maternal_2]++;
			num_allele_2[paternal_2]++;
			joint_alleles[maternal_1][maternal_2]++;
			joint_alleles[paternal_1][paternal_2]++;
		}//end for fff < PopulationSize

		int major_allele_1 = 0;
		int major_allele_2 = 0;
		for (al = 0; al < NumAlleles; al++)
		{
			freq_allele_1[al] = (num_allele_1[al]) / (2 * dNadult);
			freq_allele_2[al] = (num_allele_2[al]) / (2 * dNadult);
			if (freq_allele_1[al] > freq_allele_1[major_allele_1])
				major_allele_1 = al;
			if (freq_allele_2[al] > freq_allele_2[major_allele_2])
				major_allele_2 = al;
			for (f = 0; f < NumAlleles; f++)
				joint_alleles[al][f] = joint_alleles[al][f] / (2 * dNadult);
		}
		PopPrimFreq[chromosome_1].LociOnChrom[allele_1] = freq_allele_1[major_allele_1];
		PopMajorAllele[chromosome_1].LociOnChrom[allele_1] = major_allele_1;
		PopPrimFreq[chromosome_2].LociOnChrom[allele_2] = freq_allele_2[major_allele_2];
		PopMajorAllele[chromosome_2].LociOnChrom[allele_2] = major_allele_2;


		if (freq_allele_1[major_allele_1] != 1 && freq_allele_2[major_allele_1] != 1)//if it's polymorphic
		{
			d_allele_avgs = 0;
			count = 0;
			for (f = 0; f < NumAlleles; f++)
			{
				for (ff = 0; ff < NumAlleles; ff++)
				{
					if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
					{
						D[f][ff] = joint_alleles[f][ff] - freq_allele_1[f] * freq_allele_2[ff];
						d_allele_avgs = d_allele_avgs + fabs(D[f][ff]);
						count++;
						if (D[f][ff] < 0)
							Dmax[f][ff] = min(freq_allele_1[f] * freq_allele_2[ff], (1 - freq_allele_1[f])*(1 - freq_allele_2[ff]));
						else
							Dmax[f][ff] = min((1 - freq_allele_1[f])*freq_allele_2[ff], freq_allele_1[f] * (1 - freq_allele_2[ff]));
					}
				}
			}
			double dcount = count;
			d_allele_avgs = d_allele_avgs / dcount;

			Dprime = 0;
			bool decentDmax = true;
			for (f = 0; f < NumAlleles; f++)
			{
				for (ff = 0; ff < NumAlleles; ff++)
				{
					if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
					{
						if (Dmax[f][ff] > 0)
							Dprime = Dprime + freq_allele_1[f] * freq_allele_2[ff] * fabs(D[f][ff]) / Dmax[f][ff];
						else
							decentDmax = false;

					}
				}
			}
			if (!decentDmax)
				Dprime = -5;

		}
		PopLD_a[chromosome_1].LociOnChrom[allele_1] = Dprime;

		result.d = d_allele_avgs;
		result.dprime = Dprime;



		for (f = 0; f < NumAlleles; f++)
			delete[] joint_alleles[f];
		for (f = 0; f < NumAlleles; f++)
			delete[] D[f];
		for (f = 0; f < NumAlleles; f++)
			delete[] Dmax[f];
		delete[] joint_alleles;
		delete[] D;
		delete[] Dmax;

		return result;
	}//end Adult Pop LD

	void SamplePop()
	{
		//Pull a random sample of adults and a sample of offspring
		//Sort adults into males and females
		//Calculate allele frequencies for males, females, and offspring
		//Compare using Fst approaches
		if (ProgSample == 0)
			ProgSample = ProgenyNum;
		SampledAdults = new int[pAdultSample];
		SampledOffspring = new int[ProgSample];
		bool* PopTaken;
		bool* OffTaken;
		PopTaken = new bool[PopulationSize];
		OffTaken = new bool[ProgenyNum];
		int s, ss, t;
		int rand1, rand2;
		double fhetratio, mhetratio, ohetratio;
		fhetratio = 0;
		mhetratio = 0;
		ohetratio = 0;
		double pfem, pmal, poff, padult;
		pfem = 0;
		pmal = 0;
		poff = 0;
		padult = 0;
		int femsamp, malesamp;
		femsamp = 0;
		malesamp = 0;
		int Nt = (pAdultSample + ProgSample);
		double adultN = pAdultSample;
		double femN = 0;
		double malN = 0;
		double offN = ProgSample;

		if (pAdultSample > PopulationSize)
			pAdultSample = PopulationSize;
		if (ProgSample > ProgenyNum)
			ProgSample = ProgenyNum;

		for (t = 0; t < PopulationSize; t++)
			PopTaken[t] = false;
		for (t = 0; t < ProgenyNum; t++)
			OffTaken[t] = false;
		for (s = 0; s < pAdultSample; s++)
		{//Loop through each adult sample and choose which individual it will be
			rand1 = randnum(PopulationSize);//need to sample without replacement
			if (PopTaken[rand1] == false){
				SampledAdults[s] = rand1;
				if (Adult[SampledAdults[s]].IsFemale && Adult[SampledAdults[s]].Alive)
					femsamp++;
				if (Adult[SampledAdults[s]].IsFemale == false && Adult[SampledAdults[s]].Alive)
					malesamp++;
				PopTaken[rand1] = true;
			}
			else{
				while (PopTaken[rand1] == true){
					rand1 = randnum(PopulationSize);
				}
				if (PopTaken[rand1] == false){
					SampledAdults[s] = rand1;
					if (Adult[SampledAdults[s]].IsFemale && Adult[SampledAdults[s]].Alive)
						femsamp++;
					if (Adult[SampledAdults[s]].IsFemale == false && Adult[SampledAdults[s]].Alive)
						malesamp++;
					PopTaken[rand1] = true;
				}
			}//end else
		}//end for(s)

		FemSample = femsamp;
		femN = femsamp;
		MaleSample = malesamp;
		malN = malesamp;

		for (ss = 0; ss < ProgSample; ss++){
			rand2 = randnum(ProgenyNum);//need to sample without replacement
			if (OffTaken[rand2] == false){
				SampledOffspring[ss] = rand2;
				OffTaken[rand2] = true;
			}
			else{
				while (OffTaken[rand2] == true){
					rand2 = randnum(ProgenyNum);
				}
				if (OffTaken[rand2] == false){
					SampledOffspring[ss] = rand2;
					OffTaken[rand2] = true;
				}
			}//end else
		}//end of for(ss)

		offN = ProgSample;

		delete[] OffTaken;
		delete[] PopTaken;
	}//end SamplePop

	void SampleFsts()
	{
		int f, ff, fff, al;
		int numGenotypes;
		int Apoly, Opoly, Mpoly;
		int locus;
		int numhets, mnumhets;
		int iNmales = 0;
		for (int m = 0; m < pAdultSample; m++)
		{
			if (!Adult[m].IsFemale)
				iNmales++;
		}
		double freqhets, Ofreqhets, dNGen, Mfreqhets;
		double dNadult = pAdultSample;
		double dNoff = ProgSample;
		double dNmale = iNmales;
		double dN = pAdultSample + ProgSample;
		int* numEach = new int[NumAlleles];
		int* TnumEach = new int[NumAlleles];
		int* mtnumEach = new int[NumAlleles];
		int* MnumEach = new int[NumAlleles];
		int* OnumEach = new int[NumAlleles];
		int* numHom = new int[NumAlleles];
		int* OnumHom = new int[NumAlleles];
		int* MnumHom = new int[NumAlleles];
		double* OfreqEach = new double[NumAlleles];
		double* freqEach = new double[NumAlleles];
		double* TfreqEach = new double[NumAlleles];
		double* mtfreqEach = new double[NumAlleles];
		double* MfreqEach = new double[NumAlleles];
		dNGen = (NumAlleles*(NumAlleles + 1)) / 2;
		numGenotypes = dNGen;
		double* freqHom = new double[NumAlleles];
		double* OfreqHom = new double[NumAlleles];//these track the freq of each homozygote
		double* MfreqHom = new double[NumAlleles];
		Apoly = 0;
		Opoly = 0;
		Mpoly = 0;
		int count_a = 0;
		int count_m = 0;
		int count_o = 0;
		//figure out which loci are polymorphic

		for (f = 0; f < NumChrom; f++)
		{
			locus = 0; //keep track of all which sampled locus we're on to reduce memory load.
			for (ff = 0; ff < NumMarkers; ff++)//treat this as if it's a real sample (don't know markers)
			{
				if (SampledLoci[f].LociOnChrom[ff] == 1 && locus < NumSampledLoci)
				{
					numhets = 0;
					mnumhets = 0;
					for (al = 0; al < NumAlleles; al++)
					{
						numEach[al] = 0;
						MnumEach[al] = 0;
						TnumEach[al] = 0;
						mtnumEach[al] = 0;
						mtfreqEach[al] = 0;
						TfreqEach[al] = 0;
						freqEach[al] = 0;
						MfreqEach[al] = 0;
						numHom[al] = 0;
						MnumHom[al] = 0;
					}

					for (fff = 0; fff < pAdultSample; fff++)
					{
						if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] != Adult[SampledAdults[fff]].paternal[f].LociArray[ff])
							numhets++;
						for (al = 0; al < NumAlleles; al++)
						{
							if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] == al){
								numEach[al]++;
								TnumEach[al]++;
							}
							if (Adult[SampledAdults[fff]].paternal[f].LociArray[ff] == al){
								numEach[al]++;
								TnumEach[al]++;
							}
							if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] == al &&
								Adult[SampledAdults[fff]].paternal[f].LociArray[ff] == al)
								numHom[al]++;
						}
						if (!Adult[fff].IsFemale)
						{
							if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] != Adult[SampledAdults[fff]].paternal[f].LociArray[ff])
								mnumhets++;
							for (al = 0; al < NumAlleles; al++)
							{
								if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] == al){
									MnumEach[al]++;
									mtnumEach[al]++;
								}
								if (Adult[SampledAdults[fff]].paternal[f].LociArray[ff] == al){
									MnumEach[al]++;
									mtnumEach[al]++;
								}
								if (Adult[SampledAdults[fff]].maternal[f].LociArray[ff] == al &&
									Adult[SampledAdults[fff]].paternal[f].LociArray[ff] == al)
									MnumHom[al]++;
							}
						}
					}//end for fff < pAdultsample


					int offnumhets = 0;
					for (al = 0; al < NumAlleles; al++)
					{
						OnumEach[al] = 0;
						OfreqEach[al] = 0;
						OnumHom[al] = 0;
					}

					for (fff = 0; fff < ProgSample; fff++)
					{
						if (Progeny[SampledOffspring[fff]].maternal[f].LociArray[ff] != Progeny[SampledOffspring[fff]].paternal[f].LociArray[ff])
							offnumhets++;
						for (al = 0; al < NumAlleles; al++)
						{
							if (Progeny[SampledOffspring[fff]].maternal[f].LociArray[ff] == al){
								OnumEach[al]++;
								TnumEach[al]++;
								mtnumEach[al]++;
							}
							if (Progeny[SampledOffspring[fff]].paternal[f].LociArray[ff] == al){
								OnumEach[al]++;
								TnumEach[al]++;
								mtnumEach[al]++;
							}
							if (Progeny[SampledOffspring[fff]].maternal[f].LociArray[ff] == al &&
								Progeny[SampledOffspring[fff]].paternal[f].LociArray[ff] == al)
								OnumHom[al]++;
						}
					}//end of 2nd fff
					int adultMajor = 0;
					int offMajor = 0;
					int malMajor = 0;
					count_a = 0;
					count_m = 0;
					count_o = 0;
					for (al = 0; al < NumAlleles; al++)
					{
						MfreqEach[al] = MnumEach[al] / (2 * dNmale);
						freqEach[al] = (numEach[al]) / (2 * dNadult);
						freqHom[al] = (numHom[al] + OnumHom[al]) / (dNadult + dNoff);
						MfreqHom[al] = MnumHom[al] / dNmale;
						OfreqHom[al] = OnumHom[al] / dNoff;
						OfreqEach[al] = OnumEach[al] / (2 * dNoff);
						TfreqEach[al] = TnumEach[al] / (2 * dN);
						mtfreqEach[al] = mtnumEach[al] / (2 * (dNmale + dNoff));
						if (freqEach[al] > freqEach[adultMajor])
							adultMajor = al;
						if (OfreqEach[al] > OfreqEach[offMajor])
							offMajor = al;
						if (MfreqEach[al] > MfreqEach[malMajor])
							malMajor = al;
						if (numEach[al] != 0)
							count_a++;
						if (mtnumEach[al] != 0)
							count_m++;
						if (OnumEach[al] != 0)
							count_o++;
					}

					if (freqEach[adultMajor] != 1 || QTLtracker[f].LociOnChrom[ff] == 1)
						Polymorphisms[f].LociOnChrom[locus] = 1;
					else
						Polymorphisms[f].LociOnChrom[locus] = 0;
					if (OfreqEach[offMajor] != 1 || QTLtracker[f].LociOnChrom[ff] == 1)
						OPolymorphisms[f].LociOnChrom[locus] = 1;
					else
						OPolymorphisms[f].LociOnChrom[locus] = 0;

					AFreq1[f].LociOnChrom[locus] = freqEach[adultMajor];
					MFreq1[f].LociOnChrom[locus] = MfreqEach[malMajor];
					OFreq1[f].LociOnChrom[locus] = OfreqEach[offMajor];
					AMajor[f].LociOnChrom[locus] = adultMajor;
					MMajor[f].LociOnChrom[locus] = malMajor;
					OMajor[f].LociOnChrom[locus] = offMajor;
					A_NumAlleles[f].LociOnChrom[locus] = count_a;
					M_NumAlleles[f].LociOnChrom[locus] = count_m;
					O_NumAlleles[f].LociOnChrom[locus] = count_o;

					MalHet[f].LociOnChrom[locus] = 0;
					TotHet[f].LociOnChrom[locus] = 0;
					OffHet[f].LociOnChrom[locus] = 0;
					AdultHet[f].LociOnChrom[locus] = 0;
					for (al = 0; al < NumAlleles; al++)
					{
						//calculate expected heterozygosity
						TotHet[f].LociOnChrom[locus] = TotHet[f].LociOnChrom[locus] + (TfreqEach[al] * TfreqEach[al]);
						OffHet[f].LociOnChrom[locus] = OffHet[f].LociOnChrom[locus] + (OfreqEach[al] * OfreqEach[al]);
						AdultHet[f].LociOnChrom[locus] = AdultHet[f].LociOnChrom[locus] + (freqEach[al] * freqEach[al]);
						MalHet[f].LociOnChrom[locus] = MalHet[f].LociOnChrom[locus] + (MfreqEach[al] * MfreqEach[al]);
					}
					TotHet[f].LociOnChrom[locus] = 1 - TotHet[f].LociOnChrom[locus];
					OffHet[f].LociOnChrom[locus] = 1 - OffHet[f].LociOnChrom[locus];
					AdultHet[f].LociOnChrom[locus] = 1 - AdultHet[f].LociOnChrom[locus];
					MalHet[f].LociOnChrom[locus] = 1 - MalHet[f].LociOnChrom[locus];
					HetS[f].LociOnChrom[locus] = ((AdultHet[f].LociOnChrom[locus] * dNadult) + (OffHet[f].LociOnChrom[locus] * dNoff)) / dN;
					MaleH[f].LociOnChrom[locus] = ((MalHet[f].LociOnChrom[locus] * dNmale) + (OffHet[f].LociOnChrom[locus] * dNoff)) / (dNmale + dNoff);
					Ofreqhets = offnumhets / (dNGen*dNoff);//observed het
					freqhets = numhets / (dNGen*dNadult);
					Mfreqhets = mnumhets / (dNGen*dNmale);
					locus++;
				}//if Sampled
			}//end ff NumMarkers
		}//end for f chrom
		//Calculating two types of Fsts
		//One is the classical Fst = (Ht-Hs)/Ht
		Opoly = 0;
		Mpoly = 0;

		for (f = 0; f < NumChrom; f++)
		{
			locus = 0;
			for (ff = 0; ff < NumMarkers; ff++)
			{
				if (SampledLoci[f].LociOnChrom[ff] == 1)
				{
					F_Marker.Locus[f].LociOnChrom[locus] = ff;
					if (Polymorphisms[f].LociOnChrom[locus] == 1 || OPolymorphisms[f].LociOnChrom[locus] == 1)
					{
						if (TotHet[f].LociOnChrom[locus] != 0)
						{
							F_Marker.FstValue[f].LociOnChrom[locus] = (TotHet[f].LociOnChrom[locus] - HetS[f].LociOnChrom[locus]) / TotHet[f].LociOnChrom[locus];
							FstTracker[f].LociOnChrom[locus] = 1;
							Opoly++;
						}
						else
							F_Marker.FstValue[f].LociOnChrom[locus] = 0;
					}
					locus++;
				}
			}
		}
		NumPolymorphicLoci = Opoly;
		MOnumPolyLoci = Mpoly;
		delete[] numEach;
		delete[] OnumEach;
		delete[] MnumEach;
		delete[] MfreqEach;
		delete[] MfreqHom;
		delete[] freqHom;
		delete[] OfreqHom;
		delete[] TfreqEach;
		delete[] freqEach;
		delete[] OfreqEach;
		delete[] TnumEach;
		delete[] numHom;
		delete[] OnumHom;
		delete[] MnumHom;
		delete[] mtnumEach;
		delete[] mtfreqEach;
	}//end SampleFsts

	void CalcFstChiSq(Fst fst)
	{
		int whichchromosome, whichmarker, locus;
		// We need the number of alleles per locus

		double actualnumberofalleles;
		double degreesoffreedom;
		double pvalue;
		double chisqrteststat;
		double Ntotalsamplesize;
		Ntotalsamplesize = pAdultSample + ProgSample;
		for (whichchromosome = 0; whichchromosome < NumChrom; whichchromosome++)
		{
			locus = 0;
			for (whichmarker = 0; whichmarker < NumMarkers; whichmarker++)
			{
				if (SampledLoci[whichchromosome].LociOnChrom[whichmarker] == 1 && locus < NumSampledLoci)
				{
					actualnumberofalleles = 0;
					//This could be improved--what if adults and progeny have different alleles?
					if (A_NumAlleles[whichchromosome].LociOnChrom[locus] > 0)
					{
						if (O_NumAlleles[whichchromosome].LociOnChrom[locus] > 0)
						{
							if (A_NumAlleles[whichchromosome].LociOnChrom[locus] >= O_NumAlleles[whichchromosome].LociOnChrom[locus])
								actualnumberofalleles = A_NumAlleles[whichchromosome].LociOnChrom[locus];
							else
								actualnumberofalleles = O_NumAlleles[whichchromosome].LociOnChrom[locus];
						}
					}
					if (fst.FstValue[whichchromosome].LociOnChrom[locus] >= 0 &&
						fst.FstValue[whichchromosome].LociOnChrom[locus] <= 1)
					{
						chisqrteststat = 2 * Ntotalsamplesize*fst.FstValue[whichchromosome].LociOnChrom[locus] * (actualnumberofalleles - 1);
						degreesoffreedom = (actualnumberofalleles - 1);
						pvalue = 1 - chisqr(chisqrteststat, degreesoffreedom);
					}
					else
						pvalue = 1.01;

					fst.PValue[whichchromosome].LociOnChrom[locus] = pvalue;
					Psig[whichchromosome].LociOnChrom[locus] = pvalue;
					locus++;
				}
			}//whichmarker
		}
	}//end CalcFstChiSq

	void BenjaminiHochbergFDR(LociTracker* pvalues)
	{
		int p;
		double k;
		LociTracker* swapHolder = new LociTracker[NumChrom];
		double* CriticalValues = new double[NumChrom];
		for (int i = 0; i < NumChrom; i++)
			swapHolder[i].LociOnChrom = new double[1];
		int last = NumSampledLoci - 1;
		double m = NumSampledLoci;
		double FDRstat;
		//Sort p-values from smallest to largest
		k = 0;
		for (int chrom = 0; chrom < NumChrom; chrom++)
		{
			for (int counter = NumSampledLoci - 1; counter >= 0; counter--)
			{
				for (int index = 0; index < last; index++)
				{
					if (pvalues[chrom].LociOnChrom[index] > pvalues[chrom].LociOnChrom[index + 1])
					{//basic swap
						swapHolder[chrom].LociOnChrom[0] = pvalues[chrom].LociOnChrom[index + 1];
						//move index to index + 1
						pvalues[chrom].LociOnChrom[index + 1] = pvalues[chrom].LociOnChrom[index];
						//move swapholder to index
						pvalues[chrom].LociOnChrom[index] = swapHolder[chrom].LociOnChrom[0];
					}
				}//end index
				last--;
			}//end counter

			for (p = 0; p < NumSampledLoci; p++)
			{
				FDRstat = (p / m)*SigValue;
				if (pvalues[chrom].LociOnChrom[p] <= FDRstat)
					k = pvalues[chrom].LociOnChrom[p];
			}
			CriticalValues[chrom] = k;
			//cout<<k<<'\n';
		}//end chrom
		for (p = 0; p < NumChrom - 1; p++)
		{
			if (CriticalValues[p] <= CriticalValues[p + 1])
				CriticalValue = CriticalValues[p];
			else
				CriticalValue = CriticalValues[p + 1];
		}
		FDR95 = CriticalValue;
		for (int l = 0; l < NumChrom; l++)
			delete[] swapHolder[l].LociOnChrom;
		delete[] swapHolder;
		delete[] CriticalValues;
	}//end FDR

	void CalcLD()
	{
		//To do this, need to compare allele frequencies. 
		//D=x11-p1q1
		//x11 = observed freq of major allele in locus A and major allele in locus B
		//p1 = observed freq of major allele in locus A
		//q1 = observed freq of major allele in locus B
		//if D > 0, Dmax = min(p1q1, p2q2)
		//if D < 0, Dmax = min(p1q2, p2q1)
		//D' = D/Dmax
		double D, Dprime, Dmax;
		int numA1B1 = 0;
		int A1, B1, Alocus, Blocus, sec_locus;
		double x11, p1, q1, p2, q2;
		double dNa = pAdultSample;
		double dNo = ProgSample;
		int c, l, a, o;
		int numMales = 0;
		for (o = 0; o < pAdultSample; o++)
		{
			if (!Adult[SampledAdults[o]].IsFemale)
				numMales++;
		}
		double dNm = numMales;
		bool located;
		for (c = 0; c < NumChrom; c++)
		{
			Alocus = 0;
			for (l = 0; l < NumMarkers; l++)
			{
				if (SampledLoci[c].LociOnChrom[l] == 1)
				{
					located = false; //find the next sampled locus!
					for (int locate = l; locate < NumMarkers; locate++)
					{
						if (located == false)
						{
							if (SampledLoci[c].LociOnChrom[locate] == 1)
							{
								sec_locus = locate;
								located = true;
							}
						}
					}
					if (Alocus == NumSampledLoci - 1)
					{
						ALD[c].LociOnChrom[Alocus] = 0;
						OLD[c].LociOnChrom[Alocus] = 0;
						MLD[c].LociOnChrom[Alocus] = 0;
					}
					else
					{
						Blocus = Alocus + 1;
						if (Polymorphisms[c].LociOnChrom[Alocus] == 1)
						{
							A1 = AMajor[c].LociOnChrom[Alocus];
							B1 = AMajor[c].LociOnChrom[Blocus];
							p1 = AFreq1[c].LociOnChrom[Alocus];
							p2 = 1 - p1;
							q1 = AFreq1[c].LociOnChrom[Blocus];
							q2 = 1 - q1;
							numA1B1 = 0;
							for (a = 0; a < pAdultSample; a++)
							{
								if (Adult[SampledAdults[a]].maternal[c].LociArray[l] == A1
									&& Adult[SampledAdults[a]].maternal[c].LociArray[sec_locus] == B1)
									numA1B1++;
								if (Adult[SampledAdults[a]].paternal[c].LociArray[l] == A1
									&& Adult[SampledAdults[a]].paternal[c].LociArray[sec_locus] == B1)
									numA1B1++;
							}
							x11 = numA1B1 / (dNa * 2);
							D = x11 - (p1*q1);
							if (D < 0)
							{
								if (-(p1*q1) > -(p2*q2))
									Dmax = -p1*q1;
								else
									Dmax = -p2*q2;
							}
							if (D > 0)
							{
								if ((p1*q2) < (p2*q1))
									Dmax = p1*q2;
								else
									Dmax = p2*q1;
							}
							if (D == 0)
								Dmax = 0;
							if (Dmax != 0 && D != 0)
								Dprime = D / Dmax;
							else
								Dprime = 0;
							ALD[c].LociOnChrom[Alocus] = Dprime;
						}
						else
							ALD[c].LociOnChrom[Alocus] = -5;
						//Now for Offspring:
						if (OPolymorphisms[c].LociOnChrom[Alocus] == 1)
						{
							A1 = OMajor[c].LociOnChrom[Alocus];
							B1 = OMajor[c].LociOnChrom[Blocus];
							p1 = OFreq1[c].LociOnChrom[Alocus];
							p2 = 1 - p1;
							q1 = OFreq1[c].LociOnChrom[Blocus];
							q2 = 1 - q1;
							numA1B1 = 0;
							for (a = 0; a < ProgSample; a++)
							{
								if (Progeny[SampledOffspring[a]].maternal[c].LociArray[l] == A1
									&& Progeny[SampledOffspring[a]].maternal[c].LociArray[sec_locus] == B1)
									numA1B1++;
								if (Progeny[SampledOffspring[a]].paternal[c].LociArray[l] == A1
									&& Progeny[SampledOffspring[a]].paternal[c].LociArray[sec_locus] == B1)
									numA1B1++;
							}
							x11 = numA1B1 / (dNo * 2);
							D = x11 - (p1*q1);
							if (D < 0)
							{
								if ((p1*q1) < (p2*q2))
									Dmax = p1*q1;
								else
									Dmax = p2*q2;
							}
							if (D > 0)
							{
								if ((p1*q2) < (p2*q1))
									Dmax = p1*q2;
								else
									Dmax = p2*q1;
							}
							if (D == 0)
								Dmax = 0;
							if (Dmax != 0 && D != 0)
								Dprime = D / Dmax;
							else
								Dprime = 0;
							OLD[c].LociOnChrom[Alocus] = Dprime;
						}
						else
							OLD[c].LociOnChrom[Alocus] = -5;
					}//end of if/else
					Alocus++;
				}//sampled
			}//end loci
		}//end chrom
	}//end CalcLD

	void LongDistLD()
	{
		//actual population, not the sample
		//every possible pairwise comparison for the NumLD loci
		//first, calculate the allele freqs
		double dNadult = PopulationSize;
		int* numEach = new int[NumAlleles];
		double* freqEach = new double[NumAlleles];
		LociTracker* Adult_Poly = new LociTracker[NumChrom];
		for (int l = 0; l < NumChrom; l++)
			Adult_Poly[l].LociOnChrom = new double[NumMarkers];
		int count_a = 0;
		int f, ff, fff, al, count;
		int rand_locus1, rand_chrom1, rand_locus2, rand_chrom2;
		double D, Dprime, Dmax;
		int numA1B1 = 0;
		int A1, B1;
		double x11, p1, q1, p2, q2;
		int a;
		//figure out which loci are polymorphic
		for (f = 0; f < NumChrom; f++)
		{
			for (ff = 0; ff < NumMarkers; ff++)
			{
				for (al = 0; al < NumAlleles; al++)
				{
					numEach[al] = 0;
					freqEach[al] = 0;
				}

				for (fff = 0; fff < PopulationSize; fff++)
				{
					for (al = 0; al < NumAlleles; al++)
					{
						if (Adult[fff].maternal[f].LociArray[ff] == al)
							numEach[al]++;
						if (Adult[fff].paternal[f].LociArray[ff] == al)
							numEach[al]++;
					}

				}//end for fff < PopulationSize

				int adultMajor = 0;
				for (al = 0; al < NumAlleles; al++)
				{
					freqEach[al] = (numEach[al]) / (2 * dNadult);
					if (freqEach[al] > freqEach[adultMajor])
						adultMajor = al;
				}
				PopPrimFreq[f].LociOnChrom[ff] = freqEach[adultMajor];
				PopMajorAllele[f].LociOnChrom[ff] = adultMajor;
				if (freqEach[adultMajor] != 1)
					Adult_Poly[f].LociOnChrom[ff] = 1;
				else
					Adult_Poly[f].LociOnChrom[ff] = 0;

			}//end NumMarkers
		}//end NumChrom

		count = 0;
		meanLD = 0;
		avg_longdist_d = 0;
		while (count < NumLD)
		{
			Dprime = 0;
			rand_chrom1 = randnum(NumChrom);
			rand_locus1 = randnum(NumMarkers);
			rand_chrom2 = randnum(NumChrom);
			rand_locus2 = randnum(NumMarkers);
			while (rand_chrom1 == rand_chrom2 && rand_locus1 == rand_locus2)
			{
				rand_chrom2 = randnum(NumChrom);
				rand_locus2 = randnum(NumMarkers);
			}

			if (Adult_Poly[rand_chrom1].LociOnChrom[rand_locus1] == 1 && Adult_Poly[rand_chrom2].LociOnChrom[rand_locus2] == 1)
			{
				A1 = PopMajorAllele[rand_chrom1].LociOnChrom[rand_locus1];
				B1 = PopMajorAllele[rand_chrom2].LociOnChrom[rand_locus2];
				p1 = PopPrimFreq[rand_chrom1].LociOnChrom[rand_locus1];
				p2 = 1 - p1;
				q1 = PopPrimFreq[rand_chrom2].LociOnChrom[rand_locus2];
				q2 = 1 - q1;
				numA1B1 = 0;
				for (a = 0; a < PopulationSize; a++)
				{
					if (Adult[a].maternal[rand_chrom1].LociArray[rand_locus1] == A1
						&& Adult[a].maternal[rand_chrom2].LociArray[rand_chrom2] == B1)
						numA1B1++;
					if (Adult[a].paternal[rand_chrom1].LociArray[rand_locus1] == A1
						&& Adult[a].paternal[rand_chrom2].LociArray[rand_chrom2] == B1)
						numA1B1++;
				}
				x11 = numA1B1 / (dNadult * 2);
				D = x11 - (p1*q1);
				if (D < 0)
				{
					if (-(p1*q1) > -(p2*q2))
						Dmax = -p1*q1;
					else
						Dmax = -p2*q2;
				}
				if (D > 0)
				{
					if ((p1*q2) < (p2*q1))
						Dmax = p1*q2;
					else
						Dmax = p2*q1;
				}
				if (D == 0)
					Dmax = 0;
				if (Dmax != 0 && D != 0)
					Dprime = D / Dmax;
				else
					Dprime = 0;
				meanLD = meanLD + Dprime;
				avg_longdist_d = avg_longdist_d + D;
				count++;
			}//end if
		}//end NumLD
		//cout<<meanLD<<'\t'<<count<<'\n';
		meanLD = meanLD / NumLD;
		avg_longdist_d = avg_longdist_d / NumLD;

		for (int c = 0; c < NumChrom; c++)
		{
			long_dist_dprime.push_back(0);
		}

		for (int c = 0; c < NumChrom; c++)
		{
			count = 0;
			while (count < 100)
			{
				Dprime = 0;
				rand_locus1 = randnum(NumMarkers);
				rand_locus2 = randnum(NumMarkers);

				if (Adult_Poly[c].LociOnChrom[rand_locus1] == 1 && Adult_Poly[c].LociOnChrom[rand_locus2] == 1)
				{
					A1 = PopMajorAllele[c].LociOnChrom[rand_locus1];
					B1 = PopMajorAllele[c].LociOnChrom[rand_locus2];
					p1 = PopPrimFreq[c].LociOnChrom[rand_locus1];
					p2 = 1 - p1;
					q1 = PopPrimFreq[c].LociOnChrom[rand_locus2];
					q2 = 1 - q1;
					numA1B1 = 0;
					for (a = 0; a < PopulationSize; a++)
					{
						if (Adult[a].maternal[c].LociArray[rand_locus1] == A1
							&& Adult[a].maternal[c].LociArray[rand_chrom2] == B1)
							numA1B1++;
						if (Adult[a].paternal[c].LociArray[rand_locus1] == A1
							&& Adult[a].paternal[c].LociArray[rand_chrom2] == B1)
							numA1B1++;
					}
					x11 = numA1B1 / (dNadult * 2);
					D = x11 - (p1*q1);
					if (D < 0)
					{
						if (-(p1*q1) > -(p2*q2))
							Dmax = -p1*q1;
						else
							Dmax = -p2*q2;
					}
					if (D > 0)
					{
						if ((p1*q2) < (p2*q1))
							Dmax = p1*q2;
						else
							Dmax = p2*q1;
					}
					if (D == 0)
						Dmax = 0;
					if (Dmax != 0 && D != 0)
						Dprime = D / Dmax;
					else
						Dprime = 0;
					long_dist_dprime[c] = long_dist_dprime[c] + Dprime;
					count++;
				}//end if
			}//while
		}//for each chrom
		//	SampleLD.close();
		for (int c = 0; c < NumChrom; c++)
		{
			long_dist_dprime[c] = long_dist_dprime[c] / 100;
		}
		delete[] numEach;
		delete[] freqEach;
		for (int l = 0; l < NumChrom; l++)
			delete[] Adult_Poly[l].LociOnChrom;
		delete[] Adult_Poly;
	}//end long distance LD

	void QTLFstCalcs()
	{
		int f, ff, fff;
		double temp;
		int iNmales = 0;
		for (int m = 0; m < pAdultSample; m++)
		{
			if (!Adult[SampledAdults[m]].IsFemale)
				iNmales++;
		}
		double dNa = pAdultSample;
		double dNo = ProgSample;
		double dNm = iNmales;
		int sizeA = pAdultSample * 2;
		int sizeO = ProgSample * 2;
		int sizeM = iNmales * 2;
		double* Aalleles = new double[sizeA];//keep track of loci at each locus
		double* Malleles = new double[sizeM];
		double* Oalleles = new double[sizeO];
		double* Talleles = new double[2 * (ProgSample + pAdultSample)];
		double* MTalleles = new double[2 * (ProgSample + iNmales)];
		int acount, mcount, ocount, tcount, mtcount;
		double* ADiffAlls = new double[sizeA];
		double* ODiffAlls = new double[sizeO];
		double* MDiffAlls = new double[sizeM];
		double* TDiffAlls = new double[2 * (ProgSample + pAdultSample)];
		double* MTDiffAlls = new double[2 * (ProgSample + iNmales)];
		int* AnumEach = new int[sizeA];
		int* OnumEach = new int[sizeO];
		int* MnumEach = new int[sizeM];
		int* TnumEach = new int[2 * (ProgSample + pAdultSample)];
		int* MTnumEach = new int[2 * (ProgSample + iNmales)];

		for (ff = 0; ff < NumChrom; ff++)
		{
			for (fff = 0; fff < NumQTLs; fff++)
			{
				acount = 0;
				mcount = 0;
				tcount = 0;
				mtcount = 0;
				for (f = 0; f < pAdultSample; f++)
				{
					Aalleles[acount] = Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff];
					acount++;
					Aalleles[acount] = Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff];
					acount++;
					Talleles[tcount] = Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff];
					tcount++;
					Talleles[tcount] = Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff];
					tcount++;
					if (!Adult[SampledAdults[f]].IsFemale)
					{
						Malleles[mcount] = Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff];
						mcount++;
						Malleles[mcount] = Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff];
						mcount++;
						MTalleles[mtcount] = Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff];
						mtcount++;
						MTalleles[mtcount] = Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff];
						mtcount++;
					}
				}//end adults
				ocount = 0;
				for (f = 0; f < ProgSample; f++)
				{
					Oalleles[ocount] = Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff];
					ocount++;
					Oalleles[ocount] = Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff];
					ocount++;
					MTalleles[mtcount] = Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff];
					mtcount++;
					MTalleles[mtcount] = Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff];
					mtcount++;
					Talleles[tcount] = Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff];
					tcount++;
					Talleles[tcount] = Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff];
					tcount++;
				}
				for (int clear = tcount; clear < 2 * (ProgSample + pAdultSample); clear++)
					Talleles[clear] = 0;
				for (int clear = mtcount; clear < 2 * (ProgSample + iNmales); clear++)
					MTalleles[clear] = 0;
				//Then sort the three 'alleles' arrays and count the number of times there's a change.
				for (int lower = 0; lower < sizeA - 1; lower++)
				{
					for (int upper = lower + 1; upper < sizeA; upper++)
					{
						if (Aalleles[lower] > Aalleles[upper])
						{//basic swap
							temp = Aalleles[upper];
							//move index to index + 1
							Aalleles[upper] = Aalleles[lower];
							//move swapholder to index
							Aalleles[lower] = temp;
						}
					}//end upper
				}//end lower--Adults

				for (int lower = 0; lower < sizeO - 1; lower++)
				{
					for (int upper = lower + 1; upper < sizeO; upper++)
					{
						if (Oalleles[lower] > Oalleles[upper])
						{//basic swap
							temp = Oalleles[upper];
							//move index to index + 1
							Oalleles[upper] = Oalleles[lower];
							//move swapholder to index
							Oalleles[lower] = temp;
						}
					}//end upper
				}//end lower--Prog

				for (int lower = 0; lower < sizeM; lower++)
				{
					for (int upper = lower + 1; upper < sizeM; upper++)
					{
						if (Malleles[lower] > Malleles[upper])
						{//basic swap
							temp = Malleles[upper];
							//move index to index + 1
							Malleles[upper] = Malleles[lower];
							//move swapholder to index
							Malleles[lower] = temp;
						}
					}//end upper
				}//end lower--males
				for (int lower = 0; lower < tcount - 1; lower++)
				{
					for (int upper = lower + 1; upper < tcount; upper++)
					{
						if (Talleles[lower] > Talleles[upper])
						{//basic swap
							temp = Talleles[upper];
							//move index to index + 1
							Talleles[upper] = Talleles[lower];
							//move swapholder to index
							Talleles[lower] = temp;
						}
					}//end upper
				}//end lower--tot
				for (int lower = 0; lower < mtcount - 1; lower++)
				{
					for (int upper = lower + 1; upper < mtcount; upper++)
					{
						if (MTalleles[lower] > MTalleles[upper])
						{//basic swap
							temp = MTalleles[upper];
							//move index to index + 1
							MTalleles[upper] = MTalleles[lower];
							//move swapholder to index
							MTalleles[lower] = temp;
						}
					}//end upper
				}//end lower--tot

				int Anum = 0;
				int Mnum = 0;
				int Onum = 0;
				int tnum = 0;
				int mtnum = 0;
				//				cout<<"Chrom"<<ff<<"Locus: "<<fff<<'\n'<<"Adults: ";
				ADiffAlls[Anum] = Aalleles[0];
				for (int tally = 0; tally < sizeA - 1; tally++)
				{
					if (Aalleles[tally] != Aalleles[tally + 1]){
						Anum++;
						ADiffAlls[Anum] = Aalleles[tally + 1];
					}
					//					cout<<Aalleles[tally]<<'\t';
				}
				//				cout<<'\n'<<"Offspring: ";
				ODiffAlls[Onum] = Oalleles[0];
				for (int tally = 0; tally < sizeO - 1; tally++)
				{
					if (Oalleles[tally] != Oalleles[tally + 1]){
						Onum++;
						ODiffAlls[Onum] = Oalleles[tally + 1];
					}
					//					cout<<Oalleles[tally]<<'\t';
				}
				//				cout<<"\nAdult Males: ";
				MDiffAlls[Mnum] = Malleles[0];
				for (int tally = 0; tally < sizeM - 1; tally++)
				{
					if (Malleles[tally] != Malleles[tally + 1]){
						Mnum++;
						MDiffAlls[Mnum] = Malleles[tally + 1];
					}
					//					cout<<Malleles[tally]<<'\t';
				}
				TDiffAlls[tnum] = Talleles[0];
				for (int tally = 0; tally < tcount - 1; tally++)
				{
					if (Talleles[tally] != Talleles[tally + 1] && Talleles[tally] != 0){
						tnum++;
						TDiffAlls[tnum] = Talleles[tally + 1];
					}
				}
				MTDiffAlls[mtnum] = MTalleles[0];
				for (int tally = 0; tally < mtcount - 1; tally++)
				{
					if (MTalleles[tally] != MTalleles[tally + 1] && MTalleles[tally] != 0){
						mtnum++;
						MTDiffAlls[mtnum] = MTalleles[tally + 1];
					}
				}
				//				cout<<'\n';
				ANumQTLAlleles[ff].LociOnChrom[fff] = Anum + 1;
				ONumQTLAlleles[ff].LociOnChrom[fff] = Onum + 1;
				MNumQTLAlleles[ff].LociOnChrom[fff] = Mnum + 1;
				TNumQTLAlleles[ff].LociOnChrom[fff] = tnum + 1;
				MTNumQTLAlleles[ff].LociOnChrom[fff] = mtnum + 1;

				//go through all the alleles that are different and count up how many there are. 
				int ttrack = 0;
				int mttrack = 0;
				for (int al = 0; al < sizeA; al++)
					AnumEach[al] = 0;
				for (int al = 0; al < sizeM; al++)
					MnumEach[al] = 0;
				for (int al = 0; al < sizeO; al++)
					OnumEach[al] = 0;
				for (int al = 0; al < (2 * (ProgSample + iNmales)); al++)
					MTnumEach[al] = 0;
				for (int al = 0; al < (2 * (pAdultSample + ProgSample)); al++)
					TnumEach[al] = 0;
				for (f = 0; f < pAdultSample; f++)
				{
					for (int al = 0; al < Anum + 1; al++)
					{
						if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] == ADiffAlls[al])
							AnumEach[al]++;
						if (Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff] == ADiffAlls[al])
							AnumEach[al]++;
					}
					if (!Adult[SampledAdults[f]].IsFemale)
					{
						for (int al = 0; al < Mnum + 1; al++)
						{
							if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] == MDiffAlls[al])
								MnumEach[al]++;
							if (Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff] == MDiffAlls[al])
								MnumEach[al]++;
						}
						for (int al = 0; al < mtnum + 1; al++)
						{
							if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] == MTDiffAlls[al]){
								MTnumEach[mttrack]++;
								mttrack++;
							}
							if (Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff] == MTDiffAlls[al]){
								MTnumEach[mttrack]++;
								mttrack++;
							}
						}
					}
					for (int al = 0; al < tnum + 1; al++)
					{
						if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] == TDiffAlls[al]){
							TnumEach[ttrack]++;
							ttrack++;
						}
						if (Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff] == TDiffAlls[al]){
							TnumEach[ttrack]++;
							ttrack++;
						}
					}
				}
				for (f = 0; f < ProgSample; f++)
				{
					for (int al = 0; al < Onum + 1; al++)
					{
						if (Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff] == ODiffAlls[al])
							OnumEach[al]++;
						if (Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff] == ODiffAlls[al])
							OnumEach[al]++;
					}
					for (int al = 0; al < tnum + 1; al++)
					{
						if (Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff] == TDiffAlls[al]){
							TnumEach[ttrack]++;
							ttrack++;
						}
						if (Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff] == TDiffAlls[al]){
							TnumEach[ttrack]++;
							ttrack++;
						}
					}
					for (int al = 0; al < mtnum + 1; al++)
					{
						if (Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff] == MTDiffAlls[al]){
							MTnumEach[mttrack]++;
							mttrack++;
						}
						if (Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff] == MTDiffAlls[al]){
							MTnumEach[mttrack]++;
							mttrack++;
						}
					}
				}
			}//end NumQTL
		}//end NumChrom
		//Now that we know how many alleles there are in each subpop for each locus, 
		//we can calculate allele frequencies and expected heterozygosity and eventually Fsts
		int AdultMaxNo = ANumQTLAlleles[0].LociOnChrom[0];
		int OffMaxNo = ONumQTLAlleles[0].LociOnChrom[0];
		int MaleMaxNo = MNumQTLAlleles[0].LociOnChrom[0];
		int TotMaxNo = TNumQTLAlleles[0].LociOnChrom[0];
		int MTotMaxNo = MTNumQTLAlleles[0].LociOnChrom[0];
		for (f = 0; f < NumChrom; f++)
		{
			for (ff = 0; ff < NumQTLs; ff++)
			{
				if (ANumQTLAlleles[f].LociOnChrom[ff] > AdultMaxNo)
					AdultMaxNo = ANumQTLAlleles[f].LociOnChrom[ff];
				if (ONumQTLAlleles[f].LociOnChrom[ff] > OffMaxNo)
					OffMaxNo = ONumQTLAlleles[f].LociOnChrom[ff];
				if (MNumQTLAlleles[f].LociOnChrom[ff] > MaleMaxNo)
					MaleMaxNo = MNumQTLAlleles[f].LociOnChrom[ff];
				if (TNumQTLAlleles[f].LociOnChrom[ff] > TotMaxNo)
					TotMaxNo = TNumQTLAlleles[f].LociOnChrom[ff];
				if (MTNumQTLAlleles[f].LociOnChrom[ff] > MTotMaxNo)
					MTotMaxNo = MTNumQTLAlleles[f].LociOnChrom[ff];
			}
		}

		double* OfreqEach = new double[OffMaxNo];
		double* AfreqEach = new double[AdultMaxNo];
		double* TfreqEach = new double[TotMaxNo];
		double* MTfreqEach = new double[MTotMaxNo];
		double* MfreqEach = new double[MaleMaxNo];

		int* AnumHom = new int[AdultMaxNo];
		int* OnumHom = new int[OffMaxNo];
		int* MnumHom = new int[MaleMaxNo];

		double* AfreqHom = new double[AdultMaxNo];//these track the freq of each homozygote
		double* OfreqHom = new double[OffMaxNo];
		double* MfreqHom = new double[MaleMaxNo];

		//Now loop through all the chromosomes and the loci and count up how many of each allele there are in the pop
		for (ff = 0; ff < NumChrom; ff++)
		{
			for (fff = 0; fff < NumQTLs; fff++)
			{
				for (f = 0; f < pAdultSample; f++)
				{
					for (int al = 0; al < ANumQTLAlleles[ff].LociOnChrom[fff]; al++)
					{
						if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] ==
							Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff])
							AnumHom[al]++;
					}
					if (!Adult[SampledAdults[f]].IsFemale)
					{
						for (int al = 0; al < MNumQTLAlleles[ff].LociOnChrom[fff]; al++)
						{
							if (Adult[SampledAdults[f]].maternal[ff].allelicEffects[fff] ==
								Adult[SampledAdults[f]].paternal[ff].allelicEffects[fff])
								MnumHom[al]++;
						}
					}

				}
				for (f = 0; f < ProgSample; f++)
				{
					for (int al = 0; al < ONumQTLAlleles[ff].LociOnChrom[fff]; al++)
					{
						if (Progeny[SampledOffspring[f]].maternal[ff].allelicEffects[fff] ==
							Progeny[SampledOffspring[f]].paternal[ff].allelicEffects[fff])
							OnumHom[al]++;
					}
				}

				//calculate frequencies
				for (int al = 0; al < AdultMaxNo; al++)
				{
					if (al < ANumQTLAlleles[ff].LociOnChrom[fff])
					{
						AfreqEach[al] = AnumEach[al] / (2 * dNa);
						AfreqHom[al] = AnumHom[al] / dNa;
					}
				}
				for (int al = 0; al < OffMaxNo; al++)
				{
					if (al < ONumQTLAlleles[ff].LociOnChrom[fff])
					{
						OfreqEach[al] = OnumEach[al] / (2 * dNo);
						OfreqHom[al] = OnumHom[al] / dNo;
					}
				}
				for (int al = 0; al < MaleMaxNo; al++)
				{
					if (al < MNumQTLAlleles[ff].LociOnChrom[fff])
					{
						MfreqEach[al] = MnumEach[al] / (2 * dNm);
						MfreqHom[al] = MnumHom[al] / dNm;
					}
				}
				for (int al = 0; al < TotMaxNo; al++)
				{
					if (al < TNumQTLAlleles[ff].LociOnChrom[fff])
						TfreqEach[al] = TnumEach[al] / (2 * (ProgSample + pAdultSample));
				}
				for (int al = 0; al < MTotMaxNo; al++)
				{
					if (al < MTNumQTLAlleles[ff].LociOnChrom[fff])
						MTfreqEach[al] = MTnumEach[al] / (2 * (ProgSample + iNmales));
				}
				//now calculate expected heterozygosities
				QMalHet[ff].LociOnChrom[fff] = 0;
				QOffHet[ff].LociOnChrom[fff] = 0;
				QAdultHet[ff].LociOnChrom[fff] = 0;
				QTotHet[ff].LociOnChrom[fff] = 0;
				for (int al = 0; al < ANumQTLAlleles[ff].LociOnChrom[fff]; al++)
					QAdultHet[ff].LociOnChrom[fff] = QAdultHet[ff].LociOnChrom[fff] + (AfreqEach[al] * AfreqEach[al]);
				for (int al = 0; al < MNumQTLAlleles[ff].LociOnChrom[fff]; al++)
					QMalHet[ff].LociOnChrom[fff] = QMalHet[ff].LociOnChrom[fff] + (MfreqEach[al] * MfreqEach[al]);
				for (int al = 0; al < ONumQTLAlleles[ff].LociOnChrom[fff]; al++)
					QOffHet[ff].LociOnChrom[fff] = QOffHet[ff].LociOnChrom[fff] + (OfreqEach[al] * OfreqEach[al]);
				for (int al = 0; al < TNumQTLAlleles[ff].LociOnChrom[fff]; al++)
					QTotHet[ff].LociOnChrom[fff] = QTotHet[ff].LociOnChrom[fff] + (TfreqEach[al] * TfreqEach[al]);
				QTotHet[ff].LociOnChrom[fff] = 1 - QTotHet[ff].LociOnChrom[fff];
				QOffHet[ff].LociOnChrom[fff] = 1 - QOffHet[ff].LociOnChrom[fff];
				QAdultHet[ff].LociOnChrom[fff] = 1 - QAdultHet[ff].LociOnChrom[fff];
				QMalHet[ff].LociOnChrom[fff] = 1 - QMalHet[ff].LociOnChrom[fff];
				QHs[ff].LociOnChrom[fff] = ((QAdultHet[ff].LociOnChrom[fff] * dNa) + (OffHet[ff].LociOnChrom[fff] * dNo)) / (dNa + dNo);
			}//fff NumQtls
		}//ff NumChrom

		//Calculate Fsts for adults->off and adult males->off
		for (f = 0; f < NumChrom; f++)
		{
			for (ff = 0; ff < NumQTLs; ff++)
			{
				if (QTotHet[f].LociOnChrom[ff] != 0)
					F_QTL.FstValue[f].LociOnChrom[ff] = (QTotHet[f].LociOnChrom[ff] - QHs[f].LociOnChrom[ff]) / QTotHet[f].LociOnChrom[ff];
				else
					F_QTL.FstValue[f].LociOnChrom[ff] = 0;
			}
		}

		delete[] Aalleles;
		delete[] Malleles;
		delete[] Oalleles;
		delete[] Talleles;
		delete[] MTalleles;
		delete[] AnumEach;
		delete[] MnumEach;
		delete[] OnumEach;
		delete[] TnumEach;
		delete[] MTnumEach;
		delete[] OfreqEach;
		delete[] AfreqEach;
		delete[] MfreqEach;
		delete[] TfreqEach;
		delete[] MTfreqEach;
		delete[] AnumHom;
		delete[] OnumHom;
		delete[] MnumHom;
		delete[] AfreqHom;
		delete[] OfreqHom;
		delete[] MfreqHom;
		delete[] ADiffAlls;
		delete[] ODiffAlls;
		delete[] MDiffAlls;
		delete[] TDiffAlls;
		delete[] MTDiffAlls;
	}//end QTL fst calcs

	void WeightedFst(Fst fsts, int sampletype)
	{
		/***************************************************************************************************************
		Code to weight Fsts to generate smooth genome-wide distributions
		Uses kernel-smoothing average, following the approach in Hohenlohe et al. (2010) and Catchen et al. (2013; Stacks)
		for each region at nucleotide c, contribution of statistic at position p to the region average is
		weighted by Gaussian function:
		exp((-(p-c)*(p-c))/(2*sigma*sigma)), where sigma = 150 kb (chosen arbitrarily)
		the distribution is truncated at 3sigma.
		Shifted the moving average by a step of 100 kb
		Further weight each statistic at each nucletide position by (nk-1) where nk = # alleles sampled at site k
		***************************************************************************************************************/
		//First: Create an index of loci that differ between the two pops
		//Sub1size = adult sample size, sub2size = prog sample size
		int w, x, y;
		double p, c;
		int numToavg = 0;
		int aPItoavg = 0, oPItoavg = 0;
		double dN, odN, adN, avgweights = 0;
		int start, end;
		double coef = 0;
		int locus;
		bool calc = false;
		double weightedFst, aweightPI, oweightPI, AavgPI = 0, OavgPI = 0;

		CalculateExpectedHeterozygosity();
		locus = 0;
		for (w = 0; w < NumChrom; w++)
		{
			for (x = 0; x < NumSampledLoci; x++)
			{
				if (fsts.FstValue[w].LociOnChrom[x] >= 0)
				{
					//set up weighting parameters
					avgweights = 0;
					AavgPI = 0;
					OavgPI = 0;
					aPItoavg = 0;
					oPItoavg = 0;
					start = x - 3 * sigma;
					if (start < 0)
						start = 0;
					end = x + 3 * sigma;
					if (end > NumSampledLoci)
						end = NumSampledLoci;
					//Calculate weighted heterozygosity at every nt for both pops
					for (y = start; y <= end; y++)
					{
						if (FstTracker[w].LociOnChrom[y] == 1)
						{
							p = y;
							c = x;
							if (sampletype == 0)
								aweightPI = AdultHet[w].LociOnChrom[y] * exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
							if (sampletype == 1)
								aweightPI = MalHet[w].LociOnChrom[y] * exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
							AavgPI = AavgPI + aweightPI;
							aPItoavg++;
							coef = coef + exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
							oweightPI = OffHet[w].LociOnChrom[y] * exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
							OavgPI = OavgPI + oweightPI;
							oPItoavg++;
						}
					}
					odN = oPItoavg;
					adN = aPItoavg;
					AavgPI = AavgPI / coef;
					OavgPI = OavgPI / coef;

					//weighted Fst
					numToavg = 0;
					coef = 0;
					for (y = start; y <= end; y++)
					{
						if (fsts.FstValue[w].LociOnChrom[y] >= 0)
						{
							p = y;
							c = x;
							weightedFst = fsts.FstValue[w].LociOnChrom[y] * exp(double(-1 * ((p - c)*(p - c))) / (double(2 * sigma*sigma)));
							avgweights = avgweights + weightedFst;
							coef = coef + exp(double((-1 * ((p - c)*(p - c)))) / (double(2 * sigma*sigma)));
							numToavg++;
						}
					}
					dN = numToavg;
					weightedFst = avgweights / coef;
					fsts.WeightedFst[w].LociOnChrom[x] = weightedFst;
				}//end weighting Fsts (if Polymorphisms == 1)
				else
					fsts.WeightedFst[w].LociOnChrom[x] = 0;
			}//end x
		}//end w
	}//end weighted fst

	void CalculateFstCIs(Fst fsts)
	{
		double meanFst, stdevFst, varFst;
		int m, mm, count;
		meanFst = 0;
		count = 0;
		for (m = 0; m < NumChrom; m++)
		{
			for (mm = 0; mm < NumSampledLoci; mm++)
			{
				if (fsts.FstValue[m].LociOnChrom[mm] >= 0)
				{
					meanFst = meanFst + fsts.FstValue[m].LociOnChrom[mm];
					count++;
				}
			}
		}
		meanFst = meanFst / count;

		varFst = 0;
		count = 0;
		for (m = 0; m < NumChrom; m++)
		{
			for (mm = 0; mm < NumSampledLoci; mm++)
			{
				varFst = varFst + (meanFst - fsts.FstValue[m].LociOnChrom[mm])*(meanFst - fsts.FstValue[m].LociOnChrom[mm]);
				count++;
			}
		}
		varFst = varFst / count;
		stdevFst = sqrt(varFst);

		SmoothedCritValue99 = meanFst + 2.57583 * stdevFst;
		SmoothedCritValue98 = meanFst + 2.32635 * stdevFst;
		SmoothedCritValue95 = meanFst + 1.95996 * stdevFst;
		SmoothedCritValue90 = meanFst + 1.64485 * stdevFst;
		SmoothedCritValue80 = meanFst + 1.28155 * stdevFst;

	}//end CalculateFstCIs

	void bootstrap_resampling(Fst fsts)
	{
		if (bootstrap == true)
		{
			int b, r;
			double avgweights, dN;
			vector <double> null_fsts;
			double p, c, weighted_fst;
			double num_to_avg, coef_sum;
			int nloci, count, sigmaloc;
			int rand_locus1, rand_chrom1;
			nloci = (3 * sigma) + (3 * sigma);
			c = nloci / 2;
			//cout<<bs_replicates;
			for (r = 0; r < bs_replicates; r++)
			{
				count = 0;
				avgweights = 0;
				num_to_avg = 0;
				coef_sum = 0;
				sigmaloc = 0;
				while (count < nloci)
				{
					rand_locus1 = randnum(NumSampledLoci);
					rand_chrom1 = randnum(NumChrom);
					if (F_Marker.FstValue[rand_chrom1].LociOnChrom[rand_locus1] > 0)
					{
						p = sigmaloc;
						weighted_fst =
							fsts.WeightedFst[rand_chrom1].LociOnChrom[rand_locus1] *
							exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
						avgweights = avgweights + weighted_fst;
						coef_sum = coef_sum + exp((-1 * ((p - c)*(p - c))) / (2 * sigma*sigma));
						num_to_avg++;
						count++;
					}
				}
				dN = num_to_avg;
				weighted_fst = avgweights / coef_sum;
				null_fsts.push_back(weighted_fst);
			}//end r
			sort(null_fsts.begin(), null_fsts.end());

			//once it's sorted, the upper confidence limit
			double meanFst, stdevFst, varFst;
			meanFst = 0;
			count = 0;
			for (b = 0; b < null_fsts.size(); b++)
			{
				if (null_fsts[b] >= 0)
				{
					meanFst = meanFst + null_fsts[b];
					count++;
				}
			}
			meanFst = meanFst / count;

			varFst = 0;
			count = 0;
			for (b = 0; b < null_fsts.size(); b++)
			{
				varFst = varFst + (meanFst - null_fsts[b])*(meanFst - null_fsts[b]);
				count++;
			}
			varFst = varFst / count;
			stdevFst = sqrt(varFst);


			bsCritValue99 = meanFst + 2.57583 * stdevFst;
			bsCritValue98 = meanFst + 2.32635 * stdevFst;
			bsCritValue95 = meanFst + 1.95996 * stdevFst;
			bsCritValue90 = meanFst + 1.64485 * stdevFst;
			bsCritValue80 = meanFst + 1.28155 * stdevFst;

		}
	}//end bootstrap

	void PeakDetector()
	{
		int p, pp, ppp;
		bs_spurious = bs_detected = bs_detected95 = bs_spurious95 = 0;
		ci_spurious = ci_detected = ci_spurious95 = ci_detected95 = 0;
		fdr_spur = fdr_det = 0;
		int low, high;
		bool peak;

		for (p = 0; p < NumChrom; p++)
		{
			for (pp = 0; pp < NumSampledLoci; pp++)
			{
				if (F_Marker.WeightedFst[p].LociOnChrom[pp] > F_Marker.WeightedFst[p].LociOnChrom[pp - 1] &&
					F_Marker.WeightedFst[p].LociOnChrom[pp] > F_Marker.WeightedFst[p].LociOnChrom[pp + 1])
				{//then it's a peak, so is it spurious or not?
					low = pp - peak_window;
					high = pp + peak_window;
					peak = false;
					for (ppp = low; ppp < high; ppp++)
					{
						if (QTLtracker[p].LociOnChrom[ppp] == 1)
							peak = true;//then it's a real peak
					}//ppp
					if (F_Marker.WeightedFst[p].LociOnChrom[pp] >= bsCritValue99)
					{
						if (peak == true)
							bs_detected++;
						else
							bs_spurious++;
					}
					if (F_Marker.WeightedFst[p].LociOnChrom[pp] >= bsCritValue95)
					{
						if (peak == true)
							bs_detected95++;
						else
							bs_spurious95++;
					}//bs95
					if (F_Marker.WeightedFst[p].LociOnChrom[pp] >= SmoothedCritValue99)
					{
						if (peak == true)
							ci_detected++;
						else
							ci_spurious++;
					}
					if (F_Marker.WeightedFst[p].LociOnChrom[pp] >= SmoothedCritValue95)
					{
						if (peak == true)
							ci_detected95++;
						else
							ci_spurious95++;
					}
					if (F_Marker.WeightedFst[p].LociOnChrom[pp] >= FDR95)
					{
						if (peak == true)
							fdr_det++;
						else
							fdr_spur++;

					}
				}//end of its peak
			}//end Num Loci Sampled
		}
	}//end peak detector

	//from Adam's (trying to make mine work!!)
	void CalculateAdultMarkerAlleleFreqs(int whichchromosome, int marker, int ref, bool includeDead)
	{
		int caa;
		double *tempallelefreq = new double[NumAlleles];
		for (caa = 0; caa < NumAlleles; caa++)
			tempallelefreq[caa] = 0;

		double allelecounter = 0;
		for (caa = 0; caa < pAdultSample; caa++)
		{
			if (Adult[SampledAdults[caa]].Alive || includeDead)
			{
				tempallelefreq[Adult[SampledAdults[caa]].maternal[whichchromosome].LociArray[marker]]++;
				tempallelefreq[Adult[SampledAdults[caa]].paternal[whichchromosome].LociArray[marker]]++;
				allelecounter = allelecounter + 2;
			}
		}
		for (caa = 0; caa < NumAlleles; caa++)
			tempallelefreq[caa] = tempallelefreq[caa] / allelecounter;
		int allelesperlocus = 0;
		for (caa = 0; caa < NumAlleles; caa++)
		{
			F_Marker.Adult_AlleleFreq[whichchromosome].allele_freq[ref].LociOnChrom[caa] = tempallelefreq[caa];
			if (tempallelefreq[caa] != 0)
				allelesperlocus++;
		}
		A_NumAlleles[whichchromosome].LociOnChrom[ref] = allelesperlocus;
		delete[] tempallelefreq;
	}//end CalculateAdultMarkerAlleleFreqs

	void CalculateProgenyMarkerAlleleFreqs(int whichchromosome, int marker, int ref, bool includeDead)
	{
		int caa;
		double *tempallelefreq = new double[NumAlleles];
		for (caa = 0; caa < NumAlleles; caa++)
			tempallelefreq[caa] = 0;

		double allelecounter = 0;
		for (caa = 0; caa < ProgSample; caa++)
		{
			if (Progeny[SampledOffspring[caa]].Alive || includeDead)
			{
				tempallelefreq[Progeny[SampledOffspring[caa]].maternal[whichchromosome].LociArray[marker]]++;
				tempallelefreq[Progeny[SampledOffspring[caa]].paternal[whichchromosome].LociArray[marker]]++;
				allelecounter = allelecounter + 2;
			}
		}
		for (caa = 0; caa < NumAlleles; caa++)
			tempallelefreq[caa] = tempallelefreq[caa] / allelecounter;
		int allelesperlocus = 0;
		for (caa = 0; caa < NumAlleles; caa++)
		{
			F_Marker.Progeny_AlleleFreq[whichchromosome].allele_freq[ref].LociOnChrom[caa] = tempallelefreq[caa];
			if (tempallelefreq[caa] != 0)
				allelesperlocus++;
		}
		O_NumAlleles[whichchromosome].LociOnChrom[ref] = allelesperlocus;
		delete[] tempallelefreq;
	}

	void POFst() //doing the per-chromosome calcs here, not in the .cpp file
	{
		double* averageallelefreqs;
		int locus;
		averageallelefreqs = new double[NumAlleles];
		double HsProgeny, HsAdults, Ht, Fst;
		int whichallele, whichchromosome, marker;
		for (whichchromosome = 0; whichchromosome < NumChrom; whichchromosome++)
		{
			locus = 0;
			for (marker = 0; marker < NumMarkers; marker++)
			{
				if (SampledLoci[whichchromosome].LociOnChrom[marker] == 1 && locus < NumSampledLoci)
				{
					F_Marker.Locus[whichchromosome].LociOnChrom[locus] = marker;
					CalculateAdultMarkerAlleleFreqs(whichchromosome, marker, locus, false);
					CalculateProgenyMarkerAlleleFreqs(whichchromosome, marker, locus, false);

					for (whichallele = 0; whichallele < NumAlleles; whichallele++)
					{
						averageallelefreqs[whichallele] = (F_Marker.Progeny_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele] +
							F_Marker.Adult_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele]) / 2;
					}
					HsProgeny = 1;
					HsAdults = 1;
					Ht = 1;

					for (whichallele = 0; whichallele < NumAlleles; whichallele++)
					{
						HsProgeny = HsProgeny - (F_Marker.Progeny_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele] *
							F_Marker.Progeny_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele]);
						HsAdults = HsAdults - (F_Marker.Adult_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele] *
							F_Marker.Adult_AlleleFreq[whichchromosome].allele_freq[locus].LociOnChrom[whichallele]);
						Ht = Ht - (averageallelefreqs[whichallele] * averageallelefreqs[whichallele]);
					}
					if (Ht > 0)
						Fst = 1 - (HsProgeny + HsAdults) / (2 * Ht);
					else
						Fst = -1.0;

					F_Marker.FstValue[whichchromosome].LociOnChrom[locus] = Fst;
					TotHet[whichchromosome].LociOnChrom[locus] = Ht;
					HetS[whichchromosome].LociOnChrom[locus] = (HsProgeny + HsAdults) / 2;
					AdultHet[whichchromosome].LociOnChrom[locus] = HsAdults;
					OffHet[whichchromosome].LociOnChrom[locus] = HsProgeny;
					locus++;
				}//if it's sampled
			}//marker
		}//chromosome
		delete[] averageallelefreqs;
	}//end POFst

	void CalculateExpectedHeterozygosity()
	{
		int f, ff, al, locus;
		for (f = 0; f < NumChrom; f++)
		{
			locus = 0;
			for (ff = 0; ff < NumMarkers; ff++)
			{
				if (SampledLoci[f].LociOnChrom[ff] == 1 && locus < NumSampledLoci)
				{
					CalculateAdultMarkerAlleleFreqs(f, ff, locus, false);
					CalculateProgenyMarkerAlleleFreqs(f, ff, locus, false);
					OffHet[f].LociOnChrom[locus] = 0;
					AdultHet[f].LociOnChrom[locus] = 0;
					for (al = 0; al < NumAlleles; al++)
					{
						//calculate expected heterozygosity
						AdultHet[f].LociOnChrom[locus] = AdultHet[f].LociOnChrom[locus] + (F_Marker.Adult_AlleleFreq[f].allele_freq[locus].LociOnChrom[al] * F_Marker.Adult_AlleleFreq[f].allele_freq[locus].LociOnChrom[al]);
						OffHet[f].LociOnChrom[locus] = OffHet[f].LociOnChrom[locus] + (F_Marker.Progeny_AlleleFreq[f].allele_freq[locus].LociOnChrom[al] * F_Marker.Progeny_AlleleFreq[f].allele_freq[locus].LociOnChrom[al]);
					}
					OffHet[f].LociOnChrom[locus] = 1 - OffHet[f].LociOnChrom[locus];
					AdultHet[f].LociOnChrom[locus] = 1 - AdultHet[f].LociOnChrom[locus];
					locus++;
				}
			}
		}
	}//end calc heterozygosity

	int DeterminePopSize()
	{
		int count, p;
		count = 0;
		for (p = 0; p < CarryingCapacity; p++)
		{
			if (Adult[p].Alive)
				count++;
		}
		return count;
	}//end DeterminePopSize

	void GeneticVariance(void)
	{
		int w, x, QTL;
		double dN;
		//Calculate genetic variances per locus
		for (w = 0; w < NumChrom; w++)
		{
			for (x = 0; x < NumMarkers; x++)
			{
				if (QTLtracker[w].LociOnChrom[x] == 1)
				{
					for (int q = 0; q < NumQTLs; q++)
					{
						if (Locations[w].LociOnChrom[q] == x)
						{
							QTL = q;
							GenetMeans[w].LociOnChrom[QTL] = 0;
							//Start with adults
							for (int a = 0; a < PopulationSize; a++)
								GenetMeans[w].LociOnChrom[QTL] = GenetMeans[w].LociOnChrom[QTL] + Adult[a].maternal[w].allelicEffects[QTL]
								+ Adult[a].paternal[w].allelicEffects[QTL];
							dN = PopulationSize;
							GenetMeans[w].LociOnChrom[QTL] = GenetMeans[w].LociOnChrom[QTL] / dN;
							GenetVars[w].LociOnChrom[QTL] = 0;
							for (int a = 0; a < PopulationSize; a++)
								GenetVars[w].LociOnChrom[QTL] = GenetVars[w].LociOnChrom[QTL] +
								((Adult[a].maternal[w].allelicEffects[QTL] + Adult[a].paternal[w].allelicEffects[QTL]) - GenetMeans[w].LociOnChrom[QTL])*
								((Adult[a].maternal[w].allelicEffects[QTL] + Adult[a].paternal[w].allelicEffects[QTL]) - GenetMeans[w].LociOnChrom[QTL]);
							GenetVars[w].LociOnChrom[QTL] = GenetVars[w].LociOnChrom[QTL] / dN;
							//Now do progeny
							double dNumOff = ProgenyNum;
							for (int o = 0; o < ProgenyNum; o++)
								GenetMeans[w].LociOnChrom[QTL] = GenetMeans[w].LociOnChrom[QTL] + Progeny[o].maternal[w].allelicEffects[QTL]
								+ Progeny[o].paternal[w].allelicEffects[QTL];
							GenetMeans[w].LociOnChrom[QTL] = GenetMeans[w].LociOnChrom[QTL] / dNumOff;
							GenetVars[w].LociOnChrom[QTL] = 0;
							for (int o = 0; o < PopulationSize; o++)
								GenetVars[w].LociOnChrom[QTL] = GenetVars[w].LociOnChrom[QTL] +
								((Progeny[o].maternal[w].allelicEffects[QTL] + Progeny[o].paternal[w].allelicEffects[QTL]) - GenetMeans[w].LociOnChrom[QTL])*
								((Progeny[o].maternal[w].allelicEffects[QTL] + Progeny[o].paternal[w].allelicEffects[QTL]) - GenetMeans[w].LociOnChrom[QTL]);
							GenetVars[w].LociOnChrom[QTL] = GenetVars[w].LociOnChrom[QTL] / dNumOff;
						}
					}
				}
			}
		}
	}

	double AdamsLD(int whichchromosome, int markerA, int markerB)
	{
		int iii, jjj;
		double *allelefreqA = new double[NumAlleles];
		double *allelefreqB = new double[NumAlleles];
		double **jointAB;
		jointAB = new double *[NumAlleles];
		for (iii = 0; iii < NumAlleles; iii++)
			jointAB[iii] = new double[NumAlleles];
		double count;
		double **Dij;
		Dij = new double *[NumAlleles];
		for (iii = 0; iii < NumAlleles; iii++)
			Dij[iii] = new double[NumAlleles];
		double **Dmax;
		Dmax = new double *[NumAlleles];
		for (iii = 0; iii < NumAlleles; iii++)
			Dmax[iii] = new double[NumAlleles];


		double Dprime;

		for (iii = 0; iii < NumAlleles; iii++)
		{
			allelefreqA[iii] = 0;
			allelefreqB[iii] = 0;
			for (jjj = 0; jjj < NumAlleles; jjj++)
			{
				jointAB[iii][jjj] = 0;
			}
		}

		count = 0;
		int maternalalleleA, maternalalleleB, paternalalleleA, paternalalleleB;
		for (iii = 0; iii < PopulationSize; iii++)
		{
			maternalalleleA = Adult[iii].maternal[whichchromosome].LociArray[markerA];
			maternalalleleB = Adult[iii].maternal[whichchromosome].LociArray[markerB];
			paternalalleleA = Adult[iii].paternal[whichchromosome].LociArray[markerA];
			paternalalleleB = Adult[iii].paternal[whichchromosome].LociArray[markerB];

			allelefreqA[maternalalleleA]++;
			allelefreqA[paternalalleleA]++;
			allelefreqB[maternalalleleB]++;
			allelefreqB[paternalalleleB]++;
			jointAB[maternalalleleA][maternalalleleB]++;
			jointAB[paternalalleleA][paternalalleleB]++;

			count = count + 2;
		}

		for (iii = 0; iii < NumAlleles; iii++)
		{
			allelefreqA[iii] = allelefreqA[iii] / count;
			allelefreqB[iii] = allelefreqB[iii] / count;
			for (jjj = 0; jjj < NumAlleles; jjj++)
				jointAB[iii][jjj] = jointAB[iii][jjj] / count;
		}

		for (iii = 0; iii < NumAlleles; iii++)
		{
			for (jjj = 0; jjj < NumAlleles; jjj++)
			{
				if (allelefreqA[iii] > 0 && allelefreqB[jjj] > 0)
				{
					Dij[iii][jjj] = jointAB[iii][jjj] - allelefreqA[iii] * allelefreqB[jjj];
					if (Dij[iii][jjj] < 0)
						Dmax[iii][jjj] = min(allelefreqA[iii] * allelefreqB[jjj], (1 - allelefreqA[iii])*(1 - allelefreqB[jjj]));
					else
						Dmax[iii][jjj] = min((1 - allelefreqA[iii])*allelefreqB[jjj], allelefreqA[iii] * (1 - allelefreqB[jjj]));

				}
			}
		}

		Dprime = 0;
		bool decentDmax = true;
		for (iii = 0; iii < NumAlleles; iii++)
		{
			for (jjj = 0; jjj < NumAlleles; jjj++)
			{
				if (allelefreqA[iii] > 0 && allelefreqB[jjj] > 0)
				{
					if (Dmax[iii][jjj] > 0)
						Dprime = Dprime + allelefreqA[iii] * allelefreqB[jjj] * fabs(Dij[iii][jjj]) / Dmax[iii][jjj];
					else
						decentDmax = false;

				}
			}
		}
		if (!decentDmax)
			Dprime = -5;

		for (iii = 0; iii < NumAlleles; iii++)
			delete[] jointAB[iii];
		for (iii = 0; iii < NumAlleles; iii++)
			delete[] Dij[iii];
		for (iii = 0; iii < NumAlleles; iii++)
			delete[] Dmax[iii];
		delete[] jointAB;
		delete[] Dij;
		delete[] Dmax;
		delete[] allelefreqA;
		delete[] allelefreqB;


		return Dprime;
	}
};
