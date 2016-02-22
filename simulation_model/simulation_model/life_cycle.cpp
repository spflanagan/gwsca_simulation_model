//Author: Sarah P. Flanagan
//Date: 26 March 2015
//This simulation model to test ability to detect selection within a single population.
//That population selects mates, produces offspring (complete with mutation & recombination),
//and then undergoes viability selection before randomly achieving adulthood.


#include "rand_nums.h"
#include "chi_square.h"
#include "populations.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <array>

using namespace std;

int main()
{
	int NumReps = 1;//this is one parameter to mess with here
	string basefilename = "null_1rep_440a_160j_6alleles";
	int vs_initial = 500;
	int vs_sampled = 500;
	int ss_initial = 0;
	int ss_sampled = 0;
	int end, g, reps, ld_count;
	vector <double> ci95_det, ci95_spu, ci99_det, ci99_spu;//for sexual selection
	vector <double> vs_ci95_det, vs_ci95_spu, vs_ci99_det, vs_ci99_spu;//vs = viability selection
	vector <double> bs95_det, bs95_spu, bs99_det, bs99_spu, fdr_spu, fdr_det;
	double ci95d_SE, ci95s_SE, ci99d_SE, ci99s_SE;//standard errors of spurious and detected peaks
	double ci95d_m, ci95s_m, ci99d_m, ci99s_m;//means of spurious and detected peaks
	double bs95d_SE, bs95s_SE, bs99d_SE, bs99s_SE, fdr_SE, fdrs_SE;//standard errors of spurious and detected peaks
	double bs95d_m, bs95s_m, bs99d_m, bs99s_m, fdr_m, fdrs_m;//means of spurious and detected peaks
	int peak_n = 0, total_QTL, vs_peak_n = 0;
	Population Pop;
	ofstream summ;
	ofstream Fsts;
	ofstream plot_Fsts;
	ofstream parameters;
	ofstream QTLs;
	ofstream peak_summ;
	ofstream ld, sig_thresh;
	string ldname = basefilename + "_LD.txt";
	string QTLname = basefilename + "_QTLs.txt";
	string parametername = basefilename + "_Parameters.txt";
	string summfilename = basefilename + "_SummStats.txt";
	string adultfilename = basefilename + "_AdultSummStats.txt";
	string fstfilename = basefilename + "_Fsts.txt";
	string fdrfilename = basefilename + "_plot_Fsts.txt";
	string peaksummname = basefilename + "_peakSumm.txt";
	string sig_threshold_name = basefilename + "_sig_thresholds.txt";
	ld.open(ldname.c_str());
	summ.open(summfilename.c_str());
	Fsts.open(fstfilename.c_str());
	plot_Fsts.open(fdrfilename.c_str());
	parameters.open(parametername.c_str());
	QTLs.open(QTLname.c_str());
	sig_thresh.open(sig_threshold_name.c_str());
	sig_thresh << "rep\tgen\tnum.fdr.real\tnum.fdr.spur\tnum.99.real\tnum.99.spur\tnum.95.real\tnum.95.spur\tnum.bs99.real\tnum.bs99.spur\tnum.bs95.real\tnum.bs95.spur";

	summ << "rep\tgen\tTimepoint\twhich\tMeanLD\t" << "GenVar\t" << "PhenVar\t" << "Heritability\t" << "PopulationSize\t" << "NumMales\t"
		<< "NumFem\t" << "SexRatio\t" << "Mean Trait Value\t" << "MaleTraitMean\t" << "FemTraitMean\t" << "Mating Diff\t"
		<< "Male m'\t" << "Fem m'" << "\tbs99_real\tbs99_spurious\tci99_real\tci99spurious" << '\n';
	parameters << "NumReps\tCarryingCapacity\tMaxFecundity\tMaxEncounters\tAdultSampleSize\tOffSampleSize\tNumGen\tNumSelGen\tNumLoci"
		<< "\tNumChrom\tNumQTLs\tMaxNumAlleles\tRecombRate\tMutRate\tMutVar\tEnvVar\tSigma\tbs_reps\tNumforLD\n";
	ci95s_m = ci95d_m = ci99d_m = ci99s_m = 0;
	bs95d_m = bs95s_m = bs99d_m = bs99s_m = fdr_m = fdrs_m = 0;
	bs95d_SE = bs95s_SE = bs99d_SE = bs99s_SE = fdr_SE = fdrs_SE = 0;
	vector <double> new_within_chrom_dprime;
	double mean_dp_ldistchrom = 0;
	//vector <double> adams_avg_dprime;
	ld_info returned_data;

	for (reps = 0; reps < NumReps; reps++)
	{
		Pop.SetParameters();
		Pop.Initialize(Pop.prior_qtl);

		cout << "Started rep " << reps + 1 << '\n';
		if (reps == 0)
		{
			ld << "rep\tgen\t" << "Avg_perchrom_ldDp\t" << "PairwiseD\t" << "Dprime\t"
				<< "long_dist_d\t" << "long_distDprime";
			total_QTL = Pop.NumQTLs*Pop.NumChrom;
			QTLs << "reps";
			for (int c = 0; c < Pop.NumChrom; c++)
			{
				new_within_chrom_dprime.push_back(0);
				//	adams_avg_dprime.push_back(0);
				ld << "\tchrom_" << c << "_D" << "\tchrom_" << c << "_ldDp";
				for (int cc = 0; cc < Pop.NumQTLs; cc++)
					QTLs << '\t' << c << "." << cc;
			}
			QTLs << '\n';
			ld << '\n';
			parameters << NumReps << '\t' << Pop.CarryingCapacity << '\t' << Pop.MaximumFecundity << '\t' << Pop.MaximumEncounters << '\t'
				<< Pop.pAdultSample << '\t' << Pop.ProgSample << '\t' << Pop.Generations << '\t' << Pop.NumSelGen << '\t' << Pop.NumMarkers << '\t'
				<< Pop.NumChrom << '\t' << Pop.NumQTLs << '\t' << Pop.NumAlleles << '\t' << Pop.RecombRate << '\t' << Pop.MutationRate << '\t'
				<< Pop.MutationalVariance << '\t' << Pop.EnvironmentalVariance << '\t' << Pop.sigma << '\t' << Pop.bs_replicates << '\t' << Pop.NumLD << '\n';
		}
		QTLs << reps;
		for (int c = 0; c < Pop.NumChrom; c++)
		{
			for (int cc = 0; cc < Pop.NumQTLs; cc++)
				QTLs << '\t' << Pop.Locations[c].LociOnChrom[cc];
		}
		QTLs << '\n';

		for (g = 0; g < Pop.Generations; g++)
		{
			Pop.PopulationSize = Pop.DeterminePopSize();
			Pop.Mating(ss_initial);//random mating if 0
			if (Pop.popExtinct == true)
			{
				cout << "Extinct Generation " << g + 1 << '\n';
				break;
			}

			Pop.Mutation();
			Pop.SelectionOnPhenotypes(vs_initial);//none if 0
			Pop.DensityRegulation();
			Pop.PopulationSize = Pop.DeterminePopSize();
			Pop.CalcMeanTraitValues();
			Pop.AdultSumm();
			double mean_new_longdist_d, mean_new_longdist_dprime;
			mean_new_longdist_d = mean_new_longdist_dprime = 0;
			Pop.AvgLD = 0;
			Pop.avg_pairwise_d = 0;
			Pop.avg_longdist_d = 0;
			mean_dp_ldistchrom = 0;
			ld_count = 0;
			while (ld_count < Pop.NumLD)
			{
				returned_data = Pop.AdultPopLD(randnum(Pop.NumChrom), randnum(Pop.NumChrom), randnum(Pop.NumMarkers), randnum(Pop.NumMarkers));
				if (returned_data.dprime != -5)
				{
					mean_new_longdist_d = mean_new_longdist_d + returned_data.d;
					mean_new_longdist_dprime = mean_new_longdist_dprime + returned_data.dprime;
					ld_count++;
				}
			}
			mean_new_longdist_d = mean_new_longdist_d / Pop.NumLD;
			mean_new_longdist_dprime = mean_new_longdist_dprime / Pop.NumLD;
			for (int c = 0; c < Pop.NumChrom; c++)
			{
				for (int cc = 0; cc < Pop.NumMarkers - 1; cc++)
				{//Pairwise
					returned_data = Pop.AdultPopLD(c, c, cc, cc + 1);
					Pop.avg_pairwise_d = Pop.avg_pairwise_d + returned_data.d;
					Pop.AvgLD = Pop.AvgLD + returned_data.dprime;
				}
				//adams_avg_dprime[c] = 0;
				Pop.avg_d[c] = 0;
				ld_count = 0;
				new_within_chrom_dprime[c] = 0;
				while (ld_count < 100)
				{//random within-chromosomes
					returned_data = Pop.AdultPopLD(c, c, randnum(Pop.NumMarkers), randnum(Pop.NumMarkers));
					if (returned_data.dprime != -5)
					{
						Pop.avg_d[c] = Pop.avg_d[c] + returned_data.d;
						new_within_chrom_dprime[c] = new_within_chrom_dprime[c] + returned_data.dprime;
						mean_dp_ldistchrom = mean_dp_ldistchrom + returned_data.dprime;
						ld_count++;
					}
				}
				Pop.avg_d[c] = Pop.avg_d[c] / 100;
				new_within_chrom_dprime[c] = new_within_chrom_dprime[c] / 100;
			}
			mean_dp_ldistchrom = mean_dp_ldistchrom / (Pop.NumChrom * 100);
			Pop.AvgLD = Pop.AvgLD / (Pop.NumChrom * Pop.NumMarkers);
			Pop.avg_pairwise_d = Pop.avg_pairwise_d / (Pop.NumChrom * Pop.NumMarkers);
			ld << reps << '\t' << g << '\t' << mean_dp_ldistchrom << '\t' << Pop.avg_pairwise_d << '\t' << Pop.AvgLD
				<< '\t' << mean_new_longdist_d << '\t' << mean_new_longdist_dprime;
			for (int c = 0; c < Pop.NumChrom; c++)
			{
				ld << '\t' << Pop.avg_d[c] << '\t' << new_within_chrom_dprime[c];
			}
			ld << '\n';
			if ((g + 1) % 100 == 0)
				cout << "Generation " << g + 1 << " complete. ";
		}//end of generations

		cout << '\n';

		if (!Pop.popExtinct)
		{
			int Q;
			for (g = 0; g < Pop.NumSelGen; g++)
			{
				Pop.PopulationSize = Pop.DeterminePopSize();
				if (Pop.PopulationSize == 0)
					Pop.popExtinct = true;
				else
				{
					cout << "Sample Gen " << g + 1 << ", ";
					Pop.StandardizeGenotypes();
					Pop.Mating(ss_sampled);//mate choice
					Pop.Mutation();
					//Sample the new population
					Pop.SamplePop();
					//calculate stats
					//Pop.SampleFsts();
					Pop.POFst();
					Pop.CalcFstChiSq(Pop.F_Marker);
					Pop.BenjaminiHochbergFDR(Pop.Psig);
					Pop.CalcLD();
					for (int c = 0; c < Pop.NumChrom; c++)
					{
						for (int cc = 0; cc < Pop.NumMarkers - 1; cc++)
						{//Pairwise
							returned_data = Pop.AdultPopLD(c, c, cc, cc + 1);
							Pop.avg_pairwise_d = Pop.avg_pairwise_d + returned_data.d;
							Pop.AvgLD = Pop.AvgLD + returned_data.dprime;
						}
					}
					Pop.AvgLD = Pop.AvgLD / (Pop.NumChrom * Pop.NumMarkers);
					//Pop.LongDistLD();
					Pop.QTLFstCalcs();
					Pop.GeneticVariance();
					Pop.WeightedFst(Pop.F_Marker, 0);
					Pop.CalculateFstCIs(Pop.F_Marker);
					Pop.bootstrap_resampling(Pop.F_Marker);
					Pop.PeakDetector();
					ci95d_m = ci95d_m + Pop.ci_detected95;
					ci95s_m = ci95s_m + Pop.ci_spurious95;
					ci99d_m = ci99d_m + Pop.ci_detected;
					ci99s_m = ci99s_m + Pop.ci_spurious;
					bs95d_m = bs95d_m + Pop.bs_detected95;
					bs95s_m = bs95s_m + Pop.bs_spurious95;
					bs99d_m = bs99d_m + Pop.bs_detected;
					bs99s_m = bs99s_m + Pop.bs_spurious;
					fdr_m = fdr_m + Pop.fdr_det;
					fdrs_m = fdrs_m + Pop.fdr_spur;
					ci95_det.push_back(Pop.ci_detected95);
					ci95_spu.push_back(Pop.ci_spurious95);
					ci99_det.push_back(Pop.ci_detected);
					ci99_spu.push_back(Pop.ci_spurious);
					bs95_det.push_back(Pop.bs_detected95);
					bs95_spu.push_back(Pop.bs_spurious95);
					bs99_det.push_back(Pop.bs_spurious);
					bs99_spu.push_back(Pop.bs_detected);
					fdr_det.push_back(Pop.fdr_det);
					fdr_spu.push_back(Pop.fdr_spur);
					sig_thresh << '\n' << reps << '\t' << g << '\t' << Pop.fdr_det << '\t' << Pop.fdr_spur << '\t' << Pop.ci_detected << '\t' << Pop.ci_spurious << '\t'
						<< Pop.ci_detected95 << '\t' << Pop.ci_spurious95 << '\t' << Pop.bs_detected << '\t' << Pop.bs_spurious << '\t' << Pop.bs_detected95 << '\t'
						<< Pop.bs_spurious95;
					peak_n++;
					if (g == 0 && reps == 0)
					{
						Fsts << "TotalNumPolymorphicLoci:" << Pop.NumPolymorphicLoci << '\n';
						Fsts << "NumSampledAdults:" << Pop.pAdultSample << '\n';
						Fsts << "NumSampledOffspring:" << Pop.ProgSample << '\n';
						Fsts << "NumSampledGenerations:" << Pop.NumSelGen << '\n' << '\n';
						Fsts << "Label\tGeneration\t" << "Chrom\t" << "Locus\t" << "QTL?\t" << "GenVar\t" << "Fst\t" << "WeightedFst\t" << "ChiSqP\t"
							<< "ExpAdultHet\t" << "ExpOffHet\t" << "Ht\t" << "Hs\t" << "A.p1\t" << "O.p1\t" << "A.LD\t" << "O.LD\t"
							<< "Fst.Q\t" << "QAHet\t" << "QOffHet\t" << "QHt\t" << "QHs\n";
						plot_Fsts << "Rep\tGen\t" << "FDRCrit\t" << "Smooth99\t" << "Smooth98\t" << "Smooth95\t" << "Smooth90\t" << "Smooth80\t"
							<< "bs99\tbs98\tbs95\tbs90\tbs80";
						for (int l = 0; l < Pop.NumChrom; l++)
						{
							for (int ll = 0; ll < Pop.NumMarkers; ll++)
							{
								if (Pop.SampledLoci[l].LociOnChrom[ll] == 1)
									plot_Fsts << '\t' << l << "." << ll;
							}
						}
						plot_Fsts << '\n';
					}
					plot_Fsts << reps << '\t' << g << "a\t" << Pop.CriticalValue << '\t' << Pop.SmoothedCritValue99 << '\t' << Pop.SmoothedCritValue98 << '\t'
						<< Pop.SmoothedCritValue95 << '\t' << Pop.SmoothedCritValue90 << '\t' << Pop.SmoothedCritValue80
						<< '\t' << Pop.bsCritValue99 << '\t' << Pop.bsCritValue98 << '\t' << Pop.bsCritValue95 << '\t' << Pop.bsCritValue90
						<< '\t' << Pop.bsCritValue80;

					for (int w = 0; w < Pop.NumChrom; w++)
					{
						for (int x = 0; x < Pop.NumSampledLoci; x++)
						{
							plot_Fsts << '\t' << Pop.F_Marker.WeightedFst[w].LociOnChrom[x];
							if (Pop.QTLtracker[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] == 1)
							{
								for (int q = 0; q < Pop.NumQTLs; q++)
								{
									if (Pop.Locations[w].LociOnChrom[q] == Pop.F_Marker.Locus[w].LociOnChrom[x])
										Q = q;
								}
								Fsts << g << "_A-O\t" << w << '\t';
								Fsts << Pop.F_Marker.Locus[w].LociOnChrom[x] << '\t';
								Fsts << "YES\t" << Pop.GenetVars[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.F_Marker.FstValue[w].LociOnChrom[x] << '\t';
								Fsts << Pop.F_Marker.WeightedFst[w].LociOnChrom[x] << '\t';
								Fsts << Pop.F_Marker.PValue[w].LociOnChrom[x] << '\t';
								Fsts << Pop.AdultHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.OffHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.TotHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.HetS[w].LociOnChrom[x] << '\t';
								Fsts << Pop.AFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.MFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.OFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.ALD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
								Fsts << Pop.OLD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
								Fsts << Pop.F_QTL.FstValue[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QAdultHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QOffHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QTotHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QHs[w].LociOnChrom[Q] << '\t';
								Fsts << '\n';
							}
							else{
								if (Pop.F_Marker.FstValue[w].LociOnChrom[x] > 0)
								{
									Fsts << g << "_A-O\t" << w << '\t';
									Fsts << Pop.F_Marker.Locus[w].LociOnChrom[x] << '\t';
									Fsts << "NO\t" << 0 << '\t';
									Fsts << Pop.F_Marker.FstValue[w].LociOnChrom[x] << '\t';
									Fsts << Pop.F_Marker.WeightedFst[w].LociOnChrom[x] << '\t';
									Fsts << Pop.F_Marker.PValue[w].LociOnChrom[x] << '\t';
									Fsts << Pop.AdultHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.OffHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.TotHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.HetS[w].LociOnChrom[x] << '\t';
									Fsts << Pop.AFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.MFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.OFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.ALD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
									Fsts << Pop.OLD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
									for (int out = 0; out < 5; out++)
										Fsts << "0\t";
									Fsts << '\n';
								}
							}
						}
					}
					plot_Fsts << '\n';

					Pop.AdultSumm();
					summ << reps << '\t' << g << "\tAO\tAdults\t" << Pop.meanLD << '\t' << Pop.GenVar << '\t' << Pop.PhenVar << '\t' << Pop.Heritability << '\t'
						<< Pop.PopulationSize << '\t' << Pop.malN << '\t' << Pop.femN << '\t' << Pop.adultRatio << '\t'
						<< Pop.allphenavg << '\t' << Pop.malphenavg << '\t' << Pop.femphenavg << '\t' << Pop.allDiff << '\t'
						<< Pop.maleDiff << '\t' << Pop.femDiff << '\t' << Pop.bs_detected << '\t' << Pop.bs_spurious << '\t'
						<< Pop.ci_detected << '\t' << Pop.ci_spurious << '\n';
					Pop.SelectionOnPhenotypes(vs_sampled);
					Pop.OffSumm();
					summ << reps << '\t' << g << "\tAO\tOff\t" << '\t' << Pop.GenVar << '\t' << Pop.PhenVar << '\t' << Pop.Heritability << '\t'
						<< Pop.PopulationSize << '\t' << Pop.malN << '\t' << Pop.femN << '\t' << Pop.offRatio << '\t'
						<< Pop.allphenavg << '\t' << Pop.malphenavg << '\t' << Pop.femphenavg << '\t' << Pop.allDiff << '\t'
						<< Pop.maleDiff << '\t' << Pop.femDiff << '\t' << '\t' << '\t' << '\t' << '\n';
					//population continues

					Pop.DensityRegulation();
					Pop.CalcMeanTraitValues();
					Pop.PopulationSize = Pop.DeterminePopSize();
					//compare current adults to previous progeny
					Pop.SamplePop();
					Pop.POFst();
					//Pop.SampleFsts();
					Pop.CalcFstChiSq(Pop.F_Marker);
					Pop.BenjaminiHochbergFDR(Pop.Psig);
					Pop.CalcLD();
					for (int c = 0; c < Pop.NumChrom; c++)
					{
						for (int cc = 0; cc < Pop.NumMarkers - 1; cc++)
						{//Pairwise
							returned_data = Pop.AdultPopLD(c, c, cc, cc + 1);
							Pop.avg_pairwise_d = Pop.avg_pairwise_d + returned_data.d;
							Pop.AvgLD = Pop.AvgLD + returned_data.dprime;
						}
					}
					Pop.AvgLD = Pop.AvgLD / (Pop.NumChrom * Pop.NumMarkers);
					//Pop.LongDistLD();
					Pop.QTLFstCalcs();
					Pop.GeneticVariance();
					Pop.WeightedFst(Pop.F_Marker, 0);
					Pop.CalculateFstCIs(Pop.F_Marker);
					Pop.bootstrap_resampling(Pop.F_Marker);
					Pop.PeakDetector();
			
					for (int w = 0; w < Pop.NumChrom; w++)
					{
						for (int x = 0; x < Pop.NumSampledLoci; x++)
						{
							if (Pop.QTLtracker[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] == 1)
							{
								for (int q = 0; q < Pop.NumQTLs; q++)
								{
									if (Pop.Locations[w].LociOnChrom[q] == Pop.F_Marker.Locus[w].LociOnChrom[x])
										Q = q;
								}

								Fsts << g << "_O-A\t" << w << '\t';
								Fsts << Pop.F_Marker.Locus[w].LociOnChrom[x] << '\t';
								Fsts << "YES\t" << Pop.GenetVars[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.F_Marker.FstValue[w].LociOnChrom[x] << '\t';
								Fsts << Pop.F_Marker.WeightedFst[w].LociOnChrom[x] << '\t';
								Fsts << Pop.F_Marker.PValue[w].LociOnChrom[x] << '\t';
								Fsts << Pop.AdultHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.OffHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.TotHet[w].LociOnChrom[x] << '\t';
								Fsts << Pop.HetS[w].LociOnChrom[x] << '\t';
								Fsts << Pop.AFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.MFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.OFreq1[w].LociOnChrom[x] << '\t';
								Fsts << Pop.ALD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
								Fsts << Pop.OLD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
								Fsts << Pop.F_QTL.FstValue[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QAdultHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QOffHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QTotHet[w].LociOnChrom[Q] << '\t';
								Fsts << Pop.QHs[w].LociOnChrom[Q] << '\t';
								Fsts << '\n';
							}
							else{
								if (Pop.F_Marker.FstValue[w].LociOnChrom[x] > 0)
								{
									Fsts << g << "_O-A\t" << w << '\t' << Pop.SampledLoci[w].LociOnChrom[x] << '\t';
									Fsts << "NO\t" << 0 << '\t';
									Fsts << Pop.F_Marker.FstValue[w].LociOnChrom[x] << '\t';
									Fsts << Pop.F_Marker.WeightedFst[w].LociOnChrom[x] << '\t';
									Fsts << Pop.F_Marker.PValue[w].LociOnChrom[x] << '\t';
									Fsts << Pop.AdultHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.OffHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.TotHet[w].LociOnChrom[x] << '\t';
									Fsts << Pop.HetS[w].LociOnChrom[x] << '\t';
									Fsts << Pop.AFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.MFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.OFreq1[w].LociOnChrom[x] << '\t';
									Fsts << Pop.ALD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
									Fsts << Pop.OLD[w].LociOnChrom[int(Pop.F_Marker.Locus[w].LociOnChrom[x])] << '\t';
									for (int out = 0; out < 5; out++)
										Fsts << "0\t";
									Fsts << '\n';
								}
							}
						}
					}
					
					Pop.AdultSumm();
					summ << reps << '\t' << g << "\tOA\tAdult\t" << Pop.meanLD << '\t' << Pop.GenVar << '\t' << Pop.PhenVar << '\t' << Pop.Heritability << '\t'
						<< Pop.PopulationSize << '\t' << Pop.malN << '\t' << Pop.femN << '\t' << Pop.adultRatio << '\t'
						<< Pop.allphenavg << '\t' << Pop.malphenavg << '\t' << Pop.femphenavg << '\t' << Pop.allDiff << '\t'
						<< Pop.maleDiff << '\t' << Pop.femDiff << '\t' << Pop.bs_detected << '\t' << Pop.bs_spurious << '\t'
						<< Pop.ci_detected << '\t' << Pop.ci_spurious << '\n';
					Pop.OffSumm();
					summ << reps << '\t' << g << "\tOA\tOff\t" << '\t' << Pop.GenVar << '\t' << Pop.PhenVar << '\t' << Pop.Heritability << '\t'
						<< Pop.PopulationSize << '\t' << Pop.malN << '\t' << Pop.femN << '\t' << Pop.offRatio << '\t'
						<< Pop.allphenavg << '\t' << Pop.malphenavg << '\t' << Pop.femphenavg << '\t' << Pop.allDiff << '\t'
						<< Pop.maleDiff << '\t' << Pop.femDiff << '\t' << '\t' << '\t' << '\t' << '\n';
				}//run five times with MC and track it each time
			}
		}
		cout << "Done rep " << reps + 1 << '\n';
		Pop.DeInitialize();
	}//end of reps

	//calculate the summary stats regarding the peaks
	peak_summ.open(peaksummname.c_str());
	peak_summ << "Name\t" << "which\t" << "AvgRatioReal_99\t" << "RealSE_99\t" << "AvgRatioSpur_99\t" << "SpuriousSE_99\t"
		<< "AvgRatioReal_95\t" << "Real_SE_95\t" << "AvgRatioSpur_95\t" << "SpuriousSE_95\t" << "TotalNumQTLs\t" << "TotalN\n";
	//calculate summary stats for sexual selection
	ci95d_m = (ci95d_m / peak_n); //peak_n is the number of reps & generations that cis are calculated for.
	ci95s_m = (ci95s_m / peak_n);
	ci99d_m = (ci99d_m / peak_n);
	ci99s_m = (ci99s_m / peak_n);
	ci95d_SE = 0;
	ci95s_SE = 0;
	ci99d_SE = 0;
	ci99s_SE = 0;
	for (int c = 0; c < peak_n; c++)
	{
		ci95d_SE = ci95d_SE + ((ci95_det[c]) - ci95d_m)*((ci95_det[c]) - ci95d_m);
		ci95s_SE = ci95s_SE + ((ci95_spu[c]) - ci95s_m)*((ci95_spu[c]) - ci95s_m);
		ci99d_SE = ci99d_SE + ((ci99_det[c]) - ci99d_m)*((ci99_det[c]) - ci99d_m);
		ci99s_SE = ci99s_SE + ((ci99_spu[c]) - ci99s_m)*((ci99_spu[c]) - ci99s_m);
	}
	ci95d_SE = sqrt(ci95d_SE / peak_n);
	ci95s_SE = sqrt(ci95s_SE / peak_n);
	ci99s_SE = sqrt(ci99s_SE / peak_n);
	ci99d_SE = sqrt(ci99d_SE / peak_n);
	peak_summ << basefilename << "\tSmooth\t" << ci99d_m << '\t' << ci99d_SE << '\t' << ci99s_m << '\t' << ci99s_SE
		<< '\t' << ci95d_m << '\t' << ci95d_SE << '\t' << ci95s_m << '\t' << ci95s_SE << '\t'
		<< total_QTL << '\t' << peak_n << '\n';
	//calculate summary stats for sexual selection
	bs95d_m = (bs95d_m / peak_n) / total_QTL; //peak_n is the number of reps & generations that cis are calculated for.
	bs95s_m = (bs95s_m / peak_n) / total_QTL;
	bs99d_m = (bs99d_m / peak_n) / total_QTL;
	bs99s_m = (bs99s_m / peak_n) / total_QTL;
	fdrs_m = (fdrs_m / peak_n) / total_QTL;
	fdr_m = (fdr_m / peak_n) / total_QTL;
	bs95d_SE = 0;
	bs95s_SE = 0;
	bs99d_SE = 0;
	bs99s_SE = 0;
	fdrs_SE = 0;
	fdr_SE = 0;
	for (int c = 0; c < peak_n; c++)
	{
		bs95d_SE = bs95d_SE + ((bs95_det[c] / total_QTL) - bs95d_m)*((bs95_det[c] / total_QTL) - bs95d_m);
		bs95s_SE = bs95s_SE + ((bs95_spu[c] / total_QTL) - bs95s_m)*((bs95_spu[c] / total_QTL) - bs95s_m);
		bs99d_SE = bs99d_SE + ((bs99_det[c] / total_QTL) - bs99d_m)*((bs99_det[c] / total_QTL) - bs99d_m);
		bs99s_SE = bs99s_SE + ((bs99_spu[c] / total_QTL) - bs99s_m)*((bs99_spu[c] / total_QTL) - bs99s_m);
		fdr_SE = fdr_SE + ((fdr_det[c] / total_QTL) - fdr_m)*((fdr_det[c] / total_QTL) - fdr_m);
		fdrs_SE = fdrs_SE + ((fdr_spu[c] / total_QTL) - fdrs_m)*((fdr_spu[c] / total_QTL) - fdrs_m);
	}
	bs95d_SE = sqrt(bs95d_SE / peak_n);
	bs95s_SE = sqrt(bs95s_SE / peak_n);
	bs99s_SE = sqrt(bs99s_SE / peak_n);
	bs99d_SE = sqrt(bs99d_SE / peak_n);
	fdr_SE = sqrt(fdr_SE / peak_n);
	fdrs_SE = sqrt(fdrs_SE / peak_n);
	peak_summ << basefilename << "\tStacksBS\t" << bs99d_m << '\t' << bs99d_SE << '\t' << bs99s_m << '\t' << bs99s_SE
		<< '\t' << bs95d_m << '\t' << bs95d_SE << '\t' << bs95s_m << '\t' << bs95s_SE << '\n'
		<< "\tFDR\t\t\t\t\t" << fdr_m << '\t' << fdr_SE << '\t' << fdrs_m << '\t' << fdrs_SE << '\n';

	peak_summ.close();
	summ.close();
	Fsts.close();
	QTLs.close();
	parameters.close();
	plot_Fsts.close();

	cout << "Finished! Input integer to quit\n";
	cin >> end;
	return 0;
}

