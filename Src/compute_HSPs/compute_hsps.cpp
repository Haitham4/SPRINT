//============================================================================
// Name        : compute_hsps.cpp
// Author      : Yiwei Li
// Date		   : Date: May 2016
//============================================================================
#include "global_parameters.h"
#include "hash_table.h"
#include "PtoHSP.h"

int main(int argc, char * argv[]) {
	cout<<"-------------------------------------------------------------------\n";
	string error_msg = "In order to run SPRINT-compute HSPs, type compute_HSPs and the following options: \n -p <protein_file> (required)\n -h <hsp_output_file_name> (required)\n -Thit <an integer, the threshold Thit > (optional, default 15) \n -Tsim <an integer, the threshold Tsim > (optional, default 35) \n ";
	cout<<error_msg;
	cout<<"-------------------------------------------------------------------\n";
	for(int a = 0; a < argc; a ++){
		if(!strcmp(argv[a], "-p")){
			PROTEIN_FN = argv[a+1];
		}
		if(!strcmp(argv[a], "-h")){
			HSP_FN = argv[a+1];
		}
		if(!strcmp(argv[a], "-Thit")){
			Thit = atoi(argv[a+1]);
		}
		if(!strcmp(argv[a], "-Tsim")){
			T_kmer = atoi(argv[a+1]);
		}
	}
	cout<<"PROTEIN_FN: "<<PROTEIN_FN<<endl;
	cout<<"HSP_FN: "<<HSP_FN<<endl;
	cout<<"Thit: "<<Thit<<endl;
	cout<<"Tsim: "<<T_kmer<<endl;


	cout<<"-----------compute_hsps starts--------"<<endl;
	load_protein(PROTEIN_FN);
	cout<<"Number of Proteins: "<<num_protein<<endl;
	load_BLOSUM_convert(BLOSUM_convert);

	cout<<"-------initialize the hash tables-----"<<endl;
	HASH_TABLE ht0 = HASH_TABLE();
	HASH_TABLE ht1 = HASH_TABLE();
	HASH_TABLE ht2 = HASH_TABLE();
	HASH_TABLE ht3 = HASH_TABLE();
	#ifdef PARAL
	#pragma omp parallel 
	#endif
	{	
#ifdef PARAL	
#pragma omp sections
#endif
		{
#ifdef PARAL
#pragma omp section
#endif			
			{	
				ht0.creat_hash_table(seed_orig[0]);
			}
#ifdef PARAL
#pragma omp section
#endif			
			{
				ht1.creat_hash_table(seed_orig[1]);
			}
#ifdef PARAL			
#pragma omp section
#endif			
			{
				ht2.creat_hash_table(seed_orig[2]);
			}
#ifdef PARAL			
#pragma omp section
#endif			
			{
				ht3.creat_hash_table(seed_orig[3]);
			}
		}
	}

	cout<<"-----------build the HSP table--------"<<endl;
	//build the HSP table
	PtoHSP hsp = PtoHSP();

	clock_t t = clock();
	cout<<"-----------building the first hashtalbe (total: 4)"<<endl;
	hsp.load_hash_table(ht0);
	t = clock() - t;
	printf ("HSP table 0 took me (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);			
	
	t = clock();
	cout<<"-----------building the second hashtalbe (total: 4)"<<endl;
	hsp.load_hash_table(ht1);
	t = clock() - t;
	printf ("HSP table 1 took me (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);			
	
	t = clock();
	cout<<"-----------building the third hashtalbe (total: 4)"<<endl;
	hsp.load_hash_table(ht2);
	t = clock() - t;
	printf ("HSP table 2 took me (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);			

	t = clock();
	cout<<"-----------building the fourth hashtalbe (total: 4)"<<endl;
	hsp.load_hash_table(ht3);
	t = clock() - t;
	printf ("HSP table 3 took me (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);			

	hsp.print_HSP();
	cout<<"-----------compute_hsps finished--------"<<endl;

	return 0;
}
