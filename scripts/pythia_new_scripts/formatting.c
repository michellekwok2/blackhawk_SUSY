// Script to format multiple PYTHIA hadronization tables

#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <dirent.h>

int nb_initial_particles = 14;
int nb_final_particles = 6;
int nb_initial_energies = 250;
int nb_final_energies = 500;
double Emin_init = 5.;
double Emax_init = 100000.;
double Emin_fin = 1.e-6;
double Emax_fin = 100000.;

void save_up(){
	FILE *in = fopen("save_up.txt","r");
	if(!in){
		printf("\n Could not open save_up.txt");
		fflush(stdout);
		return;
	};
	FILE *out = fopen("table_up.txt","w+");
	
	fprintf(out,"%15s%15s%15s%15s%15s%15s%15s%15s\n","init_energy","fin_energy","photon","electron","nu_e","nu_mu","nu_tau","proton");
	
	char dummy[30];
	double temp;
	for(int i = 0;i<13;i++){
		fscanf(in,"%s",dummy);
	};
	for(int i = 0;i<nb_final_energies*nb_initial_energies;i++){
		fscanf(in,"%lf",&temp);
		fprintf(out,"%15.5e",temp);
		fscanf(in,"%lf",&temp);
		fprintf(out,"%15.5e",temp);
		for(int j = 0;j<11;j++){
			switch(j){
				case 0: case 1: case 3: case 4: case 5: case 9:{
					fscanf(in,"%lf",&temp);
					fprintf(out,"%15.5e",temp);
					break;
				};
				default:{
					fscanf(in,"%s",dummy);
					break;
				};
			};
		};
		fprintf(out,"\n");
	};
	fclose(in);
	fclose(out);
	return;
};

void read_tables(char *path,double ***result){
	FILE *in = fopen(path,"r");
	if(!in){
		printf("\n Could not open file %s",path);
		return;
	};
	rewind(in);
	char dummy[20];
	for(int i = 0;i<8;i++){
		fscanf(in,"%s",dummy);
	};
	for(int i = 0;i<nb_initial_energies;i++){
		for(int j = 0;j<nb_final_energies;j++){
			fscanf(in,"%s",dummy);
			fscanf(in,"%s",dummy);
			for(int k = 0;k<nb_final_particles;k++){
				fscanf(in,"%lf",&result[i][j][k]);
			};
		};
	};
	fclose(in);
	return;
};

void write_tables(double ****tables,double *initial_energies,double *final_energies){
	char **table_names = (char **)malloc(nb_final_particles*sizeof(char *));
	for(int i = 0;i<nb_final_particles;i++){
		table_names[i] = (char *)malloc(32*sizeof(char));
	};
	sprintf(table_names[0],"photon.txt");
	sprintf(table_names[1],"electron.txt");
	sprintf(table_names[2],"nue.txt");
	sprintf(table_names[3],"numu.txt");
	sprintf(table_names[4],"nutau.txt");
	sprintf(table_names[5],"proton.txt");
	
	for(int i = 0;i<nb_final_particles;i++){
		FILE *out = fopen(table_names[i],"w");
		fprintf(out,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","initial_energy","final_energy","photon","gluons","higgs","Wpm","Z0","electron","muon","tau","up","down","charm","strange","top","bottom");
		for(int j = 0;j<nb_initial_energies;j++){
			for(int k = 0;k<nb_final_energies;k++){
				fprintf(out,"%15.5e%15.5e",initial_energies[j],final_energies[k]);
				for(int l = 0;l<nb_initial_particles;l++){
					if(l==2){ // in this case it is e+ e- -> Higgs, no need to divide by 2
						fprintf(out,"%15.5e",tables[i][j][k][l]);
					}
					else{
						fprintf(out,"%15.5e",tables[i][j][k][l]/2.);
					};
				};
				fprintf(out,"\n");
			};
		};
		fclose(out);
	};
	return;
};

int main(){
	
	//save_up();
	
	double *initial_energies = (double *)malloc(nb_initial_energies*sizeof(double));
	for(int i = 0;i<nb_initial_energies;i++){
		initial_energies[i] = pow(10.,log10(Emin_init) + (log10(Emax_init) - log10(Emin_init))/(nb_initial_energies-1)*i);
	};
	double *final_energies = (double *)malloc(nb_final_energies*sizeof(double));
	for(int i = 0;i<nb_final_energies;i++){
		final_energies[i] = pow(10.,log10(Emin_fin) + (log10(Emax_fin) - log10(Emin_fin))/(nb_final_energies-1)*i);
	};
	double ****result = (double ****)malloc(nb_initial_particles*sizeof(double ***));
	for(int i = 0;i<nb_initial_particles;i++){
		result[i] = (double ***)malloc(nb_initial_energies*sizeof(double **));
		for(int j = 0;j<nb_initial_energies;j++){
			result[i][j] = (double **)malloc(nb_final_energies*sizeof(double *));
			for(int k = 0;k<nb_final_energies;k++){
				result[i][j][k] = (double *)malloc(nb_final_particles*sizeof(double));
			};
		};
	};
	
	char **names = (char **)malloc(nb_initial_particles*sizeof(char *));
	for(int i = 0;i<nb_initial_particles;i++){
		names[i] = (char *)malloc(32*sizeof(char));
	};
	sprintf(names[0],"table_photon.txt");
	sprintf(names[1],"table_gluon.txt");
	sprintf(names[2],"table_higgs.txt");
	sprintf(names[3],"table_Wpm.txt");
	sprintf(names[4],"table_Z0.txt");
	sprintf(names[5],"table_electron.txt");
	sprintf(names[6],"table_muon.txt");
	sprintf(names[7],"table_tau.txt");
	sprintf(names[8],"table_up.txt");
	sprintf(names[9],"table_down.txt");
	sprintf(names[10],"table_charm.txt");
	sprintf(names[11],"table_strange.txt");
	sprintf(names[12],"table_top.txt");
	sprintf(names[13],"table_bottom.txt");

	// Reading the computed tables
	
	printf("[main] : Reading computed tables...\n");

	for(int i = 0;i<nb_initial_particles;i++){
		printf("\t%i",i);
		read_tables(names[i],result[i]);
		fflush(stdout);
	};

	printf("\n[main] : Reading complete\n");

	// Format conversion

	printf("[main] : Formating tables...\n");

	double ****tables = (double ****)malloc(nb_final_particles*sizeof(double ***));
	for(int i = 0;i<nb_final_particles;i++){
		tables[i] = (double ***)malloc(nb_initial_energies*sizeof(double **));
		for(int j = 0;j<nb_initial_energies;j++){
			tables[i][j] = (double **)malloc(nb_final_energies*sizeof(double *));
			for(int k = 0;k<nb_final_energies;k++){
				tables[i][j][k] = (double *)malloc(nb_initial_particles*sizeof(double));
				for(int l = 0;l<nb_initial_particles;l++){
					tables[i][j][k][l] = result[l][j][k][i]; // format conversion
				};
			};
		};
	};
	for(int l = 0;l<nb_final_particles;l++){
		for(int k = 0;k<nb_initial_particles;k++){
			for(int i = 0;i<nb_initial_energies;i++){
				for(int j = 0;j<nb_final_energies;j++){
					if(j != 0 && j != nb_final_energies-1 && tables[l][i][j][k] != 0. && ((tables[l][i][j-1][k] + tables[l][i][j+1][k])/2./tables[l][i][j][k] > 5. || (tables[l][i][j-1][k] + tables[l][i][j+1][k])/2./tables[l][i][j][k] < 0.2)){
						tables[l][i][j][k] = (tables[l][i][j-1][k] + tables[l][i][j+1][k])/2.; // avoiding spurious peaks
						j+=2;
					};
				};
			};
		};
	};
	
	printf("[main] : Formating complete\n");

	// Writing table files

	printf("[main] : Writing into files 'particles.txt'...\n");

	write_tables(tables,initial_energies,final_energies);

	printf("[main] : Writing complete\n");
	
	return 1;
};














