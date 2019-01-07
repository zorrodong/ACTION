#include "starter.h"
#include "io.h"
#include "kernel.h"
#include "action.h"
#include "rsvd.h"


#define EPS 2.2204e-16

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: runACTION\n ");
    fprintf (stream,
		 "\t-i	--input FILE_NAME\tfullpath of the input expression file (MANDATORY)\n"
		 "\t-o	--output DIR_PATH\tfullpath of the folder to store output files (default = ./results)\n"
		 "\t-t	--type MM/TAB/CSV\tType of the input file (default = TAB)\n"
		 "\t-n	--normalize 0/1\t\twhether or not to normalize the expression matrix (default = 1)\n"
		 "\t-d	--dim INT\t\tNumber of principle components to use (default = 30)\n"
		 "\t-k	--k_min INT\t\tMinimum number of archetypes [<= DIM] to try in order to find optimal k (default = 2)\n"
		 "\t-K	--k_max INT\t\tMaximum number of archetypes [<= DIM] to try in order to find optimal k (default = 15)\n"
		 "\t-m	--min_markers INT\t\tMinimum number of markers required for each cell type (default = 30)\n"
		 "\t-v	--pval_threashold INT\t\tp-value threshold to decleare differentially expressed genes (default = 0.05)\n"
		 "\t-p	--thread_no INT\t\tNumber of threads (default = -1)\n"
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
	);	 
		
    exit (exit_code);
}




int main(int argc, char ** argv)
{
		
  int next_option;
  const char *const short_options = "i:o:t:d:k:K:p:v:m:h";
  const struct option long_options[] = {
		{"help",     0, NULL, 'h'},
		{"input",  	 1, NULL, 'i'},
		{"output", 	 1, NULL, 'o'},
		{"type", 	 1, NULL, 't'},
		{"dim", 	 1, NULL, 'd'},
		{"k_min", 	 1, NULL, 'k'},
		{"k_max", 	 1, NULL, 'K'},
		{"pval_threshold", 	 1, NULL, 'v'},
		{"min_markers", 	 1, NULL, 'm'},
		{"thread_no", 	 1, NULL, 'p'},
		{NULL,       0, NULL,  0 }		
	};

	if(argc == 1) {
		print_usage (stdout, -1);
	}
	
	
	char input_file[1024] = "", output_path[1024] = "./results", type[1024] = "TAB";
	int PCA_dim = 30;
	int k_min = 2; int k_max = 15;
	int thread_no = -1;
	int normalize = 1;
	double pval_threshold = 0.05;
	int min_markers = 30;
	
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_file, optarg);
	    		break;

			case 'o':
				strcpy(output_path, optarg);
	    		break;
	    		
			case 't':
				strcpy(type, optarg);
	    		break;

			case 'n':
				normalize = atoi(optarg);
				break;
				
	    		
			case 'p':
				thread_no = atoi(optarg);
	    		break;			
	    		
			case 'd':
				PCA_dim = atoi(optarg);
	    		break;			
	
			case 'k':
				k_min = atoi(optarg);
	    		break;			
	    		
			case 'K':
				k_max = atoi(optarg);
	    		break;						    		

			case 'v':
				pval_threshold = atof(optarg);
	    		break;						    		

			case 'm':
				min_markers = atoi(optarg);
	    		break;						    		
		}
    } while (next_option != -1);

	if(!strcmp(input_file, "")) {
		fprintf(stderr, "full path to the expression file is mandatory.\n");
		print_usage(stderr, -1);
		return -1;
	}



	sp_mat expression;
	if(!strcmp(type, "TAB")) {
		expression = read_from_table(input_file, '\t');	
	}
	else if(!strcmp(type, "CSV")) {
		expression = read_from_table(input_file, ',');	
	}
	else if(!strcmp(type, "MM")) {
		expression = read_from_mm(input_file);	
	}
	else {
		fprintf(stderr, "Unsupported file type %s\n", type);
	}
	
	
	if(expression.n_elem == 0) {
		fprintf(stderr, "Error reading file %s\n", input_file);
		return -1;
	}
	else {
		printf("Expression file imported successfully: %d genes and %d cells\n", expression.n_rows, expression.n_cols);
	}


	// values are not log-normalized
	if(expression.max() > 100) {
		printf("log-normalizing expression matrix\n");

		sp_mat::iterator it     = expression.begin();
		sp_mat::iterator it_end = expression.end();	
		for(; it != it_end; ++it) {		
			(*it) = std::log2((*it) + 1);
		}		
	}
	
	/*
	
	Projection projection = reduceGeneExpression(expression, PCA_dim, 1, 10);
	mat cell_signatures = projection.cell_signatures;
	
	ACTION_results trace = runACTION(cell_signatures, k_min, k_max, thread_no);
	vector<mat> Archs = computeMetaCells(expression, trace.C);	
	uvec cell_types = discretizeCellProfiles(expression, trace.H, pval_threshold, min_markers, 3, thread_no);
	


	// Export Results
	char fname[1024];
	FILE *fd;
		
	
	// Exporting cell types
	sprintf(fname, "%s/Celltypes.txt", output_path);		
	fd = fopen(fname, "w");			
	for(int i = 0; i < cell_types.n_elem; i++) {
		fprintf(fd, "%d\n", cell_types(i)+1);
	}		
	fclose(fd);		
	

	// Export trace of archs (Z*C)
	for (int k = k_min; k <= k_max; k++) {
		mat arch = Archs[k];

		sprintf(fname, "%s/MetaCells_k=%d.txt", output_path, k);
		arch.save(fname, csv_ascii);

	}
	

	// Export trace of cell profiles (H)
	for (int k = k_min; k <= k_max; k++) {
		mat H = trace.H[k];
		
		sprintf(fname, "%s/CellProfiles_k=%d.txt", output_path, k);		
		H.save(fname, csv_ascii);
	}	
	*/
	
	return 0;
}
