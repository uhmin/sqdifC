#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "filename.h"

#define MAXLETTER 256
#define ALL_CODON 64
#define ALPHABET 26
#define CODON_LENGTH 3
#define SEQ_NUM 10
#define true 1
#define false 0

int CONSIDER_ROUTE_IN_COUNTING_SYNONIMOUS_SITES=1;
/*
  char CODONTABLE[MAXLETTER]="/usr/bio/src/sqdifC/universal.codon";
  char SUBSTITUTION_TABLE[MAXLETTER]="/usr/bio/src/sqdifC/MiyataYasunagaTable_real.csv";
*/
char BASES[]="atgc";

typedef struct SWITCHES {
  int weighted;  /* default=1 */
  int debug;     /* default=0 */
  int ELIMINATE_TERMINAL_CODON; /* default=0 */
  int CONSIDER_ROUTE_IN_COUNTING_SYNONIMOUS_SITES; /* default=1 */
  char *infile;
  char *outfile;
  char *codonTable; 
  char *substitutionTable;
  char format; /* defaulf='T' */
  int width; /* default=0 */
  int slide; /* default=0 */
  char allCombination; /* default='T' */
} SWITCHES_;

typedef struct SEQUENCE {
  char *name;
  char *sequence;
  int length;
} SEQUENCE_;

typedef struct SEQ_DATA {
  struct SEQUENCE *sequence;
  int number;
} SEQ_DATA_;

typedef struct SUBSTITUTION_TABLE{
  double **synonimous, **nonsynonimous, **synonimousSiteNum;
} sUBtABLE;

typedef struct ROUTE{
  char *synonimous, *nonsynonimous;
  char route[ALL_CODON][MAXLETTER];
  double *weight;
  double *synonimousSiteNum;
  int *depth;
  int number;
} rOUTE;

typedef struct BASIC_TABLES{
  double **substitutionMatrix;
  char *codonTable;
  double *tableS;
  sUBtABLE substitutionTable;
} BASIC_TABLES;

SWITCHES_  read_and_check_options(int argc, char *argv[]);
SWITCHES_ read_options(int argc, char *argv[]);
SWITCHES_ loadDefaultSwitch(void);
SWITCHES_ check_width_and_slide(SWITCHES_ switches);
BASIC_TABLES  create_tables(SWITCHES_ switches);
double **load_substitution_matrix(SWITCHES_ switches);
double **make_null_substitution_matrix(void);
char **split(char *str, char *substr);
int *input_alphabet_order(char **data);
double **convert_to_substitution_matrix(char ***data, SWITCHES_ switches);

char *load_codon_table(SWITCHES_ switches);
int encode_codon(char *line);
char *decode_codon(int ID);
double *create_synonimous_table(BASIC_TABLES tables, SWITCHES_ switches);
sUBtABLE create_pathway_table(BASIC_TABLES tables, SWITCHES_ switches);
double *average_synonimous_change(struct ROUTE route, SWITCHES_ switches);
struct ROUTE trace_route(struct ROUTE route, char synonims_org, char nonsynonims_org, 
  		 char *Codon1, char *Codon2, BASIC_TABLES tables,
			 SWITCHES_ switches);
struct SEQ_DATA read_and_check_multifasta_sequence(SWITCHES_ switches);
struct SEQ_DATA read_input_multifasta_sequence(SWITCHES_ switches);
struct SEQUENCE read_one_sequence(FILE *fin);
int check_sequences(struct SEQ_DATA Sequence, SWITCHES_ switches);
int  print_header(SWITCHES_ switches);
int cue(FILE *fin);
int compare_sequences(struct SEQUENCE sequence1, struct SEQUENCE sequence2,
		      int from, int to, BASIC_TABLES tables, SWITCHES_ switches);
char **add_alignment_data_1(char **alignments, int q, char baseA, char baseB);
char **add_alignment_data_2(char **alignments, int q, char aminoA, char aminoB);
int print_text_result(struct SEQUENCE sequence1, struct SEQUENCE sequence2, 
		      int from, int to, char *alignment_buffer, int nuc_change,
		      int effectiveLength, int amino_change, int numOfSites,
		      double S, double Ts, double Ta);
int print_header_of_table_result(void);
int print_table_result(struct SEQUENCE sequence1, struct SEQUENCE sequence2,
		       int from, int to, char *alignment_buffer, int nuc_change,
		       int effectiveLength, int amino_change, int numOfSites,
		       double S, double Ts, double Ta);
char *add_alignment(char *alignment_buffer, int *buffer_length, 
		    char **line);
void err_in_codon_table(void);
int free_safe(char **P);
char *chomp(char *sequence);
void help(void);
void *malloc_(int line, size_t size);

int main(int argc, char *argv[]){
  BASIC_TABLES tables;
  struct SEQ_DATA Sequence;
  SWITCHES_ switches;
  int i, j;
  int from, to;
  int iteration;
  char SlideWindow=1;

  switches=read_and_check_options(argc, argv);
  tables=create_tables(switches);
  Sequence=read_and_check_multifasta_sequence(switches);
  print_header(switches);

  if(switches.width==0 && switches.slide==0){
    switches.width=Sequence.sequence[0].length;
    switches.slide=1;
    SlideWindow=0;
  }
  if(switches.allCombination=='F'){
    iteration=1;
  }else{
    iteration=Sequence.number;
  }

  for(i=0; i<iteration; i++){
    for(j=i+1; j<Sequence.number; j++){
      for(from=0; 
	  from+switches.width-1<Sequence.sequence[0].length; 
	  from+=switches.slide){
	to=from+switches.width;
	compare_sequences(Sequence.sequence[i], Sequence.sequence[j]
			  , from, to, tables , switches);
      }
      if(SlideWindow){
	/* show summary */
	printf("Summary\n");
	compare_sequences(Sequence.sequence[i], Sequence.sequence[j], 0, 
			  Sequence.sequence[0].length, tables, switches);
      }
      printf("\n");
    }
  }
  return 0;
}

BASIC_TABLES  create_tables(SWITCHES_ switches){
  BASIC_TABLES tables;
  /* read substitution matrix */
  if(switches.weighted){
    tables.substitutionMatrix=load_substitution_matrix(switches);
  }else{
    tables.substitutionMatrix=make_null_substitution_matrix();
  }
  /* read codon table */
  tables.codonTable=load_codon_table(switches);
  /* create table of number of synonimouse sites for each codon */
  tables.tableS=create_synonimous_table(tables, switches);
  /* create table of number of synonimous sites for each 64x64 pathways */
  tables.substitutionTable=create_pathway_table(tables, switches);
  return tables;
}

int  print_header(SWITCHES_ switches){
  /* compare sequence for all possible combination */
  if(switches.debug){
    printf("======  RESULT  ======\n");
  }
  if(switches.format=='X'){
    print_header_of_table_result();
  }
  return 0;
}

int compare_sequences(struct SEQUENCE sequence1, struct SEQUENCE sequence2, 
		      int from, int to, BASIC_TABLES tables, SWITCHES_ switches){
  double S, Ts, Ta;
  char baseA, baseB;
  int IDa, IDb;
  int p, q;
  int i;
  int numOfSites, effectiveLength;
  int nuc_change, amino_change;
  /*  char alignments[4][MAXLETTER];*/
  char **alignments;
  int width=60;
  char *alignment_buffer;
  int buffer_length=MAXLETTER*3;

  /* initialize alignment buffers */
  alignment_buffer=(char *)malloc_(__LINE__, sizeof(char)*buffer_length); alignment_buffer[0]='\0';
  alignments=(char **)malloc_(__LINE__, sizeof(char*)*4);
  for(i=0; i<4; i++) alignments[i]=(char *)malloc_(__LINE__, sizeof(char)*(width+10));

  /*
  for(p=0, q=0, S=0, Ts=0, Ta=0, amino_change=0, nuc_change=0,
	numOfSites=0, effectiveLength=0; p<sequence1.length; p++, q++){
  */
  for(p=from, q=0, S=0, Ts=0, Ta=0, amino_change=0, nuc_change=0,
	numOfSites=0, effectiveLength=0; p<to; p++, q++){
    /* get p-th base of i-th sequence */
    baseA=sequence1.sequence[p];
    /* get p-th base of j-th sequence */
    baseB=sequence2.sequence[p];

    add_alignment_data_1(alignments, q, baseA, baseB);

    if(isalpha(baseA) && isalpha(baseB)){
      effectiveLength++;
      if(baseA!=baseB){
	nuc_change++;
      }
    }
    if(p % 3 == 0){
      /* get amino acid of the codon of p-th to (p+2)-th base */
      IDa=encode_codon(sequence1.sequence+p);
      IDb=encode_codon(sequence2.sequence+p);
      if(IDa==-1 || IDb==-1 || tables.codonTable[IDa]=='\0' || tables.codonTable[IDb]=='\0'){
	/* terminal codon or invalid codon */
	continue;
      }
      add_alignment_data_2(alignments, q, 
			   tables.codonTable[IDa], tables.codonTable[IDb]);
      /* increment the number of synonimous sites */
      if(CONSIDER_ROUTE_IN_COUNTING_SYNONIMOUS_SITES){
	S+=tables.substitutionTable.synonimousSiteNum[IDa][IDb];
      }else{
	S+=(tables.tableS[IDa]+tables.tableS[IDb])/2; /* non weighted method */
      }
      /*      printf("# %s: %.3f + %s: %.3f,  %.3f\n", decode_codon(IDa), tables.tableS[IDa], decode_codon(IDb), tables.tableS[IDb], tables.substitutionTable.synonimousSiteNum[IDa][IDb]);*/
      /* increment the number of    synonimous  changes */
      Ts+=tables.substitutionTable.synonimous[IDa][IDb];
      /* increment the number of nonsynonimous  changes */
      Ta+=tables.substitutionTable.nonsynonimous[IDa][IDb];
      /* increment the amino acid changes */
      if(tables.codonTable[IDa]!=tables.codonTable[IDb]){
	amino_change++;
      }
      /* increment the number of amino acid sites */
      numOfSites+=3;
    }
    if(q==width-1){
      /* add alignment into buffer */
      alignment_buffer=add_alignment
	(alignment_buffer, &buffer_length, alignments);

      q=-1;
    }
  }

  /* add alignment into buffer (last)*/
  alignment_buffer=add_alignment
    (alignment_buffer, &buffer_length, alignments);

  if(switches.format=='T'){ /* output in text format */
    print_text_result(sequence1, sequence2, from, to, alignment_buffer, nuc_change,
		      effectiveLength, amino_change, numOfSites, S, Ts, Ta);
  }else{ /* output in tabular format */
    print_table_result(sequence1, sequence2, from, to, alignment_buffer, nuc_change,
		       effectiveLength, amino_change, numOfSites, S, Ts, Ta);
  }

  /* free alignment buffers */
  free(alignment_buffer);
  for(i=0; i<4; i++) free(alignments[i]);
  free(alignments);
  
  return 0;
}
  
char **add_alignment_data_1(char **alignments, int q, char baseA, char baseB){
  /* create alignment data */
  alignments[0][q]=baseA; alignments[0][q+1]='\0';
  if(baseA==baseB){
    alignments[1][q]='.';
  }else{
    alignments[1][q]=baseB;
  }
  alignments[1][q+1]='\0';
  alignments[2][q]=' '; alignments[2][q+1]='\0';
  alignments[3][q]=' '; alignments[3][q+1]='\0';
  return alignments;
}

char **add_alignment_data_2(char **alignments, int q, char aminoA, char aminoB){
  /* write the alignment of translated sequence */
  alignments[2][q]=aminoA; alignments[2][q+1]='\0';
  if(aminoA!=aminoB){
	alignments[3][q]=aminoB; alignments[3][q+1]='\0';
  }else{
    alignments[3][q]='.';
  }
  return alignments;
}

char *add_alignment(char *alignment_buffer, int *buffer_length, 
		    char **line){
  char *tmpAlignment;
  int tmplength;
  int i;
  for(tmplength=0, i=0; i<4; i++){
    tmplength+=strlen(line[i])+3;
  }
  tmpAlignment=(char *)malloc_(__LINE__, sizeof(char)*tmplength);
  sprintf(tmpAlignment,"%s\n%s\n %s\n %s\n\n", line[0], line[1], line[2], line[3]);
  if(*buffer_length<tmplength+strlen(alignment_buffer)){
    *buffer_length+=MAXLETTER*3;
    if((alignment_buffer=(char *)realloc(alignment_buffer, *buffer_length*sizeof(char)))==NULL){
      fprintf(stderr, "Memory err in line %d\n", __LINE__);
      exit(0);
    }
  }
  strcat(alignment_buffer, tmpAlignment);
  free(tmpAlignment);
  return alignment_buffer;
}

int print_header_of_table_result(void){
  printf("seqname 1\t");    /* sequence name 1 */
  printf("seqname 2\t");    /* sequence name 2 */
  printf("sequence length\t"); /* sequence length */
  printf("effective length\t"); /* effective length */
  printf("nucleotide changes\t");  /* nucleotide changes */
  printf("effective amino acid length\t"); /* effective amino acid length */
  printf("amino acid changes\t"); /* num of amino acid changes */
  printf("effective sites\t");   /* effective sites */
  printf("synonimous sites\t"); /* num of synonimous sites */
  printf("synonimous changes\t"); /* num of synonimous changes */
  printf("Ks\t"); /* Ks */
  printf("non-synonimous sites\t"); /* num of non-synonimous sites */
  printf("non-synonimous changes\t"); /* num of non-synonimous changes */
  printf("Ka\t"); /* Ka */
  printf("Ka/Ks\t"); /* Ka/Ks */
  printf("Jukes and Cantor's correction of Ks (CKs)\t"); /* Jukes and Cantor's correction of Ks */
  printf("CKa\t"); /* Jukes and Cantor's correction of Ka */
  printf("CKa/CKs\t"); /*Jukes and Cantor's correction of Ka/Ks */
  printf("Region");
  printf("\n");
  return 0;
}

int print_table_result(struct SEQUENCE sequence1, struct SEQUENCE sequence2
		       , int from, int to, char *alignment_buffer, int nuc_change,
		      int effectiveLength, int amino_change, int numOfSites,
		       double S, double Ts, double Ta){
  double Ksp, Kap; /* Jukes and Cantor's correction */
  double Ks, Ka;
  int seqlength;
  char *seqname;
  char *freePosition;
  int detail;
  int i;
  if(from % 3 ==0 && to % 3 == 0){
    detail=1;
  }else{
    detail=0;
  }
  seqlength=strlen(sequence1.name)+strlen(sequence2.name);
  seqname=(char *)malloc_(__LINE__, seqlength*sizeof(char));
  freePosition=seqname;
  strcpy(seqname, sequence1.name);
  strtok(seqname, " ");
  printf("%s\t",seqname);    /* sequence name 1 */
  strcpy(seqname, sequence2.name);
  strtok(seqname, " "); 
  printf("%s\t",seqname);    /* sequence name 2 */
  free_safe(&freePosition);

  printf("%d\t", sequence1.length); /* sequence length */
  printf("%d\t",effectiveLength); /* effective length */
  printf("%d\t",nuc_change);  /* nucleotide changes */
  if(detail){
    printf("%d\t",numOfSites/3); /* effective amino acid length */
    printf("%d\t",amino_change); /* num of amino acid changes */
    printf("%d\t",numOfSites);   /* effective length */
    printf("%f\t",S); /* num of synonimous sites */
    printf("%f\t",Ts); /* num of synonimous changes */
    Ks=Ts/S;
    printf("%f\t",Ks); /* Ks */
    printf("%f\t",numOfSites-S); /* num of non-synonimous sites */
    printf("%f\t",Ta); /* num of non-synonimous changes */
    Ka=Ta/(numOfSites-S);
    printf("%f\t",Ka); /* Ka */
    printf("%f\t",Ka/Ks); /* Ka/Ks */
    Ksp=-log(1-4*Ks/3)*3/4;
    printf("%f\t",Ksp); /* Jukes and Cantor's correction of Ks */
    Kap=-log(1-4*Ka/3)*3/4;
    printf("%f\t",Kap); /* Jukes and Cantor's correction of Ka */
    printf("%f\t",Kap/Ksp); /*Jukes and Cantor's correction of Ka/Ks */
  }else{
    for(i=0; i<13; i++){
      printf("-\t");
    }
  }
  printf("%d-%d", from+1, to);
  printf("\n");
  return 0;
}

int print_text_result(struct SEQUENCE sequence1, struct SEQUENCE sequence2,
		      int from, int to, char *alignment_buffer, int nuc_change,
		      int effectiveLength, int amino_change, int numOfSites,
		      double S, double Ts, double Ta){
  double Ksp, Kap; /* Jukes and Cantor's correction */
  int detail;
  if(from % 3 == 0 && to % 3 == 0){
    detail=1;
  }else{
    detail=0;
  }
  /* print sequence information */  
  printf("Input sequences:\n");
  printf("%s length= %d\n"
	 , sequence1.name, sequence1.length);
  printf("%s length= %d\n"
	 , sequence2.name, sequence2.length);
  printf("Region: %d-%d\n", from+1, to);
  if(detail){
    /* print alignment */
    printf("%s", alignment_buffer);
  }
  /* print result */
  printf("number of nucleotide changes = %d\n", nuc_change);
  printf("  effective length = %d\n", effectiveLength);
  printf(" substitution rate = %f\n",(double)nuc_change/effectiveLength);
  if(detail){
    printf("number of amino acid changes = %d\n", amino_change);
    printf("  effective length = %d\n", numOfSites/3);
    printf(" substitution rate = %f\n", (double)amino_change/numOfSites*3);
    printf("number of sites used = %d\n", numOfSites);
  }
  if(detail){
    printf("Synonimous   sites = %f\n", S);
    printf("Synonimous changes = %f\n", Ts);
    printf("Ks = %f\n", Ts/S);
    printf("non Synonimous   sites = %f\n", numOfSites-S);
    printf("non Synonimous changes = %f\n", Ta);
    printf("Ka = %f\n", Ta/(numOfSites-S));
    printf("Ka/Ks = %f\n", (Ta/(numOfSites-S))/(Ts/S));
    
    Ksp=-log(1-4*Ts/S/3)*3/4;
    Kap=-log(1-4*Ta/(numOfSites-S)/3)*3/4;
    printf("Jukes and Cantor's correction\n");
    printf("Ksp= %f\n", Ksp);
    printf("Kap= %f\n", Kap);
    printf("Kap/Ksp = %f\n", Kap/Ksp);
  }
  printf("\n");  
  return 0;
}

SWITCHES_  read_and_check_options(int argc, char *argv[]){
  SWITCHES_ switches;
  switches=read_options(argc, argv);
  switches=check_width_and_slide(switches);
  return switches;
}

SWITCHES_ read_options(int argc, char *argv[]){
  int count,infile=-1,outfile=-1, codon=-1;
  int substitution=-1, terminalCodon=-1, debug=-1;
  int outFormat=-1, synonimousRoute=-1;
  int weighted=-1, windowS=-1, slideS=-1;
  int combination=-1;
  SWITCHES_ switches;
  switches=loadDefaultSwitch();

  /* check the number of options */
  if ((argc % 2)==0) help();

  for(count=1;count<argc-1;count+=2){
    switch(argv[count][1]){
    case 'i': /* infile name */
      infile=count+1;
      break;
    case 'o': /* outfile name */
      outfile=count+1;
      break;
    case 'f': /* output format */
      outFormat=count+1;
      break;
    case 'r': /* select if consider the route in counting synonimous sites */
      synonimousRoute=count+1;
      break;
    case 'W': /* select if consider weited estimation */
      weighted=count+1;
      break;
    case 'w': /* window size */
      windowS=count+1;
      break;
    case 's': /* slide size */
      slideS=count+1;
      break;
    case 'c': /* slide size */
      combination=count+1;
      break;
    case 'C': /* name of codon table */
      codon=count+1;
      break;
    case 'S': /* name of substitution table */
      substitution=count+1;
      break;
    case 'T': /* eliminate terminal codon */
      terminalCodon=count+1;
      break;
    case 'D':
      debug=count+1;
      break;
    default:
      fprintf(stderr, "Unidentified option: %s", argv[count]);
      help();
      break;
    }
  }
  if(infile!=-1){
    switches.infile=argv[infile];
  }
  if(outfile!=-1){
    switches.outfile=argv[outfile];
  }
  if(codon!=-1){
    switches.codonTable=argv[codon];
  }
  if(combination!=-1 && argv[combination][0]=='F'){
    switches.allCombination='F';
  }
  if(substitution!=-1){
    switches.substitutionTable=argv[substitution];
  }
  if(terminalCodon!=-1 && argv[terminalCodon][0]=='T'){
    switches.ELIMINATE_TERMINAL_CODON=1;
  }
  if(outFormat!=-1 && argv[outFormat][0]=='X'){
    switches.format='X'; /* out format is tabular */
  }
  if(synonimousRoute!=-1 && argv[synonimousRoute][0]=='F'){
    switches.CONSIDER_ROUTE_IN_COUNTING_SYNONIMOUS_SITES=0;
  }
  if(weighted!=-1 && argv[weighted][0]=='F'){
    switches.weighted=0;
  }
  if(windowS!=-1){
    switches.width=atoi(argv[windowS]);
  }
  if(slideS!=-1){
    switches.slide=atoi(argv[slideS]);
  }
  if(debug!=-1 && argv[debug][0]=='T'){
    switches.debug=true;
  }
  return switches;
}

SWITCHES_ loadDefaultSwitch(void){
  struct SWITCHES switches;
  switches.infile=NULL;
  switches.outfile=NULL;
  switches.weighted=true;
  switches.debug=false;
  switches.ELIMINATE_TERMINAL_CODON=false;
  switches.CONSIDER_ROUTE_IN_COUNTING_SYNONIMOUS_SITES=true;
  switches.codonTable=CODONTABLE;
  switches.substitutionTable=SUBSTITUTION_TABLE;
  switches.format='T';
  switches.width=0;
  switches.slide=0;
  switches.allCombination='T';
  return switches;
}

SWITCHES_ check_width_and_slide(SWITCHES_ switches){
  if(switches.width==0 && switches.slide!=0){
    switches.width=51;
  }
  if(switches.width!=0 && switches.slide==0){
    switches.slide=3;
  }
  /*
  if(switches.width % 3 != 0 || switches.slide % 3 !=0){
    fprintf(stderr, "window size and slide size must be dividable by 3\n");
  }
  */
  if(switches.debug){
    printf("window=%d, slide=%d\n", switches.width, switches.slide);
  }
  return switches;
}

double **load_substitution_matrix(SWITCHES_ switches){
  FILE *fin;
  double **substitutionMatrix;
  char line[MAXLETTER];
  int i, skip=0;
  char ***data;

  data=(char ***)malloc_(__LINE__, sizeof(char **)*MAXLETTER);

  /* file open */
  if ((fin=fopen(switches.substitutionTable,"r"))==NULL){
    fprintf(stderr, "Could not open sbstitution table.[%s]\n"
	    , switches.substitutionTable);
    fprintf(stderr, "Specify it using -S option\n");
    exit(0);
  }
  skip=0;
  /* Reads the data of substitution matrix 
     and inputs to the data[][] */
  for(i=1; fgets(line, MAXLETTER, fin)!=NULL;){
    if(line[0]=='#' || skip==1) { /* skip if the first letter is # */
      /* printf("skip: %s\n", line); */
      skip=1;
      if(line[strlen(line)-1]=='\n') skip=0; /* in case if one line is very long */
    }else if(line[0]==' ' && strlen(line)>10){ /* read the order of amino acid data */
      data[0]=split(line, " "); /* data[0] contains the order of amino acid */
      i=1;
    }else{
      data[i]=split(line," ");
      data[i+1]=NULL;
      i++;
    }
  }
  fclose(fin);  
  /* sort the raw substitution matrix to alphabetical order */
  substitutionMatrix=convert_to_substitution_matrix(data, switches);
  return substitutionMatrix;
}

/* This function if like a split() function of Perl language.
   It splits the string "str" by delimiter "substr".  */
char **split(char *str, char *substr){
  char **result;
  char *p;
  int i, num_of_token;

  num_of_token=ALL_CODON;
  result=(char **)malloc_(__LINE__, sizeof(char *)*num_of_token);
  result[0]=NULL;
  p=strtok(str, substr);
  for(i=0; p!=NULL; i++){
    if(i+2>num_of_token){
      num_of_token+=ALL_CODON;
      if((result=(char **)realloc(result, sizeof(char *)*num_of_token))==NULL){
	fprintf(stderr,"Memory err on line %d\n", __LINE__);
	exit(0);
      }
    }
    if(p==NULL) break;
    result[i]=(char *)malloc_(__LINE__, (strlen(p)+1)*sizeof(char));
    strcpy(result[i], p);
    /*
    printf("%d(%s)\n",i, result[i]);
    */
    result[i+1]=NULL;
    p=strtok(NULL, substr);
  }
  return result;
}

/* sort the raw substitution matrix to alphabetical order */
double **convert_to_substitution_matrix(char ***data, SWITCHES_ switches){
  double **result;
  double value;
  int i, j;
  int s, t;

  if(switches.debug){
    printf("======  Substitution Matrix ======\n");
  }
  result=(double **)malloc_(__LINE__, sizeof(double *)*ALPHABET);
  for(i=1; data[i]!=NULL; i++){
    s=tolower(data[i][0][0])-'a';
    result[s]=(double *)malloc_(__LINE__, sizeof(double)*ALPHABET);
    for(j=1; data[i][j]!=NULL && isalpha(data[0][j-1][0]);j++){
      t=tolower(data[0][j-1][0])-'a';
      value=atof(data[i][j]);
      /*
      value=1/atof(data[i][j]);
      */
      if(switches.debug){
	printf("%c(%d) -> %c(%d): %.2f(%s)\n", (s+'a'), s, (t+'a'), t, value, data[i][j]);
      }
      result[s][t]=value;
    }
  }
  return result;
}

double **make_null_substitution_matrix(void){
  double **result;
  int i, j;
  result=(double **)malloc_(__LINE__, sizeof(double *)*ALPHABET);
  for(i=0; i<ALPHABET; i++){
    result[i]=(double *)malloc_(__LINE__, sizeof(double)*ALPHABET);
    for(j=0; j<ALPHABET; j++){
      result[i][j]=1;
    }
  }
  return result;
}

char *load_codon_table(SWITCHES_ switches){
  FILE *fin;
  char *codonTable;
  char line[MAXLETTER];
  int ID, i;

  codonTable=(char *)malloc_(__LINE__, ALL_CODON*sizeof(char));
  /* file open */
  if ((fin=fopen(switches.codonTable,"r"))==NULL){
    fprintf(stderr,"Could not open codon table.[%s]\n"
	    , switches.codonTable);
    fprintf(stderr, "Specify it using -C option\n");
    exit(0);
  }
  if(switches.debug){
    printf("======= Loading Codon Table from file %s =======\n", switches.codonTable);
  }
  while(fgets(line, MAXLETTER, fin)!=NULL){
    /* read first 3 letters and convert this codon to unique number */
    ID=encode_codon(line);
    if(switches.debug){
      printf("codon ID = %d, decoded to %s from %s", ID, decode_codon(ID), line);
    }
    codonTable[ID]='\0'; /* Default Amino Acid will be terminal */
    for(i=CODON_LENGTH; line[i]!='\0'; i++){ /* after the codon, search the Amino Acid symbol */
      if(isalpha(line[i]) || line[i]=='-' || line[i]=='*'){
	/* find alphabet or "-" or "*" */
	if(line[i]=='-' || line[i]=='*'){ /* these symbol is terminal codon */
	  line[i]='\0';
	}
	codonTable[ID]=line[i];
	break;
      }
    }
  }
  fclose(fin);
  return(codonTable);
}

/* This function reads the first 3 letters 
   and converts the codon to the unique number */
int encode_codon(char *line){
  int i, j, ID;
  int err=0;
  int legal_letter;
  int codon_number;
  char base;
  codon_number=strlen(BASES); /* count the number of bases */
  for(i=0, ID=0; i<CODON_LENGTH; i++){ /* read the first 3 lettern in 'line' */
    ID*=codon_number;
    legal_letter=0; /* reset base flag */
    for(j=0; j<codon_number; j++){
      base=tolower(line[i]);
      if(base=='u') base='t';
      if(base==BASES[j]){
	ID+=j;
	legal_letter=1; /* found legal letter */
	break;
      }
    }
    if(!legal_letter) err=1; /* could not find legal letter -> found err */
  }
  if(err) ID=-1; /* with illegal letter the codon is meaningless */
  return(ID);
}

/* This function converts number to docon */
char *decode_codon(int ID){
  char *codon;
  int i, j;
  int odd, base;
  int numOfBases;
  numOfBases=strlen(BASES);
  codon=(char *)malloc_(__LINE__, (CODON_LENGTH+1)*sizeof(char));
  
  /*  for(i=16, j=0; j<CODON_LENGTH; i/=numOfBases, j++){*/
  for(i=pow(numOfBases, CODON_LENGTH-1), j=0; j<CODON_LENGTH; i/=numOfBases, j++){
    base = ID / i;
    odd  = ID % i;
    ID-= base * i;
    if(base>=0 && base<numOfBases){
      codon[j]=BASES[base];
    }else{
      fprintf(stderr, "base = %d is not adequate!\n", base);
      exit(0);
      /*      codon[j]='N';*/
    }
  }
  codon[j]='\0';
  return(codon);
}

/* This function creates the table of 
   "number of synonimous sites in a codon" */
double *create_synonimous_table(BASIC_TABLES tables, SWITCHES_ switches){
  double *synonimous_table;
  char *codon=NULL;
  char swap_codon[4];
  char bases[]="atgc";
  int codonID, swap_codonID;
  char modifiedCodon;
  int j, k;
  char s, n;
  double synonimous;

  if(switches.debug){
    printf("====== Creating Synonimous Table ======\n");
  }
  synonimous_table=(double *)malloc_(__LINE__, sizeof(double)*ALL_CODON);
  for(codonID=0; codonID<ALL_CODON; codonID++){ /* generate all possible codon 1 */
    free_safe(&codon);
    codon=decode_codon(codonID);
    if(tables.codonTable[codonID]=='\0'){
      if(switches.debug){
	printf("%s is non coding\n", codon);
      }
      synonimous_table[codonID]=-1;
      continue;
    }
    synonimous=0;
    for(j=0; j<CODON_LENGTH; j++){   /* scan substitution position */
      strcpy(swap_codon, codon);     /* copy original codon to temporary variable */
      s=0;
      n=0;
      for(k=0; k<4; k++){            /* substitute nucleotide at j-th codon */
	swap_codon[j]=bases[k];
	swap_codonID=encode_codon(swap_codon);
	if(strcmp(swap_codon, codon)==0){
	  /* swap_codon and codon is the same */
	  continue;
	}else{
	  /* 
	     If you want to eliminate the terminal codon
	     set ELIMINATE_TERMINAL_CODON = 1 or 
	     this program will count the rout which passes
	     terminal codon as same as other route.
	     --- the following is written in Japanese ---
	     ELIMINATE_TERMINAL_CODON = 1｡｡､ﾋ､ｹ､・ﾈｽｪｻﾟ･ｳ･ﾉ･ﾌ､・ﾐﾏｩ､・
	     ｿｨ､ﾊ､､､隍ｦ､ﾋ､ｹ､・｣､筅ｷ0､ﾋ､ｹ､・ﾈ｡｢ｽｪｻﾟ･ｳ･ﾉ･ﾌ､・ﾐﾏｩ､・
	     ﾂｾ､ﾎｷﾐﾏｩ､ﾈﾆｱﾍﾍ､ﾋｰｷ､ｦ｡｣
	  */
	  if(switches.ELIMINATE_TERMINAL_CODON && tables.codonTable[swap_codonID]=='\0'){
	    if(switches.debug){
	      printf("position %d %s swapped codon is terminal codon\n", j, swap_codon);
	    }
	  }else{
	    if(switches.debug){
	      modifiedCodon=tables.codonTable[swap_codonID];
	      if(modifiedCodon=='\0') modifiedCodon='*';
		printf("position %d %s (%c) and %s (%c)"
		       , j, codon, tables.codonTable[codonID]
		       , swap_codon, modifiedCodon);
	    }
	    if(tables.codonTable[codonID]==tables.codonTable[swap_codonID]){
	      if(switches.debug){
		printf(" is synonimous\n");
	      }
	      s++;
	    }else{
	      if(switches.debug){
		printf(" is NON-synonimous\n");
	      }
	      n++;
	    }
	  }
	}
      }
      if(switches.debug){
	printf("synonimous: %d / %d\n", s, (s+n));
      }
      synonimous+=(double)s/(s+n);
    }
    if(switches.debug){
      printf("%s synonimous=%f\n", codon, synonimous);
    }
    synonimous_table[encode_codon(codon)]=synonimous;
  }
  free_safe(&codon);
  return(synonimous_table);
}

sUBtABLE create_pathway_table(BASIC_TABLES tables, SWITCHES_ switches){
  int i, j;
  char *codon=NULL, *changed_codon=NULL;
  sUBtABLE result;
  double **synonimous, **nonsynonimous, **synonimousSiteNum;
  double *tmpresult;
  struct ROUTE route;

  if(switches.debug){
    printf("====== Creating Pathway Table ======\n");
  }
  /* initiate */
  route.synonimous=(char *)malloc_(__LINE__, MAXLETTER*sizeof(char));
  route.nonsynonimous=(char *)malloc_(__LINE__, MAXLETTER*sizeof(char));
  route.weight=(double *)malloc_(__LINE__, sizeof(double)*MAXLETTER);
  route.synonimousSiteNum=(double *)malloc_(__LINE__, sizeof(double)*MAXLETTER);
  route.depth=(int *)malloc_(__LINE__, sizeof(int)*MAXLETTER);
  route.synonimous[0]=-1; route.nonsynonimous[0]=-1;
 /* table of synonimouse change etc */
  synonimous       =(double **)malloc_(__LINE__, sizeof(double *)*ALL_CODON);
  nonsynonimous    =(double **)malloc_(__LINE__, sizeof(double *)*ALL_CODON);
  synonimousSiteNum=(double **)malloc_(__LINE__, sizeof(double *)*ALL_CODON);

  for(i=0; i<ALL_CODON; i++){   /* create all types of codon */
    codon=decode_codon(i);
    synonimous[i]       =(double *)malloc_(__LINE__, sizeof(double)*ALL_CODON);
    nonsynonimous[i]    =(double *)malloc_(__LINE__, sizeof(double)*ALL_CODON);
    synonimousSiteNum[i]=(double *)malloc_(__LINE__, sizeof(double)*ALL_CODON);
    if(tables.codonTable[i]=='\0'){ /* find terminal codon */
      if(switches.debug){
	printf("----- %s is terminal codon -----\n\n", codon);
      }
      free_safe(&codon);
      continue;
    }
    for(j=0; j<ALL_CODON; j++){ /* create all types of codon */
      changed_codon=decode_codon(j);
      if(switches.debug){
	printf("----- %s(%c) to %s(%c) ------\n"
	       , codon, tables.codonTable[i]
	       , changed_codon, tables.codonTable[j]);
      }
     
      route.synonimous[0]=0;                  /* initiation of route */
      route.nonsynonimous[0]=0 ;              /* initiation of route */
      route.synonimous[1]=-1;                 /* initiation of route */
      route.nonsynonimous[1]=-1;              /* initiation of route */
      route.weight[0]=1;                      /* initiation of route */
      route.number=0; route.route[0][0]='\0'; /* initiation of route */
      route.synonimousSiteNum[0]=tables.tableS[i];   /* initiation of route */
      route.depth[0]=1;                       /* initiation of route */
      /* check all types of route form "codon" to "changed_codon" */
      route=trace_route(route, (char)0, (char)0, codon, changed_codon, 
			tables, switches);
      /* calculate synonimouse and nonsynonimous change */
      tmpresult=average_synonimous_change(route, switches);
      synonimous[i][j]       =tmpresult[0];
      nonsynonimous[i][j]    =tmpresult[1];
      synonimousSiteNum[i][j]=tmpresult[2];
      if(switches.debug){
	printf("%s to %s -> s=%.3f, a=%.3f, synonimous sites=%.3f\n\n"
	       , decode_codon(i), decode_codon(j), 
	       synonimous[i][j], nonsynonimous[i][j], synonimousSiteNum[i][j]);
      }
      free_safe(&changed_codon);
    }
    free_safe(&codon);
  }
  result.synonimous=synonimous;
  result.nonsynonimous=nonsynonimous;
  result.synonimousSiteNum=synonimousSiteNum;
  return result;
}

struct ROUTE trace_route(struct ROUTE route, char synonims_org, char nonsynonims_org, 
			 char *Codon1, char *Codon2, BASIC_TABLES tables, SWITCHES_ switches){
  char *codon1, *codon2;
  char currentRoute[MAXLETTER], RouteBKUP[MAXLETTER];
  double weightBKUP, synonimousSiteNumBKUP;
  int i, depthBKUP;
  char synonims, nonsynonims;
  int s, t;

  /* If two codons are same do nothing and return */
  if(strcmp(Codon1, Codon2)==0){
    if(switches.debug){
      printf("Same codon\n");
    }
    return route;
  }
  /* If the changed codon is terminal codon do nothing */
  if(tables.codonTable[encode_codon(Codon2)]=='\0'){
    if(switches.debug){
      printf("changed codon is terminal codon\n");
    }
    return route;
  }
  /* If current codon is terminal codon do nothing and return */
  if(tables.codonTable[encode_codon(Codon1)]=='\0'){
    if(switches.debug){
      printf("current codon is terminal codon\n");
    }
    return route;
  }

  /* create temporal data for Codon1 and Codon2 */
  codon1=(char *)malloc_(__LINE__, (strlen(Codon1)+1)*sizeof(char));
  codon2=(char *)malloc_(__LINE__, (strlen(Codon2)+1)*sizeof(char));
  strcpy(codon1, Codon1);
  strcpy(codon2, Codon2);
  /* backup original intermedial route and weight*/
  strcpy(RouteBKUP, route.route[route.number]);
  weightBKUP=route.weight[route.number];
  synonimousSiteNumBKUP=route.synonimousSiteNum[route.number];
  depthBKUP=route.depth[route.number];

  for(i=0; i<CODON_LENGTH; i++){ /* for all codon position */
    synonims   =   synonims_org; /* call original amount of    synonimous changes */
    nonsynonims=nonsynonims_org; /* call original amount of nonsynonimous changes */
    strcpy(codon1, Codon1);      /* recall original codon */
    strcpy(route.route[route.number], RouteBKUP); /* recall original intermedial route */
    route.weight[route.number]=weightBKUP;        /* recall original intermedial weight */
    route.synonimousSiteNum[route.number]=synonimousSiteNumBKUP; /* recall original */
    route.depth[route.number]=depthBKUP;          /* recall original intermedial depth */

    /* If terninal codon is found, give up this route 
       and go to the next route.  */
    if(tables.codonTable[encode_codon(codon1)]=='\0'){
      sprintf(currentRoute,"%s(term)\n", codon1);
      strcat(route.route[route.number], currentRoute);
      if(switches.debug){
	printf("%s", route.route[route.number]);
      }
      continue; /* terminal codon */
    }

    if(codon1[i]!=codon2[i]){    /* i-th codon is different */
      sprintf(currentRoute,"%s(%c) -> "
	      , codon1, tables.codonTable[encode_codon(codon1)]);
      strcat(route.route[route.number], currentRoute);
      /* substitute i-th codon of codon1 to codon2   */
      codon1[i]=codon2[i];
      /* increment synonimous site number */
      route.synonimousSiteNum[route.number]+=tables.tableS[encode_codon(codon1)];
      /* increment route depth */
      route.depth[route.number]++;
      if(tables.codonTable[encode_codon(codon1)]=='\0'){ /* the substituted codon become terminal codon */
	sprintf(currentRoute, "%s(term)\n",codon1);
	strcat(route.route[route.number], currentRoute);
	if(switches.debug){
	  printf("%s", route.route[route.number]);
	}
	continue; /* the codon become terminal codon. go to the next pathway */
      }
      s=tolower(tables.codonTable[encode_codon(Codon1)])-'a'; /* Amino Acid of Original    Codon */
      t=tolower(tables.codonTable[encode_codon(codon1)])-'a'; /* Amino Acid of Substituted Codon */
      /* multiply the weight to intermedial weight */
      route.weight[route.number]*=tables.substitutionMatrix[s][t];
      if(s==t){ /* this path is synonimous change */
	strcat(route.route[route.number], "*");
	synonims++;
      }else{    /* this path is nonsynonimous change */
	nonsynonims++;
      }
      if(strcmp(codon1, codon2)==0){ /* the substitution made the two codon to the same one */
	route.synonimous[route.number]=synonims;       /* final decision of    synonimous changes*/
	route.nonsynonimous[route.number]=nonsynonims; /* final decision of nonsynonimous changes */
	if(switches.debug){
	  printf("%s", route.route[route.number]);
	  printf("%s(%c) done! %d Syn, %d nonSyn route%d , weight=%f\n"
		 , codon1, tables.codonTable[encode_codon(codon1)]
		 , synonims, nonsynonims, route.number, route.weight[route.number]);
	}
	route.number++;                       /* preparation to the next route */
	synonims=synonims_org;                /* preparation to the next route */
	nonsynonims=nonsynonims_org;          /* preparation to the next route */
	route.synonimous[route.number]=-1;    /* initiation of the next route data */
	route.nonsynonimous[route.number]=-1; /* initiation of the next route data */
	route.weight[route.number]=1;         /* initiation of the next route data */
      }else{
	/* Could not match the codon1 and codon2. Need more substitution. */
	route=trace_route(route, synonims, nonsynonims, codon1, codon2, tables, switches);
      }
    }
  }
  return route;
}

/* calculate synonimouse and nonsynonimous change
   from the struct ROUTE */
double *average_synonimous_change(struct ROUTE route, SWITCHES_ switches){
  int i;
  double synonimous, nonsynonimous;
  double *result;
  double sumWeight;
  double synonimousSiteNum;

  result=(double *)malloc_(__LINE__, sizeof(double)*3);
  for(i=0, synonimous=0, nonsynonimous=0, sumWeight=0, synonimousSiteNum=0;
      route.synonimous[i]!=-1; i++){
    synonimous+=route.synonimous[i]*route.weight[i];
    nonsynonimous+=route.nonsynonimous[i]*route.weight[i];
    synonimousSiteNum+=route.synonimousSiteNum[i]/route.depth[i]*route.weight[i];;
    sumWeight+=route.weight[i];
  }
  if(i!=0 && sumWeight!=0){ /* escape from division by zero. Because i==0 means sumWeight==0 */
    result[0]=synonimous/sumWeight;
    result[1]=nonsynonimous/sumWeight;
    result[2]=synonimousSiteNum/sumWeight;
  }else{
    result[0]=0;
    result[1]=0;
    result[2]=0;
  }
  if(switches.debug){
    printf("num of routes=%d,  synonimous=%f,  nonsynonimous=%f,  sumWeight=%f\n"
	   ,i , synonimous, nonsynonimous, sumWeight);
    printf("   synonimous: %f / %f = %f\n", synonimous, sumWeight, result[0]);
    printf("nonsynonimous: %f / %f = %f\n", nonsynonimous, sumWeight, result[1]);
  }
  return result;
}

struct SEQ_DATA read_and_check_multifasta_sequence(SWITCHES_ switches){
  struct SEQ_DATA Sequence;
  Sequence=read_input_multifasta_sequence(switches);
  /* reject input if sequence length is not the same */
  check_sequences(Sequence, switches);
  return Sequence;
}


struct SEQ_DATA read_input_multifasta_sequence(SWITCHES_ switches){
  struct SEQ_DATA result;
  struct SEQUENCE one_sequence;
  int num_of_sequence=SEQ_NUM;
  FILE *fin;
  int i;

  result.sequence=(struct SEQUENCE *)malloc
    (sizeof(struct SEQUENCE)*num_of_sequence);
  if(switches.infile!=NULL){
    if(switches.debug){
      printf("Opening file: %s\n", switches.infile);
    }
    if ((fin=fopen(switches.infile,"r"))==NULL){
      fprintf(stderr,"Could not open infile.[%s]\n",switches.infile);
      exit(0);
    }
  }else{
    if(switches.debug){
      printf("input is stdin\n");
    }
    fin=stdin;
  }

  i=-1; result.number=i;
  result.sequence[0].sequence=NULL;
  cue(fin); /* search for header of the multifasta sequence */
  while(!feof(fin)){
    one_sequence=read_one_sequence(fin);
    if(one_sequence.name!=NULL){
      i++; result.number=i+1;
      if(num_of_sequence<i+2){ /* realloc */
	num_of_sequence+=SEQ_NUM;
	if((result.sequence=
	    (struct SEQUENCE *)realloc
	    (result.sequence, sizeof(struct SEQUENCE)*num_of_sequence))
	   ==NULL){
	  fprintf(stderr, "Memory err, on reallocing result.sequence\n");
	};
      }
      result.sequence[i]=one_sequence;
    }
  }
  fclose(fin);
  return(result);
}

struct SEQUENCE read_one_sequence(FILE *fin){
  struct SEQUENCE result;
  char *name, *sequence;
  char c;
  int i, length;

  result.name=NULL;
  result.sequence=NULL;
  result.length=0;
  length=MAXLETTER;

  /* read name */
  name=(char *)malloc_(__LINE__, length*sizeof(char));
  for(i=0;(c=fgetc(fin))!=EOF;){
    if(c=='\r' || c=='\n'){
      break;
    }
    if(length<i+2){
      length+=MAXLETTER;
      if((name=(char *)realloc(name, length*sizeof(char)))==NULL){
	fprintf(stderr,"Memory err! At line %d. Not enough memory?\n",__LINE__);
      }
    }
    name[i++]=c;
    name[i]='\0';
  }
  if(feof(fin)){
    return result;
  }
  /* read sequence */
  length=MAXLETTER;
  sequence=(char *)malloc_(__LINE__, length*sizeof(char));
  for(i=0; (c=fgetc(fin))!=EOF;){
    if(c=='>'){
      break;
    }
    if(i>length-1){
      length+=MAXLETTER;
      if((sequence=(char *)realloc(sequence, length*sizeof(char)))==NULL){
	fprintf(stderr,"Memory err! Not enough memory?\n");
      }
    }
    if(isalpha(c) || c=='-'){
      sequence[i++]=c;
      sequence[i]='\0';
    }
  }
  result.name=name;
  result.sequence=sequence;
  result.length=i;
  return result;
}

void help(void){
  fprintf(stderr, "  ======  ABOUT sqdifC =====\n");
  fprintf(stderr, " Sqdif reads nucleotide sequences of multifasta format \n");
  fprintf(stderr, " and calculates the sequence difference of all possible \n");
  fprintf(stderr, " sequence pairs based on the method of");
  fprintf(stderr, " Miyata and Yasunaga (1980). \n");
  fprintf(stderr, "Reference:\n");
  fprintf(stderr, " Miyata, T., and Yasunaga, T., 1980. Molecular evolution of mRNA: A method\n");
  fprintf(stderr, " for estimating evolutionary rates of synonimous and amino acid substitutions\n");
  fprintf(stderr, " from homologous nucleotide sequences and its application. J. Mol. Evol.\n");
  fprintf(stderr, "16: 23-36\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " Sqdif is originally written in FORTRAN in 1982 and re-written \n");
  fprintf(stderr, " in C in 2008 by Genome Information research Center (this version).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE:\n");
  fprintf(stderr, "-i Infile. stdin if default.\n");
  fprintf(stderr, "-o Outfile. stdout if default. (not yet implemented)\n");
  fprintf(stderr, "-f Output format text (T) or tabular (X). T if default.\n");
  fprintf(stderr, "-w Window size to. This value must be dividable by 3. Zero if default.\n");
  fprintf(stderr, "-s Slide size. This value must be dividable by 3. Zero if default. \n");
  fprintf(stderr, "-r Consider substitution route in counting synonimous sites.\n");
  fprintf(stderr, "     [T/F] T if default. F if Nei-Gojobori method.\n");
  fprintf(stderr, "-c [T/F] Calculate difference of all combination [T] or Compare the first sequence with others [F].  T if default.\n");
  fprintf(stderr, "-W Substitution route is weighted. [T/F] T if default. F if Nei-Gojobori method.\n");
  fprintf(stderr, "-T Eliminate route of ternilal codon when creating the table of synonimous sites. \n");
  fprintf(stderr, "     [T/F]. F if default. T if Nei-Gojobori method.\n");
  fprintf(stderr, "-C Codon table file. Currently '%s'.\n", CODONTABLE);
  fprintf(stderr, "-S Substitution table file. Currently '%s'.\n", SUBSTITUTION_TABLE);
  fprintf(stderr, "-D show debug information [T/F] T if default. \n");
  fprintf(stderr,"%s compiled on %s %s\n", __FILE__, __DATE__, __TIME__);
  /*
  fprintf(stderr, "- \n");
  fprintf(stderr, "- \n");
  fprintf(stderr, "- \n");
  fprintf(stderr, "- \n");
  */
  exit(0);
}

void err_in_codon_table(void){
  fprintf(stderr, "Codon Table is wrong: %s\n", CODONTABLE);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(0);
}

int check_sequences(struct SEQ_DATA Sequence, SWITCHES_ switches){
  int i;
  if(switches.debug){
    printf("Number of Sequence is %d\n", Sequence.number);
  }
  for(i=0; i<Sequence.number; i++){
    if(Sequence.sequence[0].length!=Sequence.sequence[i].length){
      fprintf(stderr, "Sequences must be the same lengthes\n");
      exit(0);
    }
  }
  return 0;
}

int free_safe(char **P){
  int result=0;
  if(*P!=NULL){
    free(*P);
    *P=NULL;
    result=1;
  }
  /*
  if(result==0){
    printf("Did not perform free\n");
  }
  */
  return(result);
}

int cue(FILE *fin){
  char a, b='\n';
  while((a=fgetc(fin))!=EOF){
    if(a=='>' && b=='\n'){
      break;
    }
    b=a;
  }
  return 1;
}

void *malloc_(int line, size_t size){
  void *pointer;
  pointer=malloc(size);
  if(pointer==NULL){
    fprintf(stderr, "Could not allocate memory called in line %d\n", line);
    exit(1);
  }
  return(pointer);
}
