#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <limits.h>


#define W_MAXDEL 10
#define W_MAXCOL 40
#define NELEMENT(A)  sizeof(A)/sizeof(A[0])
#define SIMPLE 231.0
#define FALSE 0
#define TRUE 1
#define MAXBASI 500000000


char *heder[]=
{
"\n\n",
" Pol3scan searches for the eukaryotic polymerase III transcription ",
"control elements B box, A box and terminator signal, discriminating",
"between tRNA genes and tRNA-like elements (e.g. SINEs).",
" \nRun program whith -h to see options.\n\n",
"\0"
};

typedef struct{
 float mat[4][W_MAXCOL];
 int   ncol;
 char  *mdel[W_MAXCOL+1];
 float probdel[W_MAXCOL+1];
 int   posdel[W_MAXDEL];
 int   ndel;
 }Weight_matrix;

typedef struct{
 int	start;
 int	end;
 int    range;
 float  vect[100];
 int 	ncol;
 }Weight_vector;

 typedef int (*fptr)(const void*,const void*);

void Look_for_option(int num,char *argv[]);
void print_message(char *mess[]);
void divide_dir(char *stringa,char *dir,char *nome);
void exiterror(const char *,const char *);
int numeric(const int *,const int *);

void Gen_comb_del(int **,const int,const int *);
void Gen_path_del(Weight_matrix *);

#define MAX(a,b)  (((a)>(b))?(a):(b))
#define limit1 "*****"
#define limit2 "*****"

int jump(FILE *,long,char,int);
int Read_gcg_matrix(char *nome,Weight_matrix *,int flag);
int Read_dist_vect(char *nome,Weight_vector *);
int Transform_log(float *vettore,int ncol,float str);
void Read_cod_tab(char *file_cod);
void stampa_matrix(Weight_matrix);




char el_nucleotidi[]="GATCN";
char terminator[]={2,2,2,2};




typedef struct{
 char nome[120];
 char *basi;
 int len;
 }sequenza;

typedef struct{
float B;
float A;
float A_B;
float T;
}Point;

typedef struct{
int pos;
int end;
float score;
char seq[W_MAXCOL+1];
}W_mot;


typedef struct{
char nome[120];
int len;
W_mot B_box;
W_mot A_box;
int A_B_dist;
int T_dist;
int T_start;
float A_B_score;
float Tot_score;
char strand;
char crt; /*  (* or +)  */

int tstart,tend;
char sequence[400];
char cod[4];
char aminoacid[4];
int stem_start;
int stem_end;

}Out;

/* new output */


void Look_for_format(FILE *,char *f);
int Read_the_sequence(FILE *,sequenza *p,char *f);
void Convert_the_sequence_into_number( char *base,int n);
void Invert_the_sequence(char base[],int n);
void Complement_the_sequence(char base[],int n);
void stampa_seq(sequenza s);

float Del_matrix_score( char *basi,Weight_matrix m,char *del,int ndel,float soglia);
float Weight_matrix_score( char *basi,Weight_matrix m,int ncol,float soglia);
float Dist_score(Weight_vector w,int dist);
int Find_terminator(sequenza seq,Point *s,Weight_vector Tim,int t);
int Find_selenocystein_Abox(sequenza seq,int t,Weight_matrix A_box,Point *sel);
float Aminoacyl_stem_score(sequenza seq,int posb,int posa);
float Stem_score(const char *a,const char *b,char sym_app[],int n);
void print_heder(FILE *fp);
void memo_out(Out *o,sequenza s,int pos_B,char *bv,int  pos_A,char *av,int T_pos,int strand,Point p,char c);
void new_stampa_out(FILE *fp,Out t);
void stampa_out(FILE *fp,Out t);


/* codon_table */
char gCODTAB[64][4];
char gAMINOTAB[64][4];
/* codon table */

/* Cutoff values */
#define gdefCT " [-31.25]"
#define gdefCA " [-28.00]"
#define gdefCB " [-14.14]"
#define gdefOUTF " [outfile]"
#define gdefINF  " [infile]"
#define gdefSM   " [1]"
#define gdefSV   " [1]"

float B_soglia=-14.14;
float A_soglia=-28.00;
float A_B_soglia=-31.25;
float Tot_soglia=-31.8;
const float Media_tim=19.70563;
const float Sdev_tim=13.78815;
float       stricym=1.0;
float	    stricyv=1.0;
int   gfold=FALSE; /* print or not the folding of aa stem */
int   grprint=FALSE; /* if true print a 1 row output */
#define   SOGLIA_TIM ((int)(Media_tim+Sdev_tim*2))
/* Cutoff values */

/* Filename */
char nomeseq[200]="infile";
char fmat_A_box[]="d_a_box.m";
char fmat_B_box[]="d_b_box.m";
char fvect_A_B[]="a_b_dist.v";
char fvect_tim[]="dist_tim.v";
char file_cod[]="cod.tab";
char filout[200]="outfile";
char no_filout[200]="n_";
/* Filename */


int main(int argc,char *argv[])
{
FILE *fpseq,*fo,*fno;

sequenza seq,*pseq;

Weight_matrix A_box,*pA,B_box,*pB;

Weight_vector A_B,*pA_B,Tim,*pTim;

Point score,*sc;

Out   tRNA,*ptRNA;

char nome[120];
char dir[100];
char format[15]="";

int t,i,d,strand;
int count=0;
int found=0;
int like_found=0;
int num_seq=0;
double num_basi=0;
int dist;
int A_len;
int A_start;
int A_end;
int forbidden; /* marker for the position of a tRNA*/
int A_del;
int A_pos;int T_pos;

float partial_A_score;
float p_A_soglia;
float best_A_B_score;
float amm_score;




pseq=&seq;
pA=&A_box; pB=&B_box; pA_B=&A_B ;ptRNA=&tRNA; pTim=&Tim;sc=&score;

score.A=0;score.B=0;score.A_B=0;

print_message(heder);

if(argc>1)
	Look_for_option(argc,argv);

if (strcmp(nomeseq,"infile")==0)
	{
	fprintf(stderr,"Insert a name for the input file : " );
	scanf("%s",nomeseq);
	}



if((fpseq=fopen64(nomeseq,"r"))==NULL)
       exiterror("\ncan't open file:\n",nomeseq);

if (strcmp(filout,"outfile")==0)
	{
	fprintf(stderr,"Insert a name for the output file : " );
	scanf("%s",filout);
	}
if((fo=fopen(filout,"w"))==NULL)
       exiterror("\ncan't open file:\n",filout);
if (grprint)
   {
   divide_dir(filout,dir,nome);
   sprintf(no_filout,"%s%s%s",dir,"n_",nome);

   if((fno=fopen(no_filout,"w"))==NULL)
       exiterror("\ncan't open file:\n",no_filout);
   }
Look_for_format(fpseq,format);

/* read all the sequence*/
while(Read_the_sequence(fpseq,pseq,format)>0)
	{
	num_seq++;
	num_basi+=(double)pseq->len+1;
	free(pseq->basi);
	}
rewind(fpseq);




/*read the B-box weight matrix   */
Read_gcg_matrix(fmat_B_box,pB,0);
Gen_path_del(pB);
for(t=0;t<4;t++)
	Transform_log(B_box.mat[t],B_box.ncol,stricym/SIMPLE);


/*read the A-box weight matrix*/
Read_gcg_matrix(fmat_A_box,pA,0);
Gen_path_del(pA);
for(t=0;t<4;t++)
	Transform_log(A_box.mat[t],A_box.ncol,stricym/SIMPLE);

Read_dist_vect(fvect_A_B,pA_B);
Transform_log(A_B.vect,A_B.ncol,stricyv/SIMPLE);

Read_dist_vect(fvect_tim,pTim);
Transform_log(Tim.vect,Tim.ncol,stricyv/SIMPLE);

Read_cod_tab(file_cod);

printf("\nSearching on: %s ...\n\n",nomeseq);

print_heder(fo);

while(Read_the_sequence(fpseq,pseq,format))
{
Convert_the_sequence_into_number(seq.basi,seq.len);

    for(strand=1;strand<=2;strand++)
	{
	forbidden=0;
	for(t=40+forbidden;t<seq.len-10;t++)
	      {
	      if ((score.B=Weight_matrix_score(seq.basi+t,B_box,11,B_soglia))>B_soglia)
		  {
		   if(score.B>-3.65 && score.B<-2.16)
			  if( (A_pos=(Find_selenocystein_Abox(seq,t,A_box,sc)))!=0)
			     {
			     score.A_B=score.B+score.A+Dist_score(A_B,t-A_pos-20);
			     T_pos=Find_terminator(seq,sc,Tim,t);
			     amm_score=Aminoacyl_stem_score(seq,t,A_pos);
			     memo_out(ptRNA,seq,t,B_box.mdel[0],A_pos,A_box.mdel[3],T_pos,strand,score,'*');
			     if (grprint) stampa_out(fo,tRNA);
			     else new_stampa_out(fo,tRNA);
			     forbidden=t+22;t=forbidden;found++;
			     continue;
			     }

		    best_A_B_score=INT_MIN;
		    A_len=strlen((char *)(A_box.mdel[0]+1))+1;
		    A_start=t-A_B.start-A_len;
		    if(forbidden>(A_end=t-A_B.end-A_len))
			A_end=forbidden;  /*dont look */
		    for(i=A_start;i>=A_end;i--) /* search for the A_box */
			{
			p_A_soglia=A_B_soglia-score.B+5;
			if((partial_A_score=Weight_matrix_score(seq.basi+i,A_box,10,p_A_soglia))>p_A_soglia)
				for(d=0;d<pow(2,A_box.ndel);d++) /* for all the deletion matrix*/
				   {
				     if ( (dist=t-i-strlen((char *)(A_box.mdel[d]+1))-1) >=A_B.start)
					if ((score.A=partial_A_score+Del_matrix_score(seq.basi+i+10,A_box,A_box.mdel[d]+10,d,A_soglia-partial_A_score))>A_soglia)
					{
					/*printf("\n A box : %s\t%d A_score %f",seq.nome,t+i,A_score);*/
					if((score.A_B=(score.B+score.A+Dist_score(A_B,dist)))>A_B_soglia)
						if(score.A_B>best_A_B_score)
							{
							best_A_B_score=score.A_B;
							A_del=d;A_pos=i;
							}
					}
				 } /*next d*/
			  } /*next i*/
			if(best_A_B_score>A_B_soglia)
				{
				score.A_B=best_A_B_score;tRNA.Tot_score=-99;
				T_pos=Find_terminator(seq,sc,Tim,t);
				amm_score=Aminoacyl_stem_score(seq,t,A_pos);
				memo_out(ptRNA,seq,t,B_box.mdel[0],A_pos,A_box.mdel[A_del],T_pos,strand,score,'+');
				if ((amm_score>=5 || amm_score==-1) &&((tRNA.Tot_score>Tot_soglia) || (tRNA.Tot_score==-99 && seq.len-3-t-11<SOGLIA_TIM)) )
					{
					if (grprint) stampa_out(fo,tRNA);
					else new_stampa_out(fo,tRNA);
					forbidden=t+22;
					t=forbidden;
					found++;
					/*Define the tRNA*/

					}
				else
					{
					if (amm_score<5) tRNA.crt='f';
					if(tRNA.Tot_score<=Tot_soglia && !(tRNA.Tot_score==-99 && seq.len-3-t-11<SOGLIA_TIM)) tRNA.crt='t';
					if(amm_score>=5) tRNA.crt='?';
					if(!strstr(tRNA.A_box.seq,"GGCCTACCAT"))
					       {if (grprint) stampa_out(fo,tRNA);
						else new_stampa_out(fo,tRNA);
						like_found++;}
					}

				}
		   } /* termin B-box*/

	      }/*next t*/
	 Invert_the_sequence(seq.basi,seq.len);
	 Complement_the_sequence(seq.basi,seq.len);
	 } /* next strand */

fprintf(stderr,"%-13s%6d/%d\ttRNA: %d\ttRNA-like: %d\r",seq.nome,++count,num_seq,found,like_found);
free(pseq->basi);  /*libera la memoria riservata alla sequenza */
}
fprintf(fo,"\n\nSequences searched: %d ( %.0f bp. both strands)\nPredicted tRNAs: %d\nPredicted tRNA-like elements: %d\n",num_seq,num_basi,found,like_found);
fclose(fpseq);fclose(fo);
if (grprint) fclose(fno);
fprintf(stderr,"\n\nResults written in %s file\n",filout);
return(0);
}

void print_message(char *mess[])
{
int i=0;
	while(*mess[i++]) fprintf(stderr,"\n%s",mess[i]);
return;
}
void print_heder(FILE *fp)
{



fprintf(fp,"________________________________________________________________________________\n");

fprintf(fp,"\n  Pol3scan");
fprintf(fp,"\n              Input file: %s",nomeseq);
fprintf(fp,"\n	       Cut-off value (total): %.2f", Tot_soglia);
fprintf(fp,"\n	       Cut-off value (A box): %.2f",A_soglia);
fprintf(fp,"\n	       Cut-off value (B box): %.2f", B_soglia);
/*fprintf(fp,"\n	       Stringency    (matrix): %.1f", stricym);
fprintf(fp,"\n               Stringency    (vector): %.1f",stricyv);  */
fprintf(fp,"\n\nRef:\nPavesi, A., Conterio, F., Bolchi, A., Dieci, G. and Ottonello, S., (1994)\n");
fprintf(fp," Identification of new eucariotic tRNA genes in genomic DNA databases by a\n multistep weight matrix");
fprintf(fp," analysis of transcriptional control regions,\n Nucleic Acids Res. 22, 1247-1256\n\n");

fprintf(fp,"________________________________________________________________________________\n");

return;
}
void Look_for_option(int num,char *argv[])
{
char *option[]=
{
" \n\n",
"  Pol3scan -o outfilename 	", gdefOUTF,
" \n\t-i infilename	", gdefINF,
" \n\tallowed formats: Genbank,EMBL,IG,fasta",
" \n\t-ct Cut-off value (total)", gdefCT,
" \n\t-ca Cut-off value (A box)", gdefCA,
" \n\t-cb Cut-off value (B box)", gdefCB,
" \n\t-sm stringency    (matrix)", gdefSM,
" \n\t-sv stringency    (vector)",gdefSV,
"\n\n",
"\0"
};

int i;
char c;
float s;
	for(i=1;i<num;i++)
		{
		if(strncmp(argv[i],"-",1)==0)
		   {
		   c=argv[i][1];
		   switch (c)
			{
			case 'r': /* 1 row output */
				grprint=TRUE;
				break;
			case 'f':
				gfold=TRUE;
				break;
			case 'i':
				strcpy(nomeseq,argv[i+1]);
				i+=1;
				break;
			case 'o':
				strcpy(filout,argv[i+1]);
				i+=1;
				break;
			case 's':
				if( (s=atof(argv[i+1]))==0)
					exiterror("stringency mast be a number"," >0");
				if(argv[i][2]=='v'|| argv[i][2]=='m' )
					(argv[i][2]=='m')?(stricym=s):(stricyv=s);
				else
					exiterror("invalid option: ",argv[i]);
				i+=1;
				break;
			case 'c':
				if( (s=atof(argv[i+1]))==0 || s>0)
					exiterror("A cut-off value mast be a number"," <0");
				argv[i][2]=toupper(argv[i][2]);
					switch(argv[i][2])
						{
						case 'B':
							B_soglia=s;
							break;
						case 'A':
							A_soglia=s;
							break;
						case 'T':
							A_B_soglia=s;
							break;
						 default : exiterror("invalid option: ",argv[i]);
						 }
				 i+=1;
				break;
			   case 'h': {printf("\r\r\r");print_message(option);exit(1);}

			default: {print_message(option);exit(1);}
		       }
		   }
		}
return;
}

void divide_dir(char *stringa,char *dir,char *nome)
{
   char tok[] = "\\/";
   char tempstr[40],*p;

   strcpy(tempstr,stringa);
   strcpy(dir,"");strcpy(nome,stringa);
   /* strtok places a NULL terminator
   in front of the token, if found */
   p = strtok(tempstr, tok);
   if (p)
	{
   /* A second call to strtok using a NULL
   as the first parameter returns a pointer
   to the character following the token  */
	while((p = strtok(NULL, tok))!=0)
		strcpy(nome,p);

	strncpy(dir,stringa,strlen(stringa)-strlen(nome));
	dir[strlen(stringa)-strlen(nome)]=0; /*termiate string*/
	}
return;
}
void Look_for_format(FILE *fp,char *format)
{
char stringa[200];

   format[0]=0;

	do{
	if(fscanf(fp,"%s",stringa)==EOF)
		exiterror("input file not in a known"," format");
	if (stringa[0]=='>') strcpy(format,"B");
	if (stringa[0]==';') strcpy(format,"I");
	if (strncmp(stringa,"ID",2)==0) strcpy(format,"EMBL");
	if (strncmp(stringa,"LOCUS",2)==0) strcpy(format,"GB");
	}while (!format[0]);
rewind(fp);
return;
}



void Gen_comb_del(int *cdel[],const int n,const int d[])
{
static int i=0;
static int p=0;
static int del[128];
int cont=0;
int t;

	while(cont<2)
	{
	       del[i]=d[cont];
			if (i<n-1)
				{
				i++;
				Gen_comb_del(cdel,n,d);
				}
			else
				{
					for (t=0;t<n;t++)
						cdel[p][t]=del[t];
				p++;
				}
	cont++;
	}
i--;
return;
}

int numeric(const int *p1,const int *p2)
{
	return(*p1-*p2);
}


/* This function returns the deletions matrix indicating the positions
whithout deletions; matdel is a sparse matrix of integer. The value of
1 means that the first position of the weight matrix mast be computed */

void Gen_path_del(Weight_matrix *pm)
{

const int sdel[]={0,1};

int *p;
int *cdel[128];
int t,i,j,k;
float sumdel;

	 if(!pm->ndel)     /*not a deletion matrix*/
		{
		pm->mdel[0]=(char *)malloc(sizeof(char)*((pm->ncol)+1));
		for(t=0;t<pm->ncol;t++)
			pm->mdel[0][t]=t;
		pm->mdel[0][t]=0;
		return;
		}
	 for (t=0;t<pow(2,pm->ndel);t++) /*reserve memory for cdel*/
	     cdel[t]=(int *)malloc(sizeof(int)*pm->ndel);

	 Gen_comb_del(cdel,pm->ndel,sdel);

	 for (t=0;t<pow(2,pm->ndel);t++)
	 {
	 j=-1;
	 pm->mdel[t]=(char *)malloc(sizeof(char)*((pm->ncol)+1));
	 pm->probdel[t]=0;
		for(i=0;i<pm->ncol; i++)
		{
			if ((p=(int *)bsearch(&i,pm->posdel,pm->ndel,sizeof(int),(fptr)numeric))!=0)
				if (cdel[t][p-(int *)&pm->posdel[0]]==0)
				{
				sumdel=0; /* deletion probability*/
				for(k=0;k<4;k++)
					sumdel+=pm->mat[k][i];
				pm->probdel[t]+=(float)log(1-sumdel);
				continue; /* deletion  */
				}
		j++;
		pm->mdel[t][j]=i;
		}
	      pm->mdel[t][++j]=0;
	  }
for (t=0;t<pow(2,pm->ndel);t++)
	free(cdel[t]);
return;
}

void exiterror(const char *a,const char *b)
{
fprintf(stderr,"%s%s\n",a,b);
exit(-1);
}

int jump(FILE *fp,long beg,char ch,int times)
{
int t;
char stringa[100];
fseek(fp,beg,0);

	for(t=0;t<=times;t++)
	  do{
	    if(fscanf(fp,"%s",stringa)==EOF)
		exiterror("\nnon c'e' il secondo separatore :\n",limit2);
	    if(strcmp(stringa,limit2)==0)
		return(0); /* stop */
	    }while(stringa[strlen(stringa)-1]!=ch);
return(1);
}

void Read_cod_tab(char *fname)
{
FILE *fp;
int i;
char cod[3];
char dummy[30];

if((fp=fopen(fname,"r"))==NULL)
       exiterror("\ncan't open file:\n",fname);
fgets(dummy,30,fp); /* read first line */

for (i=0;i<64;i++)
	{
	fscanf(fp,"%s",gAMINOTAB[i]);
	fscanf(fp,"%s",cod);
	Convert_the_sequence_into_number(cod,3);
	Invert_the_sequence(cod,3);
	Complement_the_sequence(cod,3);
	gCODTAB[i][0]=el_nucleotidi[cod[0]];
	gCODTAB[i][1]=el_nucleotidi[cod[1]];
	gCODTAB[i][2]=el_nucleotidi[cod[2]];
	gCODTAB[i][3]=0;
	}
return;
}


int Read_gcg_matrix(char *nf,Weight_matrix *p,int percent)
{
int t,i,volte;
FILE *fp;
char stringa[100];
long begin;
int data[W_MAXCOL];
int totcount[W_MAXCOL];
int maxcount=INT_MIN;

p->ndel=0;

if((fp=fopen(nf,"r"))==NULL)
       exiterror("\ncan't open file:\n",nf);

       do{
       if(fscanf(fp,"%s",stringa)==EOF)
	 exiterror("\nerrore nel formato del file:\n",nf);
       }while(strcmp(stringa,limit1)!=0);


begin=ftell(fp);   /*remember start point*/
volte=0;

	i=0;
	while(jump(fp,begin,'l',volte)!=0) /*finds the lines containing total*/
		{
		volte++;
		while(fscanf(fp,"%d",&totcount[i]))
			{
			maxcount=MAX(maxcount,totcount[i]);
			i++;
			}
		}
	p->ncol=i;
	totcount[i]=0;

	       for(i=0;totcount[i]!=0;i++)     /* deletion */
			if (totcount[i]<maxcount)
			{
			p->posdel[p->ndel]=i;
			p->ndel++;
			}

       for(t=0;t<4;t++)
		{
		volte=0;i=0;
		while(jump(fp,begin,el_nucleotidi[t],volte)!=0)
			{
			volte++;
			while(fscanf(fp,"%d",&data[i]))
				{
				p->mat[t][i]=data[i];
				if(p->mat[t][i]!=0)
					if(percent) /* percentual data ?*/
					{
					p->mat[t][i]*=totcount[i];
					p->mat[t][i]/=(100*maxcount);
					}
					else
					p->mat[t][i]/=maxcount;
				i++;
				}
			 p->mat[t][i]=0;
			 }
		}
fclose(fp);
return (0);
}

int Read_dist_vect(char *nf,Weight_vector *w)
{

FILE *fp;
int dummy;
char stringa[100];
int i,t,first,sec;
int tot;



if((fp=fopen(nf,"r"))==NULL)
       exiterror("\ncan't open file:\n",nf);

       do{
       if(fscanf(fp,"%s",stringa)==EOF)
	 exiterror("\nerrore nel formato del file:\n",nf);
       }while(strcmp(stringa,limit1)!=0); /*find 1ø marker*/

	if(fscanf(fp,"%d%d%f",&w->start,&sec,&w->vect[0])==EOF)
	  exiterror("\nerrore nel formato del file:\n",nf);

	w->range=sec-w->start;
	tot=(int)w->vect[0];
	i=0;
	do{
	  dummy=fgetc(fp);
	  if(dummy==EOF)
	    exiterror("\nerrore nel formato del file:\n",nf);
	     if(isdigit(dummy))
		{
		ungetc(dummy,fp);
		i++;
		if(fscanf(fp,"%d%d%f",&first,&sec,&w->vect[i])==EOF || sec-first!=w->range)
		       exiterror("\nerrore nel formato del file:\n",nf);
		tot+=(int)w->vect[i];
		}
	  }while(dummy!='*');

       w->end=sec;
       w->ncol=i+1;
       for(t=0;t<=i;t++)
	     if(w->vect[t]!=0)
		w->vect[t]/=(float)tot;




return(0);
}



int Transform_log(float *p,int n,float lim)
{
int i;
		for(i=0;i<n;i++)
			if(p[i]==0)
				 p[i]=log(lim);
			else
				 p[i]=log(p[i]);
return(1);
}

float Dist_score(Weight_vector w,int s)
{
int intrv;
	  if( (intrv=((int)floor((s-w.start-.1)/w.range))) <0)
		intrv=0;
	  return(w.vect[intrv]);
}

float Weight_matrix_score( char *basi,Weight_matrix m,int n,float soglia)
{
float score=0;
int i=0;

	   do{
	      score+=(*(basi+i)!=4)?(m.mat[*(basi+i)][i]):((float)log(.25));
	      i++;
	      } while(i<n && score>soglia);

return(score);
}

float Del_matrix_score( char *basi,Weight_matrix m,char *del,int ndel,float soglia)
{
float score;
int i=0;

score=m.probdel[ndel]; /*probability score for the matrix*/

	   do{
	      score+=(*(basi+i)!=4)?(m.mat[*(basi+i)][del[i]]):((float)log(.25));
	      i++;
	      } while(del[i]!=0 && score>soglia);

return(score);
}


int Find_selenocystein_Abox(sequenza seq,int t,Weight_matrix A_box,Point *sel)
{
const char s_motiv1[]={0,0,2,3,2,0,0,0,0,2};
const char s_motiv2[]={0,0,2,3,2,0,2,0,0,2};
const char s_motiv3[]={0,0,2,3,3,0,0,0,0,2};

int distsec1=30;
int distsec2=40;
int i;


for(i=t-distsec1-10;i>t-distsec2-10;i--)
	 if ((memcmp(seq.basi+i,s_motiv1,10)==0)||(memcmp(seq.basi+i,s_motiv2,10)==0)||(memcmp(seq.basi+i,s_motiv3,10)==0))
	     {
	     sel->A=Del_matrix_score(seq.basi+i-10,A_box,A_box.mdel[3],3,INT_MIN);
	     return(i-10);/*found !*/
	     }

return(0);
}



int Find_terminator(sequenza seq,Point *s,Weight_vector Tim,int t)
{
int i,T_end;

T_end=( (t+11+Tim.end) <(seq.len-4) )?(t+11+Tim.end):(seq.len);
	for(i=t+11+Tim.start;i<T_end-3;i++)
		if(memcmp(seq.basi+i,terminator,NELEMENT(terminator))==0)
			{
			s->T=s->A_B+Dist_score(Tim,i-t-11);
			return(i);   /* terminator position*/
			}
return(0);
}
float Stem_score(const char *a,const char *b,char sym_app[],int n)
{
int i;
float s_score=0;
	for(i=0;i<=n;i++)
		{
		sym_app[i]=' ';
		if(a[i]==4 || b[n-i]==4)
			{sym_app[i]='?';s_score+=(1.5)/4.0;}
		else
			{
			if(a[i]+b[n-i]==3) /* G/C or A/T */
				{s_score++;sym_app[i]='|';}
		/*	if(abs(a[i]-b[n-i])==2) */ /* G/T or A/C */
			if(abs(a[i]-b[n-i])==2 && a[i]+b[n-i]==2)  /* G/T */ /*modificato 06/06/2005 RIC*/
				{s_score+=0.5;sym_app[i]='.';}
			}
		}
sym_app[i]=0;
return s_score;
}




float Aminoacyl_stem_score(sequenza s,int posb,int posa)
{
int start_B_stem=posb+14;
int start_A_stem=posa-6;
int i,t,pos;
char sym_app[3][20],pred[20],event[20]="";
float d_score[3]={0.0,0.0,0.0};
float s_score=0;



		if((start_B_stem+6+1)>s.len)
			{if(gfold) printf("\nfew basis in the sequence %s\n",s.nome);return -1.0;}
       s_score=Stem_score(s.basi+start_A_stem,s.basi+start_B_stem,sym_app[1],6);
       pos=0;
       if (s_score<5) /* loock for deletion or insertion */
	 {
	  for(t=-1;t<2;t+=1)
		{
		d_score[1+t]=Stem_score(s.basi+start_A_stem,s.basi+start_B_stem+t,sym_app[1+t],6);
		 if(d_score[1+t]>s_score)
			{
			s_score=d_score[1+t];
			pos=t;
			}
		 }
	 }

	strcpy(pred,(s_score>=5)?("tRNA"):("not a tRNA"));

	if(pos==-1) strcpy(event,"Inserted base");
	if(pos==+1) strcpy(event,"Deleted base");

if (gfold)
 {
 printf("\n%s %f  \t\t%s  %s\n\t\t",s.nome,s_score,pred,event);
	for(i=0;i<8;i++) /* print last base too*/
		printf("%c",el_nucleotidi[s.basi[start_B_stem+7-i+pos]]);
 printf("\n\t\t %s\n\t\t ",sym_app[pos+1]);
	for(i=0;i<7;i++)
		printf("%c",el_nucleotidi[s.basi[start_A_stem+i]]);
 }


return(s_score);
}

void Invert_the_sequence(char base[],int n)
{
int i;
char t;
int halfbase;

	halfbase=floor(n/2);
	for(i=0;i<halfbase;i++)
	{
		t=base[i];
		base[i]=base[n-i-1];
		base[n-i-1]=t;
	}
}
void Complement_the_sequence(char base[],int n)
{
int t;
	for(t=0;t<n;t++)
	    if (base[t]!=4)
		base[t]=abs(base[t]-3);
}
void Convert_the_sequence_into_number( char *base,int n)
{
int i;
 char *p;


	for (i=0;i<n;i++)
	{
		if ((p=(char*)strchr(el_nucleotidi,*base))!=NULL)
			 *base=(p-el_nucleotidi);
		else    /* base N */
			 *base=4;

	 base++;
	 }
return;
}

int Read_the_sequence(FILE *fps,sequenza *pseq,char *format)
{

int number=0; /* len of the sequence */
int dummy;
char riga[2001]; /*modified 25/02/2004 */
char separetor[10];
char startchar[10]=";IL>";



	do{		   /*read until next sequence or EOF*/
	 dummy=fgetc(fps);
	       if(dummy==EOF)
		 return(0); /* all the file was read */
	  }while(strchr(startchar,(char)dummy)==NULL);
ungetc(dummy,fps);



	     switch (format[0])
		{
		case 'E':
			strcpy(separetor,"/");
			strcpy(startchar,"I");
			fgets(riga,2000,fps); /* read first line */
			if (strncmp(riga,"ID",2)!=0)
				{ printf("\nBad line in input file:\n%s\n",riga); return(-1);}
			if(!sscanf(riga,"%*s%s",pseq->nome))
				{ printf("\nBad line in input file:\n%s\n",riga); return(-1);}
			do{
			   if(fgets(riga,2000,fps)==NULL) { printf("\nBad line in input file\n"); return(-1);}
			  }while(strncmp(riga,"SQ",2)!=0);

		       break;
		case 'B':
			strcpy(separetor,">");
			strcpy(startchar,">");
			do{
			  if(fgets(riga,2000,fps)==NULL) exiterror("BAD LINE IN INPUT FILE","\n");
			  }while ( riga[0]!='>');
			if(!sscanf(riga+1,"%s",pseq->nome))
			     { printf("\nBad line in input file:\n%s\n",riga); return(-1);}
			break;
		case 'I':
			strcpy(separetor,"12");
			strcpy(startchar,";");
			do{
			  if(fgets(riga,2000,fps)==NULL) exiterror("BAD LINE IN INPUT FILE","\n");
			  }while ( riga[0]==';');
			if(!sscanf(riga,"%s",pseq->nome))
			     { printf("\nBad line in input file:\n%s\n"); return(-1);}
			break;
		case 'G':
			strcpy(separetor,"/L");
			strcpy(startchar,"L");
			 if(fgets(riga,2000,fps)==NULL)  exiterror("BAD LINE IN INPUT FILE","\n");
					if(!sscanf(riga,"%*s%s",pseq->nome))
				{ printf("\nBad line in input file:\n%s\n",riga); return(-1);}
			do{
			   if(fgets(riga,2000,fps)==NULL)  exiterror("BAD LINE IN INPUT FILE","\n");
			   }while (strstr(riga,"ORIGIN")==NULL);

			do{
			   dummy=fgetc(fps);
					if(dummy==EOF) break;

			   if (isalnum(dummy)&&(char)dummy!='1') if (fgets(riga,2000,fps)==NULL) exiterror("BAD LINE IN INPUT FILE","\n");
			   if ((char)dummy=='1') break;
			   } while(1);

			break;
		 }

if (NULL==(pseq->basi=(char *)malloc(sizeof(char)*MAXBASI+1))) exiterror(" ----- "," Out of memory");

	do{       /* read the nucleotides */
	  dummy=fgetc(fps);
		if(dummy==EOF) break;
	      /*	{ printf("\nSequence not in a correct format:\n%s\n",pseq->nome); return(-1);}*/
	  if(isalpha(dummy))
		{
		pseq->basi[number]=toupper((char)dummy);
		number++;
		 if(number>(MAXBASI-1))
			{fprintf(stderr,"%s is too long. The upper limit is %d bp. \n",pseq->nome,MAXBASI);exit(-1);}
		}
	   }while(strchr(separetor,(char)dummy)==NULL);
	   ungetc(dummy,fps);

pseq->len=number; /*-1; Modificato il 9/2/2004 RIC*/
pseq->basi[number]=0;
return(1);
}
void memo_out(Out *tRNA,sequenza s,int posb,char  bv[],int posa,char av[],int post,int strand,Point score,char c)
{
int i,d;
int no_term_score=-99;
int lenA=0;

     strcpy(tRNA->nome,s.nome);
     tRNA->len=s.len;
     /* B box */

     tRNA->B_box.score=score.B;
	 i=0;
	 do{
	   tRNA->B_box.seq[i]=*(s.basi+posb+bv[i]);
	   tRNA->B_box.seq[i]=el_nucleotidi[tRNA->B_box.seq[i]];
	   i++;
	   }while(bv[i]!=0);
	   tRNA->B_box.seq[i]=0;

     tRNA->B_box.pos=(strand>1)?(s.len-posb-1):(posb);
     tRNA->B_box.end=(strand>1)?(tRNA->B_box.pos-10):(tRNA->B_box.pos+10);

     /* A box */

     tRNA->A_box.score=score.A_B-score.B;
	i=0;d=0;
	 do{
	   if(av[i-d]!=i)
		{
		tRNA->A_box.seq[i]='-'; /*deleted position*/
		d++;
		}
	   else
		{

		tRNA->A_box.seq[i]=*(s.basi+posa+i-d);
		tRNA->A_box.seq[i]=el_nucleotidi[tRNA->A_box.seq[i]];
		}
	   i++;
	   }while(av[i-d]!=0);
	 tRNA->A_box.seq[i]=0;
      tRNA->A_B_score=score.A_B;
      lenA=i-d; /*i-d is the lenght of the A box*/

      tRNA->A_box.pos=(strand>1)?(s.len-posa-1):(posa);
      tRNA->A_box.end=(strand>1)?(tRNA->A_box.pos-lenA+1):(tRNA->A_box.pos+lenA-1);

      tRNA->A_B_dist=posb-(posa+lenA);

      /* terminator signal */
      if (post) /* there is a terminator */
	{
	tRNA->Tot_score=score.T;
	tRNA->T_start=(strand>1)?(s.len-post-1):(post);
	tRNA->T_dist=post-posb-11;
	}
      else
	{
	tRNA->Tot_score=no_term_score;
	tRNA->T_start=0;
	}

      tRNA->strand=(strand>1)?'C':'D';

      tRNA->crt=c;

      for (i=0;i<3;i++)
	{
	tRNA->cod[i]=*(s.basi+posa+lenA+8+i);
	tRNA->cod[i]=el_nucleotidi[tRNA->cod[i]];
	}
      tRNA->cod[3]=0;

	  for(i=0;i<64;i++)
	      if(strncmp(tRNA->cod,gCODTAB[i],3)==0)
		{
		strcpy(tRNA->aminoacid,gAMINOTAB[i]);
		break;
		}


      for (i=posa-6;i<=posb+21;i++){
      	   /*Bug in Sel tRNA sequence Modificato RIC*/
	   tRNA->sequence[i-(posa-6)]=el_nucleotidi[s.basi[i]];
	}
       tRNA->sequence[i-(posa-6)]=0;

      if (tRNA->strand=='D')
	{
	tRNA->tstart=tRNA->A_box.pos+1-6;
	tRNA->tend=tRNA->B_box.pos+1+21;
	}
	else
	{
	tRNA->tend=tRNA->A_box.pos+1+6;
	tRNA->tstart=tRNA->B_box.pos+1-21;
	}
	 if (tRNA->tstart<=1) tRNA->tstart=1; /*to compensate for a Bug Modificato RIC*/
	 if (tRNA->tend>=s.len) tRNA->tend=s.len; /*to compensate for a Bug Modificato RIC*/

return;
}




void stampa_out(FILE *fp,Out t)
{
	fprintf(fp,"%s, %u  %d %s,%.6f ",t.nome,t.len,t.B_box.pos+1,t.B_box.seq,t.B_box.score);
	if(t.Tot_score==-99)
		fprintf(fp," %d %s,%.5f %.0f %d %d %c,%c \n",t.A_box.pos+1,t.A_box.seq,t.A_B_score,t.Tot_score,t.A_B_dist,t.T_dist,t.strand,t.crt);
	else
		fprintf(fp," %d %s,%.5f %.5f  %d  %d %c,%c \n",t.A_box.pos+1,t.A_box.seq,t.A_B_score,t.Tot_score,t.A_B_dist,t.T_dist,t.strand,t.crt);
}

char *strlwr(char *s)
{
int i;
 for(i=0;s[i]!=0;i++)
 if(s[i]<='Z' && s[i]>='A') s[i]=s[i]-'A'+'a';
return(s);
}



void new_stampa_out(FILE *fp,Out t)
{
static char oldname[120];
char reason[100];

fprintf(fp,"\n\nsequence name: %s",t.nome);
fprintf(fp,"\nstrand: %s\n",(t.strand=='D'?"direct":"complementary"));
fprintf(fp,"B box signal:%4d %4d ; %-22s ;%14s %6.2f\n",t.B_box.pos+1,t.B_box.end+1,strlwr(t.B_box.seq),"score:",t.B_box.score);
fprintf(fp,"A box signal:%4d %4d ; %-22s ;%14s %6.2f %s\n",t.A_box.pos+1,t.A_box.end+1,strlwr(t.A_box.seq),"score:",t.A_box.score,(t.crt=='*')?"(selenocysteine motif)":"");
if (t.crt=='+' ||t.crt=='*') /* tRNA*/
	{
	if (t.T_start>0)
	fprintf(fp,"Term. signal:%4d %4d ; %-22s ;%14s %6.2f\n",t.T_start+1,(t.strand=='D'?t.T_start+1+3:t.T_start+1-3),"tttt","score (Tot):",t.Tot_score);
	else
	fprintf(fp,"Term. signal: %-33s ;%14s %6.2f\n", "not enough bases on 3' end","score (A-B):", t.A_B_score);
	fprintf(fp,"  Prediction: ---- %s ----\n","tRNA");
	fprintf(fp,"%13s %d %d\n","bounds:",t.tstart,t.tend);
	fprintf(fp,"%13s %s\n","sequence:",strlwr(t.sequence));
	if (t.crt=='*') fprintf(fp,"%13s %s  (tRNA-%s) \n","anticodon:","tca","SeC(e)");
		else fprintf(fp,"%13s %s  (tRNA-%s) \n","anticodon:",strlwr(t.cod),t.aminoacid);
	}
else /*tRNA-like*/
	{
	if (t.T_start>0)
	fprintf(fp,"Term. signal:%4d %4d ; %-22s ;%14s %6.2f\n",t.T_start+1,(t.strand=='D'?t.T_start+4:t.T_start-3),"tttt","score (Tot):",t.Tot_score);
	else
	fprintf(fp,"Term. signal: %-33s ;%14s %6.2f\n", "not found","score (A-B):", t.A_B_score);

	if(t.crt=='t'||t.crt=='?') strcpy(reason,(t.T_start>0)?"Termination signal is too distant":"Termination signal not found");
	else strcpy(reason,"Insufficient base-pairing in the aminoacyl stem");
	fprintf(fp,"  Prediction: ---- %s %s ----\n","tRNA-like element",reason);
	fprintf(fp,"%13s %d %d\n","bounds:",t.tstart,t.tend);
	fprintf(fp,"%13s %s\n","sequence:",strlwr(t.sequence));
        fprintf(fp,"%13s %s  (tRNA-%s) \n","anticodon:","???","Pseudo");
        }
strcpy(oldname,t.nome);
return;
}
