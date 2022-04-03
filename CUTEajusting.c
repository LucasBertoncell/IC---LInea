#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"
// teste a CUTE 100 pontos

//  --------------------------io.c

typedef struct {//----------struct containing all possibly necessary max and min and nbins
  double dim1_max;
  double dim1_min;
  double dim1_minTheta; //------------------------------change, add int dim1_minTheta

  int dim1_nbin;
  double dim2_max;


  int dim2_nbin;
  double dim3_min;
  double dim3_max;
  int dim3_nbin;
  int logbin;

} Binner;

void process_binner(Binner binner) // binning function
{
  //////
  // Check that binning options make sense   logarthmic binnning doesnt have zero as min 

  //checking for correct input
  if(binner.logbin<0) {
    fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
    exit(1);
  }
  if(binner.dim1_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim1 %d\n",binner.dim1_nbin);
    exit(1);
  }
  if(binner.dim2_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim2 %d\n",binner.dim2_nbin);
    exit(1);
  }
  if(binner.dim3_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim3 %d\n",binner.dim3_nbin);
    exit(1);
  }
  if(binner.dim1_max<=0) {
    fprintf(stderr,"CUTE: wrong dim1_max %lf\n",binner.dim1_max);
    exit(1);
  }



  if(binner.logbin) {
    if((binner.dim1_min<=0) || (binner.dim1_min>=binner.dim1_max)) {
      fprintf(stderr,"CUTE: wrong lower limit for logarithmic binning %lf\n",binner.dim1_min);
      exit(1);
    }
  }


  if((binner.dim3_max<=0)||(binner.dim3_min<0)||
     (binner.dim3_max<=binner.dim3_min)) {
    fprintf(stderr,"CUTE: wrong boundaries for dim3 (%lf , %lf)\n",
	    binner.dim3_min,binner.dim3_max);
    exit(1);
  }

  logbin=binner.logbin;
  if(logbin)
    n_logint=binner.dim1_nbin/log10(binner.dim1_max/binner.dim1_min);
  if(corr_type==0) {
    nb_dz=binner.dim1_nbin;
    i_dz_max=1./binner.dim1_max;
  }
  else  if(corr_type==1) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
  }


  //-----------------------------------------------------

  else if(corr_type==2) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./(binner.dim1_max - binner.dim1_minTheta); //-------------------------change
    log_r_max=log10(binner.dim1_max);
  }

  /// acf  above  -----------------------------------------
  

  else if(corr_type==3) {
    nb_rl=binner.dim2_nbin;
    i_rl_max=1./binner.dim2_max;
    nb_rt=binner.dim1_nbin;
    i_rt_max=1./binner.dim1_max;
    log_rt_max=log10(binner.dim1_max);
  }

  else if(corr_type==4) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./binner.dim1_max;
    log_r_max=log10(binner.dim1_max);
    nb_mu=binner.dim2_nbin;


  else if(corr_type==5) {
    nb_theta=binner.dim1_nbin;

    i_theta_max=1./(DTORAD*binner.dim1_max-);
    log_th_max=log10(DTORAD*binner.dim1_max);
    nb_dz=binner.dim2_nbin;
    i_dz_max=1./binner.dim2_max;
    nb_red=binner.dim3_nbin;
    i_red_interval=1./(binner.dim3_max - binner.dim3_min);
    red_0=binner.dim3_min;
  }
  else {
    fprintf(stderr,"WTF!?\n"); ///??RUDE
    exit(1);
  }
}
// definicao do theta min = 0, delta = 
void read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  Binner binner;

  binner.dim1_max=-1;
  binner.dim1_min=-1;
  binner.dim1_minTheta=-1;//-------------------------change

  binner.dim2_max=-1;
  binner.dim3_min=-1;
  binner.dim3_max=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.logbin=-1;

  print_info("*** Reading run parameters \n");

  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"data_filename="))
      sprintf(fnameData1,"%s",s2);
    else if(!strcmp(s1,"data_filename_2="))
      sprintf(fnameData2,"%s",s2);
    else if(!strcmp(s1,"random_filename="))
      sprintf(fnameRandom1,"%s",s2);
    else if(!strcmp(s1,"random_filename_2="))
      sprintf(fnameRandom2,"%s",s2);
    else if(!strcmp(s1,"RR_filename="))
      sprintf(fnameRR,"%s",s2);
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"corr_type=")) {
      if(!strcmp(s2,"radial")) corr_type=0;
      else if(!strcmp(s2,"angular")) corr_type=1;
      else if(!strcmp(s2,"monopole")) corr_type=2;
      else if(!strcmp(s2,"3D_ps")) corr_type=3;
      else if(!strcmp(s2,"3D_rm")) corr_type=4;
      else if(!strcmp(s2,"full")) corr_type=5;
      else {
  fprintf(stderr,"CUTE: wrong corr type %s.",s2);
  fprintf(stderr," Possible types are \"radial\", \"angular\", \"full\",");
  fprintf(stderr," \"monopole\", \"3D_ps\" and \"3D_rm\".\n");
      }
    }
    else if(!strcmp(s1,"omega_M="))
      omega_M=atof(s2);
    else if(!strcmp(s1,"omega_L="))
      omega_L=atof(s2);
    else if(!strcmp(s1,"w="))
      weos=atof(s2);
    else if(!strcmp(s1,"radial_aperture="))
      aperture_los=atof(s2)*DTORAD;

    else if(!strcmp(s1,"dim1_max="))
      binner.dim1_max=atof(s2);


    else if(!strcmp(s1,"dim1_min=")) // ---------------------change
      binner.dim1_minTheta=atof(s2);


    else if(!strcmp(s1,"dim1_min_logbin="))
      binner.dim1_min=atof(s2);
    else if(!strcmp(s1,"dim2_max="))
      binner.dim2_max=atof(s2);
    else if(!strcmp(s1,"dim3_max="))
      binner.dim3_max=atof(s2);
    else if(!strcmp(s1,"dim3_min="))
      binner.dim3_min=atof(s2);
    else if(!strcmp(s1,"dim1_nbin="))
      binner.dim1_nbin=atoi(s2);
    else if(!strcmp(s1,"dim2_nbin="))
      binner.dim2_nbin=atoi(s2);
    else if(!strcmp(s1,"dim3_nbin="))
      binner.dim3_nbin=atoi(s2);
    else if(!strcmp(s1,"log_bin="))
      binner.logbin=atoi(s2);
    else if(!strcmp(s1,"use_pm="))
      use_pm=atoi(s2);
    else if(!strcmp(s1,"n_pix_sph=")) {
      n_side_cth=atoi(s2);
      n_side_phi=2*n_side_cth;
    }
    else
      fprintf(stderr,"CUTE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  process_binner(binner);

  if(strcmp(fnameData2,"file_none") || strcmp(fnameRandom2,"file_none"))
    use_two_catalogs=1;

  check_params();

  print_info("\n");
}

Catalog *read_catalog(char *fname,np_t *sum_w,np_t *sum_w2)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog *cat=my_malloc(sizeof(Catalog));

  print_info("*** Reading catalog ");
#ifdef _VERBOSE
  print_info("from file %s",fname);
#endif
  print_info("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  ng=linecount(fd);
  rewind(fd);
  print_info("  %d lines in the catalog\n",ng);

  //Allocate catalog memory
  cat->np=ng;
  cat->red=(double *)my_malloc(cat->np*sizeof(double));
  cat->cth=(double *)my_malloc(cat->np*sizeof(double));
  cat->phi=(double *)my_malloc(cat->np*sizeof(double));
#ifdef _WITH_WEIGHTS
  cat->weight=(double *)my_malloc(cat->np*sizeof(double));
#endif //_WITH_WEIGHTS

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  *sum_w=0;
  *sum_w2=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,weight;
    int st=read_line(fd,&zz,&cth,&phi,&weight);

    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    if(zz<0) {
      fprintf(stderr,"Wrong redshift = %lf %d\n",zz,ii+1);
      exit(1);
    }
    if((cth>1)||(cth<-1)) {
      fprintf(stderr,"Wrong cos(theta) = %lf %d\n",cth,ii+1);
      exit(1);
    }
    phi=wrap_phi(phi);

    cat->red[i_dat]=zz;
    cat->cth[i_dat]=cth;
    cat->phi[i_dat]=phi;
#ifdef _WITH_WEIGHTS
    cat->weight[i_dat]=weight;
    (*sum_w)+=weight;
    (*sum_w2)+=weight*weight;
#else //_WITH_WEIGHTS
    (*sum_w)++;
    (*sum_w2)++;
#endif //_WITH_WEIGHTS
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  print_info("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

#ifdef _WITH_WEIGHTS
  print_info("  Effective n. of particles: %lf\n",(*sum_w));
#else //_WITH_WEIGHTS
  print_info("  Total n. of particles read: %d\n",(*sum_w));
#endif //_WITH_WEIGHTS

  print_info("\n");
  return cat;
}

Catalog_f read_catalog_f(char *fname,int *np)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog_f cat;

  print_info("*** Reading catalog ");
#ifdef _VERBOSE
  print_info("from file %s",fname);
#endif
  print_info("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  ng=linecount(fd);
  *np=ng;
  rewind(fd);

  //Allocate catalog memory
  cat.np=ng;
  cat.pos=(float *)my_malloc(3*cat.np*sizeof(float));

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,rr,sth,dum_weight;
    int st=read_line(fd,&zz,&cth,&phi,&dum_weight);
    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    sth=sqrt(1-cth*cth);
    if(corr_type!=1)
      rr=z2r(zz);
    else
      rr=1;

    cat.pos[3*i_dat]=(float)(rr*sth*cos(phi));
    cat.pos[3*i_dat+1]=(float)(rr*sth*sin(phi));
    cat.pos[3*i_dat+2]=(float)(rr*cth);
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  print_info("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

  print_info("\n");
  return cat;
}
