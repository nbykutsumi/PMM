/* tc array for ATMS */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void tc_atms(float *Tc_S1_buf, float *Tc_S2_buf, float *Tc_S3_buf, float *Tc_S4_buf,
          int nl_S1, int ns_S1, int nc_S1, 
          int nl_S2, int ns_S2, int nc_S2, 
          int nl_S3, int ns_S3, int nc_S3, 
          int nl_S4, int ns_S4, int nc_S4, 
          float *lat_S1_buf, float *lon_S1_buf,
          char *qual_S1_buf, char *qual_S2_buf, char *qual_S3_buf, char *qual_S4_buf, 
          float *incidenceAngle_S1_buf,
          float *tc)
{

  int nc_out= nc_S1 + nc_S2 + nc_S3 + nc_S4;
  int i,j,k,idx,idx1, qual;
  float tb, tb2[nc_out];
  float lat1, lon1, inc;

  for (i=0; i< nl_S1*ns_S1*nc_out; i++) tc[i]= -999.0;
  for (i=0; i< nc_out; i++) tb2[i]= -999.0;

  for (i=0; i< nl_S1; i++) {
    for (j=0; j<ns_S1; j++) {
      idx= i*ns_S1+j;
      lat1= lat_S1_buf[idx];
      lon1= lon_S1_buf[idx];
      inc= incidenceAngle_S1_buf[idx];
      if ( fabs(lat1) > 90 || fabs(lon1) > 180 ) continue;

      qual= qual_S1_buf[idx];
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S1; k++) {
          idx1= i*ns_S1*nc_S1 + j*nc_S1 + k;
          tb= Tc_S1_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k]= tb;
          }
        }
      }

      qual= qual_S2_buf[idx];
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S2; k++) {
          idx1= i*ns_S2*nc_S2 + j*nc_S2 + k;
          tb= Tc_S2_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1]= tb;
          }
        }
      }

      qual= qual_S3_buf[idx];
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S3; k++) {
          idx1= i*ns_S3*nc_S3 + j*nc_S3 + k;
          tb= Tc_S3_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1+nc_S2]= tb;
          }
        }
      }

      qual= qual_S4_buf[idx];
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S4; k++) {
          idx1= i*ns_S4*nc_S4 + j*nc_S4 + k;
          tb= Tc_S4_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1+nc_S2+nc_S3]= tb;
          }
        }
      }

      printf("%4d %4d %8.3f %8.3f %8.3f ", i,j,lat1,lon1,inc);
      for (k= 0; k< nc_out; k++) {
        idx1= i*ns_S1*nc_out + j*nc_out + k;
        tc[idx1]= tb2[k];
        printf("%7.2f ", tb2[k]);
      }
      printf("\n");

    }
  }

}

