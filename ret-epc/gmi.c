/* tc array for GMI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define RADEARTH 6371
#define RTD 57.29578
#define DTR 0.017453

void tc_gmi(float *Tc_S1_buf, float *Tc_S2_buf, 
          int nl_S1, int ns_S1, int nc_S1, 
          int nl_S2, int ns_S2, int nc_S2, 
          float *lat_S1_buf, float *lon_S1_buf, float *lat_S2_buf, float *lon_S2_buf,
          char *qual_S1_buf, char *qual_S2_buf, 
          float *incidenceAngle_S1_buf, float *incidenceAngle_S2_buf,
          float *tc)
{

  int nc_out= nc_S1 + nc_S2;
  int i,j,k,idx,idx1, idx2, qual, i2min, j2min, i2, j2, itest;
  float tb, tb2[nc_out];
  float lat1, lon1, lat2, lon2, km, km_min, inc, inc2;

  for (i=0; i< nl_S1*ns_S1*nc_out; i++) tc[i]= -999.0;

  for (i=0; i< nl_S1; i++) {
    for (j=0; j<ns_S1; j++) {
      for (k=0; k< nc_out; k++) tb2[k]= -999.0;
      idx= i*ns_S1+j;
      lat1= lat_S1_buf[idx];
      lon1= lon_S1_buf[idx];
      inc= incidenceAngle_S1_buf[idx];
      if ( fabs(lat1) > 90 || fabs(lon1) > 180 ) continue;
      /*printf("%4d %3d %8.3f %8.3f %8.3f ", i,j, inc, lat1, lon1);*/

      qual= qual_S1_buf[idx];
      /*printf("%4d ", qual);*/
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
      /*for (k= 0; k< nc_S1; k++) printf("%7.2f ", tb2[k]);*/

      km_min= 1000.0;
      i2min= j2min= -1;
      /* Matching S2 line is 4-12 scans before S1 scan, and S2 sample is within +/-28 samples of S1 sample */
      /* but the spacecraft can re-orient so make it +/-16 lines to be sure */
          
      /* works but only fore direction     for (i2=i; i2 >= i-16; i2--) { */
      for (itest=0; itest<= 32; itest++) {
        if ( km_min < 6.0 ) break;
        i2= ( itest % 2 == 0 ) ? i-(itest/2) : i+(itest/2)+1;
        if ( i2 < 0 || i2 >= nl_S2 ) continue;
        for (j2=j-30; j2 <= j+30; j2++) {
          if ( j2 < 0 || j2 >= ns_S2 ) continue;
          idx2= i2*ns_S2 + j2;
          qual= qual_S2_buf[idx2];
          if ( qual < 0 ) continue;
          lat2= lat_S2_buf[idx2];
          lon2= lon_S2_buf[idx2];
          if ( fabs(lat2) > 90 || fabs(lon2) > 180 ) continue;
          if ( fabs(lat1-lat2) > 1.0 ) continue;
          km= RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2));
          if ( km < km_min ) {
            km_min= km;
            i2min= i2; j2min= j2;
            if ( km_min < 6.0 ) break;
          }
        }
      }

      if ( km_min < 7.0 && i2min >= 0 && j2min >= 0 ) {
        idx2= i2min*ns_S2 + j2min;
        lat2= lat_S2_buf[idx2];
        lon2= lon_S2_buf[idx2];
        inc2= incidenceAngle_S2_buf[idx];
        qual= qual_S2_buf[idx2];
        /*printf("%4d ", qual);*/
        if ( qual >= 0 ) {
          for (k= 0; k< nc_S2; k++) {
            idx1= i2min*ns_S2*nc_S2 + j2min*nc_S2 + k;
            tb= Tc_S2_buf[idx1];
            if ( tb > 20 && tb < 400 ) {
              if ( qual > 0 ) tb*= -1.0;
              tb2[k+nc_S1]= tb;
            }
          }
        }
        /*printf("%4d %3d %8.3f %8.3f %8.3f km=%8.3f ", i2min,j2min, inc2, lat2, lon2, km_min);*/
        /*for (k= 0; k< nc_S2; k++) printf("%7.2f ", tb2[nc_S1+k]);*/
      }

      for (k= 0; k< nc_out; k++) {
        idx1= i*ns_S1*nc_out + j*nc_out + k;
        tc[idx1]= tb2[k];
        /*printf("%7.2f ", tb2[k]);*/
      }
      /*printf("\n");*/

    }
  }

}

