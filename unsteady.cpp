/*************************************************************************
   MAIN CODE:  code.cpp (C++)
   VERSION:  4.0 (1/15/14)
   AUTHOR:  Xiaofan Yang
   DESCRIPTION:  Main code: UMSM
**7***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"
#include "poisson.h"
#include "grid.h"
#include "plot.h"
#include "function.h"
#include "global.h"
#include <gsl.h>

main(void) {
  FILE *fpt, *fcp;
  int time, stime=1;
  int itime = itimeNo;
  int irestart = restartOrNot;
  int i, j, k, strtNo, cnt, cnt2;

  dx = dxValue;
  dy = dyValue;
  dz = dzValue;
  dt = 100;
  period = dt*36*24;
  fpt = fopen("result/vel", "w");
  fcp = fopen("result/pre", "w");

/*--- Initaialize the whole domain ---*/
  initialize();
  if(irestart==0) {
    soil();
  }
  else {
    readRestart(&stime);
  }
  grid();
  stag();
  time = stime;
  itime = stime+itime-1; 
  plot3dl1(0);
  while(time <= itime) { //loop1 for level1 grid
    printf("\n%d ", time);
//  start calculation
    boundary();
    caconv();
    //moving(time);
    findFc();
    bforce();
    fforce();
    caff();
    poisson2();
    boundary();
    caff();
    ca();
    boundary();
    poisson();
    boundary();
    checkCont(time);
    cleanup();
    newToOld1();
    if(time%outputInt == 0) {
      plot3dl1(time/outputInt);
    }
    time++;
  }

/**** END of the initialization ****/

/*--- Start computation ---*/

/*
    plotu(fpt);
    plotcp1(fcp);
*/
/**** End of Calation Loop ****/

/*--- Output ---*/
  fclose(fcp);
  fclose(fpt);
  restart(time);
/**** End of Output ****/
}
