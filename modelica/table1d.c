/* 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or    
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of    
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
GNU General Public License for more details.

You should have received a copy of the GNU General Public License    
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "table1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include <time.h>

int log_this(const char *msg)
{
  FILE* flog = fopen("/tmp/.log","a");
  fprintf(flog, msg);
  fclose(flog);
  return 0;
}  

int log_int(const int msg)
{
  FILE* flog = fopen("/tmp/.log","a");
  fprintf(flog, "%d", msg);
  fclose(flog);
  return 0;
}  

int initTable(const char* filename, const double offset)
{
  log_this("Init. <");log_this(filename);log_this(">\n");
  if(filename[0] == 0) return 0;

  double u,y;
  double **data = NULL;
  
  // Open file
  FILE* fd = fopen(filename,"r");

  int n = 0;
  while(fscanf(fd,"%lf,%lf",&u,&y)!=EOF){
    data = (double**) realloc(data, (n+1)*sizeof(double*));
    data[n] = (double*) calloc(2,sizeof(double));

    data[n][0] = u - offset;
    data[n][1] = y;
    
    n++;
  }

  // Close file
  fclose(fd);

  srandom(time(NULL));

  int key_info = 
    random();

  int shmid_info = 
    shmget((key_t) key_info, 2 * sizeof(int), 0644 | IPC_CREAT | IPC_EXCL);

  while((shmid_info < 0) && (errno == EEXIST)){
    key_info = random();
    shmid_info = shmget((key_t) key_info, 2 * sizeof(int), 0644 | IPC_CREAT | IPC_EXCL);
  }

  log_this("Creado shm info\n"); log_int(key_info); log_this("\n");

  int *info = 
    (int*) shmat(shmid_info, NULL, 0);
  info[0] = n;

  int key_data = 
    random();

  int shmid_data = 
    shmget((key_t) key_data, 2 * n * sizeof(double), 0644 | IPC_CREAT | IPC_EXCL);

  while((shmid_data < 0) && (errno == EEXIST)){
    key_data = random();
    shmid_data = shmget((key_t) key_data, 2 * n * sizeof(double), 0644 | IPC_CREAT | IPC_EXCL);
  }
  log_this("Creado shm data: "); log_int(key_data); log_this("\n");

  double *shdata = 
    (double*) shmat(shmid_data, NULL, 0);

  int i;
  for(i=0;i<n;i++){
    shdata[2*i] = data[i][0];
    shdata[2*i+1] = data[i][1];
  }

  shmdt(shdata);

  info[1] = key_data;
  shmdt(info);

  log_this("N: ");log_int(n);log_this("\nData: ");log_int(key_data);log_this("\n");
  log_this("Fin init\n");

  return key_info;
}

double interpolateTable(const int key, const double u)
{
  //  log_this("Interpolate <");log_int(key);log_this(">\n");
  // BUG??: I don't know why but sometimes this is called without a valid key.
  if(key == 0) return 0;

  int shmid_info = shmget((key_t) key, 2 * sizeof(int), 0644);
  if(shmid_info < 0) exit(-1);
  //  log_this("shm info");

  int *info = (int*) shmat(shmid_info, NULL, 0);
  if(info==(int*)-1) exit(-1);

  int n = info[0];
  int key_data = info[1];

  shmdt(info);

  int shmid_data = shmget((key_t) key_data, 2 * n * sizeof(double), 0644);
  if(shmid_data < 0) exit(-1);
  //  log_this("shm data");
  
  double *data = (double*) shmat(shmid_data, NULL, 0);
  if(data == (double*) -1) exit(-1);

  double val;

  if(n==1){
    //  Constant function
    val = data[1];
  }
  else if(u<=data[0]){
    // u is less than min{u}
    val = (data[3]-data[1])/(data[2]-data[0])*u + (data[1]*data[2] - data[0]*data[3])/(data[2]-data[0]) ;

  }else if(u > data[2*n-2]){
    // us is grater than max{u}
    val = ((data[2*n-1]-data[2*n-3])/(data[2*n-2]-data[2*n-4]))*u + (data[2*n-3]*data[2*n-2] - data[2*n-4]*data[2*n-1])/(data[2*n-2]-data[2*n-4]);

  } else {
    int i;
    for(i=0;i<n-1;i++){
      if(u<data[2*i+2]){
	val = ((data[2*i+3]-data[2*i+1])/(data[2*i+2]-data[2*i]))*u + (data[2*i+1]*data[2*i+2] - data[2*i]*data[2*i+3])/(data[2*i+2]-data[2*i]);
	break;
      }
    }

  }

  shmdt(data);
  return val;
}

int stopTable(const int key)
{
  if(key == 0) return 0 ;

  int shmid_info = shmget((key_t) key, 2 * sizeof(int), 0644);
  if(shmid_info < 0) exit(-1);

  int *info = (int*) shmat(shmid_info, NULL, 0);
  if(info==(int*)-1) exit(-1);

  int n = info[0];
  int key_data = info[1];

  shmdt(info);
  shmctl(shmid_info, IPC_RMID, NULL);

  int shmid_data = shmget((key_t) key_data, 2 * n * sizeof(double), 0644);
  if(shmid_data < 0) exit(-1);

  shmctl(shmid_data, IPC_RMID, NULL);  

  return 0;

}
