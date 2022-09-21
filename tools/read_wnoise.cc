#include<iostream>
#include<fstream>
#include<vector>
#include<cstdio>
#include<stdlib.h>

int main(int argc, char *argv[]){
  unsigned Ni, Nj, Nk;
  char fname[128];

  sprintf(fname, "%s.bin", argv[1]);
  int downsampling_factor = atoi(argv[2]);

	std::ifstream ifs(fname,std::ios::binary|std::ios::in);

	ifs.read( reinterpret_cast<char*> (&Ni), sizeof(unsigned) );
  ifs.read( reinterpret_cast<char*> (&Nj), sizeof(unsigned) );
  ifs.read( reinterpret_cast<char*> (&Nk), sizeof(unsigned) );

  std::cout << "Ni = " << Ni << std::endl;
  std::cout << "Nj = " << Nj << std::endl;
  std::cout << "Nk = " << Nk << std::endl;

  // for(int i = 0; i < Ni*Nj*Nk; i++){
  //   double x;
	//   ifs.read(reinterpret_cast<char*> (&x), sizeof(double) );

  //   if (i % downsampling_factor == 0){
  //     std::cout << x << std::endl;
  //   }
  // }

  int flag_i, flag_j, flag_k;
  for (int i = 0; i < Ni; i++){
    flag_i = i % downsampling_factor == 0 ? 1 : 0;
    for (int j = 0; j < Nj; j++){
      flag_j = j % downsampling_factor == 0 ? 1 : 0;
      for (int k = 0; k < Nk; k++){
        flag_k = k % downsampling_factor == 0 ? 1 : 0;

        double x;
	      ifs.read(reinterpret_cast<char*> (&x), sizeof(double) );
        if ( flag_i == 1 && flag_j == 1 && flag_k == 1){
          std::cout << x << std::endl;
        }

      }
    }
  }
	ifs.close();
}    
  
