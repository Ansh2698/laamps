#include<bits/stdc++.h>
#include <fstream>
#include <stdlib.h>
#include <string>
using namespace std;
// Reading the File water.txt
double box_size=62.00;
int number_frames=3;
double number_bins=200;
double inc=box_size/(2*number_bins);
double r_cut=box_size/2;
const double PI=3.141592;
vector<string> Read_Data(){
    fstream lammps;
    vector<string>vec;
    lammps.open("water.txt",ios::in);
    if (lammps.is_open()){   //checking whether the file is open
      string data;
      while(getline(lammps,data)){
          vec.push_back(data);
      }
    }
    lammps.close(); //close the file object.
    return vec;
}
// Defining Coordinates for the Frame
vector<vector<double>> coordinates_for_frame(vector<string>Data,int start,int end){
    vector<vector<double>>Coordinates;
    for(int i=start;i<=end;i++){
        string d=Data[i];
        stringstream is(d);
        string s;
        int count=0;
        vector<double>vec;
        while(is>>s){
            if(count==1 && s!="2"){  // Taking Only oxygen atom
                break;
            }
            if(count>=3){
                vec.push_back(stof(s.c_str()));
            }
            count++;
        }
        if(!vec.empty()){
            Coordinates.push_back(vec);
        }
    }
    return Coordinates;
}
//Distance Between the Atoms limiting with the Boundary
double Distance_bw_atoms(vector<double>a,vector<double>b){
    double x=abs(a[0]-b[0]);
    x=min(x,abs(box_size-x));
    cout<<x<<endl;
    double y=abs(a[1]-b[1]);
    y=min(y,abs(box_size-y));
    cout<<y<<endl;
    double z=abs(a[2]-b[2]);
    z=min(z,abs(box_size-z));
    cout<<z<<endl;
    return sqrt((x*x)+(y*y)+(z*z));
}
//Find the distribution of distances while looping the frames
void Distribution(vector<string>Data,double num_atom){
    vector<double>Distribution(number_bins,0);
    int n;
    for(int i=0;i<number_frames;i++){
        double st=(i*num_atom)+((i+1)*9);
        double end=st+num_atom;
        vector<vector<double>>Coord=coordinates_for_frame(Data,st,end);
        n=Coord.size();
        for(int j=0;j<=n-2;j++){
            for(int k=j+1;k<n;k++){
                double distance=Distance_bw_atoms(Coord[j],Coord[k]);
                if(distance<r_cut){
                    int a=(int)(distance/inc);
                    Distribution[a]=Distribution[a]+2;
                }
            }
        }
    }
    //Normalization of Distribution
    for(int i=0;i<number_bins;i++){
        Distribution[i]=Distribution[i]/number_frames;
    }
    double density=n/(box_size*box_size*box_size);
    for(int i=0;i<number_bins;i++){
        double r1=Distribution[i];
        double r2=r1+inc;
        double vol_bin=(4/3)*PI*((r2*r2*r2)-(r1*r1*r1));
        double n_ideal=vol_bin*density;
        Distribution[i]=Distribution[i]/n_ideal;
    }
    return;
    // No only to show the plot of Distribution
}
int main(){
    vector<string>Data=Read_Data();
    double number_of_atom=stoi(Data[3]);
    Distribution(Data,number_of_atom);
    return 0;
}