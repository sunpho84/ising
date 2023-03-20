#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>

using namespace std;

/// Custom types
using Coord=int;
using Dir=int;
using Site=int;
using Spin=int;
using Coords=array<Coord,2>;
using Neighs=array<Site,4>;

/// Input parameters
int L;
double Beta;
int nEvol;
int inputSeed;

/// Configuration
vector<Spin> conf;

/// Helpers
int V;
vector<Neighs> neighs;
mt19937 gen;
array<double,17> expTable;

using Instant=chrono::steady_clock::time_point;

/// Machine clock type
Instant now()
{
  return chrono::steady_clock::now();
}

double durationInSec(const Instant beg,const Instant end)
{
  return chrono::duration<double>(end-beg).count();
}

Coords coordsOfSite(const Site site)
{
  return {site%L,site/L};
}

Site siteOfCoords(const Coord x,const Coord y)
{
  return x+L*y;
}

Spin getRndSpin()
{
  return uniform_int_distribution<int>(0,1)(gen)*2-1;
}

double magnetization()
{
  int sumSites=0;
  
  for(Site site=0;site<V;site++)
    sumSites+=conf[site];
  
  return (double)sumSites/V;
 }

double energy()
{
  int energy=0;
  
  for(Site site=0;site<V;site++)
    for(Dir dir=0;dir<2;dir++)
      energy-=conf[site]*conf[neighs[site][dir]];
  
  return energy;
}

double energyOfSite(const Site site)
{
  int energy=0;
  
  for(Dir dir=0;dir<4;dir++)
      energy-=conf[site]*conf[neighs[site][dir]];
  
  return energy;
}

void setup()
{
  /// Compute total volume
  V=L*L;
  
  /// Reset the number generator
  gen.seed(inputSeed);
  
  /// Resize conf and lookup table for neighbors
  conf.resize(V,+1);
  neighs.resize(V);
  
  /// Draw a conf
  for(Site site=0;site<V;site++)
    conf[site]=getRndSpin();
  
  /// Define site neighbors
  for(int site=0;site<V;site++)
    {
      const Coords coords=coordsOfSite(site);
      const Coord y=coords[0],x=coords[1];
      
      neighs[site][0]=siteOfCoords((x+1)%L,y);
      neighs[site][1]=siteOfCoords(x,(y+1)%L);
      neighs[site][2]=siteOfCoords((x+L-1)%L,y);
      neighs[site][3]=siteOfCoords(x,(y+L-1)%L);
    }
  
  for(int i=0;i<=16;i++)
    expTable[i]=exp(-Beta*(i-8));
}

void updateSite(const Site site)
{
  // const double oldEn=energy();
  const double oldSiteEn=energyOfSite(site);
  const Spin oldSpin=conf[site];
  
  conf[site]=getRndSpin();
  // const double newEn=energy();
  const double newSiteEn=energyOfSite(site);
  
  // cout<<newEn-oldEn<<" "<<newSiteEn-oldSiteEn<<endl;
  
  //const double pAcc=exp(-Beta*(newSiteEn-oldSiteEn));
  // const double pAcc=exp(-Beta*(newEn-oldEn));
  const double pAcc=expTable[newSiteEn-oldSiteEn+8];
  //cout<<newEn-oldEn+8<<" "<<pAcc<<" "<<pAcc2<<endl;
  const double r=uniform_real_distribution<double>(0,1)(gen);
  
  if(r>=pAcc)
    conf[site]=oldSpin;
}

void updateConf()
{
  for(Site site=0;site<V;site++)
    updateSite(site);
}

void storeConf()
{
  ofstream confFile("conf.dat");
  for(Site site=0;site<V;site++)
    if(conf[site]==+1)
      {
	const Coords c=coordsOfSite(site);
	confFile<<c[0]<<" "<<c[1]<<endl;
      }
}

int main()
{
  cout<<"L? "<<endl;
  cin>>L;
  cout<<"beta? "<<endl;
  cin>>Beta;
  cout<<"nEvol? "<<endl;
  cin>>nEvol;
  cout<<"seed? "<<endl;
  cin>>inputSeed;
  
  setup();
  
  ofstream out("measure.dat");
  auto beg=now();
  for(int iEvol=0;iEvol<nEvol;iEvol++)
    {
      out<<iEvol<<" "<<magnetization()<<endl;
      updateConf();
    }
  auto end=now();
  cout<<durationInSec(beg,end)<<" seconds"<<endl;
  
  storeConf();
  
  return 0;
}
