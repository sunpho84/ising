#include <src/config.hpp>

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

/// Measurements
int cachedEnergy;
int cachedSpinSum;
vector<int> energy;
vector<double> magnetization;

/// Configuration
vector<Spin> conf;
int nAcc;

/// Helpers
int V;
vector<Neighs> neighs;
vector<mt19937> gen;
array<double,17> expTable;

///////////////////////// Time measurements ////////////////////////////////////////

/// Alias for moment in time
using Instant=chrono::steady_clock::time_point;

/// Get current moment
Instant now()
{
  return chrono::steady_clock::now();
}

/// Takes the difference between two instants as second
double durationInSec(const Instant beg,const Instant end)
{
  return chrono::duration<double>(end-beg).count();
}

///////////////////////////// Geometry over the lattice ////////////////////////////////////

/// Computes the coordinate of a given site
Coords coordsOfSite(const Site site)
{
  return {site%L,site/L};
}

/// Computes the site given the coordinates
Site siteOfCoords(const Coord x,const Coord y)
{
  return x+L*y;
}

//////////////////////////// Random number generation /////////////////////////////////////

/// Returns +/-1
Spin drawRndSpin(mt19937& gen)
{
  return uniform_int_distribution<int>(0,1)(gen)*2-1;
}

/// Returns uniform number distributed between [0;1)
double drawUniformNumber(mt19937& gen)
{
  return uniform_real_distribution<double>(0,1)(gen);
}

///////////////////////// Measurements ////////////////////////////////////////

/// Computes the sum of spin
Spin measureSpinSum()
{
  Spin sumSpin=0;
  
  for(Site site=0;site<V;site++)
    sumSpin+=conf[site];
  
  return sumSpin;
}

/// Compute the magnetization
double measureMagnetization()
{
  return measureSpinSum()/(double)V;
}

/// Compute the cached magnetization
double getCachedMagnetization()
{
  return cachedSpinSum/(double)V;
}

/// Computes the total energy
int measureEnergy()
{
  int energy=0;
  
  for(Site site=0;site<V;site++)
    for(Dir dir=0;dir<2;dir++)
      energy-=conf[site]*conf[neighs[site][dir]];
  
  return energy;
}

/// Computes the energy relative to a given site
int energyOfSite(const Site site)
{
  int energy=0;
  
  for(Dir dir=0;dir<4;dir++)
      energy-=conf[site]*conf[neighs[site][dir]];
  
  return energy;
}

/////////////////////////////////////////////////////////////////

/// Setup the simulation
void setup()
{
  nAcc=0;
  
  magnetization.reserve(nEvol);
  energy.reserve(nEvol);
  
  /// Compute total volume
  V=L*L;
  
  /// Reset the number generator
  gen.resize(V);
  for(Site site=0;site<V;site++)
    gen[site].seed(inputSeed+site);
  
  /// Resize conf and lookup table for neighbors
  conf.resize(V,+1);
  neighs.resize(V);
  
  /// Draw a conf
  for(Site site=0;site<V;site++)
    conf[site]=drawRndSpin(gen[site]);
  
  /// Define site neighbors
  for(int site=0;site<V;site++)
    {
      const Coords coords=coordsOfSite(site);
      const Coord x=coords[0],y=coords[1];
      
      neighs[site][0]=siteOfCoords((x+1)%L,y);
      neighs[site][1]=siteOfCoords(x,(y+1)%L);
      neighs[site][2]=siteOfCoords((x+L-1)%L,y);
      neighs[site][3]=siteOfCoords(x,(y+L-1)%L);
    }
  
  for(int i=0;i<=16;i++)
    expTable[i]=exp(-Beta*(i-8));
  
  cachedEnergy=measureEnergy();
  cachedSpinSum=measureSpinSum();
}

/////////////////////////////////////////////////////////////////

/// Gets probability to accept/rejecta
double getPacc(const int dE)
{
  if(useLookupTableAcceptReject)
    return expTable[dE+8];
  else
    return exp(-Beta*dE);
}

/// Holds the results of the update of a site
struct SiteUpdRes
{
  bool acc;
  Spin dM;
  int dE;
};

/// Accept/reject
SiteUpdRes acceptReject(const int dE,const Site site,const Spin oldSpin)
{
  const double pAcc=getPacc(dE);
  
  const double r=drawUniformNumber(gen[site]);
  const bool acc=(r<pAcc);
  
  //nAcc+=acc;
  
  if(not acc)
    {
      conf[site]=oldSpin;
      return {0,0,0};
    }
  else
    return {1,conf[site]-oldSpin,dE};
}

/////////////////////////////////////////////////////////////////

/// Updates a single site
SiteUpdRes updateSite(const Site site)
{
  const Spin oldSpin=conf[site];
  int oldEn;
  if(useLocalEnergyChange)
    oldEn=energyOfSite(site);
  else
    oldEn=measureEnergy();
  
  conf[site]=drawRndSpin(gen[site]);
  
  int newEn;
  if(useLocalEnergyChange)
    newEn=energyOfSite(site);
  else
    newEn=measureEnergy();
  
  return acceptReject(newEn-oldEn,site,oldSpin);
}

/////////////////////////////////////////////////////////////////

/// Fully sweep a configuration
void updateConf()
{
  for(int par=0;par<2;par++)
    {
#pragma omp parallel for reduction(+:cachedEnergy,cachedSpinSum,nAcc)
      for(Site redSite=0;redSite<V/2;redSite++)
	{
	  const Coord y=redSite/L;
	  const Coord x=2*redSite%(L/2)+(par^(y%2));
	  const Site site=siteOfCoords(x,y);
	  const SiteUpdRes res=updateSite(site);
	  nAcc+=res.acc;
	  if(useMeasurementCache)
	    {
	      cachedEnergy+=res.dE;
	      cachedSpinSum+=res.dM;
	    }
	}
    }
}

/// Performs the measurements
void makeMeasurements()
{
  if(useMeasurementCache)
    {
      magnetization.push_back(getCachedMagnetization());
      energy.push_back(cachedEnergy);
    }
  else
    {
      magnetization.push_back(measureMagnetization());
      energy.push_back(measureEnergy());
    }
}

/// Store the configuration
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
  
  auto beg=now();
  for(int iEvol=0;iEvol<nEvol;iEvol++)
    {
      makeMeasurements();
      updateConf();
    }
  auto end=now();
  cout<<durationInSec(beg,end)<<" seconds"<<endl;
  cout<<"Pacc: "<<(double)nAcc/(nEvol*V)<<endl;
  
  storeConf();
  
  /////////////////////////////////////////////////////////////////
  
  ofstream magnetizationFile("magnetization.dat");
  ofstream energyFile("energy.dat");
  for(int iEvol=0;iEvol<nEvol;iEvol++)
    {
      magnetizationFile<<iEvol<<" "<<magnetization[iEvol]<<endl;
      energyFile<<iEvol<<" "<<energy[iEvol]<<endl;
    }
  
  return 0;
}
