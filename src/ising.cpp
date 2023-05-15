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
vector<bool> parOfSite;
vector<int> eoSiteOfSite;
vector<Site> siteOfEoSite;
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

  for(int par=0;par<2;par++)
    for(Site eoSite=0;eoSite<V;eoSite++)
      sumSpin+=conf[eoSite];
  
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
  
  for(Site eoSite=0;eoSite<V;eoSite++)
    for(Dir dir=0;dir<2;dir++)
      energy-=conf[eoSite]*conf[neighs[eoSite][dir]];
  
  return energy;
}

/// Computes the energy relative to a given site
int energyOfSite(const Site eoSite)
{
  int energy=0;
  
  for(Dir dir=0;dir<4;dir++)
      energy-=conf[eoSite]*conf[neighs[eoSite][dir]];
  
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
  
  /// Resize conf and lookup table for neighbors
  conf.resize(V,+1);
  neighs.resize(V);
  
  parOfSite.resize(V);
  eoSiteOfSite.resize(V);
  for(Site site=0;site<V;site++)
    {
      const Coords coords=coordsOfSite(site);
      parOfSite[site]=(coords[0]+coords[1])%2;
    }
  
  siteOfEoSite.resize(V);
  
  for(Site eoSite=0;eoSite<V;eoSite++)
    {
      const int par=eoSite/(V/2);
      const int redSite=eoSite%(V/2);
      
      const Coord y=redSite/(L/2);
      const Coord x=2*(redSite%(L/2))+(par^(y%2));
      const Site site=siteOfCoords(x,y);
	
      siteOfEoSite[eoSite]=site;
      eoSiteOfSite[site]=eoSite;
    }
  
  /// Reset the number generator
  gen.resize(V);
  for(Site site=0;site<V;site++)
    gen[site].seed(inputSeed+site);
  
  /// Draw a conf
  for(Site eoSite=0;eoSite<V;eoSite++)
    conf[eoSite]=drawRndSpin(gen[eoSite]);
  
  /// Define site neighbors
  for(Site eoSite=0;eoSite<V;eoSite++)
    {
      const Site site=siteOfEoSite[eoSite];
      const Coords coords=coordsOfSite(site);
      const Coord x=coords[0],y=coords[1];
	
	neighs[eoSite][0]=eoSiteOfSite[siteOfCoords((x+1)%L,y)];
	neighs[eoSite][1]=eoSiteOfSite[siteOfCoords(x,(y+1)%L)];
	neighs[eoSite][2]=eoSiteOfSite[siteOfCoords((x+L-1)%L,y)];
	neighs[eoSite][3]=eoSiteOfSite[siteOfCoords(x,(y+L-1)%L)];
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
SiteUpdRes acceptReject(const int dE,const Site eoSite,const Spin oldSpin)
{
  const double pAcc=getPacc(dE);
  
  const double r=drawUniformNumber(gen[eoSite]);
  const bool acc=(r<pAcc);
  
  //nAcc+=acc;
  
  if(not acc)
    {
      conf[eoSite]=oldSpin;
      return {0,0,0};
    }
  else
    return {1,conf[eoSite]-oldSpin,dE};
}

/////////////////////////////////////////////////////////////////

/// Updates a single site
SiteUpdRes updateSite(const Site eoSite)
{
  const Spin oldSpin=conf[eoSite];
  int oldEn;
  if(useLocalEnergyChange)
    oldEn=energyOfSite(eoSite);
  else
    oldEn=measureEnergy();
  
  conf[eoSite]=drawRndSpin(gen[eoSite]);
  
  int newEn;
  if(useLocalEnergyChange)
    newEn=energyOfSite(eoSite);
  else
    newEn=measureEnergy();
  
  return acceptReject(newEn-oldEn,eoSite,oldSpin);
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
	  // const Site site=siteoOfSiteo[par][redSite];
// 	  const Coord y=redSite/(L/2);
// 	  const Coord x=2*(redSite%(L/2))+(par^(y%2));

      // #pragma omp parallel for reduction(+:cachedEnergy,cachedSpinSum,nAcc)
      // for(Site site=0;site<V;site++)
      // 	if(parOfSite[site]){
	  // const Coords c=coordsOfSite(site);
	  // if((c[0]+c[1])%2==par)
	    {
	  
	      const SiteUpdRes res=updateSite(redSite);
	  nAcc+=res.acc;
	  if(useMeasurementCache)
	    {
	      cachedEnergy+=res.dE;
	      cachedSpinSum+=res.dM;
	    }
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
    if(conf[siteOfEoSite[site]]==+1)
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
