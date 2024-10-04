/*
  A tool calculates the equiverlance ratio given a pltfile
  (currently only works for Hydrogen or Carbon Fuels)
*/



#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_GpuLaunch.H>

#include <PelePhysics.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    std::string outfile(getFileRoot(plotFileName) + "_phi");
    pp.query("outfile",outfile);


    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idYin = -1;
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    auto eos = pele::physics::PhysicsType::eos();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
    }
    if (idYin<0)
      Print() << "Cannot find required data in pltfile" << std::endl;

    const int nCompIn = NUM_SPECIES;
    const int nCompOut = 1;

    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Xs start here
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    outNames[0] = "phi";

    // get number of each element in each species
    int ncf[NUM_SPECIES * NUM_ELEMENTS];
    CKNCF(ncf);
    // element names
    amrex::Vector<std::string> ename;
    CKSYME_STR(ename);
    // id for each element
    int Hid(-1), Oid(-1), Nid(-1), Cid(-1);
    for (int n=0; n<NUM_ELEMENTS; n++) {
      if (ename[n]=="H") Hid=n;
      if (ename[n]=="O") Oid=n;
      if (ename[n]=="N") Nid=n;
      if (ename[n]=="C") Cid=n;
    }

    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    const int nGrow = 0;

    for (int lev=0; lev<Nlev; ++lev)
    {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      outdata[lev].reset(new MultiFab(ba,dm,nCompOut,nGrow));
      MultiFab indata(ba,dm,nCompIn,nGrow);
      int coord = 0;
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb, coord, &(is_per[0]));
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	  const Box& bx = mfi.tilebox();
	  Array4<Real> const& Y = indata.array(mfi);
	  Array4<Real> const& phi = (*outdata[lev]).array(mfi);

	  AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
          {
	      // convert mass fraction to mole fraction
	      Real Yl[NUM_SPECIES];
	      Real Xl[NUM_SPECIES];
	      for (int n=0; n<NUM_SPECIES; ++n) {
		  Yl[n] = Y(i,j,k,idYlocal+n);
	      }
	      eos.Y2X(Yl,Xl);
	      
	      // sum atoms
	      Real sumH(0.);
	      Real sumO(0.);
	      Real sumC(0.);
	      for (int n=0; n<NUM_SPECIES; ++n) {
		  sumH += Xl[n] * ncf[n*NUM_ELEMENTS+Hid];
		  sumO += Xl[n] * ncf[n*NUM_ELEMENTS+Oid];
		  sumC += Xl[n] * ncf[n*NUM_ELEMENTS+Cid];
	      }

              // calculate equivalence ratio
	      phi(i,j,k,0) = (sumH+4.0*sumC)/(2.0*sumO);
          });
      }
      Print() << "Derive finished for level " << lev << "\n";
    }

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
