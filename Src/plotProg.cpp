/*
  --- plotProg.cpp ---

  A tool calculates the progress variable of a given number of variables from a pltfile.
  The tool can also normalise the a given variables by a specific pltfile.

*/



#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=plt????? vars= ?? [options] \n";
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

    // pass outfile name
    std::string outfile(getFileRoot(plotFileName) + "_prog");
    pp.query("outfile",outfile);

    // get progress variable names
    int nComp(pp.countval("vars"));
    if (nComp <= 0)  Abort("No progress varbiable names given");
    Vector<std::string> inVarNames;
    inVarNames.resize(nComp);
    pp.getarr("vars",inVarNames);
    
    // normalise by specific plotfile
    std::string normalisePltName; pp.query("normalisePltName",normalisePltName);
    bool normalisePlt=1;
    if (normalisePltName.empty())
	normalisePlt=0;
    
    Print() << "progress variable names: ";
    for (int comp = 0; comp < nComp; comp++)
	Print() << inVarNames[comp] << " ";
    Print() << "\n";

    // output data names
    Vector<std::string> outVarNames = inVarNames;
    for (int comp = 0; comp < nComp; comp++)
	outVarNames[comp] = inVarNames[comp]+"_prog";
    
    // data used to normalise
    DataServices normaliseDataServices(normalisePltName, fileType);
    if( ! normaliseDataServices.AmrDataOk() && 1 == normalisePlt) {
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& normaliseAmrData = normaliseDataServices.AmrDataRef(); // this does send a runtime warning error if not normalising

    // data to process
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }

    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);

    const int nGrow = 0;

    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i) destFillComps[i] = i;

    // find min/max to be normalised
    Vector<Real> prog_burnt(nComp);
    Vector<Real> prog_unburnt(nComp);

    Print() << "nComp = " << nComp << std::endl;
    for (int comp=0; comp < nComp; comp++) {
	if (normalisePlt) {
	    normaliseAmrData.MinMax(normaliseAmrData.ProbDomain()[0],inVarNames[comp],0,prog_burnt[comp], prog_unburnt[comp]);    
	}
	else {
	    amrData.MinMax(amrData.ProbDomain()[0],inVarNames[comp],0,prog_burnt[comp], prog_unburnt[comp]);
	}
    }
   
    for (int lev = 0; lev < Nlev; ++lev) {      
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba,dm,nComp,nGrow);
      MultiFab indata(ba,dm,nComp,nGrow);

      int coord = 0;
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb, coord, &(is_per[0]));            
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inVarNames,destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;
      
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  auto const& out_a = outdata[lev].array(mfi);
	  auto const& in_a = indata.array(mfi);
	  amrex::ParallelFor(bx, [=]
			     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	      {
		  for (int n = 0; n < nComp; n++) {
		      out_a(i,j,k,n) = (in_a(i,j,k,n) - prog_unburnt[n]) / (prog_burnt[n] - prog_unburnt[n]);
		  }
	      });
	}
    }
    
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outVarNames,
                                   geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
