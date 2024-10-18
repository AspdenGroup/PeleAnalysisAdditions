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
      AmrData::SetVerbose(false);

    std::string plotFileName; pp.get("infile",plotFileName);
    std::string fuelName1="H2";  pp.get("fuelName1",fuelName1);
    std::string fuelName2="CH4"; pp.get("fuelName2",fuelName2);
    Real chiX(1.0); pp.get("chiX",chiX); // blend fraction chi=1 => 100% fuel1 (based on mole-fraction)
    Real M1(2.0);  pp.get("M1",M1); // molar mass fuel 1
    Real M2(16.0); pp.get("M2",M2); // molar mass fuel 2
    Real chiY=M1*chiX/(M1*chiX+M2*(1.-chiX)); // convert X-based chi to Y-based chi using molar masses
    
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    bool initNormalise=false;

    DataServices dataServicesInit("plt00000", fileType);
    if (!dataServicesInit.AmrDataOk()) {
      std::cout << "Cannot find initial condition - normalising using current plotfile" << std::endl;
    } else {
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Normalising using initial condition" << std::endl;
      initNormalise=true;
    }

    AmrData& initAmrData = dataServicesInit.AmrDataRef();
    
    //init_mech();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    int idTin = -1;
    int idY1in = -1;
    int idY2in = -1;
    int idI1in = -1;
    int idI2in = -1;
    //Vector<std::string> spec_names = GetSpecNames();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string TName   = "temp";
    const std::string spName1 = "Y("+fuelName1+")";
    const std::string spName2 = "Y("+fuelName2+")";
    const std::string prName1 = "I_R("+fuelName1+")";
    const std::string prName2 = "I_R("+fuelName2+")";
    //std::cout << spName << std::endl;
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == TName)   idTin  = i;
      if (plotVarNames[i] == spName1) idY1in = i;
      if (plotVarNames[i] == spName2) idY2in = i;
      if (plotVarNames[i] == prName1) idI1in = i;
      if (plotVarNames[i] == prName2) idI2in = i;
    }
    if (idY1in<0 || idY2in<0 || idTin<0 || idI1in<0 || idI2in<0) {
      Print() << "Cannot find required data in pltfile" << std::endl;
      Print() << "idTin  = " << idTin << std::endl;
      Print() << "idY1in = " << idY1in << std::endl;
      Print() << "idY2in = " << idY2in << std::endl;
      Print() << "idI1in = " << idI1in << std::endl;
      Print() << "idI2in = " << idI2in << std::endl;
      Abort();
    }
    const int nCompIn  = 5;
    //const int idPhiout = 0;
    const int nCompOut = 5;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    const int idTlocal  = 0; // T  in/out
    const int idY1local = 1; // Y1 in/out
    const int idY2local = 2; // Y2 in/out
    const int idI1local = 3; // I1 in
    const int idI2local = 4; // I2 in

    const int idBlocal  = 3; // blend out
    const int idIBlocal = 4; // IB out

    destFillComps[idTlocal]  = idTlocal;
    destFillComps[idY1local] = idY1local;
    destFillComps[idY2local] = idY2local;
    destFillComps[idI1local] = idI1local;
    destFillComps[idI2local] = idI2local;
    
    inNames[idTlocal]   = "temp";
    inNames[idY1local]  = "Y("+fuelName1+")";
    inNames[idY2local]  = "Y("+fuelName2+")";
    inNames[idI1local]  = "I_R("+fuelName1+")";
    inNames[idI2local]  = "I_R("+fuelName2+")";

    outNames[idTlocal]  = "prog_temp";
    outNames[idY1local] = "prog_"+fuelName1;
    outNames[idY2local] = "prog_"+fuelName2;
    outNames[idBlocal]  = "prog_blend";
    outNames[idIBlocal] = "I_R(blend)";

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    
    Vector<std::unique_ptr<MultiFab>> outdata(Nlev);
    const int nGrow = 0;
    int b[3] = {1, 1, 1};
    Real Y_fuel1_u, Y_fuel1_b, Y_fuel2_u, Y_fuel2_b, T_u, T_b;
    if(initNormalise) {
      if(initAmrData.MinMax(initAmrData.ProbDomain()[0],"Y("+fuelName1+")",0,Y_fuel1_b,Y_fuel1_u) &&
	 initAmrData.MinMax(initAmrData.ProbDomain()[0],"Y("+fuelName2+")",0,Y_fuel2_b,Y_fuel2_u) &&
	 initAmrData.MinMax(initAmrData.ProbDomain()[0],"temp",0,T_u,T_b)) {
	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Found min/max" << std::endl;
      } else {
	std::cout << "Could not find min/max" << std::endl;
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
      }
    } else {
       if(amrData.MinMax(amrData.ProbDomain()[0],"Y("+fuelName1+")",0,Y_fuel1_b,Y_fuel1_u) &&
	  amrData.MinMax(amrData.ProbDomain()[0],"Y("+fuelName2+")",0,Y_fuel2_b,Y_fuel2_u) &&
	  amrData.MinMax(amrData.ProbDomain()[0],"temp",0,T_u,T_b)){
	 if (ParallelDescriptor::IOProcessor())
	   std::cout << "Found min/max" << std::endl;
       } else {
	 std::cout << "Could not find min/max" << std::endl;
	 DataServices::Dispatch(DataServices::ExitRequest, NULL);
       }
    }
   

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
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	
        const Box& bx = mfi.tilebox();
	Array4<Real> const& varIn  = indata.array(mfi);
        Array4<Real> const& varOut = (*outdata[lev]).array(mfi);

        AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
        {
	  varOut(i,j,k,idTlocal)  = (varIn(i,j,k,idTlocal)-T_u)/(T_b-T_u);
          varOut(i,j,k,idY1local) = 1-varIn(i,j,k,idY1local)/Y_fuel1_u;
          varOut(i,j,k,idY2local) = 1-varIn(i,j,k,idY2local)/Y_fuel2_u;
	  varOut(i,j,k,idBlocal)  = chiY*varOut(i,j,k,idY1local) + (1-chiY)*varOut(i,j,k,idY2local);
	  varOut(i,j,k,idIBlocal) = chiY*varIn(i,j,k,idI1local)  + (1-chiY)*varIn(i,j,k,idI2local);
        });
      }

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_prog");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
