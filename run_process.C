// script runs DANSS analysis job
// the main task is to load proper objects (libs & classes) used by it
// before start of the main job.

// loads libs required for application run
void load_libs_and_classes()
{

  {
	 // LIBS

	 // loaded libs, MUST finished with empty terminating string
	 //	 TString lnames[] = {"config", "bz2", ""};
	 TString lnames[] = {"DmtpcCore", "DmtpcAnalysis", ""};

	 // status of loaded libs (desribed in TSystem::Load method doc):
	 TString load_status[] = {"ok", "already loaded", "failed to load", "version_mismatch"};
  
	 int i=-1;
	 while (lnames[++i] != "" ) {
		if (!i) cout << "Load libs:\n";
		TString lname=(lnames[i].BeginsWith("lib")) ? lnames[i] : Form("lib%s",lnames[i].Data());
		cout << Form("  %s... ",lname.Data());
		int ierr=gSystem->Load(lname.Data());
		if (ierr<0) ierr=1-ierr; // convert negative code to the positive number = indexes of load_status[]
		if (ierr>3)
		  Fatal("load_libs_and_classes",Form("Unknown error code=%d! Check doc of probably updated TSystem::Load()!",-ierr+1));
		
		cout << load_status[ierr] << "\n";
		if (ierr<2) continue;

		Fatal("load_libs_and_classes","Please, resolve lib loading problem and try again!");
	 }
  }

  {
	 // CLASSES

	 // loaded classes, MUST finished with empty terminating string
	 // cnames[j].C files expected for class implepentation code
	 TString cmacro,cnames[] = {""};

	 int i=-1;
	 while (cnames[++i] != "" ) {
		if (!i) cout << "Load classes:\n";
		cmacro.Form("%s.C+",cnames[i].Data());
		cout << Form("  %s... ",cmacro.Data());
		if (gROOT->LoadMacro(cmacro))
		  Fatal("load_libs_and_classes","Error to load macro='%s'",cmacro.Data());
		cout << "ok\n";
	 }
  }  

} // load_libs_and_classes

// set of macro command pars MUST BE IDENTICAL to pars of $mac script
// belo as they are all passed to it w/o processing
void run_process
(
 int run =1
 ,int cam = 0
 , const char * dir = "/scratch3/wparker2/dmtpc2/data/2017/10/raw/"
 )
{  
  Info("run_process","Started");
  gSystem->AddIncludePath("-I$DMTPC_HOME/DmtpcCore/include");
  load_libs_and_classes();

  TString macro_file="process.C+g";
  Info("run_process","Loading main macro='%s'... ",macro_file.Data());
  
  if (gROOT->LoadMacro(macro_file))
    Fatal("run_process","Can't load macro='%s'",macro_file.Data());

  // MAIN JOB 
  Info("run_process","It seems like ok, starting main job now... ");

  process(run,cam,dir);
  
  Info("run_process","Finished");
} // run_process
