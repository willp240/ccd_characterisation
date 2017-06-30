#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TColor.h>

#include <CameraInfo.hh>
#include <Event.hh>
#include <DatasetHeader.hh>
#include <Dataset.hh>

void calc_rms(int run=1, int cam = 0, const char * dir = "/scratch3/wparker2/dmtpc2/data/2017/06/raw/")
{

  bool dbg=true;

  int portion=30, ncur=0;

  double stops[5] = {0,0.34,0.61,0.84,1.0};
  gStyle->SetOptStat(0); 
  double red[5] = {0.0,0.0,0.87,1.0,0.51};
  double green[5] = {0.0,0.81,1.00,0.2,0.0};
  double blue[5] = {0.51,1.0,0.12,0.0,0.0};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TH2 *image, *image_bgr;
  
  TString fout=Form("%s/hptpc_test_R%05d.calc_rms.root",dir,run);
  TFile * foutFile = new TFile(fout,"RECREATE");

  TTree t("pix_t","Pixel Tree");
  int nbins = 1000;  
  TH1F* intensity=new TH1F("intensity","intensity",nbins,-100,100);

  double rms = 0;
  double mean = 0;
  double chi = 0;
  int runnum = 0;
  double ccdTemp;
  int nevents = 0;

  dmtpc::core::Dataset *d = new dmtpc::core::Dataset;
  d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run)); 

  int n = d->nevents();
  runnum = d->event()->run();
  std::cout << n << std::endl;
  int ncam= d->getHeader()->ncameras;

  if (cam<0 || cam>=ncam)
    Fatal("calc_rms","Wrong camera ID=%d out of allowed range=[0,%d]",cam,ncam-1);

  //t.Branch("intensity",&intensity);
  t.Branch("mean",&mean);
  t.Branch("rms",&rms);
  // t.Branch("runnum",&runnum);
  t.Branch("nevents",&nevents);
  t.Branch("chi",&chi);
  t.Branch("ccdTemp",&ccdTemp);

  Info("calc_rms","Process %d frames of run #%d: %d active camera(s)",n,runnum,ncam);
  for (int i = 0; i < n; i++) {

    if (++ncur>portion) {
      Info("calc_rms","Reopen file at event=%d",i);

      d->close();
      delete d;

      d = new dmtpc::core::Dataset;
      d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run)); 
      ncur=0;
    }

    Info("calc_rms","Process event %d",i);
    d->getEvent(i);
    ccdTemp = d->event()->ccdConfig(cam)->ccdTemp;
    
    if (dbg)
      Info("calc_rms","  Cam %d: temperature: %g degrees",cam,ccdTemp);

    image = (TH2*) d->event()->getImage(cam)->getVisible()->Clone();
    image_bgr = (TH2*) d->biasAvg(cam)->Clone();
    image->Add(image_bgr,-1);

    if (dbg)
      Info("calc_rms","images are ok");

    intensity->Reset(); // YS: should we reset?

    for (int m=0; m<image->GetNbinsX(); m++){
      for (int j=0; j<image->GetNbinsY(); j++) {
	float pixval=image->GetBinContent(m,j);
		  
	intensity->Fill(pixval);
      }
    }

    if (dbg)	
      Info("calc_rms","intensity histo is ok");

    TF1 * gFit = new TF1("gFit","gaus"); 
    intensity->Fit(gFit,"Q");
    rms = intensity->GetRMS(); //gFit->GetParameter(2);
    mean = gFit->GetParameter(1);
    chi = gFit->GetChisquare();

    if (dbg)
      Info("calc_rms","fit is ok");

    t.Fill();

    if (dbg)
      Info("calc_rms","tree is ok");


    delete image, image_bgr;
    delete gFit;

  } 
  Info("calc_rms","end of cycle");
  
  intensity->SetStats(kTRUE);
  gStyle->SetOptStat(1111);
  intensity->SetXTitle("ADU"); 
  intensity->SetYTitle("Number of Entries");
  
  t.Write();

  //TCanvas *c2 = new TCanvas("c2","",0,0,600,600);
  //c2->cd();
  //intensity->Draw();
  //intensity->Write();
  //c2->Update();
  // c2->SetLogy();
  
  foutFile->Write();
  foutFile->Close();

  Info("calc_rms","Output data has been saved in file='%s'",fout.Data());
}


