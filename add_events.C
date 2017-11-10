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

void add_events(int run=1, int cam = 0, const char * dir = "/scratch3/wparker2/dmtpc2/data/2017/08/raw")
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

  TH2 *image_sum0, *image_sum1,*image_sum2, *image_sum3;
  TH2 *image0, *image1, *image2, *image3;
  

  ///////////Load up the data
  TString fout=Form("/scratch3/wparker2/dmtpc2/data/2017/hptpc_test_R%05d.add_events.root",run);

  TFile * foutFile = new TFile(fout,"RECREATE");

  TTree t("pix_t","Pixel Tree");
  int nbins = 1000;  
  
  TH1F* TotalIntensity0=new TH1F("TotalIntensity0","Total Intensity0",nbins,-100,100);
  TH1F* TotalIntensity1=new TH1F("TotalIntensity1","Total Intensity1",nbins,-100,100);
  TH1F* TotalIntensity2=new TH1F("TotalIntensity2","Total Intensity2",nbins,-100,100);
  TH1F* TotalIntensity3=new TH1F("TotalIntensity3","Total Intensity3",nbins,-100,100);

  double rms = 0, totalrms=0;
  double mean = 0, totalmean=0;
  double chi = 0;
  int runnum = 0;
  double ccdTemp;
  int nevents = 0;

  dmtpc::core::Dataset *d = new dmtpc::core::Dataset;
  d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run)); 

  int n = d->nevents();
  runnum = d->event()->run();
  int ncam= d->getHeader()->ncameras;

  if (cam<0 || cam>=ncam)
    Fatal("add_events","Wrong camera ID=%d out of allowed range=[0,%d]",cam,ncam-1);

  //t.Branch("intensity",&intensity);
  t.Branch("mean",&mean);
  t.Branch("rms",&rms);
  // t.Branch("runnum",&runnum);
  t.Branch("nevents",&nevents);
  t.Branch("chi",&chi);
  t.Branch("ccdTemp",&ccdTemp);

  ///////// Initially asign sum  Histo to one image to get correct dimensions, then set to 0
  d->getEvent(0);
  image_sum0 = (TH2D*) d->event()->getImage(0)->getVisible()->Clone();
  image_sum1 = (TH2D*) d->event()->getImage(1)->getVisible()->Clone();
  image_sum2 = (TH2D*) d->event()->getImage(2)->getVisible()->Clone();
  image_sum3 = (TH2D*) d->event()->getImage(3)->getVisible()->Clone();
  for (int m=1; m<=image_sum0->GetNbinsX(); m++){
    for (int j=1; j<=image_sum0->GetNbinsY(); j++) {

      image_sum0->SetBinContent(m,j,0);
      image_sum1->SetBinContent(m,j,0);
      image_sum2->SetBinContent(m,j,0);
      image_sum3->SetBinContent(m,j,0);

    }
  }

  Info("add_events","Process %d frames of run #%d: %d active camera(s)",n,runnum,ncam);
  for (int i = 0; i < n; i++) {
    
    if (++ncur>portion) {
      Info("add_events","Reopen file at event=%d",i);

      d->close();
      delete d;

      d = new dmtpc::core::Dataset;
      d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run)); 
      ncur=0;
    }

    Info("add_events","Process event %d",i);
    d->getEvent(i);
    ccdTemp = d->event()->ccdConfig(0)->ccdTemp;
    
    if (dbg)
      Info("add_events","  Cam %d: temperature: %g degrees",0,ccdTemp);

    image0 = (TH2*) d->event()->getImage(0)->getVisible()->Clone();
    image1 = (TH2*) d->event()->getImage(1)->getVisible()->Clone();
    image2 = (TH2*) d->event()->getImage(2)->getVisible()->Clone();
    image3 = (TH2*) d->event()->getImage(3)->getVisible()->Clone();
    
    if (dbg)
      Info("add_events","images are ok");


    for (int m=1; m<=image0->GetNbinsX(); m++){
      for (int j=1; j<=image0->GetNbinsY(); j++) {
	
	float pixval0=image0->GetBinContent(m,j);
	float temp0=image_sum0->GetBinContent(m,j);
	image_sum0->SetBinContent(m,j, pixval0+temp0);
	
	//	std::cout << "image0: " << pixval0 << " temp0: " << temp0 << " image_sum0: "<< image_sum0->GetBinContent(m,j) << std::endl;

	float pixval1=image1->GetBinContent(m,j);
        float temp1=image_sum1->GetBinContent(m,j);
        image_sum1->SetBinContent(m,j, pixval1+temp1);

	float pixval2=image2->GetBinContent(m,j);
        float temp2=image_sum2->GetBinContent(m,j);
        image_sum2->SetBinContent(m,j, pixval2+temp2);

	float pixval3=image3->GetBinContent(m,j);
        float temp3=image_sum3->GetBinContent(m,j);
        image_sum3->SetBinContent(m,j, pixval3+temp3);
	
      }
    }

    if (dbg)	
      Info("add_events","intensity histo is ok");

    if (dbg)
      Info("add_events","fit is ok");

    t.Fill();

    if (dbg)
      Info("add_events","tree is ok");


    delete image0;
    delete image1;
    delete image2;
    delete image3;

  } 

  Info("add_events","end of cycle");
   
  t.Write();
  
  for (int m=1; m<=image_sum0->GetNbinsX(); m++){
    for (int j=1; j<=image_sum0->GetNbinsY(); j++) {
  
      float pixval0b = image_sum0->GetBinContent(m,j);
      TotalIntensity0->Fill(pixval0b/n);

      float pixval1b = image_sum0->GetBinContent(m,j);
      TotalIntensity1->Fill(pixval1b/n);

      float pixval2b = image_sum0->GetBinContent(m,j);
      TotalIntensity2->Fill(pixval2b/n);

      float pixval3b = image_sum0->GetBinContent(m,j);
      TotalIntensity3->Fill(pixval3b/n);

    }
  }
  TF1 * gFit0 = new TF1("gFit0","gaus");
  TotalIntensity0->Fit(gFit0,"Q");
  TotalIntensity0->SetStats(kTRUE);
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(56);
  TotalIntensity0->SetXTitle("ADU");
  TotalIntensity0->SetYTitle("Number of Entries");
  image_sum0->Write();

TF1 * gFit1 = new TF1("gFit1","gaus");
TotalIntensity1->Fit(gFit1,"Q");
TotalIntensity1->SetStats(kTRUE);
gStyle->SetOptStat(1111);
gStyle->SetPalette(56);
TotalIntensity1->SetXTitle("ADU");
TotalIntensity1->SetYTitle("Number of Entries");
image_sum1->Write();

TF1 * gFit2 = new TF1("gFit2","gaus");
TotalIntensity2->Fit(gFit2,"Q");
TotalIntensity2->SetStats(kTRUE);
gStyle->SetOptStat(1111);
gStyle->SetPalette(56);
TotalIntensity2->SetXTitle("ADU");
TotalIntensity2->SetYTitle("Number of Entries");
image_sum2->Write();

TF1 * gFit3 = new TF1("gFit3","gaus");
TotalIntensity3->Fit(gFit3,"Q");
TotalIntensity3->SetStats(kTRUE);
gStyle->SetOptStat(1111);
gStyle->SetPalette(56);
TotalIntensity3->SetXTitle("ADU");
TotalIntensity3->SetYTitle("Number of Entries");
image_sum3->Write();


  //TCanvas *c2 = new TCanvas("c2","",0,0,600,600);
  //c2->cd();
  //intensity->Draw();
  //intensity->Write();
  //c2->Update();
  // c2->SetLogy();
  
  foutFile->Write();
  foutFile->Close();

  Info("add_events","Output data has been saved in file='%s'",fout.Data());
}


