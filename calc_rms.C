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

// returns [xmi,xma] part of histogramm, xmi=0 means first bin
// xma=0 means last bin
TH1F* hist_part(TH1F* h, Double_t xmi=0, Double_t xma=0, TString hsuff = "_part") {

  if (!h)
	 Fatal("hist_part","No histogram given!");
  
  //  cout << "In hist_part\n";
  h->GetXaxis()->UnZoom();
  Double_t hmin = h->GetXaxis()->GetXmin();
  Double_t hmax = h->GetXaxis()->GetXmax();
  Double_t bw =   h->GetXaxis()->GetBinWidth(2);
  //  Double_t ne = h->GetEntries();
  Double_t r0 = h->Integral();
  //  Double_t emean = h->GetMean();
  TString title= h->GetTitle(), name=h->GetName();

  //  cout << "Histo: name='" << name << "', title='" <<
  //	 title << "', xmin=" << hmin << ", xmax=" << hmax << 
  //	 "\n  Integral=" << r0 << ", Entries=" << ne << 
  //	 ", mean energy=" << emean << endl;
  
  if (r0<=0) {
	 Error("hist_part","Empty histo!");
	 return 0;
  }

  Double_t xmin = ( (xmi>hmin) && (xmi<hmax) ) ? xmi : hmin;
  Double_t xmax = ( (xma>xmin) && (xma<hmax) ) ? xma : hmax;
  
  Int_t binMin = h->FindBin(xmin);
  Int_t binMax = h->FindBin(xmax);

  //  Double_t r = h->Integral(binMin,binMax);
  //  Double_t rr = r*100./ne;  
  //  cout << "Integral(" << xmin << "," << xmax << ")=" << 
  //	 r << " (" << rr << " % of total)\n";
  
  Int_t nbins = binMax-binMin+1;
  if (nbins<2) {
	 Error("hist_part","Bad range in histo='%s': <2 bin in [xmin, xmax]=[%g,%g]",
			 name.Data(),xmin,xmax);
	 return 0;
  }
  TString h2name = h->GetName(); h2name += hsuff;
  TString h2title = h->GetTitle(); h2title += Form(": part[%d,%d]",binMin,binMax);
  
  // cout << "Created hist='"<< h2name << \n"', title '"<< h2title << \n"', number of bins=" <<
  // 	 nbins << ", range=[" << xmin-bw/2 << "," << xmax+bw/2 << "]\n";
  
  TH1F *h2 = new TH1F(h2name,h2title,nbins,xmin-bw/2,xmax+bw/2);
  h2->SetDirectory(0);

  Int_t j=0;
  for (Int_t i=binMin; i<=binMax; i++)
	 h2->SetBinContent(++j,h->GetBinContent(i));
  
  return h2;
} // hist_part

void calc_rms(int run=1, int cam = 0, const char * dir = "/scratch3/wparker2/dmtpc2/data/2017/10/raw")
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

  TH2 *image, *image_bgr, *image_ave;
  
  //  TString fout=Form("%s/hptpc_test_R%05d.calc_rms.root",dir,run);
  TString fout=Form("/scratch3/wparker2/dmtpc2/data/2017/10/hptpc_test_R%05d_cam%d.calc_rms2.root",run,cam);

  TFile * foutFile = new TFile(fout,"RECREATE");

  TTree t("pix_t","Pixel Tree");
  int nbins = 1000;  
  TH1F* intensity=new TH1F("intensity","intensity",nbins,-100,100);
  TH1F* TotalIntensity=new TH1F("TotalIntensity","Total Intensity",nbins,-100,100);

  double rms = 0, totalrms=0;
  double mean = 0, totalmean=0, sig;
  double chi = 0;
  int runnum = 0;
  double ccdTemp;
  int nevents = 0;

  dmtpc::core::Dataset *d = new dmtpc::core::Dataset;
  d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run)); 

  int n = d->nevents();
  runnum = d->event()->run();
  std::cout << "Number of events: " << n << std::endl;
  int ncam= d->getHeader()->ncameras;
  std::cout << "Number of cameras: " << ncam << std::endl;

  if (cam<0 || cam>=ncam)
    Fatal("calc_rms","Wrong camera ID=%d out of allowed range=[0,%d]",cam,ncam-1);

  //t.Branch("intensity",&intensity);
  t.Branch("mean",&mean);
  t.Branch("rms",&rms);
  // t.Branch("runnum",&runnum);
  t.Branch("nevents",&nevents);
  t.Branch("chi",&chi);
  t.Branch("ccdTemp",&ccdTemp);

  ///////// Initially asign average Histo to one image to get correct dimensions, then set to 0
  d->getEvent(0);
  image_ave = (TH2*) d->event()->getImage(cam)->getVisible()->Clone();
  for (int m=1; m<=image_ave->GetNbinsX(); m++){
    for (int j=1; j<=image_ave->GetNbinsY(); j++) {

      image_ave->SetBinContent(m,j,0);

    }
  }


  Info("calc_rms","Process %d frames of run #%d: %d active camera(s)",n,runnum,ncam);
  for (int i = 0; i < n; i++) {
    
    if (++ncur>portion) {//Had crashes when ran over too many images at once so close file and reopen after #portion events
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
    std::cout << image->GetNbinsX() << " " << image->GetNbinsY() << std::endl;
    if (dbg)
      Info("calc_rms","images are ok");

    intensity->Reset(); // YS: should we reset?

    for (int m=1; m<=image->GetNbinsX(); m++){
      for (int j=1; j<=image->GetNbinsY(); j++) {
	float pixval=image->GetBinContent(m,j);
	
	float temp = image_ave->GetBinContent(m,j);
	image_ave->SetBinContent(m,j, pixval+temp);
	
	intensity->Fill(pixval);
      }
    }

    if (dbg)	
      Info("calc_rms","intensity histo is ok");

    TF1 * gFit = new TF1("gFit","gaus"); 
    intensity->Fit(gFit,"Q");

    mean = gFit->GetParameter(1);
    sig  = gFit->GetParameter(2);
    chi = gFit->GetChisquare();

    double nsig=4;
    double xmin=mean-nsig*sig,
      xmax=mean+nsig*sig;

    //TH1F *intensity_short=hist_part(intensity,xmin,xmax);

    //    intensity_short->Draw();
    //return;

        rms = intensity->GetRMS(); //gFit->GetParameter(2);
	//rms = intensity_short->GetRMS(); //gFit->GetParameter(2);

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

  for (int m=1; m<=image_ave->GetNbinsX(); m++){
    for (int j=1; j<=image_ave->GetNbinsY(); j++) {
      float temp=image_ave->GetBinContent(m,j);
      image_ave->SetBinContent(m,j, temp/n);
      float pixval2 = image_ave->GetBinContent(m,j);
      
      TotalIntensity->Fill(pixval2);
    }
  }
  TF1 * gFit2 = new TF1("gFit2","gaus");
  TotalIntensity->Fit(gFit2,"Q");
  TotalIntensity->SetStats(kTRUE);
  gStyle->SetOptStat(1111);
  TotalIntensity->SetXTitle("ADU");
  TotalIntensity->SetYTitle("Number of Entries");
  image_ave->Write();
  
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


