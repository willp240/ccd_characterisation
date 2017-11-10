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
#include <TGraphErrors.h>

#include <CameraInfo.hh>
#include <Event.hh>
#include <DatasetHeader.hh>
#include <Dataset.hh>


void save_pictures(TString fpref, TString cname="cv") {

  TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(cname);
  if (!c2)
    Fatal("save_pictures","Can't find canvas='%s'",cname.Data());
  //  c2->SetGrid();
  gStyle->SetPalette(51);
  //  TString fout,sext[]={"png", "pdf", "C", ""};
  TString fout,sext[]={"png", "C", ""};

  int j=-1; while (sext[++j] != "") {  // save
    fout=Form("%s.%s",fpref.Data(),sext[j].Data());
    c2->SaveAs(fout.Data());
    //    Info("save_pictures","saved file='%s'",fout.Data());
  }
} // save_pictures


void make_freqHistos(TH2 *hsub, TString pdir, int runNum, int camNum, int eventNum, double mean[][4],double rms[][4],bool save,TString hname){

  int xmin = hsub->GetXaxis()->GetFirst();
  int ymin = hsub->GetYaxis()->GetFirst();
  int xmax = hsub->GetXaxis()->GetLast();
  int ymax = hsub->GetYaxis()->GetLast();

  double mu, delta;

  vector<double> vals;
  double min = 1e99;
  double max = -1e99;
  for (int x = xmin; x < xmax; x++)
    {
      for (int y = ymin; y < ymax; y++)
        {
          double val = hsub->GetBinContent(x,y);

          if (val < min) min = val;
          if (val > max) max = val;
          vals.push_back(val);
        }
    }

  int nbins = int(max - min + 0.5);

  TH1F ph("ph","ph",nbins,min,max);

  for (unsigned i = 0; i < vals.size(); i++)
    {
      ph.Fill(vals[i]);
    }

  TCanvas *pixhist_c = new TCanvas("pixhist_c", "Pixel Histogram",600,600);

  pixhist_c->cd();
  ph.SetLineColor(1);
  ph.SetFillColor(1);
  //ph.Draw();
  //pixhist_c->Update();
  mu=ph.GetMean();
  delta=ph.GetRMS();

  TString  fout_mask=Form("%s/%08d/pixel_hist_r%08d_e%03d_c%d",pdir.Data(),runNum,runNum,eventNum,camNum);
  if(save){
    ph.Draw();
    pixhist_c->Update();
    save_pictures(fout_mask,"pixhist_c");
    ph.Write(hname);
  }
  mean[eventNum+1][camNum+1]=mu;
  rms[eventNum+1][camNum+1]=delta;

}//make_freqHistos


void process(int run=1, int cam = 0, const char * dir = "/scratch3/wparker2/dmtpc2/data/2017/09/raw")
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
  

  ///////////Load up the data hard coded
  TString fout=Form("/scratch3/wparker2/dmtpc2/data/2017/hptpc_test_R%05d.process.root",run);
 
  //  TFile * foutFile = new TFile(fout,"RECREATE");

  //Hardcode where save images
  TString pdir="pict_oct", fout_mask;
  gSystem->mkdir(pdir);
  //name of directory where to save
  TString idir="/home/wparker2/s3/dmtpc2/data/2017/10"
    ,fmask="hptpc_test"
    ,hname
    ;
  TString gname = "";
  //open data
  dmtpc::core::Dataset *d = new dmtpc::core::Dataset;
  d->open(TString::Format("%s/hptpc_test_R%05d.raw.root",dir,run));
 
  int n = d->nevents();
  int rnum = d->event()->run();
  int ncam= d->getHeader()->ncameras;

  TString rdir = Form("%s/%08d/",pdir.Data(),rnum);
  gSystem->mkdir(rdir);

  double mean[n][4];
  double rms[n][4];
  double ccdtemp[n][4];
  double Nev[n];
  double nevEr[n];
  bool save_pict = false;
  bool make_freqHist = true;
  //initialise arrays  
for (int i=0; i<n; i++)
    {
      for ( int j=0; j <ncam ; j++)
        {
	  mean[i][j]=0;
	  rms[i][j]=0;
	  Nev[i]=i;
	  nevEr[i]=0;
        }
    }

  double norm_cf = (n>0) ? 1./n : 1;

  if (cam<0 || cam>=ncam)
    Fatal("process","Wrong camera ID=%d out of allowed range=[0,%d]",cam,ncam-1);

  TH2D* hsub_all[4]={0,0,0,0};
  TH2D* hraw_all[4]={0,0,0,0};

  Info("process","Process %d frames of run #%d: %d active camera(s)",n,rnum,ncam);
  for (int i = 0; i < n; i++) {
    //over all events
    d->getEvent(i);
    for (int j=0; j<ncam; j++) {
      //over all cams
      ccdtemp[i][j] = d->event()->ccdConfig(j)->ccdTemp;
    
      if (dbg)
	Info("process","  Cam %d: temperature: %g degrees Event %d",j,ccdtemp[i][j],i);
      //Get image and bias
      TH2D * hraw =     (TH2D*) d->event()->ccdData(j)->Clone();
      TH2D * hbias_avg= (TH2D*) d->biasAvg(j)->Clone();

      if (!hraw) {
        Warning("my_first","Ev %d, camera %d: can't take camera image, skip it!",i,j);
        continue;
      }
      //subtract bias
      hname = hraw->GetName();
      hname += "_sub";
      TH2D * hsub =     (TH2D*) hraw->Clone(hname);
      hsub->Add(hbias_avg,-1);
      //add to total sum of all bias subtractr
      if (!hsub_all[j]) {
	hname += "_all";
	hsub_all[j] = (TH2D*) hsub->Clone(hname);
	hsub_all[j]->Reset();
      }
      else
	hsub_all[j]->Add(hsub,norm_cf);
      //add to total sum of all raw images
      if (!hraw_all[j]) {
	hname = hraw->GetName();
	hname += "_raw_all";
	hraw_all[j] = (TH2D*) hraw->Clone(hname);
	hraw_all[j]->Reset();
      }
      else
	hraw_all[j]->Add(hraw,norm_cf);

      //hsub->GetZaxis()->SetRangeUser(-200,250);
      gStyle->SetPalette(51);
      //hraw->Draw("colz");
      //gPad->Update();

      if (save_pict) {
	fout_mask=Form("%s/%08d/hraw_r%08d_e%03d_c%d",pdir.Data(),rnum,rnum,i,j);
	save_pictures(fout_mask,"cmy");
      }

      if(make_freqHist)
	{
	  make_freqHistos(hsub,pdir.Data(),rnum,j,i,mean,rms,0,"");
	}

      delete hsub;
      delete hraw;
      delete hbias_avg;
    }//Matches for(int j=0; j<ncam; j++)
    
  } // Matches   for(int i=0; i < n; i++) {

  if (dbg)
    Info("process","tree is ok");
   
  double meanVec[n];
  double rmsVec[n];
  double ccdtempvec[n];
  double aveMean[2][4];
  double aveRMS[2][4];
  for (int i=0; i<ncam; i++)
    {
      for (int j=0; j<2; j++)
        {
	  aveMean[j][i]=0;
	  aveRMS[j][i]=0;
        }
    }
  TFile * foutFile = new TFile(fout,"RECREATE");
  TCanvas *cmy;
  cmy = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cmy");
  if (!cmy)
    cmy = new TCanvas("cmy","cmy",10,10,600,600); 

  //Make plots of raw and bias subtracted pixel intensity histogram, and draw raw and bias subbed images 
  for (int j=0; j<ncam; j++) {
    hname = Form("Intensity_sub_%d",j);
    make_freqHistos(hsub_all[j],pdir.Data(),rnum,j,0,aveMean,aveRMS,1,hname);
    hname = Form("Intensity_raw_%d",j);
    make_freqHistos(hraw_all[j],pdir.Data(),rnum,j,1,aveMean,aveRMS,1,hname);
    //hsub_all[j]->GetZaxis()->SetRangeUser(-200,250);
    gStyle->SetPalette(56);
    hsub_all[j]->Draw("colz");
    gPad->Update();
    fout_mask=Form("%s/%08d/hsub_all_r%08d_c%d",pdir.Data(),rnum,rnum,j);
    //save_pictures(fout_mask,"cmy");
    hname = Form("hsub_%d",j);
    hsub_all[j]->Write(hname);
    
    hraw_all[j]->Draw("colz");
    gPad->Update();
    fout_mask=Form("%s/%08d/hraw_all_r%08d_c%d",pdir.Data(),rnum,rnum,j);
    //save_pictures(fout_mask,"cmy");
    hname = Form("hraw_%d",j);
    hraw_all[j]->Write(hname);
    //Plot mean pixel value, rms , and temp across a run
    for (int i=0; i < n; i++)
      {
	meanVec[i]=mean[i][j];
	rmsVec[i]=rms[i][j];
	ccdtempvec[i]=ccdtemp[i][j];
      }

    TGraphErrors *g = new TGraphErrors(n,Nev,rmsVec,nevEr,nevEr);
    gname = Form("RMS_%d",j);
    g->SetTitle(gname);
    g->SetMarkerColor(4);
    g->SetMarkerStyle(21);
    g->Draw("ap");
    g->Write(gname);

    //fout_mask=Form("%s/%08d/noiseGraph_hsub_r%08d_c%d",pdir.Data(),rnum,rnum,j);
    //save_pictures(fout_mask,"cmy");

    TGraphErrors *g2 = new TGraphErrors(n,Nev,meanVec,nevEr,nevEr);
    gname = Form("Mean_%d",j);
    g2->SetTitle(gname);
    g2->SetMarkerColor(4);
    g2->SetMarkerStyle(21);
    g2->Draw("ap");
    g2->Write(gname);

    //fout_mask=Form("%s/%08d/meanGraph_hsub_r%08d_c%d",pdir.Data(),rnum,rnum,j);
    //save_pictures(fout_mask,"cmy");

    TGraphErrors *g3 = new TGraphErrors(n,Nev,ccdtempvec,nevEr,nevEr);
    gname = Form("ccdTemp_%d",j);
    g3->SetTitle(gname);
    g3->SetMarkerColor(4);
    g3->SetMarkerStyle(21);
    g3->Draw("ap");
    g3->Write(gname);

    //fout_mask=Form("%s/%08d/tempGraph_hsub_r%08d_c%d",pdir.Data(),rnum,rnum,j);
    //save_pictures(fout_mask,"cmy");

  }
  
  foutFile->Write();
  foutFile->Close();

  Info("process","Output data has been saved in file='%s'",fout.Data());
}


