#include <assert.h>
#include "TH2.h"
#include "TH2I.h"
#include "TLine.h"
#include "TList.h"
#include "TFile.h"
#include "TText.h"
#include "TCanvas.h"
#include <TTimeStamp.h>
#include <math.h>
#include "Dataset.hh"
#include "Event.hh"
#include "Image.hh"
#include "CameraInfo.hh"

#include "CcdPedMaker.hh" 

//=======================================
//=======================================
void dmtpc::skim::CcdPedMaker::initDims(const dmtpc::core::Dataset *ds , int camid,  int camsrc) {

  camId=camid;
  camDaq=camsrc;

  printf(" CcdPedMaker()::initDims  camId=%d  camDaq=%d  \n", camId,camDaq);
  coreName=Form("cam%d",camId);

  const dmtpc::core::Event * eve=ds->event();
  const TTimeStamp * eveTime=eve->getTimestamp();
  const dmtpc::core::Image *image=eve->getImage(camDaq); assert(image);
  const dmtpc::core::CameraInfo *info=image->getInfo(); assert(info);
  camSN=info->serialNumber; 
  int eveId=eve->ev();
  int ncam_in =eve->nccd();
 
  printf(" First frame=%d image camDaq=%d  SN=%s   temp/C=%.1f\n", eveId,camDaq, info->serialNumber.Data(), info->ccdTemp);
  eveTime->Print();
  printf(" timeStamp sec=%d , exposure length =%d (msec)\n", eveTime->GetSec(),info->exposureTime);
  assert(camDaq>=0); 
  assert(camDaq<ncam_in);
 

  //  fullNumPix=3084; // hardcoded 4-shooter image size, add if/else for 1-shooter

  runId=eve->run();

  TH2I * h2raw=(TH2I * )eve->ccdData(camDaq); assert(h2raw);
  //....  further rebin image if requested
  if(par_userRebin!=1) {
    h2raw->Rebin2D(par_userRebin,par_userRebin);
  }
  
  const int nbinsx =  h2raw->GetNbinsX();
  const int nbinsy =  h2raw->GetNbinsY();
  const double xmin = h2raw->GetXaxis()->GetXmin();
  const double xmax = h2raw->GetXaxis()->GetXmax();
  const double ymin = h2raw->GetYaxis()->GetXmin();
  const double ymax = h2raw->GetYaxis()->GetXmax();
  int nReb=(int)(xmax/nbinsx+0.1);

  printf("CcdPedMaker::initDims, Event header info: ncam_in=%d,nbinsx=%d,nbinsy=%d, xmin/max=%.1f/%.1f, ymin/max=%.f/%.1f guess rb=%d\n", ncam_in, nbinsx,nbinsy, xmin,xmax,ymin,ymax,nReb);

  //  if(fullNumPix!=xmax) printf(" WARN, for this run=%d  the raw images range has  pixel/histo-bin  mismatch of %.6f\n",eve->run(),1.*xmax/fullNumPix);

  // I never tested code for  a non-square CCD image format
  assert(nbinsx==nbinsy);
  assert(xmin==ymin);
  assert(xmax==ymax);

 //--------------- save values needed elsewere ------------
  nchX=nbinsx;
  nchY=nbinsy;
  agregReb=nReb;

  // find  gross range of pedestal, including masked channel for now
  h1Dped=new TH1D(coreName+"_aduLo","raw ADU spectrum lower end, all input; ADU",1000,0,50000);
  for(int jchX=1; jchX<=nchX; jchX++)  // range 1,N
    for(int jchY=1; jchY<=nchY; jchY++) { // range 1,N
      int mBin=h2raw->GetBin(jchX, jchY); //  pick any 2D histo for this
      float rawVal=h2raw->GetBinContent(mBin);    
      h1Dped->Fill(rawVal);
    }  
  int peakBin=h1Dped->GetMaximumBin();
  float peakAbscia= h1Dped->GetBinCenter(peakBin);
  printf("%s-frame0  MPV ADU=%.0f\n",coreName.Data(),peakAbscia);
  assert(peakAbscia>500); // sth went worng with range of ADU
  setPar(peakAbscia-150*sqrt(agregReb),peakAbscia+350); // aduLo, pedAvrHi

  // pedestal determination
  par_pedWindowHalf=60*sqrt(agregReb); // in ADU
  
  // define sensitive pixel limits, note Hbook convention for the output
  if (nchX==1542) { // 2x2 bining
    actLim.loX=actLim.loY=11;  
    actLim.hiX=actLim.hiY=1538;     
  } else if (nchX==1028) { // 3x3 bining
    actLim.loX=actLim.loY=7;  
    actLim.hiX=actLim.hiY=1026;     
  } else if (nchX==771) { // 4x4 bining
    actLim.loX=actLim.loY=6;  
    actLim.hiX=actLim.hiY=769;     
  } else if (nchX==385) { // 8x8 bining
    actLim.loX=actLim.loY=4;  
    actLim.hiX=actLim.hiY=384;     
  } else {
    assert(3==98); // implement other binnings
  }

  maxMbin=h2raw->GetBin(nchX, nchY); // the largest 1D bin index of interest

  // create containers for output calibration
  hPedStat=new TH2S(coreName+"_pedStat",coreName+Form(" pedestal status ,SN=%s ,run=%d",camSN.Data(),runId),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  hPedAvr=new TH2F(coreName+"_pedAvr",coreName+" pedestal AVR err=RMS ,SN="+camSN,nbinsx,xmin,xmax,nbinsy,ymin,ymax);

  hPedRms=new TH2F(coreName+"_pedRms",coreName+" pedestal RMS ,SN="+camSN,nbinsx,xmin,xmax,nbinsy,ymin,ymax);  hPedRms->SetMinimum(7);

  // for convenience this monitoring hist is created here as well
  hPedPeakSum=new TH2F(coreName+"_pedPeakSum",coreName+" ped peak ampl (frames), accepted so far, "+camSN,nbinsx,xmin,xmax,nbinsy,ymin,ymax);


  printf( "CcdPedMaker()::initDims end, pedWindowHalf=%.1f (ADU) , CUTS:  peakSumFrac=%.2f,  rms=[%.1f,%.1f] ADU, valH=%.1f ADU, can=%p  userRbin=%d\n", par_pedWindowHalf,cut_peakSum, cut_pedRmsL, cut_pedRmsH, cut_pedValH,can,par_userRebin);

  
  initHisto();
  initMask();
}


//=======================================
//=======================================
void dmtpc::skim::CcdPedMaker::saveHisto(TFile *fd, int flag) {

  printf(" %s ped-histos  are written  to '%s' , flag=%d" ,coreName.Data(),fd->GetName(),flag);
  for(int i=0;i<mxH;i++) {
    if(hA[i]) hA[i]->Write();
  }

  h1Dped->Write();
  if (flag==1) {
    printf(" saved also big histos\n");
    hBigPed->Write(); // tmp
  }
  printf(" , save Ok\n");
}






//=======================================
//=======================================
void dmtpc::skim::CcdPedMaker::ingest_pedSpecBig(TString fName) {

  TFile *fd=new TFile(fName);
  assert(fd->IsOpen());
  
  hBigPed=(TH2S*) fd->Get(Form("cam%d_bigPed",camId));
  if(hBigPed==0) fd->ls();
  assert(hBigPed);

  printf(" ingested histo %s from '%s'  with %.1f entries\n" ,hBigPed->GetName(),fd->GetName(),hBigPed->GetEntries());
  // verify at least the X,Y dimensions
  // printf("aa %d %d\n",hBigPed->GetNbinsX() , maxMbin);
  assert(hBigPed->GetNbinsX() == maxMbin);

  TH1D *hA0=(TH1D *)fd->Get(Form("cam%d_case",camId));
  assert(hA0);
  int kBin=hA0->GetXaxis()->FindBin("FRAME");
  assert(kBin>1);

  float nFrame=hA0->GetBinContent(kBin);
  printf(" seen %.0f frames kBin=%d\n",nFrame,kBin);
  assert(nFrame>50); // some minimal body of data is needed to find pedestals
  NtotEve=nFrame;

  hA[0]->Fill("FRAME",NtotEve); // propagate this to final histo

}

//=======================================
//=======================================
void dmtpc::skim::CcdPedMaker::adjustCuts(){ 
  cut_peakSum*=NtotEve;

  TH1 *h;
  h=hA[11];
  { TList *tl=h->GetListOfFunctions(); assert(tl);
    TLine *ln=new TLine(cut_peakSum, 0, cut_peakSum, 1e6);    tl->Add(ln);
    ln->SetLineColor(kMagenta); ln->SetLineStyle(2);
  }

  printf( " adjust CUTS    peakSum=%.1f (frames) , totFrames=%d\n",cut_peakSum,NtotEve);

}

//=======================================
//=======================================
void dmtpc::skim::CcdPedMaker::initHisto() {

  memset(hA,0,sizeof(hA)); // input spectra, unweighted
  TH1* h;
  int mxFrames=500; // no impact on performance

  int nCase=mxPedStat;
  hA[0]=h=new TH1D(coreName+"_case","QA of pedestals; ; channel count",nCase,0,nCase);  h->GetXaxis()->SetTitleOffset(0.4);  h->GetXaxis()->SetLabelSize(0.06);  h->GetXaxis()->SetTitleSize(0.05); h->SetMinimum(0.8);  h->SetLineColor(kBlue); h->SetFillColor(42);h->SetLineWidth(2);  h->SetMarkerSize(2);//<-- large text // no-W
    
  TText mesg(0.3,0.5,"fixMe17");
  TString nameA[mxPedStat]={"active","lowSum","pedRmsLo","pedRmsHi","pedAvrHi","good","FRAME"};

  // prescale random rejection factor for different types of error
  int  xpdfFract[mxPedStat]={ -1, 1, 1,1,100,10000,-1};

  for(int i=0;i<mxPedStat;i++) {
    hA[0]->Fill(nameA[i],0.);
    if(i==0) continue;
    if(i>=6) continue;
    if(can==0) continue; // no plotting of individual pixels
    can->Clear();
    pdfName[i]=Form("R%d_pedQA_%s.pdf",runId,nameA[i].Data());
    //  open output PSDF files
    pdfFract[i]=xpdfFract[i];
    TString aa=Form("CCD pixels type=%s  for %s run R%d,  selection: 1/%d",nameA[i].Data(),coreName.Data(),runId,pdfFract[i]);
    printf("aa=%s=\n",aa.Data());
    mesg.SetText(0.05,0.5,aa);
    mesg.Draw();
    can->Print( pdfName[i]+"(","pdf");
   
  }

  int aduHi=par_aduLo+400*sqrt(agregReb), nAduBins=60; // defines range of pedestals spectra considered as reasonable

  // histos 1...9 
  hA[1]=0;

  hA[9]=0;
  // 1D histos 10 ..19 :  QA of pedestal determination
  hA[10]=0;
  
  hA[11]=h=new TH1D(coreName+"_peakSum_cut","profile of ped peak sum, accepted so far; frame count; pixel count", mxFrames,-0.5 ,mxFrames-0.5);
 

  hA[12]=h=new TH1D(coreName+"_peakAmpl","profile of ped peak amplitude, accepted so far; frame count; pixel count", mxFrames/2,-0.5,mxFrames/2.-0.5);
  

  hA[13]=h=new TH1D(coreName+"_pedRms_cut","profile of pedestal width, accepted so far; pedestal RMS (ADU); pixel count", 100,0,100);
  { TList *tl=h->GetListOfFunctions(); assert(tl);
    TLine *ln=new TLine(cut_pedRmsH, 0, cut_pedRmsH, 1e6);    tl->Add(ln);
    ln->SetLineColor(kMagenta); ln->SetLineStyle(2);
    ln=new TLine(cut_pedRmsL, 0, cut_pedRmsL, 1e6);    tl->Add(ln);
    ln->SetLineColor(kMagenta); ln->SetLineStyle(2);
  }

  hA[14]=h=new TH1D(coreName+"_pedVal_cut","profile of pedestal average, accepted so far; pedestal value(ADU); pixel count", 100,par_aduLo,aduHi);
  { TList *tl=h->GetListOfFunctions(); assert(tl);
    TLine *ln=new TLine( cut_pedAvrHi, 0, cut_pedAvrHi, 1e6);    tl->Add(ln);
    ln->SetLineColor(kMagenta); ln->SetLineStyle(2);
  }
  
  // 2D histos 20 ..29 :  QA of pedestal determination

  hA[20]=hPedStat;
  h=hA[21]=hPedAvr;
  hA[22]=hPedRms;
  hA[23]=hPedPeakSum;
  
  h->SetMinimum(par_aduLo); 
  h->SetMaximum(aduHi); 

  //TH1S : histograms with one short per channel. Maximum bin content = 32767
  hBigPed=new TH2S(coreName+"_bigPed",coreName+Form(" pedestal spectra run=%d ; unified mBin maps to (iX,iY*NX) ; ADU",runId), maxMbin,0.5,maxMbin+0.5,nAduBins,par_aduLo,aduHi);

  printf("pedMaker::initHisto(run=%d)  works w/ %d pixels, aduLo=%.0f pedAvrHi=%.0f\n",runId, maxMbin,par_aduLo, cut_pedAvrHi);

}


