// plot pedestal of selected pixel for multiple runs.

enum {mxF=26};
int runA[mxF]={1022003, 1022004, 1022005, 1022006, 1022007, 1022008, 1022009,
	       1022010, 1022011, 1022012,  1022013, 1022014, 1022015, 1022016, 1022017, 1022018, 1022019,
	       1022021,   1022022,   1022023,   1022024,  1022025,   1022026,  1022028,   1023001,   1023002};

TFile *fdA[mxF];


void plPixPedRun(  int xPix=800, int yPix=1800, int nReb=16){
  TString inpPath="~/calib/doneCalib/";
  //TString inpPath="./";

  // open files 
  for(int j=0; j<mxF;j++) {
    TString fname=inpPath+Form("m3_Michael_R%d.m3ped.root",runA[j]);
    TFile *fd=new TFile(fname); 
    assert(fd->IsOpen());
    fdA[j]=fd;
    //printf("%d opened %s\n",j,fd->GetName());
    //break;
  }

  TString coreName=Form("pix_x%d_y%d",xPix,yPix);
  TH1F *h1=new TH1F(coreName, coreName+Form(" avr pedestal vs. run, nReb=%d; run index (day22+); (ADU)",nReb), mxF, 0.5, mxF+0.5);
  h1->SetMarkerStyle(4);
  
  c=new TCanvas(coreName,coreName);
  populatePixel(xPix, yPix, h1,nReb);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  h1->Draw();
  h1->Fit("pol0","","R",10,22.4);
  gPad->SetGrid();
  gPad->SetLeftmarin(0.15)
}

//-----------------------------------
void populatePixel( int xPix, int yPix, TH1F *h1, int nReb){  
  TString avrName="cam3_pedAvr";
  float par_nFrame=200;
  float fact=nReb*nReb;

  for(int j=0; j<mxF;j++) {
    TH2I *h2avr=fdA[j]->Get(avrName); assert(h2avr);
    h2avr->Rebin2D(nReb,nReb);
 
    int xBin=h2avr->GetXaxis()->FindBin(xPix);
    int yBin=h2avr->GetYaxis()->FindBin(yPix);

    double val=h2avr->GetBinContent(xBin, yBin)/fact;
    double err=h2avr->GetBinError(xBin, yBin)/fact/sqrt(par_nFrame);

    h1->SetBinContent(j+1,val);
    h1->SetBinError(j+1,err);

    printf("j=%d  %d/%d, val=%.1f sig=%.3f \n",j,xBin, yBin,val,err);
    //
     
  } 

}
