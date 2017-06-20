{

  //  noise_in_ADU = sqrt[ (read_noise_in_e*ADU_per_e_for_this_camera)**2 + (dark_rate_in_e_per_pixel_per_s xposure_length_in_s * pixels * ADU_per_e_for_this_camera)**2 ]                                                  

  // const int ncam = 3;
  const int ncam = 1;
  /* std::map<TString, int> camera;
    cameras["PL0141514"] = 0;
    cameras["PL0251514"] = 1;    
    cameras["PL0261514"] = 2;*/

  double pixels = 3056*3056;
  
  
  double read_noise[ncam];
  double adu_pe[ncam];
  double dark[ncam];
  
  //read_noise[0] = 11.2;
  //read_noise[1] = 10.2;
  read_noise[0] = 11.3; //9.6;                                                                                  
  
  //adu_pe[0] = 1.55;
  //adu_pe[1] = 1.52;
  adu_pe[0] = 1.53; //1.52;                                                                                     
  
  
  //dark should be:  dark[0] = 0.6/(100*adu_pe[1]*pixels); but going to multiply by pixels and gain in calculation, so just omitting here
  //dark[0] = 2.5/(100);
  //dark[1] = 0.6/(100);
  dark[0] = 40.8/(100); //0.7/(100);                                                                            
  
  double exp[14] = {1.0,2.0,5.0,10.0,15.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,120.0};
 
  double noise[14];
  //TGraph * graphs[ncam];
  //TCanvas * c1[ncam];
  
  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  TMultiGraph *mg = new TMultiGraph();
  //TGraph *gr1[ncam]; //= new TGraph(); 
  //gr1->SetLineColor(kBlue);
  //TGraph *gr2 = new TGraph(); 
 //gr2->SetLineColor(kRed);
  
 //for(int c = 0; c < ncam; c++){
  for(int i = 0; i < 14; i++){
    noise[i] = (TMath::Sqrt(read_noise[0]*adu_pe[0]*read_noise[0]*adu_pe[0] + dark[0]*exp[i]*dark[0]*exp[i]))/1.41;
    //}
   
    gr1 = new TGraph(14,exp,noise);
    gr1->SetLineColor(kRed);
    /*graphs[c] = new TGraph(5,exp,noise);
      c1[c] = new TCanvas();
      c1[c]->cd();
      graphs[c]->SetTitle(Form("Predicted Noise vs. Exposure Cam %d",c));
      graphs[c]->GetXaxis()->SetTitle("Exposure Time (s)");
      graphs[c]->GetYaxis()->SetTitle("Predicted Noise (ADU)");
      graphs[c]->Draw("alp");*/
    gr1->SetLineWidth(4);

   }

  Double_t exp_time[3] = {2.0,15.0,30.0};
  Double_t noise_meassured[3] = {12.7,12.7,12.7};
  Double_t ex[3] = {.1,.1,.5};
  Double_t ey[3] = {.8,.7,2.};
  
  TGraphErrors *gr2 = new TGraphErrors(3,exp_time,noise_meassured,ex,ey);
  // gr2->SetLineColor(kRed);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerColor(kBlue);
  
  mg->Add(gr1);
  mg->Add(gr2);
  
  c1->cd(); 
  mg->Draw("alp");
  mg->SetTitle(Form("Predicted Noise vs. Exposure Time"));                                          
  mg->GetXaxis()->SetTitle("Exposure Time (s)");                                                        
  mg->GetYaxis()->SetTitle("Predicted Noise (ADU)");                                                    
  //  mg->Draw("alp");
  c1->Update();
  
}
