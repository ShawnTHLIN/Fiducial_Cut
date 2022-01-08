/////////////////////////////////////////////////////////
//**#include "analyze.h"
//**#include "HistoClass.h"
//#include "TrigEJ309.h"
//#include "VerifyEJ309.h"
//#include "MainEJ290.h"
//#include "constants.h"
//#include "main.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>
//#include <TCutG.h>

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <time.h>
#include <vector>
#include <thread>
//#include <mpi.h>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"


using namespace std;



// File 
const string MainFileName =     "Data_CH6@DT5730B_1173_20211224_Cs137_5cm_800V.bin";


const int Recordlength = 992;
const int bits = 192;
const int numberOfBytesOfEvent = 312;       // bits devided by 8 (8bits/byte)

int64_t neutron_count =0 ;


//Variable

int64_t fATrigNumofRead;
int64_t fBTrigNumofRead;
int64_t main_detector_signal_counter;
int64_t A_detector_signal_counter;
int64_t B_detector_signal_counter;

//int waveCounter;
//const int numOfWavesToLookAt = 100;

//Number of threads
int world_rank = 0;


//Main EJ-290
ifstream ROIinFile;
ifstream fROIIn;
ifstream MatrixinFile;


int16_t fROIboard;
int16_t fROIchannel;
int64_t fROItimestamp;
int16_t fROIenergy;
int16_t fROIshortEnergy;
int32_t fROIflags;
int16_t fROInumOfWaveSamplesToBeRead; //record length in samples
int16_t sample;
double  fROIlistPsp;
vector<double> fROIcurrentWave; //vector of current wave
double fROIlongGate, fROIshortGate, fROIbaseline, fROIpsp;
//int64_t fROITOF;

//File output

ofstream TrigoutFile, Trigwave;
ofstream BTrigoutFile, BTrigwave;
ofstream ROIoutFile, ROIwave;
ofstream BROIoutFile;
ofstream MatrixoutFile;

string OutTrigFileName, OutTrigWaveFileName;
string OutBTrigFileName, OutBTrigWaveFileName;
string OutROIFileName, OutROIFileWaveName, OutMatrixFileName;
string OutBROIFileName;



//Function
void ROIopenBinaryFile(string inputFileName);
void ROIreadBinaryWave(int64_t start_pos, int64_t end_pos, int world_rank);
void ROIcloseFile();
void WriteToFile();
void WriteArrayToFile();
void Write2DArrayToFile();

void ATrigopenBinaryFile(string inputFileName);
void ATrigreadBinaryWave();
void ATrigcloseFile();

void BTrigopenBinaryFile(string inputFileName);
void BTrigreadBinaryWave();
void BTrigcloseFile();

//TH1F *TrigwaveHisto[numOfWavesToLookAt];
//TH1F *MainwaveHisto[numOfWavesToLookAt];

void ROOTcreateWaveHistos();
void ROOTwriteToFile();
int64_t numberOfEvent;

double fROIPSP;


int Channel[16384]= { };
int Array_2D[100][16384/4]= { }; // PSD VS Channel 




int main(int argc, char** argv){
	
	// read the file
	time_t time_start, time_end;
	time_start = time(NULL);

	ROIopenBinaryFile(MainFileName);


	// open output file for write //
	stringstream ss;
	ss << world_rank;
	string NumofFile = ss.str();

	OutROIFileName =  "20211227_gamma_events.txt";  //need change
	/*
	Rule of output file name
	D 3" EJ309 to source
	d 2" EJ309 to source
	N: neutron ; g1,g2 gamma trigger detector
	wo,w : without, with
	*/
	OutMatrixFileName = "20220108_gamma_PSD_Channel.txt";
	//OutROIFileName =  "EJ290@70-140_EN" + NumofFile + ".txt";
	//OutBROIFileName = "EJ290@50-200_EN" + NumofFile + ".txt";
	
	ROIoutFile.open(OutROIFileName, std::ios_base::app);
	//MatrixoutFile.open(OutMatrixFileName, std::ios_base::app);
	//BROIoutFile.open(OutBROIFileName,std::ios_base::app);
	
	
	//ROOTcreateWaveHistos();
	int64_t start_pos = 0;
	int64_t end_pos = 0;
	ROIreadBinaryWave(start_pos, end_pos, world_rank);


	//ROIinFile.seekg(0);

	//ATrigcloseFile();
	//BTrigcloseFile();
	ROIcloseFile();

	//cout<<"End of analysis"<<endl;


	//WriteArrayToFile();
	cout<<"before"<< endl;
	fgetc(stdin);
	Write2DArrayToFile();
	cout<<"after"<< endl;
	fgetc(stdin);

	time_end = time(NULL);
	cout << time_end - time_start << " second" << endl;

	cout << neutron_count << " neutrons" << endl;
	
	//ROOTwriteToFile();
	

	// MPI_Finalize();
	return 0;
}


/////////////////EJ309-gamma//////////////////////////

void ROIopenBinaryFile(string inputFileName){
	//cout << inputFileName << endl;
	ROIinFile.open(inputFileName.c_str(), ios::binary | ios::in);
	if (!ROIinFile.good()) cout << "MainEJ309 file is not open" << '\n'; //check file is open

	//initize some things when you open the file
}



void ROIcloseFile(){
	ROIinFile.close();
}

void ROIreadBinaryWave(int64_t start_pos, int64_t end_pos, int nThread){
	

	
	fROIIn.open(MainFileName.c_str(), ios::binary | ios::in);
	fROIIn.seekg(start_pos, ios::beg);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
	//fROIcurrentWave.clear();

	while (1){
		
		//Aneutron = false;
		
		/*
		if (waveCounter==100){
			break;
		}
		*/
		fROIcurrentWave.clear();
		fROIIn.read((char*)&fROIboard, sizeof(int16_t));
		fROIIn.read((char*)&fROIchannel, sizeof(int16_t));
		fROIIn.read((char*)&fROItimestamp, sizeof(int64_t)); //pico seconds
		fROIIn.read((char*)&fROIenergy, sizeof(int16_t));
		fROIIn.read((char*)&fROIshortEnergy, sizeof(int16_t));
		fROIIn.read((char*)&fROIflags, sizeof(int32_t));
		fROIIn.read((char*)&fROInumOfWaveSamplesToBeRead, sizeof(int16_t));
		fROIIn.read((char*)&sample, sizeof(int16_t));
		int16_t sample;
		for (int i = 0; i < fROInumOfWaveSamplesToBeRead; i++)
		{
			fROIIn.read((char*)&sample, sizeof(int16_t));
			fROIcurrentWave.push_back(sample);
		}
		
		
		
		//double fROIPSP;
		
		if (fROIenergy == 0) fROIPSP = 0;

		else { fROIPSP = (fROIenergy - fROIshortEnergy) / ((double)fROIenergy *1.0); }


		
		// display term 
		/*
		cout << "__________________mainEJ290______________" << endl;
		cout << "This is for the first event of the file" << endl;
		cout << "board: " << fROIboard << endl;
		cout << "ch: " << fROIchannel << endl;
		cout << "ts: " << fROItimestamp << endl;
		cout << "E: " << fROIenergy << endl;
		cout << "sE: " << fROIshortEnergy << endl;
		cout << "PSP:" << fROIPSP << endl;
		cout << "rec. length: " << fROInumOfWaveSamplesToBeRead << endl;
		cout << "_________________________________________" << endl << endl;
		fgetc(stdin);
		*/

		// For fiducial cut
		/*	
		if (( fROIPSP > (nd[fROIenergy-1]/1000) )&&( fROIPSP < (nu[fROIenergy-1]/1000))){
			neutron_count++;
		//if (( fROIPSP > 0.31 ) && ( fROIPSP < 0.44 )){
		//if ( fROIPSP > 0 ) {
			//cout << "okayyyyyyy"<< endl;
			ATrigreadBinaryWave();
			
		}
		*/
		
		if ( fROIPSP < 0 ){
			continue;
		}


		//int ch = fROIenergy; // for sorting
		
	


		//cout <<" before "<< Channel[ch]<<"   "<< ch << endl;
		
		////// !!! change 
		//if (( fROIPSP > 0.2 ) && ( fROIPSP < 0.4 )){
		//Channel[ch]++;

		//WriteToFile();
		//}
		int Ch_2D = fROIenergy/4;
		int PSD_2D = 100*fROIPSP;
		
		if (Ch_2D<10000 && Ch_2D>0 && PSD_2D>0 && PSD_2D < 100){
				Array_2D[PSD_2D][Ch_2D]++;
		}


		//cout<<"Channel  "<< Ch_2D << " PSD  "<< PSD_2D << endl;
		//fgetc(stdin);


		
		
		/*
		if ( ch >= 16382 ){
			cout <<fROIenergy << endl;
			fgetc(stdin);
		}
		*/
		//cout << Channel[ch]<<"   "<< ch << endl;


		//if (( fROIPSP > 0.3 ) && ( fROIPSP < 0.5 )){
		neutron_count++;
		//ATrigreadBinaryWave();
		
		//WriteToFile();


		//}




		main_detector_signal_counter++;
		if (main_detector_signal_counter % 10000 == 0){
			cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
		
		}
		/*
		if (main_detector_signal_counter == 100){
			cout<<"end of the run at "<< main_detector_signal_counter <<endl;
			break;
		}
		*/

		if(!fROIIn){
			cout<<"break in main loop"<< endl;
			cout << "Neutron counts: "<<neutron_count << endl<<endl;;
			cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
			fROIIn.close();
			break;
		}		

		/*
		if (fROIIn.tellg() >= end_pos || fROIIn.eof()){
			fROIIn.close();
			fATrigIn.close();
			fBTrigIn.close();
			cout << "End of thread: " << nThread << endl;
			return;
		}
		*/
		
		}
		if(!fROIIn){
		cout<<"break in outside loop"<< endl;
		cout << "Neutron counts: "<<neutron_count << endl<<endl;;
		cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
		fROIIn.close();
	}		
	return;
}

void WriteToFile(){

		ROIoutFile  <<  ((fROItimestamp*1.0) / 1000.0) <<"  "<< fROIPSP <<"  " << fROIenergy << endl;

}

void WriteArrayToFile(){
		for(int count = 0; count < 16384; count ++){
			MatrixoutFile << Channel[count] << endl ;
			//cout<< Channel[count]<<endl;
	}
}

void Write2DArrayToFile(){
	MatrixoutFile.open(OutMatrixFileName, std::ios_base::app);
	if(MatrixoutFile.is_open())
	{
	for(int j=0; j<100; j++) {
    	for(int i=0; i<16384/4; i++){
        	MatrixoutFile  << Array_2D[j][i] <<" ";

			}
		 MatrixoutFile << endl;
		}
	
	}
	MatrixoutFile.close();
}
