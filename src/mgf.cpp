
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List parseMgf(String filename) {
  FILE* fp; // File being read
  char* line; // Current line

  // We will parse only
  // BEGIN IONS, TITLE, SCANS, RTINSECONDS, CHARGE, PEPMASS
  // into spectrumInfo data frame.

  // The ions per each spectrum will be stored as sub data frames
  // This is theoretically a bit wasteful (the data.frame metadata will be repeated)
  // But in practice we process spectrum at a time, so it is what it is

  size_t bufferLength = 1024;
  line = (char*)malloc(bufferLength);
  if (!line) {
    stop("Cannot allocate memory");
  }

  fp = fopen(filename.get_cstring(), "r");
  if (!fp) {
    free(line);
    stop("Cannot open file");
  }

  fseek(fp, 0, SEEK_END);
  unsigned long fileSize = (unsigned long)ftell(fp);
  fseek(fp, 0, SEEK_SET);

  int state = 0;
  size_t len = 0;
  long lineNum = 0;

  // Kick it off by reading first line
  len=getline(&line, &bufferLength, fp);
  lineNum++;

  // Buffers for spectrum mz/intensities
  const int maxFragments = 102400;
  double mzVals[maxFragments];
  double intVals[maxFragments];

  int spectrumNumber = -1;
  int fragment = 0; // After last populated fragment

  std::vector<std::string> charge;
  std::vector<std::string> scans;
  std::vector<double> pepmass;
  std::vector<std::string> title;
  std::vector<double> rtInSeconds;
  std::vector<double> mz;
  std::vector<double> intensity;
  std::vector<int> firstEntry;
  std::vector<int> lastEntry;

  mz.reserve(fileSize / 16);
  intensity.reserve(fileSize / 16);

  while(!feof(fp)) {
    // This loop would alternate between various state loops that
    // do just one job for that specific state.
    // Each block, when it exits, would leave the next line that made it quit
    // in the buffer, so next state can process it

    // Waiting for begin ions
    while (len != -1 && state==0) {
      if(strncmp(line, "BEGIN IONS\n", len)==0) {
        printf("Begin ions on line %ld\n", lineNum);
        spectrumNumber++; // New spectrum started

        title.push_back("");
        rtInSeconds.push_back(NAN);
        charge.push_back("");
        pepmass.push_back(NAN);
        scans.push_back("");
        int currentEntry = mz.size()+1;
        firstEntry.push_back(currentEntry);
        lastEntry.push_back(currentEntry-1);
        fragment = 0;

        state = 1; // In begin ions
      }
      len=getline(&line, &bufferLength, fp);
      lineNum++;
    }

    while (len != -1 && state==1) {
      // After begin ions
      if(strncmp(line, "TITLE=", 6)==0) { // 6==strlen(TITLE=)
        title[spectrumNumber] = std::string(line, 6, len-6-1);
      }
      else if(strncmp(line, "RTINSECONDS=", 12)==0) { // 12==strlen(RTINSECONDS=)
        double num;
        if(sscanf(line+12, "%lf", &num)==-1) {
          fclose(fp);
          if(line) free(line);
          fprintf(stderr, "Error on row %ld", lineNum);
          stop("Malformed RTINSECONDS entry");
        }
        // rtInSeconds = std::string(line, 12, len-12-1);
        rtInSeconds[spectrumNumber] = num;
      }
      else if(strncmp(line, "CHARGE=", 7)==0) { // 7==strlen(CHARGE=)
        charge[spectrumNumber] = std::string(line, 7, len-7-1);
      }
      else if(strncmp(line, "SCANS=", 6)==0) { // 6==strlen(SCANS=)
        scans[spectrumNumber] = std::string(line, 6, len-6-1);
      }
      else if(strncmp(line, "PEPMASS=", 8)==0) { // 8==strlen(PEPMASS=)
        double num;
        if(sscanf(line+8, "%lf", &num)==-1) {
          fclose(fp);
          if(line) free(line);
          fprintf(stderr, "Error on row %ld", lineNum);
          stop("Malformed PEPMASS entry");
        }
        pepmass[spectrumNumber] = num;
      }
      else if(line[0] >= '0' && line[0] <= '9') {
        // We are reading mz/intensity pairs
        state=2;
        fragment=0;
        break; // Leave the line in buffer
      }
      else if(strncmp(line, "END IONS\n", len)==0) {
        printf("End ions on line %ld\n", lineNum);

        if (fragment>0) {
          mz.reserve(mz.size()+fragment);
          intensity.reserve(intensity.size()+fragment);

          mz.insert(mz.end(), mzVals, mzVals+fragment);
          intensity.insert(intensity.end(), intVals, intVals+fragment);

          lastEntry[spectrumNumber] = firstEntry[spectrumNumber] + fragment - 1;

          fragment = 0;
        }

        state = 0; // Wait for begin ions
      }
      len=getline(&line, &bufferLength, fp);
      lineNum++;
    }

    // Numbers
    while(len != -1 && state == 2) {
      if (line[0] >= '0' && line[0] <= '9') {
        if(sscanf(line, "%lf %lf", &(mzVals[fragment]), &(intVals[fragment])) == -1) {
          fclose(fp);
          if(line) free(line);
          fprintf(stderr, "Error on row %ld", lineNum);
          stop("Malformed mz/intensity pair");
        }
        fragment++;
        if(fragment >= maxFragments) {
          fclose(fp);
          if(line) free(line);
          fprintf(stderr, "Error on row %ld", lineNum);
          stop("Too many fragments in a spectrum");
        }
        len=getline(&line, &bufferLength, fp);
        lineNum++;
      } else {
        // Still after BEGIN IONS
        state=1;
      }
    }
  }

  fclose(fp);

  if (line) {
    free(line);
  }

  return List::create(
    Named("spectra")=DataFrame::create(
      Named("title")=StringVector(title.begin(), title.end()),
      Named("rtInSeconds")=NumericVector(rtInSeconds.begin(), rtInSeconds.end()),
      Named("pepmass")=NumericVector(pepmass.begin(), pepmass.end()),
      Named("charge")=StringVector(charge.begin(), charge.end()),
      Named("scans")=StringVector(scans.begin(), scans.end()),
      Named("firstEntry")=NumericVector(firstEntry.begin(), firstEntry.end()),
      Named("lastEntry")=NumericVector(lastEntry.begin(), lastEntry.end())
    ),
    Named("fragments")=DataFrame::create(
      Named("mz")=NumericVector(mz.begin(), mz.end()),
      Named("intensity")=NumericVector(intensity.begin(), intensity.end())
    )
  );
}
