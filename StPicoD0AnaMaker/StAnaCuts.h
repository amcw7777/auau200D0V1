#ifndef CHARMMAKERCUTS_H
#define CHARMMAKERCUTS_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */


namespace anaCuts
{
	// path to lists of triggers prescales
	// lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
	//std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";

	//event
	int const centrality = 3; // StRefMultCorr::getCentralityBin9() < centrality; 0-3 bins are for 40-80%
	unsigned int const triggers[5] = { 
		530003,    // vpdmb-5-p-nobsmd-hlt 
		430903,    // vpdmb-5-p-nobsmd-hlt 
		430904,    // vpdmb-5-p-nobsmd 
		430906,    // vpdmb-5-p-nobsmd 
		430907
	};    // vpdmb-5-p-nobsmd 
	// event
	float const vz = 6.0;
	float const vzVpdVz = 3.0;
	static const int NCENT = 4;
	static const int NQUAL = 6;
	float const centCuts[NCENT] = {220, 70, 30, 0};

	//tracking
	int  const nHitsFit = 20;
	bool const requireHFT = true;
	float const minDca = 0.0050;
	float const minPt = 0.2;
	float const eta   = 1.0;

	// PID
	float const nSigmaPion = 3.0;
	float const nSigmaKaon = 2.0;
	float const nSigmaProton = 3.0;

	//
	float const cosTheta = 0.95; // minimum
	float const dcaDaughters = 0.0100; // maximum
	float const decayLength = 0.0050; // minimum
	float const minD0Mass = 1.6;
	float const maxD0Mass = 2.2;
	float const minKPiXMass = 1.6;
	float const maxKPiXMass = 2.6;

	// histograms kaonPion pair cuts
	int   const nPtBins = 5;
	float const PtBinsEdge[nPtBins+1] = {0., 1., 2., 3., 5., 15.};//this is for optimaized cut
	float const qaNHitsFit = 20;
	float const qaNSigmaKaon = 2.0;
	float const qaRapidityCut = 1.0;
	float const qaDcaV0ToPv[nPtBins] = {0.0061, 0.0049, 0.0038, 0.0038, 0.0040};
	float const qaDecayLength[nPtBins] = {0.0145, 0.0181, 0.0212, 0.0247, 0.0259};
	float const qaCosTheta[nPtBins] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};//0.995
	float const qaDcaDaughters[nPtBins] = {0.0084, 0.0066, 0.0057, 0.0050, 0.0060}; //0.0050;
	float const qaKDca[nPtBins] = {0.0103, 0.0091, 0.0095, 0.0079, 0.0058};//0.008, // minimum
	float const qaPDca[nPtBins] = {0.0110, 0.0111, 0.0086, 0.0081, 0.0062};//0.008
}
#endif
