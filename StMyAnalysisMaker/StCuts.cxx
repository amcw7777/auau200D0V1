/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

// #include "StCuts.h"

#include "Rtypes.h"
namespace mycuts
{
   // path to lists of triggers prescales
   // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
   std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
   //event
   UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html
   float const vz = 6.0;// cm.
   float const vzVpdVz = 3.0; // 3 cm.

   //tracking
   int const nHitsFit = 20;
   bool const requireHFT = true;
   float const minPt = 0.6;

   //pions
   float const nSigmaPion = 3.0;

   //kaons
   float const nSigmaKaon = 2.0;
   float const kTofBetaDiff = 0.03;

   // tree kaonPion pair cuts
   // float const cosTheta = 0.995; // minimum
   // float const dcaDaughters = 0.0050; // maximum
   // float const decayLength = 0.0030; // minimum
   float const minMass = 1.6;
   float const maxMass = 2.1;

   // histograms kaonPion pair cuts
   float const qaNHitsFit = 20;
   float const qaNSigmaKaon = 2.0;
   // float const kDca = 0.008; // minimum
   // float const pDca = 0.008;
   // Hadron cuts
   float const hadronPtMin= 0.2;
   float const hadronPtMax = 2.0;
   float const corDetaMin = 0.75;
   float const corDetaMax = 10;
    
    //============================== new cuts: 5 cent bin X 6 pt bin
    const int nCent = 5;
    const float CentEdge[nCent+1] = { -0.5, 1.5, 3.5, 5.5, 6.5, 8.5 };
    const char nameCent[nCent][100] = {"60-80%", "40-60%", "20-40%", "10-20%", "0-10%"};
    
    const int nPtBins = 6;
    const float PtEdge[nPtBins+1] = {0, 0.5, 1., 2., 3., 5., 10.};
    
    // default
    float const cosTheta = 0.95;
    float const kDca[nCent][nPtBins] = {
        { 0.0113, 0.0103, 0.0081, 0.0066, 0.0046, 0.0038},  //60%-80%
        { 0.0110, 0.0112, 0.0081, 0.0063, 0.0064, 0.0044},  //40%-60%
        { 0.0098, 0.0089, 0.0074, 0.0085, 0.0063, 0.0049},  //20%-40%
        { 0.0111, 0.0099, 0.0091, 0.0097, 0.0062, 0.0061},  //10%-20%
        { 0.0104, 0.0099, 0.0073, 0.0086, 0.0067, 0.0056}   //0-10%
    };
    float const pDca[nCent][nPtBins] = {
        { 0.0100, 0.0096, 0.0093, 0.0094, 0.0059, 0.0050},  //60%-80%
        { 0.0107, 0.0106, 0.0097, 0.0078, 0.0063, 0.0056},  //40%-60%
        { 0.0117, 0.0106, 0.0097, 0.0066, 0.0064, 0.0056},  //20%-40%
        { 0.0098, 0.0110, 0.0101, 0.0086, 0.0070, 0.0061},  //10%-20%
        { 0.0109, 0.0106, 0.0081, 0.0092, 0.0080, 0.0057}   //0-10%
    };
    float const dcaV0ToPv[nCent][nPtBins] = {
        { 0.0075, 0.0066, 0.0064, 0.0050, 0.0058, 0.0038},  //60%-80%
        { 0.0065, 0.0064, 0.0046, 0.0049, 0.0054, 0.0057},  //40%-60%
        { 0.0045, 0.0048, 0.0042, 0.0043, 0.0052, 0.0055},  //20%-40%
        { 0.0052, 0.0049, 0.0043, 0.0042, 0.0043, 0.0050},  //10%-20%
        { 0.0061, 0.0046, 0.0042, 0.0041, 0.0037, 0.0048}   //0-10%
    };
    float const dcaDaughters[nCent][nPtBins] = {
        { 0.0073, 0.0088, 0.0092, 0.0082, 0.0083, 0.0104},  //60%-80%
        { 0.0076, 0.0087, 0.0090, 0.0082, 0.0101, 0.0093},  //40%-60%
        { 0.0078, 0.0067, 0.0069, 0.0066, 0.0073, 0.0099},  //20%-40%
        { 0.0070, 0.0070, 0.0062, 0.0067, 0.0076, 0.0085},  //10%-20%
        { 0.0066, 0.0079, 0.0061, 0.0063, 0.0076, 0.0068}   //0-10%
    };
    float const decayLength[nCent][nPtBins] = {
        { 0.0150, 0.0107, 0.0175, 0.0187, 0.0164, 0.0175},  //60%-80%
        { 0.0140, 0.0133, 0.0190, 0.0201, 0.0215, 0.0219},  //40%-60%
        { 0.0149, 0.0170, 0.0205, 0.0236, 0.0234, 0.0237},  //20%-40%
        { 0.0151, 0.0173, 0.0204, 0.0240, 0.0237, 0.0231},  //10%-20%
        { 0.0128, 0.0163, 0.0222, 0.0214, 0.0241, 0.0253}   //0-10%
    };
    
    // tight
    float const cosTheta1 = 0.95;
    float const kDca1[nCent][nPtBins] = {
        { 0.0160, 0.0142, 0.0046, 0.0088, 0.0079, 0.0071},  //60%-80%
        { 0.0120, 0.0072, 0.0093, 0.0081, 0.0071, 0.0059},  //40%-60%
        { 0.0151, 0.0089, 0.0074, 0.0093, 0.0082, 0.0071},  //20%-40%
        { 0.0116, 0.0099, 0.0091, 0.0143, 0.0109, 0.0095},  //10%-20%
        { 0.0117, 0.0117, 0.0073, 0.0119, 0.0098, 0.0086}   //0-10%
    };
    float const pDca1[nCent][nPtBins] = {
        { 0.0072, 0.0121, 0.0104, 0.0141, 0.0096, 0.0070},  //60%-80%
        { 0.0163, 0.0132, 0.0114, 0.0125, 0.0140, 0.0136},  //40%-60%
        { 0.0117, 0.0107, 0.0083, 0.0151, 0.0132, 0.0135},  //20%-40%
        { 0.0098, 0.0141, 0.0101, 0.0123, 0.0106, 0.0130},  //10%-20%
        { 0.0152, 0.0106, 0.0081, 0.0118, 0.0109, 0.0133}   //0-10%
    };
    float const dcaV0ToPv1[nCent][nPtBins] = {
        { 0.0057, 0.0057, 0.0038, 0.0032, 0.0031, 0.0031},  //60%-80%
        { 0.0057, 0.0046, 0.0036, 0.0038, 0.0036, 0.0041},  //40%-60%
        { 0.0040, 0.0039, 0.0031, 0.0030, 0.0038, 0.0052},  //20%-40%
        { 0.0051, 0.0038, 0.0034, 0.0030, 0.0033, 0.0051},  //10%-20%
        { 0.0040, 0.0037, 0.0030, 0.0034, 0.0033, 0.0038}   //0-10%
    };
    float const dcaDaughters1[nCent][nPtBins] = {
        { 0.0081, 0.0078, 0.0067, 0.0076, 0.0069, 0.0131},  //60%-80%
        { 0.0076, 0.0050, 0.0073, 0.0055, 0.0084, 0.0138},  //40%-60%
        { 0.0066, 0.0070, 0.0057, 0.0059, 0.0054, 0.0098},  //20%-40%
        { 0.0061, 0.0055, 0.0056, 0.0047, 0.0053, 0.0061},  //10%-20%
        { 0.0073, 0.0064, 0.0061, 0.0039, 0.0053, 0.0098}   //0-10%
    };
    float const decayLength1[nCent][nPtBins] = {
        { 0.0100, 0.0109, 0.0212, 0.0144, 0.0168, 0.0476},  //60%-80%
        { 0.0140, 0.0176, 0.0242, 0.0259, 0.0282, 0.0291},  //40%-60%
        { 0.0108, 0.0223, 0.0248, 0.0149, 0.0305, 0.0303},  //20%-40%
        { 0.0213, 0.0187, 0.0259, 0.0258, 0.0282, 0.0100},  //10%-20%
        { 0.0128, 0.0196, 0.0266, 0.0253, 0.0392, 0.0357}   //0-10%
    };
    
    // loose
    float const cosTheta2 = 0.95;
    float const kDca2[nCent][nPtBins] = {
        { 0.0086, 0.0072, 0.0056, 0.0050, 0.0020, 0.0020},      //60%-80%
        { 0.0114, 0.0081, 0.0081, 0.0061, 0.0040, 0.0035},      //40%-60%
        { 0.0098, 0.0089, 0.0074, 0.0062, 0.0050, 0.0035},      //20%-40%
        { 0.0094, 0.0111, 0.0068, 0.0067, 0.0050, 0.0035},      //10%-20%
        { 0.0099, 0.0091, 0.0073, 0.0070, 0.0057, 0.0035}       //0-10%
    };
    float const pDca2[nCent][nPtBins] = {
        { 0.0072, 0.0072, 0.0053, 0.0048, 0.0023, 0.0020},      //60%-80%
        { 0.0090, 0.0083, 0.0081, 0.0064, 0.0040, 0.0035},      //40%-60%
        { 0.0092, 0.0110, 0.0083, 0.0066, 0.0050, 0.0035},      //20%-40%
        { 0.0085, 0.0103, 0.0078, 0.0068, 0.0051, 0.0035},      //10%-20%
        { 0.0102, 0.0106, 0.0081, 0.0082, 0.0053, 0.0035}       //0-10%
    };
    float const dcaV0ToPv2[nCent][nPtBins] = {
        { 0.0086, 0.0084, 0.0080, 0.0079, 0.0084, 0.0090},      //60%-80%
        { 0.0064, 0.0080, 0.0064, 0.0082, 0.0077, 0.0098},      //40%-60%
        { 0.0080, 0.0080, 0.0056, 0.0048, 0.0082, 0.0091},      //20%-40%
        { 0.0078, 0.0061, 0.0055, 0.0048, 0.0050, 0.0092},      //10%-20%
        { 0.0065, 0.0060, 0.0047, 0.0048, 0.0051, 0.0095}       //0-10%
    };
    float const dcaDaughters2[nCent][nPtBins] = {
        { 0.0143, 0.0116, 0.0107, 0.0114, 0.0148, 0.0126},      //60%-80%
        { 0.0090, 0.0099, 0.0093, 0.0098, 0.0127, 0.0146},      //40%-60%
        { 0.0094, 0.0086, 0.0081, 0.0095, 0.0124, 0.0139},      //20%-40%
        { 0.0106, 0.0090, 0.0091, 0.0081, 0.0085, 0.0141},      //10%-20%
        { 0.0084, 0.0080, 0.0071, 0.0085, 0.0099, 0.0147}       //0-10%
    };
    float const decayLength2[nCent][nPtBins] = {
        { 0.0123, 0.0109, 0.0139, 0.0144, 0.0100, 0.0119},      //60%-80%
        { 0.0149, 0.0133, 0.0160, 0.0173, 0.0129, 0.0124},      //40%-60%
        { 0.0100, 0.0131, 0.0180, 0.0186, 0.0131, 0.0118},      //20%-40%
        { 0.0111, 0.0143, 0.0188, 0.0177, 0.0197, 0.0145},      //10%-20%
        { 0.0128, 0.0153, 0.0188, 0.0187, 0.0196, 0.0148}       //0-10%
    };

}

