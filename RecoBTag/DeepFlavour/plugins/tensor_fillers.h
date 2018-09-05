#ifndef RecoBTag_DeepFlavour_tensor_fillers_h
#define RecoBTag_DeepFlavour_tensor_fillers_h

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"

namespace btagbtvdeep {

  // Note on setting tensor values:
  // Instead of using the more convenient tensor.matrix (etc) methods,
  // we can exploit that in the following methods values are set along
  // the innermost (= last) axis. Those values are stored contiguously in
  // the memory, so it is most performant to get the pointer to the first
  // value and use pointer arithmetic to iterate through the next pointers.

inline  void jet_tensor_filler(tensorflow::Tensor & tensor,
                         std::size_t jet_n,
                         const btagbtvdeep::DeepFlavourFeatures & features) {

    float* ptr = &tensor.matrix<float>()(jet_n, 0);

    // jet variables
    const auto & jet_features = features.jet_features;
    *ptr     = jet_features.pt;
    *(++ptr) = jet_features.eta;
    // number of elements in different collections
    *(++ptr) = features.c_pf_features.size();
    *(++ptr) = features.n_pf_features.size();
    *(++ptr) = features.sv_features.size();
    *(++ptr) = features.npv;
    // variables from ShallowTagInfo
    const auto & tag_info_features = features.tag_info_features;
    *(++ptr) = tag_info_features.trackSumJetEtRatio;
    *(++ptr) = tag_info_features.trackSumJetDeltaR;
    *(++ptr) = tag_info_features.vertexCategory;
    *(++ptr) = tag_info_features.trackSip2dValAboveCharm;
    *(++ptr) = tag_info_features.trackSip2dSigAboveCharm;
    *(++ptr) = tag_info_features.trackSip3dValAboveCharm;
    *(++ptr) = tag_info_features.trackSip3dSigAboveCharm;
    *(++ptr) = tag_info_features.jetNSelectedTracks;
    *(++ptr) = tag_info_features.jetNTracksEtaRel;

  }

inline  void c_pf_tensor_filler(tensorflow::Tensor & tensor,
                          std::size_t jet_n,
                          std::size_t c_pf_n,
                          const btagbtvdeep::ChargedCandidateFeatures & c_pf_features) {

    float* ptr = &tensor.tensor<float, 3>()(jet_n, c_pf_n, 0);

    *ptr     = c_pf_features.btagPf_trackEtaRel;
    *(++ptr) = c_pf_features.btagPf_trackPtRel;
    *(++ptr) = c_pf_features.btagPf_trackPPar;
    *(++ptr) = c_pf_features.btagPf_trackDeltaR;
    *(++ptr) = c_pf_features.btagPf_trackPParRatio;
    *(++ptr) = c_pf_features.btagPf_trackSip2dVal;
    *(++ptr) = c_pf_features.btagPf_trackSip2dSig;
    *(++ptr) = c_pf_features.btagPf_trackSip3dVal;
    *(++ptr) = c_pf_features.btagPf_trackSip3dSig;
    *(++ptr) = c_pf_features.btagPf_trackJetDistVal;
    *(++ptr) = c_pf_features.ptrel;
    *(++ptr) = c_pf_features.drminsv;
    *(++ptr) = c_pf_features.vtx_ass;
    *(++ptr) = c_pf_features.puppiw;
    *(++ptr) = c_pf_features.chi2;
    *(++ptr) = c_pf_features.quality;

  }

inline  void n_pf_tensor_filler(tensorflow::Tensor & tensor,
                          std::size_t jet_n,
                          std::size_t n_pf_n,
                          const btagbtvdeep::NeutralCandidateFeatures & n_pf_features) {

    float* ptr = &tensor.tensor<float, 3>()(jet_n, n_pf_n, 0);

    *ptr     = n_pf_features.ptrel;
    *(++ptr) = n_pf_features.deltaR;
    *(++ptr) = n_pf_features.isGamma;
    *(++ptr) = n_pf_features.hadFrac;
    *(++ptr) = n_pf_features.drminsv;
    *(++ptr) = n_pf_features.puppiw;

  }

inline  void sv_tensor_filler(tensorflow::Tensor & tensor,
                          std::size_t jet_n,
                          std::size_t sv_n,
                          const btagbtvdeep::SecondaryVertexFeatures & sv_features) {

    float* ptr = &tensor.tensor<float, 3>()(jet_n, sv_n, 0);

    *ptr     = sv_features.pt;
    *(++ptr) = sv_features.deltaR;
    *(++ptr) = sv_features.mass;
    *(++ptr) = sv_features.ntracks;
    *(++ptr) = sv_features.chi2;
    *(++ptr) = sv_features.normchi2;
    *(++ptr) = sv_features.dxy;
    *(++ptr) = sv_features.dxysig;
    *(++ptr) = sv_features.d3d;
    *(++ptr) = sv_features.d3dsig;
    *(++ptr) = sv_features.costhetasvpv;
    *(++ptr) = sv_features.enratio;

  }
  
inline  void seed_tensor_filler(tensorflow::Tensor & tensor,
                          std::size_t jet_n,
                          std::size_t seed_n,
                          const btagbtvdeep::SeedingTrackFeatures & seed_features) {
    
    
    btagbtvdeep::SeedingTrackFeatures seed_features_mean;
    btagbtvdeep::SeedingTrackFeatures seed_features_std;

    //             svars_mean_BIS=numpy.array([0.35833012, -0.000795975586, -0.0183847683, 0.129604101, -0.00631919497, -0.0519152703, -3.65962456, 1.71627243, -4.48109269, 0.932680501, -2.70270809, 1.63019572, -1.76747277, 0.926253178, 0.112491541, 0.288080422, 0.037803297, 3.9065406, 15.1244042, -4.68685508, -1.47584773])
    //             svars_std_BIS=numpy.array([0.29620678, 1.17231085, 1.77591521, 0.03431668, 4.34793758, 4.90624525, 1.68269796, 1.3829665, 1.99893739, 1.6451007, 2.98658557, 1.48345856, 4.57733355, 1.64872808, 0.21725973, 0.37320698, 0.18225898, 1.35253554, 5.8365239, 1.94533493, 1.69962213])

    seed_features_mean.seed_pt=0.35833012; 
    seed_features_std.seed_pt=0.29620678;
    seed_features_mean.seed_eta=-0.000795975586; 
    seed_features_std.seed_eta=1.17231085;
    seed_features_mean.seed_phi=-0.0183847683; 
    seed_features_std.seed_phi=1.77591521;
    seed_features_mean.seed_mass=0.129604101; 
    seed_features_std.seed_mass=0.03431668;    
    seed_features_mean.seed_dz=-0.00631919497; 
    seed_features_std.seed_dz=4.34793758;
    seed_features_mean.seed_dxy=-0.0519152703; 
    seed_features_std.seed_dxy=4.90624525;
    seed_features_mean.seed_3D_ip=-3.65962456; 
    seed_features_std.seed_3D_ip=1.68269796;
    seed_features_mean.seed_3D_sip=1.71627243; 
    seed_features_std.seed_3D_sip=1.3829665;
    seed_features_mean.seed_2D_ip=-4.48109269; 
    seed_features_std.seed_2D_ip=1.99893739;
    seed_features_mean.seed_2D_sip=0.932680501; 
    seed_features_std.seed_2D_sip=1.6451007;    
    seed_features_mean.seed_3D_signedIp=-2.70270809;  
    seed_features_std.seed_3D_signedIp=2.98658557;
    seed_features_mean.seed_3D_signedSip=1.63019572; 
    seed_features_std.seed_3D_signedSip=1.48345856;
    seed_features_mean.seed_2D_signedIp=-1.76747277; 
    seed_features_std.seed_2D_signedIp=4.57733355;
    seed_features_mean.seed_2D_signedSip=0.926253178; 
    seed_features_std.seed_2D_signedSip=1.64872808;  
    seed_features_mean.seed_3D_TrackProbability=0.; 
    seed_features_std.seed_3D_TrackProbability=1.;
    seed_features_mean.seed_2D_TrackProbability=0.; 
    seed_features_std.seed_2D_TrackProbability=1.;
    seed_features_mean.seed_chi2reduced=0.037803297; 
    seed_features_std.seed_chi2reduced=0.18225898;
    seed_features_mean.seed_nPixelHits=3.9065406; 
    seed_features_std.seed_nPixelHits=1.35253554;
    seed_features_mean.seed_nHits=15.1244042; 
    seed_features_std.seed_nHits=5.8365239;
    seed_features_mean.seed_jetAxisDistance=-4.68685508; 
    seed_features_std.seed_jetAxisDistance=1.94533493;
    seed_features_mean.seed_jetAxisDlength=-1.47584773; 
    seed_features_std.seed_jetAxisDlength=1.69962213;


    float* ptr = &tensor.tensor<float, 3>()(jet_n, seed_n, 0);    
    
//     *ptr     = seed_features.seed_pt;
//     *(++ptr) = seed_features.seed_eta;
//     *(++ptr) = seed_features.seed_phi;
//     *(++ptr) = seed_features.seed_mass;    
//     *(++ptr) = seed_features.seed_dz;
//     *(++ptr) = seed_features.seed_dxy;
//     *(++ptr) = seed_features.seed_3D_ip;
//     *(++ptr) = seed_features.seed_3D_sip;
//     *(++ptr) = seed_features.seed_2D_ip;
//     *(++ptr) = seed_features.seed_2D_sip;    
//     *(++ptr) = seed_features.seed_3D_signedIp;
//     *(++ptr) = seed_features.seed_3D_signedSip;
//     *(++ptr) = seed_features.seed_2D_signedIp;
//     *(++ptr) = seed_features.seed_2D_signedSip;  
//     *(++ptr) = seed_features.seed_3D_TrackProbability;
//     *(++ptr) = seed_features.seed_2D_TrackProbability;
//     *(++ptr) = seed_features.seed_chi2reduced;
//     *(++ptr) = seed_features.seed_nPixelHits;
//     *(++ptr) = seed_features.seed_nHits;
//     *(++ptr) = seed_features.seed_jetAxisDistance;
//     *(++ptr) = seed_features.seed_jetAxisDlength;
    
        *ptr     = (seed_features.seed_pt!=0)*(seed_features.seed_pt-seed_features_mean.seed_pt)/seed_features_std.seed_pt;
        *(++ptr) = (seed_features.seed_eta!=0)*(seed_features.seed_eta-seed_features_mean.seed_eta)/seed_features_std.seed_eta;
        *(++ptr) = (seed_features.seed_phi!=0)*(seed_features.seed_phi-seed_features_mean.seed_phi)/seed_features_std.seed_phi;
        *(++ptr) = (seed_features.seed_mass!=0)*(seed_features.seed_mass-seed_features_mean.seed_mass)/seed_features_std.seed_mass;
        *(++ptr) = (seed_features.seed_dz!=0)*(seed_features.seed_dz-seed_features_mean.seed_dz)/seed_features_std.seed_dz;
        *(++ptr) = (seed_features.seed_dxy!=0)*(seed_features.seed_dxy-seed_features_mean.seed_dxy)/seed_features_std.seed_dxy;
        *(++ptr) = (seed_features.seed_3D_ip!=0)*(seed_features.seed_3D_ip-seed_features_mean.seed_3D_ip)/seed_features_std.seed_3D_ip;
        *(++ptr) = (seed_features.seed_3D_sip!=0)*(seed_features.seed_3D_sip-seed_features_mean.seed_3D_sip)/seed_features_std.seed_3D_sip;
        *(++ptr) = (seed_features.seed_2D_ip!=0)*(seed_features.seed_2D_ip-seed_features_mean.seed_2D_ip)/seed_features_std.seed_2D_ip;
        *(++ptr) = (seed_features.seed_2D_sip!=0)*(seed_features.seed_2D_sip-seed_features_mean.seed_2D_sip)/seed_features_std.seed_2D_sip;    
        *(++ptr) = (seed_features.seed_3D_signedIp!=0)*(seed_features.seed_3D_signedIp-seed_features_mean.seed_3D_signedIp)/seed_features_std.seed_3D_signedIp;
        *(++ptr) = (seed_features.seed_3D_signedSip!=0)*(seed_features.seed_3D_signedSip-seed_features_mean.seed_3D_signedSip)/seed_features_std.seed_3D_signedSip;
        *(++ptr) = (seed_features.seed_2D_signedIp!=0)*(seed_features.seed_2D_signedIp-seed_features_mean.seed_2D_signedIp)/seed_features_std.seed_2D_signedIp;
        *(++ptr) = (seed_features.seed_2D_signedSip!=0)*(seed_features.seed_2D_signedSip-seed_features_mean.seed_2D_signedSip)/seed_features_std.seed_2D_signedSip; 
        *(++ptr) = (seed_features.seed_3D_TrackProbability!=0)*(seed_features.seed_3D_TrackProbability-seed_features_mean.seed_3D_TrackProbability)/seed_features_std.seed_3D_TrackProbability;
        *(++ptr) = (seed_features.seed_2D_TrackProbability!=0)*(seed_features.seed_2D_TrackProbability-seed_features_mean.seed_2D_TrackProbability)/seed_features_std.seed_2D_TrackProbability;
        *(++ptr) = (seed_features.seed_chi2reduced!=0)*(seed_features.seed_chi2reduced-seed_features_mean.seed_chi2reduced)/seed_features_std.seed_chi2reduced;
        *(++ptr) = (seed_features.seed_nPixelHits!=0)*(seed_features.seed_nPixelHits-seed_features_mean.seed_nPixelHits)/seed_features_std.seed_nPixelHits;
        *(++ptr) = (seed_features.seed_nHits!=0)*(seed_features.seed_nHits-seed_features_mean.seed_nHits)/seed_features_std.seed_nHits;
        *(++ptr) = (seed_features.seed_jetAxisDistance!=0)*(seed_features.seed_jetAxisDistance-seed_features_mean.seed_jetAxisDistance)/seed_features_std.seed_jetAxisDistance;
        *(++ptr) = (seed_features.seed_jetAxisDlength!=0)*(seed_features.seed_jetAxisDlength-seed_features_mean.seed_jetAxisDlength)/seed_features_std.seed_jetAxisDlength;


  }
  
  
inline  void neighbourTracks_tensor_filler(tensorflow::Tensor & tensor,
                          std::size_t jet_n,
                          std::size_t seed_n,
                          const btagbtvdeep::SeedingTrackFeatures & seed_features) {

    
    std::vector<btagbtvdeep::TrackPairFeatures> neighbourTracks_features = seed_features.seed_nearTracks;
    
    btagbtvdeep::TrackPairFeatures tp_features_mean;
    btagbtvdeep::TrackPairFeatures tp_features_std;                
    
//     tvars_mean_BIS=numpy.array([0.469070252, 0.00920111268, 0.0106877711, 0.0307052882, -0.010131558, 0.131154315, -4.92932753, 0.332683788, -5.61626292, -0.314876937, -7.51496875, -2.63608467, 0.0848869341, 0.160081443, -1.07783469, 0.000207277479, 0.000345547684, 0.00090279539, 0.0847763365, 0.160232948, -1.07809807, 0.000270699497, 0.000152446623, 8.16506509e-05, 0.183894253, 0.196232466, 0.168856898, 0.278172289, 0.957906411, 0.958911898, -2.70696144, -2.70802778, -3.85457875, 0.751149045, -0.949525003, 1.20534152])
//     tvars_std_BIS=numpy.array([0.32821956, 1.23755055, 1.7664503, 5.3951208, 5.92355635, 0.03247655, 1.53365801, 1.14905112, 1.88311619, 1.48372824, 2.4413915, 1.85501459, 0.77155018, 0.63087902, 5.29532305, 0.0060867, 0.01456222, 0.03966832, 0.77167966, 0.63068495, 5.29534954, 0.02733926, 0.01348257, 0.00613393, 0.86610864, 0.72329339, 0.77200299, 0.6750665, 0.19763346, 0.1966655, 2.00225683, 2.00257068, 1.68567486, 0.43145278, 1.87436868, 1.08686211])

    
    tp_features_mean.nearTracks_pt=0.469070252; 
    tp_features_std.nearTracks_pt=0.32821956;
    tp_features_mean.nearTracks_eta=0.00920111268; 
    tp_features_std.nearTracks_eta= 1.23755055;
    tp_features_mean.nearTracks_phi=0.0106877711; 
    tp_features_std.nearTracks_phi=1.7664503;
    tp_features_mean.nearTracks_dz=0.0307052882; 
    tp_features_std.nearTracks_dz=5.3951208;
    tp_features_mean.nearTracks_dxy=-0.010131558; 
    tp_features_std.nearTracks_dxy=5.92355635;
    tp_features_mean.nearTracks_mass=0.131154315; 
    tp_features_std.nearTracks_mass=0.03247655;
    tp_features_mean.nearTracks_3D_ip=-4.92932753; 
    tp_features_std.nearTracks_3D_ip=1.53365801;
    tp_features_mean.nearTracks_3D_sip=0.332683788; 
    tp_features_std.nearTracks_3D_sip=1.14905112;
    tp_features_mean.nearTracks_2D_ip=-5.61626292; 
    tp_features_std.nearTracks_2D_ip=1.88311619;
    tp_features_mean.nearTracks_2D_sip=-0.314876937; 
    tp_features_std.nearTracks_2D_sip=1.48372824;
    tp_features_mean.nearTracks_PCAdist=-7.51496875; 
    tp_features_std.nearTracks_PCAdist=2.4413915;
    tp_features_mean.nearTracks_PCAdsig=-2.63608467; 
    tp_features_std.nearTracks_PCAdsig=1.85501459;
    tp_features_mean.nearTracks_PCAonSeed_x=0.0848869341; 
    tp_features_std.nearTracks_PCAonSeed_x=0.77155018;
    tp_features_mean.nearTracks_PCAonSeed_y=0.160081443; 
    tp_features_std.nearTracks_PCAonSeed_y=0.63087902;
    tp_features_mean.nearTracks_PCAonSeed_z=-1.07783469; 
    tp_features_std.nearTracks_PCAonSeed_z=5.29532305;
    tp_features_mean.nearTracks_PCAonSeed_xerr=0.000207277479; 
    tp_features_std.nearTracks_PCAonSeed_xerr=0.0060867;
    tp_features_mean.nearTracks_PCAonSeed_yerr=0.000345547684; 
    tp_features_std.nearTracks_PCAonSeed_yerr=0.01456222;
    tp_features_mean.nearTracks_PCAonSeed_zerr=0.00090279539; 
    tp_features_std.nearTracks_PCAonSeed_zerr=0.03966832;
    tp_features_mean.nearTracks_PCAonTrack_x=0.0847763365; 
    tp_features_std.nearTracks_PCAonTrack_x=0.77167966;
    tp_features_mean.nearTracks_PCAonTrack_y=0.160232948; 
    tp_features_std.nearTracks_PCAonTrack_y=0.63068495;
    tp_features_mean.nearTracks_PCAonTrack_z=-1.07809807; 
    tp_features_std.nearTracks_PCAonTrack_z=5.29534954;
    tp_features_mean.nearTracks_PCAonTrack_xerr=0.000270699497; 
    tp_features_std.nearTracks_PCAonTrack_xerr=0.02733926;
    tp_features_mean.nearTracks_PCAonTrack_yerr=0.000152446623; 
    tp_features_std.nearTracks_PCAonTrack_yerr=0.01348257;
    tp_features_mean.nearTracks_PCAonTrack_zerr=8.16506509e-05; 
    tp_features_std.nearTracks_PCAonTrack_zerr=0.00613393;
    tp_features_mean.nearTracks_dotprodTrack=0.183894253; 
    tp_features_std.nearTracks_dotprodTrack=0.86610864;
    tp_features_mean.nearTracks_dotprodSeed=0.196232466; 
    tp_features_std.nearTracks_dotprodSeed=0.72329339;
    tp_features_mean.nearTracks_dotprodTrackSeed2D=0.168856898; 
    tp_features_std.nearTracks_dotprodTrackSeed2D=0.77200299;
    tp_features_mean.nearTracks_dotprodTrackSeed3D=0.278172289; 
    tp_features_std.nearTracks_dotprodTrackSeed3D=0.6750665;
    tp_features_mean.nearTracks_dotprodTrackSeedVectors3D=0.957906411; 
    tp_features_std.nearTracks_dotprodTrackSeedVectors3D=0.19763346;
    tp_features_mean.nearTracks_dotprodTrackSeedVectors2D=0.958911898; 
    tp_features_std.nearTracks_dotprodTrackSeedVectors2D=0.1966655;
    tp_features_mean.nearTracks_PCAonSeed_pvd=-2.70696144; 
    tp_features_std.nearTracks_PCAonSeed_pvd=2.00225683;
    tp_features_mean.nearTracks_PCAonTrack_pvd=-2.70802778; 
    tp_features_std.nearTracks_PCAonTrack_pvd=2.00257068;
    tp_features_mean.nearTracks_PCAjetAxis_dist=-3.85457875; 
    tp_features_std.nearTracks_PCAjetAxis_dist=1.68567486;
    tp_features_mean.nearTracks_PCAjetMomenta_dotprod=0.751149045; 
    tp_features_std.nearTracks_PCAjetMomenta_dotprod=0.43145278;
    tp_features_mean.nearTracks_PCAjetDirs_DEta=-0.949525003; 
    tp_features_std.nearTracks_PCAjetDirs_DEta=1.87436868;
    tp_features_mean.nearTracks_PCAjetDirs_DPhi=1.20534152; 
    tp_features_std.nearTracks_PCAjetDirs_DPhi=1.08686211;
        
    
    for(unsigned int t_i=0; t_i<neighbourTracks_features.size(); t_i++){  

    float* ptr = &tensor.tensor<float, 3>()(jet_n, t_i,  0);    
        
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_pt;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_eta;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_phi;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_mass;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dz;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dxy;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_3D_ip;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_3D_sip;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_2D_ip;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_2D_sip;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAdist;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAdsig;      
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_x;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_y;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_z;      
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_xerr;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_yerr;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_zerr;      
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_x;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_y;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_z;      
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_xerr;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_yerr;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_zerr; 
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodTrack;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodSeed;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed2D;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed3D;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors2D;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors3D;      
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonSeed_pvd;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAonTrack_pvd;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAjetAxis_dist;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAjetMomenta_dotprod;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DEta;
//     *(++ptr) = neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DPhi;    
    
            
    *ptr = (neighbourTracks_features[t_i].nearTracks_pt!=0)*(neighbourTracks_features[t_i].nearTracks_pt-tp_features_mean.nearTracks_pt)/tp_features_std.nearTracks_pt;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_eta!=0)*(neighbourTracks_features[t_i].nearTracks_eta-tp_features_mean.nearTracks_eta)/tp_features_std.nearTracks_eta;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_phi!=0)*(neighbourTracks_features[t_i].nearTracks_phi-tp_features_mean.nearTracks_phi)/tp_features_std.nearTracks_phi;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dz!=0)*(neighbourTracks_features[t_i].nearTracks_dz-tp_features_mean.nearTracks_dz)/tp_features_std.nearTracks_dz;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dxy!=0)*(neighbourTracks_features[t_i].nearTracks_dxy-tp_features_mean.nearTracks_dxy)/tp_features_std.nearTracks_dxy;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_mass!=0)*(neighbourTracks_features[t_i].nearTracks_mass-tp_features_mean.nearTracks_mass)/tp_features_std.nearTracks_mass;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_3D_ip!=0)*(neighbourTracks_features[t_i].nearTracks_3D_ip-tp_features_mean.nearTracks_3D_ip)/tp_features_std.nearTracks_3D_ip;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_3D_sip!=0)*(neighbourTracks_features[t_i].nearTracks_3D_sip-tp_features_mean.nearTracks_3D_sip)/tp_features_std.nearTracks_3D_sip;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_2D_ip!=0)*(neighbourTracks_features[t_i].nearTracks_2D_ip-tp_features_mean.nearTracks_2D_ip)/tp_features_std.nearTracks_2D_ip;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_2D_sip!=0)*(neighbourTracks_features[t_i].nearTracks_2D_sip-tp_features_mean.nearTracks_2D_sip)/tp_features_std.nearTracks_2D_sip;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAdist!=0)*(neighbourTracks_features[t_i].nearTracks_PCAdist-tp_features_mean.nearTracks_PCAdist)/tp_features_std.nearTracks_PCAdist;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAdsig!=0)*(neighbourTracks_features[t_i].nearTracks_PCAdsig-tp_features_mean.nearTracks_PCAdsig)/tp_features_std.nearTracks_PCAdsig;    
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_x!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_x-tp_features_mean.nearTracks_PCAonSeed_x)/tp_features_std.nearTracks_PCAonSeed_x;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_y!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_y-tp_features_mean.nearTracks_PCAonSeed_y)/tp_features_std.nearTracks_PCAonSeed_y;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_z!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_z-tp_features_mean.nearTracks_PCAonSeed_z)/tp_features_std.nearTracks_PCAonSeed_z;    
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_xerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_xerr-tp_features_mean.nearTracks_PCAonSeed_xerr)/tp_features_std.nearTracks_PCAonSeed_xerr;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_yerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_yerr-tp_features_mean.nearTracks_PCAonSeed_yerr)/tp_features_std.nearTracks_PCAonSeed_yerr;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_zerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_zerr-tp_features_mean.nearTracks_PCAonSeed_zerr)/tp_features_std.nearTracks_PCAonSeed_zerr;   
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_x!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_x-tp_features_mean.nearTracks_PCAonTrack_x)/tp_features_std.nearTracks_PCAonTrack_x;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_y!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_y-tp_features_mean.nearTracks_PCAonTrack_y)/tp_features_std.nearTracks_PCAonTrack_y;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_z!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_z-tp_features_mean.nearTracks_PCAonTrack_z)/tp_features_std.nearTracks_PCAonTrack_z;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_xerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_xerr-tp_features_mean.nearTracks_PCAonTrack_xerr)/tp_features_std.nearTracks_PCAonTrack_xerr;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_yerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_yerr-tp_features_mean.nearTracks_PCAonTrack_yerr)/tp_features_std.nearTracks_PCAonTrack_yerr;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_zerr!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_zerr-tp_features_mean.nearTracks_PCAonTrack_zerr)/tp_features_std.nearTracks_PCAonTrack_zerr;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrack!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrack-tp_features_mean.nearTracks_dotprodTrack)/tp_features_std.nearTracks_dotprodTrack;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodSeed!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodSeed-tp_features_mean.nearTracks_dotprodSeed)/tp_features_std.nearTracks_dotprodSeed;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed2D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed2D-tp_features_mean.nearTracks_dotprodTrackSeed2D)/tp_features_std.nearTracks_dotprodTrackSeed2D;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed3D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeed3D-tp_features_mean.nearTracks_dotprodTrackSeed3D)/tp_features_std.nearTracks_dotprodTrackSeed3D;
    
    /*    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors2D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors2D-tp_features_mean.nearTracks_dotprodTrackSeedVectors2D)/tp_features_std.nearTracks_dotprodTrackSeedVectors2D;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors3D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors3D-tp_features_mean.nearTracks_dotprodTrackSeedVectors3D)/tp_features_std.nearTracks_dotprodTrackSeedVectors3D;   */ 
    
    //forse va riscambiato    
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors3D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors3D-tp_features_mean.nearTracks_dotprodTrackSeedVectors3D)/tp_features_std.nearTracks_dotprodTrackSeedVectors3D;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors2D!=0)*(neighbourTracks_features[t_i].nearTracks_dotprodTrackSeedVectors2D-tp_features_mean.nearTracks_dotprodTrackSeedVectors2D)/tp_features_std.nearTracks_dotprodTrackSeedVectors2D; 
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonSeed_pvd!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonSeed_pvd-tp_features_mean.nearTracks_PCAonSeed_pvd)/tp_features_std.nearTracks_PCAonSeed_pvd;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAonTrack_pvd!=0)*(neighbourTracks_features[t_i].nearTracks_PCAonTrack_pvd-tp_features_mean.nearTracks_PCAonTrack_pvd)/tp_features_std.nearTracks_PCAonTrack_pvd;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAjetAxis_dist!=0)*(neighbourTracks_features[t_i].nearTracks_PCAjetAxis_dist-tp_features_mean.nearTracks_PCAjetAxis_dist)/tp_features_std.nearTracks_PCAjetAxis_dist;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAjetMomenta_dotprod!=0)*(neighbourTracks_features[t_i].nearTracks_PCAjetMomenta_dotprod-tp_features_mean.nearTracks_PCAjetMomenta_dotprod)/tp_features_std.nearTracks_PCAjetMomenta_dotprod;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DEta!=0)*(neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DEta-tp_features_mean.nearTracks_PCAjetDirs_DEta)/tp_features_std.nearTracks_PCAjetDirs_DEta;
    *(++ptr) = (neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DPhi!=0)*(neighbourTracks_features[t_i].nearTracks_PCAjetDirs_DPhi-tp_features_mean.nearTracks_PCAjetDirs_DPhi)/tp_features_std.nearTracks_PCAjetDirs_DPhi;

    }


  }
  

}

#endif
