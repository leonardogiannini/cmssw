#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "RecoBTag/ONNXRuntime/interface/tensor_fillers.h"

using namespace cms::Ort;

class DeepCombinedONNXJetTagsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit DeepCombinedONNXJetTagsProducer(const edm::ParameterSet&, const ONNXRuntime*);
  ~DeepCombinedONNXJetTagsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ONNXRuntime*);

private:
  typedef std::vector<reco::DeepFlavourTagInfo> TagInfoCollection;
  typedef reco::JetTagCollection JetTagCollection;

  void produce(edm::Event&, const edm::EventSetup&) override;

  void make_inputs(unsigned i_jet, const reco::DeepFlavourTagInfo& taginfo);

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> flav_names_;
  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;

  const double min_jet_pt_;
  const double max_jet_eta_;

  enum InputIndexes {
    kGlobal = 0,
    kChargedCandidates = 1,
    kNeutralCandidates = 2,
    kVertices = 3,
    kGlobal1 = 4,
    kSeedingTracks = 5,
    kNeighbourTracks = 6
  };
  constexpr static unsigned n_features_global_ = 15;
  constexpr static unsigned n_cpf_ = 25;
  constexpr static unsigned n_features_cpf_ = 16;
  constexpr static unsigned n_npf_ = 25;
  constexpr static unsigned n_features_npf_ = 6;
  constexpr static unsigned n_sv_ = 4;
  constexpr static unsigned n_features_sv_ = 12;
  constexpr static unsigned n_features_global1_ = 4;
  constexpr static unsigned n_seed_ = 10;
  constexpr static unsigned n_features_seed_ = 21;
  constexpr static unsigned n_neighbor_ = 20;
  constexpr static unsigned n_features_neighbor_ = 36;

  const static std::vector<unsigned> input_sizes_;

  // hold the input data
  FloatArrays data_;
};

const std::vector<unsigned> DeepCombinedONNXJetTagsProducer::input_sizes_{n_features_global_,
                                                                          n_cpf_* n_features_cpf_,
                                                                          n_npf_* n_features_npf_,
                                                                          n_sv_* n_features_sv_,
                                                                          n_features_global1_,
                                                                          n_seed_* n_features_seed_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_,
                                                                          n_neighbor_* n_features_neighbor_};

DeepCombinedONNXJetTagsProducer::DeepCombinedONNXJetTagsProducer(const edm::ParameterSet& iConfig,
                                                                 const ONNXRuntime* cache)
    : src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
      input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")) {
  // get output names from flav_names
  for (const auto& flav_name : flav_names_) {
    produces<JetTagCollection>(flav_name);
  }

  assert(input_names_.size() == input_sizes_.size());
}

DeepCombinedONNXJetTagsProducer::~DeepCombinedONNXJetTagsProducer() {}

void DeepCombinedONNXJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfDeepFlavourJetTags
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfDeepFlavourTagInfos"));
  desc.add<std::vector<std::string>>("input_names",
                                     {"input_1_DFla",
                                      "input_2_DFla",
                                      "input_3_DFla",
                                      "input_4_DFla",
                                      "input_1",
                                      "input_2",
                                      "input_3",
                                      "input_4",
                                      "input_5",
                                      "input_6",
                                      "input_7",
                                      "input_8",
                                      "input_9",
                                      "input_10",
                                      "input_11",
                                      "input_12"});
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoBTag/Combined/data/DeepVertex/phase1_deepvertexcombined.onnx"));
  desc.add<std::vector<std::string>>("output_names", {"dense_13"});
  desc.add<std::vector<std::string>>("flav_names", std::vector<std::string>{"probb", "probc", "probuds", "probg"});
  desc.add<double>("min_jet_pt", 15.0);
  desc.add<double>("max_jet_eta", 2.5);

  descriptions.add("pfDeepCombinedJetTags", desc);
}

std::unique_ptr<ONNXRuntime> DeepCombinedONNXJetTagsProducer::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void DeepCombinedONNXJetTagsProducer::globalEndJob(const ONNXRuntime* cache) {}

void DeepCombinedONNXJetTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  data_.clear();

  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  if (!tag_infos->empty()) {
    unsigned int good_taginfo_count = 0;
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto& jet_ref = (*tag_infos)[jet_n].jet();
      if (jet_ref->pt() > min_jet_pt_ && std::fabs(jet_ref->eta()) < max_jet_eta_)
        good_taginfo_count++;
    }

    // init data storage w correct size
    for (const auto& len : input_sizes_) {
      data_.emplace_back(good_taginfo_count * len, 0);
    }

    // initialize output collection
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }

    // convert inputs
    unsigned int inputs_done_count = 0;
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto& jet_ref = (*tag_infos)[jet_n].jet();
      if (jet_ref->pt() > min_jet_pt_ && std::fabs(jet_ref->eta()) < max_jet_eta_) {
        const auto& taginfo = (*tag_infos)[jet_n];
        make_inputs(inputs_done_count, taginfo);
        inputs_done_count++;
      }
    }

    // run prediction
    assert(inputs_done_count == good_taginfo_count);
    const auto outputs = globalCache()->run(input_names_, data_, {}, output_names_, good_taginfo_count)[0];
    assert(outputs.size() == flav_names_.size() * good_taginfo_count);

    // get the outputs
    unsigned i_output = 0;
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto& jet_ref = (*tag_infos)[jet_n].jet();
      for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
        if (jet_ref->pt() > min_jet_pt_ && std::fabs(jet_ref->eta()) < max_jet_eta_) {
          (*(output_tags[flav_n]))[jet_ref] = outputs[i_output];
          ++i_output;
        } else
          (*(output_tags[flav_n]))[jet_ref] = -2;
      }
    }
  } else {
    // create empty output collection
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }

  // put into the event
  for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
  }
}

void DeepCombinedONNXJetTagsProducer::make_inputs(unsigned i_jet, const reco::DeepFlavourTagInfo& taginfo) {
  const auto& features = taginfo.features();
  float* ptr = nullptr;
  const float* start = nullptr;
  unsigned offset = 0;

  // jet and other global features
  offset = i_jet * input_sizes_[kGlobal];
  ptr = &data_[kGlobal][offset];
  // jet variables
  const auto& jet_features = features.jet_features;
  start = ptr;
  *ptr = jet_features.pt;
  *(++ptr) = jet_features.eta;
  // number of elements in different collections
  *(++ptr) = features.c_pf_features.size();
  *(++ptr) = features.n_pf_features.size();
  *(++ptr) = features.sv_features.size();
  *(++ptr) = features.npv;
  // variables from ShallowTagInfo
  const auto& tag_info_features = features.tag_info_features;
  *(++ptr) = tag_info_features.trackSumJetEtRatio;
  *(++ptr) = tag_info_features.trackSumJetDeltaR;
  *(++ptr) = tag_info_features.vertexCategory;
  *(++ptr) = tag_info_features.trackSip2dValAboveCharm;
  *(++ptr) = tag_info_features.trackSip2dSigAboveCharm;
  *(++ptr) = tag_info_features.trackSip3dValAboveCharm;
  *(++ptr) = tag_info_features.trackSip3dSigAboveCharm;
  *(++ptr) = tag_info_features.jetNSelectedTracks;
  *(++ptr) = tag_info_features.jetNTracksEtaRel;
  assert(start + n_features_global_ - 1 == ptr);

  // c_pf candidates
  auto max_c_pf_n = std::min(features.c_pf_features.size(), (std::size_t)n_cpf_);
  offset = i_jet * input_sizes_[kChargedCandidates];
  for (std::size_t c_pf_n = 0; c_pf_n < max_c_pf_n; c_pf_n++) {
    const auto& c_pf_features = features.c_pf_features[c_pf_n];
    ptr = &data_[kChargedCandidates][offset + c_pf_n * n_features_cpf_];
    start = ptr;
    *ptr = c_pf_features.btagPf_trackEtaRel;
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
    assert(start + n_features_cpf_ - 1 == ptr);
  }

  // n_pf candidates
  auto max_n_pf_n = std::min(features.n_pf_features.size(), (std::size_t)n_npf_);
  offset = i_jet * input_sizes_[kNeutralCandidates];
  for (std::size_t n_pf_n = 0; n_pf_n < max_n_pf_n; n_pf_n++) {
    const auto& n_pf_features = features.n_pf_features[n_pf_n];
    ptr = &data_[kNeutralCandidates][offset + n_pf_n * n_features_npf_];
    start = ptr;
    *ptr = n_pf_features.ptrel;
    *(++ptr) = n_pf_features.deltaR;
    *(++ptr) = n_pf_features.isGamma;
    *(++ptr) = n_pf_features.hadFrac;
    *(++ptr) = n_pf_features.drminsv;
    *(++ptr) = n_pf_features.puppiw;
    assert(start + n_features_npf_ - 1 == ptr);
  }

  // sv candidates
  auto max_sv_n = std::min(features.sv_features.size(), (std::size_t)n_sv_);
  offset = i_jet * input_sizes_[kVertices];
  for (std::size_t sv_n = 0; sv_n < max_sv_n; sv_n++) {
    const auto& sv_features = features.sv_features[sv_n];
    ptr = &data_[kVertices][offset + sv_n * n_features_sv_];
    start = ptr;
    *ptr = sv_features.pt;
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
    assert(start + n_features_sv_ - 1 == ptr);
  }

  // jet and other global features
  offset = i_jet * input_sizes_[kGlobal1];
  ptr = &data_[kGlobal1][offset];
  start = ptr;
  jet4vec_tensor_filler(ptr, jet_features);
  assert(start + n_features_global1_ - 1 == ptr);

  // seeds
  auto max_seed_n = std::min(features.seed_features.size(), (std::size_t)n_seed_);
  offset = i_jet * input_sizes_[kSeedingTracks];
  for (std::size_t seed_n = 0; seed_n < max_seed_n; seed_n++) {
    const auto& seed_features = features.seed_features[seed_n];
    ptr = &data_[kSeedingTracks][offset + seed_n * n_features_seed_];
    start = ptr;
    seedTrack_tensor_filler(ptr, seed_features);
    assert(start + n_features_seed_ - 1 == ptr);
  }

  // neighbours
  offset = i_jet * input_sizes_[kNeighbourTracks];
  for (std::size_t seed_n = 0; seed_n < max_seed_n; seed_n++) {
    const auto& neighbourTracks_features = features.seed_features[seed_n].nearTracks;
    auto max_neighbour_n = std::min(neighbourTracks_features.size(), (std::size_t)n_neighbor_);
    for (std::size_t neighbour_n = 0; neighbour_n < max_neighbour_n; neighbour_n++) {
      ptr = &data_[kNeighbourTracks + seed_n][offset + neighbour_n * n_features_neighbor_];
      start = ptr;
      neighbourTrack_tensor_filler(ptr, neighbourTracks_features[neighbour_n]);
      assert(start + n_features_neighbor_ - 1 == ptr);
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepCombinedONNXJetTagsProducer);
